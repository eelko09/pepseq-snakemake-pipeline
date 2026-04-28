import numpy as np
import os
import pandas as pd
import rpy2.robjects as ro
import q2_PSEA.actions.splines as splines
import qiime2
import q2_PSEA.utils as utils
import tempfile

import time
import concurrent.futures
import multiprocessing

from math import isnan, log, pow
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from q2_pepsirf.format_types import PepsirfContingencyTSVFormat
from q2_PSEA.actions.r_functions import INTERNAL


pandas2ri.activate()


def make_psea_table(
        ctx,
        scores_file,
        pairs_file,
        peptide_sets_file,
        threshold,
        epitope_map=None,
        p_val_thresh=0.05,
        nes_thresh=1,
        species_taxa_file="",
        species_color_file="",
        min_size=15,
        max_size=2000,
        permutation_num=10000,  # as per original PSEA code
        spline_type="r-smooth",
        degree=3,
        dof=None,
        table_dir="./psea_table_outdir",
        pepsirf_binary="pepsirf",
        iterative_analysis=True,
        iter_tables_dir="",
        max_workers=None,
        summary_tables_dir="./psea_ae_summary_tables",
        vis_outputs_dir=None,
        seed=149
):    
    start_time = time.perf_counter()

    volcano = ctx.get_action("ps-plot", "volcano")
    zscatter = ctx.get_action("ps-plot", "zscatter")
    aeplots = ctx.get_action("ps-plot", "aeplots")

    assert spline_type in splines.SPLINE_TYPES, \
        f"'{spline_type}' is not a valid spline method!"
    assert not os.path.exists(table_dir), \
        f"'{table_dir}' already exists! Please move or remove this directory."
    assert not os.path.exists(summary_tables_dir), \
        f"'{summary_tables_dir}' already exists! Please move or remove this directory."
    if iterative_analysis:
        assert ".gmt" in peptide_sets_file.lower(), \
            "You are running iterative analysis without a GMT peptide sets file."
    if max_workers != None:
        assert max_workers <= multiprocessing.cpu_count(), \
            f"Max workers excedes {multiprocessing.cpu_count()}, the number of CPUs on your machine."
    if vis_outputs_dir != None:
        assert not os.path.exists(vis_outputs_dir), \
            f"'{vis_outputs_dir}' already exists! Please move or remove this directory."
        os.mkdir(vis_outputs_dir)

    os.mkdir(table_dir)
    os.mkdir(summary_tables_dir)

    pairs = list()
    pair_2_title = dict()
    with open(pairs_file, "r") as fh:
        # skip header line
        fh.readline()
        for line in fh.readlines():
            line_tup = tuple(line.replace("\n", "").split("\t"))
            pair = line_tup[0:2]
            pairs.append(pair)

            if( len(line_tup) > 2):
                pair_2_title[pair] = line_tup[2]
            else:
                pair_2_title[pair] = ""

    scores = pd.read_csv(scores_file, sep="\t", index_col=0)
    processed_scores = process_scores(scores, pairs)

    scores_file_split = scores_file.rsplit("/", 1)
    if len(scores_file_split) > 1:
        processed_scores_file = f"transformed_{scores_file_split[1]}"
    else:
        processed_scores_file = f"transformed_{scores_file_split[0]}"

    # temporary directory to hold iterative analysis tables
    with tempfile.TemporaryDirectory() as temp_peptide_sets_dir:
        if iterative_analysis:
            if iter_tables_dir:
                if not os.path.exists(iter_tables_dir):
                    os.mkdir(iter_tables_dir)
                else:
                    print(
                        f"\nWarning: the directory '{iter_tables_dir}' already exists; files may"
                        " be overwritten!"
                    )

            # run iterative peptide analysis, writes to temp_peptide_sets_dir files
            pair_pep_sets_file_dict = run_iterative_peptide_analysis(
                pairs=pairs,
                processed_scores=processed_scores,
                og_peptide_sets_file=peptide_sets_file,
                species_taxa_file=species_taxa_file,
                threshold=threshold,
                permutation_num=permutation_num,
                min_size=min_size,
                max_size=max_size,
                spline_type=spline_type,
                degree=degree,
                dof=dof,
                p_val_thresh=p_val_thresh,
                nes_thresh=nes_thresh,
                peptide_sets_out_dir = temp_peptide_sets_dir,
                iter_tables_dir=iter_tables_dir,
                max_workers=max_workers,
                seed=seed
            )
        else:
            # each pair will have the same gmt file
            pair_pep_sets_file_dict = dict.fromkeys(pairs, peptide_sets_file)

        with tempfile.TemporaryDirectory() as tempdir:
            pos_nes_count_dict = dict()
            neg_nes_count_dict = dict()
            zero_nes_count_dict = dict()

            processed_scores.to_csv(processed_scores_file, sep="\t")

            taxa_access = "species_name"
            # Stream spline data to disk to avoid holding all points in memory.
            spline_path = os.path.join(tempdir, "spline_data.tsv")
            spline_file = open(spline_path, "w", encoding="utf-8")
            spline_file.write("x\ty\tpair\n")


            if not dof:
                dof = ro.NULL
            if not species_taxa_file:
                taxa_access = "ID"
            
            PSEAScores = []
            # Stream per-pair significant taxa to disk to avoid large in-memory matrices.
            pos_events_path = os.path.join(tempdir, "pos_events.tsv")
            neg_events_path = os.path.join(tempdir, "neg_events.tsv")
            zero_events_path = os.path.join(tempdir, "zero_events.tsv")
            pos_events_fh = open(pos_events_path, "w", encoding="utf-8")
            neg_events_fh = open(neg_events_path, "w", encoding="utf-8")
            zero_events_fh = open(zero_events_path, "w", encoding="utf-8")
            pos_events_fh.write("pair\ttaxa\n")
            neg_events_fh.write("pair\ttaxa\n")
            zero_events_fh.write("pair\ttaxa\n")

            with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
                pair_futures = [executor.submit(create_fgsea_table_for_pair,
                                pair,
                                processed_scores,
                                pair_pep_sets_file_dict[ pair ],
                                species_taxa_file,
                                threshold,
                                permutation_num,
                                min_size,
                                max_size,
                                spline_type,
                                degree,
                                dof,
                                p_val_thresh,
                                nes_thresh,
                                False,
                                seed,
                                table_dir
                                ) for pair in pairs]

                for future in concurrent.futures.as_completed(pair_futures):
                    result = future.result()
                    x = result[0]
                    yfit = result[1]
                    table_prefix = result[2]
                    pair = result[3]

                    # Write spline points as they complete.
                    for xi, yi in zip(x.tolist(), yfit.tolist()):
                        spline_file.write(f"{xi}\t{yi}\t{table_prefix}\n")
                    # used_pairs.append((pair[0], pair[1], pair_2_title[pair]))

                    # populate event matrix with species that are significant this pair
                    tableDf = pd.read_csv(f"{table_dir}/{table_prefix}_psea_table.tsv", sep="\t")
                    # Persist per-pair scores as soon as the table is available.
                    art = ctx.make_artifact('FeatureData[PSEAScores]', tableDf)
                    PSEAScores.append(art)
                    pos_taxa = []
                    neg_taxa = []
                    zero_taxa = []
                    for i, row in tableDf.iterrows():
                        taxa = row[taxa_access]

                        if row["p.adjust"] < p_val_thresh and np.absolute(row["NES"]) > nes_thresh:
                            if row["NES"] > 0:
                                pos_taxa.append(taxa)
                                if taxa not in pos_nes_count_dict.keys():
                                    pos_nes_count_dict[taxa] = 0
                                pos_nes_count_dict[taxa] += 1

                            elif row["NES"] < 0:
                                neg_taxa.append(taxa)
                                if taxa not in neg_nes_count_dict.keys():
                                    neg_nes_count_dict[taxa] = 0
                                neg_nes_count_dict[taxa] += 1

                            # output message if nes is 0
                            else:
                                zero_taxa.append(taxa)
                                if taxa not in zero_nes_count_dict.keys():
                                    zero_nes_count_dict[taxa] = 0
                                zero_nes_count_dict[taxa] += 1

                    pos_events_fh.write(f"{table_prefix}\t{';'.join(pos_taxa)}\n")
                    neg_events_fh.write(f"{table_prefix}\t{';'.join(neg_taxa)}\n")
                    zero_events_fh.write(f"{table_prefix}\t{';'.join(zero_taxa)}\n")

            pos_events_fh.close()
            neg_events_fh.close()
            zero_events_fh.close()

            pos_taxa_order = [k for k, _ in sorted(pos_nes_count_dict.items(), key=lambda item: item[1], reverse=True)]
            neg_taxa_order = [k for k, _ in sorted(neg_nes_count_dict.items(), key=lambda item: item[1], reverse=True)]

            write_event_matrix(
                pos_events_path,
                os.path.join(summary_tables_dir, "Positive_NES_taxa_matrix.tsv"),
                pos_taxa_order,
            )
            write_event_matrix(
                neg_events_path,
                os.path.join(summary_tables_dir, "Negative_NES_taxa_matrix.tsv"),
                neg_taxa_order,
            )

            pos_nes_count_dict = {k:v for k, v in sorted(pos_nes_count_dict.items(), key=lambda item: item[1], reverse=True)}
            neg_nes_count_dict = {k:v for k, v in sorted(neg_nes_count_dict.items(), key=lambda item: item[1], reverse=True)
            }

            # create column sums for positive and negative NES
            with open(os.path.join(summary_tables_dir, "Positive_NES_AE.tsv"), "w") as pos_file:
                pos_file.write(f"Species\tEvents\n")
                for taxa in pos_nes_count_dict.keys():
                    pos_file.write(f"{taxa}\t{pos_nes_count_dict[taxa]}\n")
            with open(os.path.join( summary_tables_dir, "Negative_NES_AE.tsv"), "w") as neg_file:
                neg_file.write(f"Species\tEvents\n")
                for taxa in neg_nes_count_dict.keys():
                    neg_file.write(f"{taxa}\t{neg_nes_count_dict[taxa]}\n")
            if len(zero_nes_count_dict) > 0:
                print("\n")
                zero_pairs_by_taxa = read_event_pairs(zero_events_path)
                for taxa in zero_nes_count_dict.keys():
                    if zero_nes_count_dict[taxa] > 1:
                        print(f"{zero_nes_count_dict[taxa]} events for {taxa}, which has an NES of 0")
                    else:
                        print(f"{zero_nes_count_dict[taxa]} event for {taxa}, which has an NES of 0")
                    print("The pairs which this occurred are: ")
                    for pair_label in zero_pairs_by_taxa.get(taxa, []):
                        print(pair_label)
            '''
            pd.DataFrame(used_pairs).to_csv(
                f"{tempdir}/used_pairs.tsv", sep="\t",
                header=False, index=False
            )
            '''
            spline_file.close()

            processed_scores_art = ctx.make_artifact(
                type="FeatureTable[Zscore]",
                view=processed_scores_file,
                view_type=PepsirfContingencyTSVFormat
            )
        
            scatter_plot, = zscatter(
                zscores=processed_scores_art,
                pairs_file=pairs_file,
                spline_file=spline_path,
                p_val_access="p.adjust",
                le_peps_access="core_enrichment",
                taxa_access=taxa_access,
                highlight_data=table_dir,
                highlight_threshold=p_val_thresh,
                colors_file=species_color_file,
                vis_outputs_dir=vis_outputs_dir
            )

            volcano_plot, = volcano(
                xy_dir=table_dir,
                xy_access=["NES", "p.adjust"],
                taxa_access=taxa_access,
                x_threshold=nes_thresh,
                y_threshold=p_val_thresh,
                xy_labels=["Enrichment score", "Adjusted p-values"],
                pairs_file=pairs_file,
                colors_file=species_color_file,
                vis_outputs_dir=vis_outputs_dir
            )

            ae_plot, = aeplots( 
                pos_nes_ae_file=os.path.join(summary_tables_dir, "Positive_NES_AE.tsv"),
                neg_nes_ae_file=os.path.join(summary_tables_dir, "Negative_NES_AE.tsv"),
                xy_access=["Events", "Species"],
                xy_labels=["Number of AEs in cohort", "Species"],
                colors_file=species_color_file,
                vis_outputs_dir=vis_outputs_dir
            )

    end_time = time.perf_counter()

    print(f"\nFinished in {round(end_time-start_time, 2)} seconds")
    
    return scatter_plot, volcano_plot, ae_plot, PSEAScores

def _compute_pair_fit_and_residuals(processed_scores, pair, spline_type, degree, dof):
    data_sorted = processed_scores.loc[:, pair].sort_values(by=pair[0])
    x = data_sorted.loc[:, pair[0]].to_numpy()
    y = data_sorted.loc[:, pair[1]].to_numpy()

    if spline_type == "py-smooth":
        yfit = splines.smooth_spline(x, y)
    elif spline_type == "py-LinearGAM":
        yfit = splines.smooth_gam(x, y)
    elif spline_type == "cubic":
        yfit = splines.R_SPLINES.cubic_spline(x, y, degree, dof)
    else:
        yfit = splines.R_SPLINES.smooth_spline(x, y)

    vals = data_sorted.loc[:, pair].to_numpy()
    maxZ = pd.Series(np.maximum(vals[:, 0], vals[:, 1]), index=data_sorted.index)
    deltaZ = pd.Series(y - yfit, index=data_sorted.index)

    return x, yfit, maxZ, deltaZ

def create_fgsea_table_for_pair(
    pair,
    processed_scores,
    pep_sets_file,
    species_taxa_file,
    threshold,
    permutation_num,
    min_size,
    max_size,
    spline_type,
    degree,
    dof,
    p_val_thresh,
    nes_thresh,
    iteration,
    seed,
    table_dir="",
    pair_fit_cache=None,   # NEW
):
    print(f"Working on pair ({pair[0]}, {pair[1]})...")

    table_prefix = f"{pair[0]}~{pair[1]}"

    # each pair has a different peptide set file
    processed_scores, peptide_sets = utils.remove_peptides(processed_scores, pep_sets_file)

    if pair_fit_cache is None:
        x, yfit, maxZ, deltaZ = _compute_pair_fit_and_residuals(
            processed_scores, pair, spline_type, degree, dof
        )
    else:
        x, yfit, maxZ_all, deltaZ_all = pair_fit_cache
        idx = processed_scores.index
        maxZ = maxZ_all.reindex(idx)
        deltaZ = deltaZ_all.reindex(idx)

    table = INTERNAL.psea(
        maxZ,
        deltaZ,
        peptide_sets,
        species_taxa_file,
        threshold,
        permutation_num,
        min_size,
        max_size,
        seed
    )
    with (ro.default_converter + pandas2ri.converter).context():
        table = ro.conversion.get_conversion().rpy2py(table)

    if iteration:
        return table

    table.to_csv(f"{table_dir}/{table_prefix}_psea_table.tsv", sep="\t", index=False)
    return x, yfit, table_prefix, pair


def process_scores(scores, pairs):
    base, offset = 2, 3
    reps = pd.unique(np.array(pairs).ravel())
    out = scores.loc[:, reps].add(base ** offset).clip(lower=1)
    return np.log2(out) - offset


def run_iterative_peptide_analysis(
    pairs,
    processed_scores,
    og_peptide_sets_file,
    species_taxa_file,
    threshold,
    permutation_num,
    min_size,
    max_size,
    spline_type,
    degree,
    dof,
    p_val_thresh,
    nes_thresh,
    peptide_sets_out_dir,
    iter_tables_dir,
    max_workers,
    seed
    ) -> dict:

    iteration_num = 1
    
    # initialize gmt dict
    gmt_dict = dict()
    with open(og_peptide_sets_file, "r") as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip().split("\t")
            # species : set of peptides
            gmt_dict[ line[0] ] = set( line[2:] )

    # each pair should have its own gmt
    pair_gmt_dict = dict()
    # keep track of which pairs are fully processed
    sig_species_found_dict = dict()
    # keep a dict for each pair and it's tested species
    tested_species_dict = dict()
    for pair in pairs:
        pair_gmt_dict[pair] = {k: set(v) for k, v in gmt_dict.items()}
        sig_species_found_dict[ pair ] = True
        tested_species_dict[ pair ] = set()

    # keep a dict for the output gmt file of each pair
    pair_sets_filename_dict = dict()
    
    if not dof:
        dof = ro.NULL
    
    pair_fit_cache = {
    pair: _compute_pair_fit_and_residuals(processed_scores, pair, spline_type, degree, dof)
    for pair in pairs
    }

    # loop until no other significant peptides were found
    while any(sig_species_found_dict.values()):

        print("\nIteration:", iteration_num)

        if iter_tables_dir:
            iter_out_dir = f"{iter_tables_dir}/Iteration_{iteration_num}"
            if not os.path.exists(iter_out_dir):
                os.mkdir(iter_out_dir)
        else:
            iter_out_dir = ""

        # -------------------------------
        # note: rpy2 is not compatible with multithreading, only multiprocessing
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            pair_futures = [executor.submit(
                            run_iterative_process_single_pair,
                            pair,
                            tested_species_dict[pair],
                            pair_gmt_dict[pair],
                            processed_scores,
                            species_taxa_file,
                            threshold,
                            permutation_num,
                            min_size,
                            max_size,
                            spline_type,
                            degree,
                            dof,
                            p_val_thresh,
                            nes_thresh,
                            peptide_sets_out_dir,
                            iter_out_dir,
                            seed,
                            pair_fit_cache[pair],   # NEW
                            ) for pair in pairs if sig_species_found_dict[pair]]

            # concurrent.futures.wait(pair_futures, timeout=None, return_when=concurrent.futures.ALL_COMPLETED)

            for future in concurrent.futures.as_completed(pair_futures):
                res_tup = future.result()
                pair = res_tup[4]
                sig_species_found_dict[pair] = res_tup[0]
                pair_sets_filename_dict[pair] = res_tup[1]
                tested_species_dict[pair] = res_tup[2]
                pair_gmt_dict[pair]= res_tup[3]        
        # -------------------------------

        iteration_num += 1

    print("\nEnd of Iterative Peptide Analysis\n")
    return pair_sets_filename_dict


def run_iterative_process_single_pair(
    pair,
    tested_species,
    gmt_dict,
    processed_scores,
    species_taxa_file,
    threshold,
    permutation_num,
    min_size,
    max_size,
    spline_type,
    degree,
    dof,
    p_val_thresh,
    nes_thresh,
    peptide_sets_out_dir,
    iter_out_dir,
    seed,
    pair_fit_cache,   # NEW
    ):

    sig_species_found = False

    # create gmt file for this pair (filtered gmt from prev iteration)
    pair_sets_filename = f"{peptide_sets_out_dir}/{pair[0]}_{pair[1]}".replace(".","-") + ".gmt"

    write_gmt_from_dict(pair_sets_filename, gmt_dict)

    table = create_fgsea_table_for_pair(
        pair=pair,
        processed_scores=processed_scores,
        pep_sets_file=pair_sets_filename,
        species_taxa_file=species_taxa_file,
        threshold=threshold,
        permutation_num=permutation_num,
        min_size=min_size,
        max_size=max_size,
        spline_type=spline_type,
        degree=degree,
        dof=dof,
        p_val_thresh=p_val_thresh,
        nes_thresh=nes_thresh,
        iteration=True,
        seed=seed,
        pair_fit_cache=pair_fit_cache,   # NEW
    )
    # sort the table by ascending p-value (lowest on top)
    table = table.sort_values(by=["p.adjust"], ascending=True)

    if iter_out_dir:
        table.to_csv(f"{iter_out_dir}/{pair}.tsv", sep="\t")

    # iterate through each row
    for index, row in table.iterrows():

        # test for significant species that has not already been used for this pair
        if row["p.adjust"] < p_val_thresh and np.absolute(row["NES"]) > nes_thresh \
                                    and row["ID"] not in tested_species:

            print(f"Found {row['species_name']} in {pair} to be significant")

            # set sig species found to true for this pair
            sig_species_found = True

            tested_species.add(row["ID"])

            # take out all tested peptides from peptide set in gmt for all other species in the gmt
            all_tested_peps = set(row["all_tested_peptides"].split("/"))
            for gmt_species in gmt_dict.keys():
                if gmt_species != row['ID']:
                    gmt_dict[ gmt_species ] = gmt_dict[ gmt_species ] - all_tested_peps

            # only get top significant species
            break

    return (sig_species_found, pair_sets_filename, tested_species, gmt_dict, pair)


def write_gmt_from_dict(outfile_name, gmt_dict)->None:
    with open( outfile_name, "w" ) as gmt_file:
        for species in gmt_dict.keys():
            gmt_file.write(f"{species}\t\t")

            for peptide in gmt_dict[ species ]:
                gmt_file.write(f"{peptide}\t")

            gmt_file.write("\n")


def write_event_matrix(events_path, output_path, taxa_order):
    """Write a pair-by-taxa event matrix from streamed event records."""
    with open(output_path, "w", encoding="utf-8") as out_f:
        out_f.write("Pair\t" + "\t".join(taxa_order) + "\n")
        with open(events_path, "r", encoding="utf-8") as in_f:
            next(in_f, None)  # header
            for line in in_f:
                line = line.rstrip("\n")
                if not line:
                    continue
                pair_label, taxa_str = line.split("\t", 1)
                taxa_set = set(filter(None, taxa_str.split(";")))
                row = ["1" if taxa in taxa_set else "0" for taxa in taxa_order]
                out_f.write(pair_label + "\t" + "\t".join(row) + "\n")


def read_event_pairs(events_path):
    """Return dict of taxa -> list of pair labels for non-empty event records."""
    taxa_to_pairs = {}
    with open(events_path, "r", encoding="utf-8") as in_f:
        next(in_f, None)  # header
        for line in in_f:
            line = line.rstrip("\n")
            if not line:
                continue
            pair_label, taxa_str = line.split("\t", 1)
            if not taxa_str:
                continue
            for taxa in taxa_str.split(";"):
                taxa_to_pairs.setdefault(taxa, []).append(pair_label)
    return taxa_to_pairs
