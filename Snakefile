import os
from pathlib import Path

# ------------------------------------------------------------
# Global paths and outputs
# ------------------------------------------------------------

OUTPUT_DIR = config.get("output_dir", "results")
QC_DIR = f"{OUTPUT_DIR}/qc"
LOG_DIR = f"{OUTPUT_DIR}/logs"
RUN_DIR = f"{OUTPUT_DIR}/run"
Q2_WORKDIR = config.get("qiime_workdir", workflow.basedir)

# ------------------------------------------------------------
# Inputs (required in config)
# ------------------------------------------------------------
RAW_COUNTS = config["raw_counts_tsv"]
SHORTENED_COUNTS = f"{QC_DIR}/00_shortened_counts.tsv"
SHORTEN_MAP = f"{QC_DIR}/00_column_name_map.tsv"
PAIRS_FILE = config["pairs_tsv"]

# ------------------------------------------------------------
# QC outputs (ordered by pipeline stage)
# ------------------------------------------------------------
RC_OUT = f"{QC_DIR}/01_rc_filtered.tsv"
RC_FAIL_OUT = f"{QC_DIR}/01_rc_failed_samples.tsv"
RC_TOTALS_OUT = f"{QC_DIR}/01_rc_total_readcounts.tsv"
# Temporary pass-list used to build RC filtered matrix via pepsirf subjoin.
RC_PASS_OUT = f"{QC_DIR}/01_rc_pass_samples.tsv"
ZERO_OUT = f"{QC_DIR}/02_zero_filtered.tsv"
ZERO_FAIL_OUT = f"{QC_DIR}/02_zero_failed_samples.tsv"
ZERO_ALL_COUNTS_OUT = f"{QC_DIR}/02_zero_all_counts.tsv"
ZERO_DOWNSAMPLED_COUNTS_OUT = f"{QC_DIR}/02_zero_downsampled_counts.tsv"
# Temporary pass-list used to build zero-count filtered matrix via pepsirf subjoin.
ZERO_PASS_OUT = f"{QC_DIR}/02_zero_pass_samples.tsv"
CORR_PAIRS_OUT = f"{QC_DIR}/03_pairs_pearson.tsv"
CORR_PAIRS_ALL_OUT = f"{QC_DIR}/03_pairs_pearson_all.tsv"
CORR_DATA_OUT = f"{QC_DIR}/03_correlation_filtered.tsv"
CORR_FAIL_OUT = f"{QC_DIR}/03_correlation_failed_samples.tsv"
# Temporary pass-list used to build correlation filtered matrix via pepsirf subjoin.
CORR_PASS_OUT = f"{QC_DIR}/03_correlation_pass_samples.tsv"
QC_FAILURE_SUMMARY_OUT = f"{QC_DIR}/04_qc_failed_samples_summary.tsv"
ZSCORE_CORR_PAIRS_OUT = f"{QC_DIR}/05_zscore_pairs_pearson.tsv"
ZSCORE_CORR_PAIRS_ALL_OUT = f"{QC_DIR}/05_zscore_pairs_pearson_all.tsv"
ZSCORE_CORR_DATA_OUT = f"{QC_DIR}/05_zscore_correlation_filtered.tsv"
ZSCORE_CORR_FAIL_OUT = f"{QC_DIR}/05_zscore_correlation_failed_samples.tsv"
# Temporary pass-list used to build Z-score correlation filtered matrix via pepsirf subjoin.
ZSCORE_CORR_PASS_OUT = f"{QC_DIR}/05_zscore_correlation_pass_samples.tsv"
RUN_TOKEN = f"{RUN_DIR}/.run_token.txt"

# ------------------------------------------------------------
# AutomateBins outputs
# ------------------------------------------------------------
AUTOBINS_PREFIX = config.get("automatebins_outprefix", "qc_filtered")
AUTOBINS_DIR = f"{OUTPUT_DIR}/bins"
AUTOBINS_FINAL = f"{AUTOBINS_DIR}/{AUTOBINS_PREFIX}_b300r1_bins.tsv"
AUTOBINS_NEG_CTRL = f"{AUTOBINS_DIR}/{AUTOBINS_PREFIX}_sblkOnly_CS_forBins.tsv"

# ------------------------------------------------------------
# AutoPepsirf / PSEA outputs and flags
# ------------------------------------------------------------
AUTOPEPSIRF_RUN_DIR = f"{OUTPUT_DIR}/autopepsirf"
RUN_AUTOPEPSIRF = config.get("run_autopepsirf", False)
RUN_PSEA = config.get("run_psea", False)
RUN_CORRELATION_QC = config.get("run_correlation_qc", True)
RUN_ZSCORE_CORRELATION_QC = config.get("run_zscore_correlation_qc", False)
FORCE_RERUN_WITH_TOKEN = config.get("force_rerun_with_token", False)
AUTOPEPSIRF_CONDA_PATH = config["autopepsirf_conda_path"]
PYTHON_SCRIPTS_DIR = config.get("python_scripts_dir", "python_scripts")
RULE_RESOURCES = config.get("rule_resources", {})

if RUN_PSEA and not RUN_AUTOPEPSIRF:
    raise ValueError("run_psea requires run_autopepsirf to be true.")
if RUN_ZSCORE_CORRELATION_QC and not RUN_AUTOPEPSIRF:
    raise ValueError("run_zscore_correlation_qc requires run_autopepsirf to be true.")

AUTOPEPSIRF_EXACT_Z_THRESH = str(config["autopepsirf_exact_z_thresh"])
AUTOPEPSIRF_EXACT_CS_THRESH = str(config["autopepsirf_exact_cs_thresh"])
AUTOPEPSIRF_HDI = config["autopepsirf_hdi"]

# Labels used to build output directory names
AUTOPEPSIRF_Z_LABEL = AUTOPEPSIRF_EXACT_Z_THRESH.replace(",", "-").replace(" ", "")
AUTOPEPSIRF_CS_LABEL = AUTOPEPSIRF_EXACT_CS_THRESH.replace(",", "-").replace(" ", "")
try:
    AUTOPEPSIRF_HDI_LABEL = str(int(round(float(AUTOPEPSIRF_HDI) * 100)))
except (TypeError, ValueError):
    AUTOPEPSIRF_HDI_LABEL = str(AUTOPEPSIRF_HDI).replace(".", "")

# AutoPepsirf outputs
AUTOPEPSIRF_DIFFENRICH_DIR = (
    f"{AUTOPEPSIRF_RUN_DIR}/"
    f"diffEnrich-aps_{AUTOPEPSIRF_Z_LABEL}Z_{AUTOPEPSIRF_CS_LABEL}CS-HDI{AUTOPEPSIRF_HDI_LABEL}"
)
AUTOPEPSIRF_PEPSIRF_DIR = f"{AUTOPEPSIRF_RUN_DIR}/" + config["autopepsirf_pepsirf_tsv_dir"]
AUTOPEPSIRF_RAW_INPUT = f"{AUTOPEPSIRF_RUN_DIR}/autopepsirf_raw_no_sblk.tsv"

# AutoPepsirf outputs used as PSEA inputs
PSEA_ZSCORE_TSV = (
    f"{AUTOPEPSIRF_PEPSIRF_DIR}/"
    f"{config['autopepsirf_tsv_base_str']}_Z-HDI{AUTOPEPSIRF_HDI_LABEL}.tsv"
)
PSEA_ZSCORE_SHORT_TSV = (
    f"{OUTPUT_DIR}/psea/"
    f"{config['autopepsirf_tsv_base_str']}_Z-HDI{AUTOPEPSIRF_HDI_LABEL}_short.tsv"
)
PSEA_PAIRS_TSV = f"{OUTPUT_DIR}/psea/pairs_filtered.tsv"

# ------------------------------------------------------------
# Optional flags and helper arguments
# ------------------------------------------------------------
ZERO_SEED_ARG = (
    f"--seed {config['zero_filter_seed']}"
    if config.get("zero_filter_seed") is not None
    else ""
)
ZERO_HIST_OUTPUT = config.get("zero_hist_output")
AUTOBINS_NOSCATT_ARG = "--no-scatterplots" if config.get("automatebins_no_scatterplots", False) else ""
AUTOBINS_MAX_PAIRS_ARG = (
    f"--max-pairs {config['automatebins_max_pairs']}"
    if config.get("automatebins_max_pairs") is not None
    else ""
)


def script_path(config_key: str, default_filename: str) -> str:
    """Resolve rule script path from an override key or shared scripts directory."""
    return config.get(config_key, f"{PYTHON_SCRIPTS_DIR}/{default_filename}")


def rule_threads(rule_name: str, default: int) -> int:
    return int(RULE_RESOURCES.get(rule_name, {}).get("threads", default))


def rule_mem_mb(rule_name: str, default: int) -> int:
    return int(RULE_RESOURCES.get(rule_name, {}).get("mem_mb", default))


def rule_runtime(rule_name: str, default: int) -> int:
    """Runtime in minutes for schedulers/profiles that map this resource to walltime."""
    return int(RULE_RESOURCES.get(rule_name, {}).get("runtime", default))


def expand_output_dir_placeholder(path_value):
    """Allow config paths like '{output_dir}/qc/file.png'."""
    if not isinstance(path_value, str):
        return path_value
    return path_value.replace("{output_dir}", OUTPUT_DIR)


def as_abs_path(path_value: str) -> str:
    """Resolve workflow-relative paths against the Snakefile directory."""
    if os.path.isabs(path_value):
        return path_value
    return str((Path(workflow.basedir) / path_value).resolve())


def resolve_output_path(path_value):
    """Expand {output_dir} placeholder and resolve to absolute path."""
    if not path_value:
        return ""
    return as_abs_path(expand_output_dir_placeholder(path_value))


def optional_run_token(_wildcards):
    """Conditionally inject run-token dependency to force full recomputation."""
    return RUN_TOKEN if FORCE_RERUN_WITH_TOKEN else []


def correlation_scatter_flags() -> str:
    # Optional scatterplot flags for correlation QC
    if not config.get("correlation_scatterplots", False):
        return ""
    flags = [
        "--scatterplots",
        f'--scatterplot-dir {config.get("correlation_scatterplot_dir", "results/qc/pairwise_scatterplots")}',
        f'--plot-format {config.get("correlation_plot_format", "png")}',
        f'--plot-dpi {config.get("correlation_plot_dpi", 300)}',
    ]
    max_pairs = config.get("correlation_max_pairs", None)
    if max_pairs is not None:
        flags.append(f"--max-pairs {max_pairs}")
    return " ".join(flags)

# ------------------------------------------------------------
# Final targets: rule all builds these
# ------------------------------------------------------------
ALL_TARGETS = [
    SHORTENED_COUNTS,
    SHORTEN_MAP,
    RC_OUT,
    RC_FAIL_OUT,
    RC_TOTALS_OUT,
    CORR_PAIRS_OUT,
    CORR_PAIRS_ALL_OUT,
    CORR_DATA_OUT,
    CORR_FAIL_OUT,
    ZERO_OUT,
    ZERO_FAIL_OUT,
    ZERO_ALL_COUNTS_OUT,
    ZERO_DOWNSAMPLED_COUNTS_OUT,
    QC_FAILURE_SUMMARY_OUT,
    AUTOBINS_FINAL,
    AUTOBINS_NEG_CTRL,
]
if RUN_AUTOPEPSIRF:
    ALL_TARGETS.append(AUTOPEPSIRF_RAW_INPUT)
    ALL_TARGETS.append(AUTOPEPSIRF_DIFFENRICH_DIR)
    ALL_TARGETS.append(AUTOPEPSIRF_PEPSIRF_DIR)
if RUN_ZSCORE_CORRELATION_QC:
    ALL_TARGETS.append(ZSCORE_CORR_PAIRS_OUT)
    ALL_TARGETS.append(ZSCORE_CORR_PAIRS_ALL_OUT)
    ALL_TARGETS.append(ZSCORE_CORR_DATA_OUT)
    ALL_TARGETS.append(ZSCORE_CORR_FAIL_OUT)
if RUN_PSEA:
    ALL_TARGETS.append(PSEA_PAIRS_TSV)
    ALL_TARGETS.append(f"{OUTPUT_DIR}/psea/psea-make-psea-table")


rule all:
    input:
        ALL_TARGETS,


# ------------------------------------------------------------
# Run token: forces re-runs via dependency on a temp file
# ------------------------------------------------------------
rule init_run_token:
    output:
        temp(RUN_TOKEN)
    threads: rule_threads("init_run_token", 1)
    resources:
        mem_mb=rule_mem_mb("init_run_token", 500),
        runtime=rule_runtime("init_run_token", 5)
    params:
        run_dir=RUN_DIR,
    shell:
        """
        mkdir -p "{params.run_dir}"
        date +%s > {output:q}
        """


# ------------------------------------------------------------
# 00: Shorten column names (keeps Sblk* full length)
# ------------------------------------------------------------
rule shorten_columns:
    input:
        RAW_COUNTS
    output:
        data=SHORTENED_COUNTS,
        mapping=SHORTEN_MAP,
    log:
        f"{LOG_DIR}/00_shorten_columns.log"
    threads: rule_threads("shorten_columns", 1)
    resources:
        mem_mb=rule_mem_mb("shorten_columns", 2000),
        runtime=rule_runtime("shorten_columns", 30)
    params:
        script=script_path("shorten_columns_script", "shorten_column_names.py"),
        qc_dir=QC_DIR,
        log_dir=LOG_DIR,
    shell:
        """
        mkdir -p "{params.qc_dir}" "{params.log_dir}"
        python3 "{params.script}" {input:q} {output.data:q} --map-output {output.mapping:q} > {log:q} 2>&1
        """


# ------------------------------------------------------------
# 01: Read count (RC) threshold filter
# ------------------------------------------------------------
rule rc_filter:
    input:
        data=SHORTENED_COUNTS,
        run_token=optional_run_token,
    output:
        pass_samples=temp(RC_PASS_OUT),
        data=RC_OUT,
        failed=RC_FAIL_OUT,
        totals=RC_TOTALS_OUT,
    log:
        f"{LOG_DIR}/01_rc_filter.log"
    threads: rule_threads("rc_filter", 1)
    resources:
        mem_mb=rule_mem_mb("rc_filter", 4000),
        runtime=rule_runtime("rc_filter", 60)
    params:
        script=script_path("rc_script", "QC_RC_thresh.py"),
        cutoff=config["rc_cutoff"],
        qc_dir=QC_DIR,
        log_dir=LOG_DIR,
        pepsirf_conda_path=AUTOPEPSIRF_CONDA_PATH,
    shell:
        """
        mkdir -p "{params.qc_dir}" "{params.log_dir}"
        python3 "{params.script}" {input.data:q} {params.cutoff} \
            -o {output.pass_samples:q} \
            --failed-samples-output {output.failed:q} \
            --all-samples-output {output.totals:q} > {log:q} 2>&1
        conda run -p "{params.pepsirf_conda_path}" pepsirf subjoin \
            -i "{input.data},{output.pass_samples}" \
            -o {output.data:q} >> {log:q} 2>&1
        """


# ------------------------------------------------------------
# 02: Downsample + zero-count distribution filter
# ------------------------------------------------------------
rule zero_count_filter:
    input:
        data=RC_OUT,
        run_token=optional_run_token,
    output:
        pass_samples=temp(ZERO_PASS_OUT),
        data=ZERO_OUT,
        failed=ZERO_FAIL_OUT,
        all_counts=ZERO_ALL_COUNTS_OUT,
        downsampled_counts=ZERO_DOWNSAMPLED_COUNTS_OUT,
    log:
        f"{LOG_DIR}/02_zero_count_filter.log"
    threads: rule_threads("zero_count_filter", 1)
    resources:
        mem_mb=rule_mem_mb("zero_count_filter", 6000),
        runtime=rule_runtime("zero_count_filter", 90)
    params:
        script=script_path("zero_filter_script", "downsample2max_zero_filter.py"),
        downsample_max=config["downsample_max"],
        min_reads=config.get("downsample_min", 0),
        sd_mult=config.get("zero_sd_mult", 3.0),
        seed_arg=ZERO_SEED_ARG,
        hist_output=(
            expand_output_dir_placeholder(ZERO_HIST_OUTPUT)
            if ZERO_HIST_OUTPUT
            else ""
        ),
        log_dir=LOG_DIR,
        pepsirf_conda_path=AUTOPEPSIRF_CONDA_PATH,
    shell:
        """
        mkdir -p "{params.log_dir}"
        if [ -n "{params.hist_output}" ]; then
            python3 "{params.script}" \
                -d {input.data:q} \
                -o {output.pass_samples:q} \
                --max {params.downsample_max} \
                --min {params.min_reads} \
                --sd-mult {params.sd_mult} \
                --failed-samples-output {output.failed:q} \
                --all-zero-counts-output {output.all_counts:q} \
                --downsampled-counts-output {output.downsampled_counts:q} \
                {params.seed_arg} \
                --hist-output {params.hist_output:q} > {log:q} 2>&1
        else
            python3 "{params.script}" \
                -d {input.data:q} \
                -o {output.pass_samples:q} \
                --max {params.downsample_max} \
                --min {params.min_reads} \
                --sd-mult {params.sd_mult} \
                --failed-samples-output {output.failed:q} \
                --all-zero-counts-output {output.all_counts:q} \
                --downsampled-counts-output {output.downsampled_counts:q} \
                {params.seed_arg} > {log:q} 2>&1
        fi
        conda run -p "{params.pepsirf_conda_path}" pepsirf subjoin \
            -i "{input.data},{output.pass_samples}" \
            -o {output.data:q} >> {log:q} 2>&1
        """


# ------------------------------------------------------------
# 03: Pairwise Pearson correlation QC (optional)
# ------------------------------------------------------------
rule correlation_filter:
    input:
        data=ZERO_OUT,
        run_token=optional_run_token,
        pairs=lambda wildcards: PAIRS_FILE if RUN_CORRELATION_QC else [],
    output:
        pairs=CORR_PAIRS_OUT,
        pairs_all=CORR_PAIRS_ALL_OUT,
        pass_samples=temp(CORR_PASS_OUT),
        data=CORR_DATA_OUT,
        failed=CORR_FAIL_OUT,
    log:
        f"{LOG_DIR}/03_correlation_filter.log"
    threads: rule_threads("correlation_filter", 1)
    resources:
        mem_mb=rule_mem_mb("correlation_filter", 8000),
        runtime=rule_runtime("correlation_filter", 120)
    params:
        script=script_path("correlation_script", "pairwise_pearson_from_pairs.py"),
        cutoff=config["pearson_cutoff"],
        enabled=str(RUN_CORRELATION_QC).lower(),
        scatter_flags=correlation_scatter_flags(),
        log_dir=LOG_DIR,
        pepsirf_conda_path=AUTOPEPSIRF_CONDA_PATH,
    shell:
        """
        mkdir -p "{params.log_dir}"
        if [ "{params.enabled}" = "true" ]; then
            python3 "{params.script}" {input.pairs:q} {input.data:q} {params.cutoff} \
                -o {output.pairs:q} \
                --all-correlations-output {output.pairs_all:q} \
                --passing-samples-output {output.pass_samples:q} \
                --failed-samples-output {output.failed:q} \
                {params.scatter_flags} > {log:q} 2>&1
            conda run -p "{params.pepsirf_conda_path}" pepsirf subjoin \
                -i "{input.data},{output.pass_samples}" \
                -o {output.data:q} >> {log:q} 2>&1
        else
            cp {input.data:q} {output.data:q}
            printf "Pair1\\tPair2\\tPearson Correlation\\tComparison\\n" > {output.pairs:q}
            printf "Pair1\\tPair2\\tPearson Correlation\\tComparison\\n" > {output.pairs_all:q}
            head -n 1 {input.data:q} | tr '\t' '\n' | tail -n +2 > {output.pass_samples:q}
            printf "Sample\\tFailed QC Step\\tFailure Value\\n" > {output.failed:q}
            printf "Correlation QC disabled; passthrough copy created.\\n" > {log:q}
        fi
        """
# ------------------------------------------------------------
# 04: Aggregate QC failures
# ------------------------------------------------------------
rule summarize_qc_failures:
    input:
        rc=RC_FAIL_OUT,
        corr=CORR_FAIL_OUT,
        zero=ZERO_FAIL_OUT,
        run_token=optional_run_token,
    output:
        QC_FAILURE_SUMMARY_OUT
    log:
        f"{LOG_DIR}/04_qc_failure_summary.log"
    threads: rule_threads("summarize_qc_failures", 1)
    resources:
        mem_mb=rule_mem_mb("summarize_qc_failures", 2000),
        runtime=rule_runtime("summarize_qc_failures", 30)
    params:
        script=script_path("qc_failure_merge_script", "merge_qc_failures.py"),
        log_dir=LOG_DIR,
    shell:
        """
        mkdir -p "{params.log_dir}"
        python3 "{params.script}" {input.rc:q} {input.corr:q} {input.zero:q} -o {output:q} > {log:q} 2>&1
        """


# ------------------------------------------------------------
# 05: Automate bins (uses final QC-filtered data)
# ------------------------------------------------------------
rule automate_bins:
    input:
        data=CORR_DATA_OUT,
        run_token=optional_run_token,
    output:
        bins=AUTOBINS_FINAL,
        neg_ctrl=AUTOBINS_NEG_CTRL,
    log:
        f"{LOG_DIR}/05_automate_bins.log"
    threads: rule_threads("automate_bins", 1)
    resources:
        mem_mb=rule_mem_mb("automate_bins", 12000),
        runtime=rule_runtime("automate_bins", 240)
    params:
        script=as_abs_path(script_path("automate_bins_script", "automateBins.py")),
        input_data=as_abs_path(CORR_DATA_OUT),
        workdir=as_abs_path(OUTPUT_DIR),
        log_file=as_abs_path(f"{LOG_DIR}/05_automate_bins.log"),
        sblk_id=config["sblk_id"],
        outprefix=AUTOBINS_PREFIX,
        top_fraction=config.get("automatebins_top_fraction", 0.20),
        min_readcount_sum=config.get("automatebins_min_readcount_sum", 0),
        noscatt_arg=AUTOBINS_NOSCATT_ARG,
        max_pairs_arg=AUTOBINS_MAX_PAIRS_ARG,
        plot_format=config.get("automatebins_plot_format", "png"),
        plot_dpi=config.get("automatebins_plot_dpi", 300),
        log_dir=LOG_DIR,
    shell:
        """
        mkdir -p "{params.log_dir}" "{params.workdir}"
        cd "{params.workdir}"
        python3 "{params.script}" \
            -i {params.input_data:q} \
            -s "{params.sblk_id}" \
            -o "{params.outprefix}" \
            --top-fraction {params.top_fraction} \
            --min-readcount-sum {params.min_readcount_sum} \
            --plot-format "{params.plot_format}" \
            --plot-dpi {params.plot_dpi} \
            {params.noscatt_arg} \
            {params.max_pairs_arg} > "{params.log_file}" 2>&1
        """


# ------------------------------------------------------------
# 06: Prepare AutoPepsirf input (remove Sblks used in bins, restore full length names)
# ------------------------------------------------------------
rule prepare_autopepsirf_raw:
    input:
        raw=CORR_DATA_OUT,
        neg_ctrl=AUTOBINS_NEG_CTRL,
        name_map=SHORTEN_MAP,
        run_token=optional_run_token,
    output:
        AUTOPEPSIRF_RAW_INPUT
    log:
        f"{LOG_DIR}/06_prepare_autopepsirf_raw.log"
    threads: rule_threads("prepare_autopepsirf_raw", 1)
    resources:
        mem_mb=rule_mem_mb("prepare_autopepsirf_raw", 4000),
        runtime=rule_runtime("prepare_autopepsirf_raw", 60)
    params:
        script=script_path("remove_controls_script", "remove_control_samples.py"),
        log_dir=LOG_DIR,
    shell:
        """
        mkdir -p "{params.log_dir}" "{AUTOPEPSIRF_RUN_DIR}"
        python3 "{params.script}" \
            --raw {input.raw:q} \
            --controls {input.neg_ctrl:q} \
            --name-map {input.name_map:q} \
            --output {output:q} > {log:q} 2>&1
        """


# ------------------------------------------------------------
# 07: Run AutoPepsirf (optional)
# ------------------------------------------------------------
rule run_autopepsirf:
    input:
        raw=AUTOPEPSIRF_RAW_INPUT,
        bins=AUTOBINS_FINAL,
        neg_ctrl=AUTOBINS_NEG_CTRL,
        run_token=optional_run_token,
    output:
        diffenrich=directory(AUTOPEPSIRF_DIFFENRICH_DIR),
        pepsirf=directory(AUTOPEPSIRF_PEPSIRF_DIR),
        zscore=PSEA_ZSCORE_TSV,
    log:
        f"{LOG_DIR}/07_autopepsirf.log"
    threads: rule_threads("run_autopepsirf", 1)
    resources:
        mem_mb=rule_mem_mb("run_autopepsirf", 32000),
        runtime=rule_runtime("run_autopepsirf", 720)
    params:
        negative_id=config["autopepsirf_negative_id"],
        exact_z_thresh=AUTOPEPSIRF_EXACT_Z_THRESH,
        exact_cs_thresh=AUTOPEPSIRF_EXACT_CS_THRESH,
        raw_constraint=config["autopepsirf_raw_constraint"],
        hdi=AUTOPEPSIRF_HDI,
        exact_zenrich_thresh=config["autopepsirf_exact_zenrich_thresh"],
        tsv_base_str=config["autopepsirf_tsv_base_str"],
        pepsirf_tsv_dir_path=AUTOPEPSIRF_PEPSIRF_DIR,
        conda_path=config["autopepsirf_conda_path"],
        outdir=AUTOPEPSIRF_DIFFENRICH_DIR,
        raw_path=as_abs_path(AUTOPEPSIRF_RAW_INPUT),
        bins_path=as_abs_path(AUTOBINS_FINAL),
        neg_ctrl_path=as_abs_path(AUTOBINS_NEG_CTRL),
        workdir=Q2_WORKDIR,
        q2_home=f"{OUTPUT_DIR}/.qiime2_home",
        tmpdir=f"{OUTPUT_DIR}/.tmp",
        npm_cache=f"{OUTPUT_DIR}/.tmp/npm-cache",
        log_dir=LOG_DIR,
    shell:
        """
        mkdir -p "{params.log_dir}" "{params.workdir}" "{params.pepsirf_tsv_dir_path}" \
            "{params.q2_home}" "{params.tmpdir}" "{params.npm_cache}"
        cd "{params.workdir}"
        export HOME="{params.q2_home}"
        export TMPDIR="{params.tmpdir}"
        export NPM_CONFIG_CACHE="{params.npm_cache}"
        conda run -p "{params.conda_path}" qiime autopepsirf diffEnrich-tsv \
            --p-raw-data-filepath {params.raw_path:q} \
            --p-bins-filepath {params.bins_path:q} \
            --p-negative-control-filepath {params.neg_ctrl_path:q} \
            --p-flexible-reps-source \
            --p-negative-id {params.negative_id} \
            --p-exact-z-thresh {params.exact_z_thresh} \
            --p-exact-cs-thresh {params.exact_cs_thresh} \
            --p-raw-constraint {params.raw_constraint} \
            --p-hdi {params.hdi} \
            --p-exact-zenrich-thresh {params.exact_zenrich_thresh} \
            --p-tsv-base-str "{params.tsv_base_str}" \
            --p-pepsirf-tsv-dir "{params.pepsirf_tsv_dir_path}" \
            --output-dir "{params.outdir}" >> {log:q} 2>&1
        """


# ------------------------------------------------------------
# 08: Shorten autopepsirf Z-score TSV for PSEA input (optional)
# ------------------------------------------------------------
rule shorten_psea_zscore:
    input:
        zscore=PSEA_ZSCORE_TSV,
        run_token=optional_run_token,
    output:
        PSEA_ZSCORE_SHORT_TSV
    log:
        f"{LOG_DIR}/08_shorten_psea_zscore.log"
    threads: rule_threads("shorten_psea_zscore", 1)
    resources:
        mem_mb=rule_mem_mb("shorten_psea_zscore", 4000),
        runtime=rule_runtime("shorten_psea_zscore", 60)
    params:
        script=script_path("shorten_columns_script", "shorten_column_names.py"),
        log_dir=LOG_DIR,
    shell:
        """
        mkdir -p "{params.log_dir}" "{OUTPUT_DIR}/psea"
        python3 "{params.script}" {input.zscore:q} {output:q} \
            --map-output /dev/null \
            --drop-prefix Sblk \
            --drop-nan > {log:q} 2>&1
        """


# ------------------------------------------------------------
# 08.5: Filter pairs to those present in shortened Z-score file
# ------------------------------------------------------------
rule filter_psea_pairs:
    input:
        pairs=PAIRS_FILE,
        zscore=PSEA_ZSCORE_SHORT_TSV,
        run_token=optional_run_token,
    output:
        PSEA_PAIRS_TSV
    log:
        f"{LOG_DIR}/08_filter_psea_pairs.log"
    threads: rule_threads("filter_psea_pairs", 1)
    resources:
        mem_mb=rule_mem_mb("filter_psea_pairs", 2000),
        runtime=rule_runtime("filter_psea_pairs", 30)
    params:
        script=script_path("filter_pairs_script", "filter_pairs_by_columns.py"),
        log_dir=LOG_DIR,
    shell:
        """
        mkdir -p "{params.log_dir}" "{OUTPUT_DIR}/psea"
        python3 "{params.script}" {input.pairs:q} {input.zscore:q} \
            -o {output:q} > {log:q} 2>&1
        """


# ------------------------------------------------------------
# 08.6: Optional Z-score correlation QC on autopepsirf output
# ------------------------------------------------------------
rule zscore_correlation_filter:
    input:
        data=PSEA_ZSCORE_SHORT_TSV,
        pairs=PAIRS_FILE,
        run_token=optional_run_token,
    output:
        pairs=ZSCORE_CORR_PAIRS_OUT,
        pairs_all=ZSCORE_CORR_PAIRS_ALL_OUT,
        pass_samples=temp(ZSCORE_CORR_PASS_OUT),
        data=ZSCORE_CORR_DATA_OUT,
        failed=ZSCORE_CORR_FAIL_OUT,
    log:
        f"{LOG_DIR}/08_zscore_correlation_filter.log"
    threads: rule_threads("zscore_correlation_filter", 1)
    resources:
        mem_mb=rule_mem_mb("zscore_correlation_filter", 8000),
        runtime=rule_runtime("zscore_correlation_filter", 120)
    params:
        script=script_path("correlation_script", "pairwise_pearson_from_pairs.py"),
        cutoff=config.get("zscore_pearson_cutoff", config["pearson_cutoff"]),
        log_cutoff_arg=(
            f"--log-pearson-cutoff {config['zscore_log_pearson_cutoff']}"
            if config.get("zscore_log_pearson_cutoff") is not None
            else ""
        ),
        log_dir=LOG_DIR,
        pepsirf_conda_path=AUTOPEPSIRF_CONDA_PATH,
    shell:
        """
        mkdir -p "{params.log_dir}" "{QC_DIR}"
        python3 "{params.script}" {input.pairs:q} {input.data:q} {params.cutoff} \
            {params.log_cutoff_arg} \
            -o {output.pairs:q} \
            --all-correlations-output {output.pairs_all:q} \
            --passing-samples-output {output.pass_samples:q} \
            --failed-samples-output {output.failed:q} > {log:q} 2>&1
        conda run -p "{params.pepsirf_conda_path}" pepsirf subjoin \
            -i "{input.data},{output.pass_samples}" \
            -o {output.data:q} >> {log:q} 2>&1
        """


# ------------------------------------------------------------
# 09: Run PSEA workflow (optional; requires AutoPepsirf)
# ------------------------------------------------------------
rule run_psea:
    input:
        run_token=optional_run_token,
        zscore_tsv=PSEA_ZSCORE_SHORT_TSV,
        pairs_tsv=PSEA_PAIRS_TSV,
        peptide_sets_gmt=config["psea_peptide_sets_gmt"],
        species_tsv=config["psea_species_tsv"],
    output:
        psea_table_dir=directory(f"{OUTPUT_DIR}/psea/psea-make-psea-table"),
    log:
        f"{LOG_DIR}/09_run_psea.log"
    threads: rule_threads("run_psea", 1)
    resources:
        mem_mb=rule_mem_mb("run_psea", 32000),
        runtime=rule_runtime("run_psea", 720)
    params:
        conda_path=config["psea_conda_path"],
        threshold=config.get("psea_threshold", 0.75),
        pval_thresh=config.get("psea_pval_thresh", 0.05),
        nes_thresh=config.get("psea_nes_thresh", 1),
        min_size=config.get("psea_min_size", 3),
        max_size=config.get("psea_max_size", 5000),
        permutation_num=config.get("psea_permutation_num", 10000),
        spline_type=config.get("psea_spline_type", "r-smooth"),
        degree=config.get("psea_degree", 3),
        vis_outputs_dir=resolve_output_path(config.get("psea_vis_outputs_dir")),
        vis_outputs_parent=(
            os.path.dirname(resolve_output_path(config.get("psea_vis_outputs_dir")))
            if config.get("psea_vis_outputs_dir")
            else ""
        ),
        vis_outputs_flag=(
            f'--p-vis-outputs-dir "{resolve_output_path(config.get("psea_vis_outputs_dir"))}"'
            if config.get("psea_vis_outputs_dir")
            else ""
        ),
        table_dir=f"{OUTPUT_DIR}/psea/psea_table_outdir",
        iterative_analysis=config.get("psea_iterative_analysis", True),
        summary_tables_dir=f"{OUTPUT_DIR}/psea/psea_ae_summary_tables",
        seed=config.get("psea_seed", 149),
        zscore_path=as_abs_path(PSEA_ZSCORE_SHORT_TSV),
        pairs_path=as_abs_path(PSEA_PAIRS_TSV),
        peptide_sets_path=as_abs_path(config["psea_peptide_sets_gmt"]),
        species_path=as_abs_path(config["psea_species_tsv"]),
        workdir=Q2_WORKDIR,
        q2_home=f"{OUTPUT_DIR}/.qiime2_home",
        log_dir=LOG_DIR,
    shell:
        """
        mkdir -p "{params.log_dir}" "{params.workdir}" "{params.q2_home}"
        if [ -n "{params.vis_outputs_parent}" ]; then
            mkdir -p "{params.vis_outputs_parent}"
        fi
        cd "{params.workdir}"
        export HOME="{params.q2_home}"

        conda run -p "{params.conda_path}" qiime psea make-psea-table \
            --p-scores-file {params.zscore_path:q} \
            --p-pairs-file {params.pairs_path:q} \
            --p-peptide-sets-file {params.peptide_sets_path:q} \
            --p-threshold {params.threshold} \
            --p-p-val-thresh {params.pval_thresh} \
            --p-nes-thresh {params.nes_thresh} \
            --p-species-taxa-file {params.species_path:q} \
            --p-min-size {params.min_size} \
            --p-max-size {params.max_size} \
            --p-permutation-num {params.permutation_num} \
            --p-spline-type {params.spline_type} \
            --p-degree {params.degree} \
            --p-table-dir {params.table_dir} \
            --p-iterative-analysis {params.iterative_analysis} \
            --p-summary-tables-dir {params.summary_tables_dir} \
            --p-seed {params.seed} \
            {params.vis_outputs_flag} \
            --verbose \
            --output-dir "{output.psea_table_dir}" >> {log:q} 2>&1

        """
