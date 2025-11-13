process RESISTANCE_REPORT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ab/abf1008542dc33f026e7c9ced0760345fbf2ad6a18ca54e644236fd6ad95f030/data' :
        'community.wave.seqera.io/library/jinja2_pandas_python:5f56d3297f91b6d9' }"

    input:
    path sierralocal_json, stageAs: "sierralocal_json/*"
    path mutation_csv    , stageAs: "mutation_tables/*"
    path resistance_csv  , stageAs: "resistance_tables/*"
    path nextclade_csv   , stageAs: "nextclade_reports/*"
    path consensus       , stageAs: "consensus/*"
    path html_template
    path css_file

    output:
    path("*.html")      , emit: mutation_csv
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:  // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    def args = task.ext.args ?: ''
    def ivar_consensus_params = task.ext.args2 ?: '-t N/A -q N/A -m N/A -n N'
    def prefix = task.ext.prefix ?: 'resistance'

    """
    resistance_report.py \\
        --sierralocal_folder ./sierralocal_json \\
        --mutation_folder ./mutation_tables \\
        --resistance_folder ./resistance_tables \\
        --nextclade_folder ./nextclade_reports \\
        --consensus_folder ./consensus \\
        --ivar_consensus_params "'${ivar_consensus_params}'" \\
        --template ${html_template} \\
        --css ${css_file} \\
        --output_html ${prefix}.html \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
