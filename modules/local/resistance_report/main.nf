process RESISTANCE_REPORT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'xx' :
        'xx' }"

    input:
    path sierralocal_json, stageAs: "sierralocal_json/*"
    path mutation_csv    , stageAs: "mutation_tables/*"
    path resistance_csv  , stageAs: "resistance_tables/*"
    path nextclade_csv   , stageAs: "nextclade_reports/*"
    path consensus       , stageAs: "consensus/*"
    path annotation      , stageAs: "gff/*"

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
        --gff_folder ./gff \\
        --ivar_consensus_params "'${ivar_consensus_params}'" \\
        --output_html ${prefix}.html \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
