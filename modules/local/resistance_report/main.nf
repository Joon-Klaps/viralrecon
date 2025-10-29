process RESISTANCE_REPORT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5e/5ee6e81aff2205d76ad8755d2181f8ea1dd747daa92aaf9dcba943b69aa9f458/data' :
        'community.wave.seqera.io/library/matplotlib_pandas_python_r-sys_pruned:23244d66110fcdf2' }"

    input:
    path(json)
    path(codfreq)

    output:
    path "*.csv"       , emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:  // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    def args = task.ext.args ?: ''
    """
    resistance_report.py \\
        --sierralocal_dir . \\
        --codfreq_dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
