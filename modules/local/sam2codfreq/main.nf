process SAM2CODFREQ {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.13.2 bioconda::pysam=0.23.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.23.3--py39hdd5828d_1' :
        'biocontainers/pysam:0.23.3--py39hdd5828d_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(profile)

    output:
    tuple val(meta), path("*.codfreq") , emit: codfreq
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:  // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    sam2codfreq.py \\
        --bam $bam \\
        --profile $profile \\
        --output ${prefix}.codfreq \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
