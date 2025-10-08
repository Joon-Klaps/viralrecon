process IVAR_VARIANTS_TO_VCF {
    tag "$meta.id"

    conda "conda-forge::python=3.13.2 conda-forge::matplotlib=3.10.1 conda-forge::pandas=2.2.3 conda-forge::r-sys=3.4.3 conda-forge::regex=2024.11.6 conda-forge::scipy=1.15.2 conda-forge::biopython=1.85"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/79e3c1e22b660e6a4f3655b1aeced00469a42bdff308be6e44910f4de0210ea0/data' :
        'community.wave.seqera.io/library/biopython_matplotlib_pandas_python_pruned:46d87e2ad1f8a063' }"

    input:
    tuple val(meta), path(tsv)
    path fasta
    path header

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:  // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ivar_variants_to_vcf.py \\
        $tsv \\
        ${prefix}.vcf \\
        --fasta $fasta \\
        $args \\
        > ${prefix}.variant_counts.log

    cat $header ${prefix}.variant_counts.log > ${prefix}.variant_counts_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
