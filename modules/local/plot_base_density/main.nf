process PLOT_BASE_DENSITY {
    tag "$fasta"
    label 'process_medium'

    conda "conda-forge::r-base=4.2 conda-forge::r-reshape2=1.4.4 conda-forge::r-optparse=1.7.5 conda-forge::r-ggplot2=3.5.1 conda-forge::r-scales=1.3.0 conda-forge::r-viridis=0.6.5 conda-forge::r-tidyverse=1.3.2 bioconda::bioconductor-biostrings=2.66.0 bioconda::bioconductor-complexheatmap=2.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5d/5d4923bc2b9c2d8c83468eddf1e32c0bfecffab8192e0e2ee69340fe267560cd/data' :
        'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-complexheatmap_r-base_r-ggplot2_pruned:6c26995a32d99713' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.pdf'), emit: pdf
    tuple val(meta), path('*.tsv'), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plot_base_density.r \\
        --fasta_files $fasta \\
        --prefixes $prefix \\
        --output_dir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
