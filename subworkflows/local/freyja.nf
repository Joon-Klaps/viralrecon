//
// Determening and deconvulion variants
//

include { FREYJA_VARIANTS   } from '../../modules/local/'
include { FREYJA_VARIANTS   } from '../../modules/local/'


workflow FREYJA {
    take:
    bam          // channel: [ val(meta), [ /path/to/bam ] ]
    fasta        // /path/to/reference.fasta


    main:

    ch_versions = Channel.empty()

    //
    // Variant calling
    // TODO: Create the freyja module for calling variants & depths
    //
    ch_freyja_variants = Channel.empty()
    ch_freyja_depth    = Channel.empty()
    FREYJA_VARIANTS(
        bam,
        fasta)
    ch_freyja_variants = FREYJA_VARIANTS.out.variants
    ch_freyja_depth    = FREYJA_VARIANTS.out.depth


    //
    // demix and define minimum variant abundances
    // TODO: make a module that demixes the variant file
    //
    ch_freyja_abundacies    = Channel.empty()
    FREYJA_DEMIX(
        ch_freyja_variants,
        ch_freyja_depth)
    ch_freyja_abundacies = FREYJA_DEMIX.out.abundancies

    //
    // Perform bootstrapping to get more accurate estimates of abundancies
    // TODO: make the freyja module for bootstrapping
    //
    FREYJA_BOOTSTRAP(
        ch_freyja_variants,
        ch_freyja_depth)


    emit:
    versions         = ch_versions         // channel: [ versions.yml ]
}
