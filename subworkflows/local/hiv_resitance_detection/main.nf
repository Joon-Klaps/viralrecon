//
// HIV resistance detecton
//

include { SIERRALOCAL  } from '../../../modules/local/sierralocal'

workflow HIV_RESISTANCE {
    take:
    consensus    // channel: [ val(meta), [ consensus ] ]

    main:
    ch_versions = Channel.empty()

    //
    // HIV resistance detecton
    //

    SIERRALOCAL (consensus)

    ch_versions      = ch_versions.mix(SIERRALOCAL.out.versions)

    emit:
    sierralocal_results = SIERRALOCAL.out.json     // channel: [ val(meta), [ json ] ]

    versions            = ch_versions              // channel: [ versions.yml ]
}
