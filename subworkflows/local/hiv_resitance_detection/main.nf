//
// HIV resistance detecton
//

include { SIERRALOCAL         } from '../../../modules/local/sierralocal'
include { SAM2CODFREQ         } from '../../../modules/local/sam2codfreq'
include { RESISTANCE_REPORT  } from '../../../modules/local/resistance_report'

workflow HIV_RESISTANCE {
    take:
    consensus    // channel: [ val(meta), [ consensus ] ]
    bam          // channel: [ val(meta), [ bam ], [bai] ]
    profile_json // path   : /path/to/HIV1.json

    main:
    ch_versions = Channel.empty()

    //
    // HIV resistance detection
    //

    SIERRALOCAL (consensus)

    ch_versions      = ch_versions.mix(SIERRALOCAL.out.versions)

    SAM2CODFREQ (
        bam,
        profile_json
    )

    RESISTANCE_REPORT(
        SIERRALOCAL.out.json.collect{ it[1] }.ifEmpty([]),
        SAM2CODFREQ.out.codfreq.collect{ it[1] }.ifEmpty([])
    )

    ch_versions      = ch_versions.mix(RESISTANCE_REPORT.out.versions)

    emit:
    sierralocal_results = SIERRALOCAL.out.json       // channel: [ val(meta), [ json ] ]
    sam2codfreq_results = SAM2CODFREQ.out.codfreq    // channel: [ val(meta), [ codfreq ] ]
    resistance_report   = RESISTANCE_REPORT.out.csv  // channel: [ val(meta), [ csv  ] ]

    versions            = ch_versions                // channel: [ versions.yml ]
}
