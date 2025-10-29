//
// HIV resistance detecton
//

include { SIERRALOCAL                                        } from '../../../modules/local/sierralocal'
include { LIFTOFF                                            } from '../../../modules/nf-core/liftoff'
include { GFF2JSON                                           } from '../../../modules/local/gff2json'
include { BAM2CODFREQ                                        } from '../../../modules/local/bam2codfreq'
include { RESISTANCE_REPORT                                  } from '../../../modules/local/resistance_report'

workflow HIV_RESISTANCE {
    take:
    consensus      // channel: [ val(meta), [ consensus ] ]
    bam            // channel: [ val(meta), [ bam ], [bai] ]
    fasta          // path   : genome.fasta
    hiv_sequence   // path   : /path/to/codfreq.fasta
    hiv_annotation // path   : /path/to/codfreq.gff

    main:
    ch_versions = Channel.empty()

    //
    // HIV resistance detection
    //

    SIERRALOCAL (consensus)

    ch_versions      = ch_versions.mix(SIERRALOCAL.out.versions)

    LIFTOFF (
        fasta.map { f -> [ [id: f.baseName], f ] },
        hiv_sequence,
        hiv_annotation,
        []
    )

    ch_versions      = ch_versions.mix(LIFTOFF.out.versions)

    GFF2JSON (
        fasta,
        LIFTOFF.out.gff3
    )

    ch_versions      = ch_versions.mix(GFF2JSON.out.versions)


    BAM2CODFREQ (
        bam,
        GFF2JSON.out.profile_json
    )

    ch_versions      = ch_versions.mix(BAM2CODFREQ.out.versions)

    RESISTANCE_REPORT(
        SIERRALOCAL.out.json.collect{ it[1] }.ifEmpty([]),
        BAM2CODFREQ.out.codfreq.collect{ it[1] }.ifEmpty([])
    )

    ch_versions      = ch_versions.mix(RESISTANCE_REPORT.out.versions)

    emit:
    sierralocal_results = SIERRALOCAL.out.json       // channel: [ val(meta), [ json ] ]
    bam2codfreq_results = BAM2CODFREQ.out.codfreq    // channel: [ val(meta), [ codfreq ] ]
    resistance_report   = RESISTANCE_REPORT.out.csv  // channel: [ val(meta), [ csv  ] ]

    versions            = ch_versions                // channel: [ versions.yml ]
}
