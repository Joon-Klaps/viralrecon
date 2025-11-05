//
// HIV resistance detecton
//

include { SIERRALOCAL                                        } from '../../../modules/local/sierralocal'
include { LIFTOFF                                            } from '../../../modules/nf-core/liftoff'
include { GFF2JSON                                           } from '../../../modules/local/gff2json'
include { BAM2CODFREQ                                        } from '../../../modules/local/bam2codfreq'
include { RESISTANCE_REPORT                                  } from '../../../modules/local/resistance_report'
include { ADDITIONAL_ANNOTATION as HIV_RESISTANCE_ANNOTATION } from '../additional_annotation'

workflow HIV_RESISTANCE {
    take:
    consensus      // channel: [ val(meta), [ consensus ] ]
    bam            // channel: [ val(meta), [ bam ], [bai] ]
    fasta          // path   : genome.fasta
    gff            // path   : genome.gff
    vcf            // channel: [ val(meta), [ vcf ] ]
    tbi            // channel: [ val(meta), [ tbi ] ]
    pangolin       // channel: [ val(meta), [ csv ] ]

    main:
    ch_versions = Channel.empty()

    //
    // HIV resistance detection
    //

    SIERRALOCAL (consensus)

    ch_versions = ch_versions.mix(SIERRALOCAL.out.versions)

    if (params.genome != 'codfreq') {

        // Not setted as param because we can't changed these files for codfreq and sierralocal to work toghether.
        // Using the ones in assets instead of the ones in test-datasets to avoid download issues when running offline
        codfreq_sequence   = file("$projectDir/assets/codfreq.fasta", checkIfExists: true)
        codfreq_annotation = file("$projectDir/assets/codfreq.gff", checkIfExists: true)
        LIFTOFF (
            fasta.map { f -> [ [id: f.baseName], f ] },
            codfreq_sequence,
            codfreq_annotation,
            []
        )
        ch_versions = ch_versions.mix(LIFTOFF.out.versions)

        GFF2JSON (
            fasta,
            LIFTOFF.out.gff3
        )

        HIV_RESISTANCE_ANNOTATION (
            vcf,
            tbi,
            fasta,
            LIFTOFF.out.gff3.map { it[1] },
            pangolin
        )
        ch_versions = ch_versions.mix(HIV_RESISTANCE_ANNOTATION.out.versions)

    } else {
        GFF2JSON (
            fasta,
            gff.map { f -> [ [id: f.baseName], f ] }
        )
    }

    ch_versions = ch_versions.mix(GFF2JSON.out.versions)

    BAM2CODFREQ (
        bam,
        GFF2JSON.out.profile_json
    )

    ch_versions = ch_versions.mix(BAM2CODFREQ.out.versions)

    RESISTANCE_REPORT(
        SIERRALOCAL.out.json.join(BAM2CODFREQ.out.codfreq, by: [0])
    )

    ch_versions = ch_versions.mix(RESISTANCE_REPORT.out.versions)

    emit:
    sierralocal_results = SIERRALOCAL.out.json                 // channel: [ val(meta), [ json ] ]
    bam2codfreq_results = BAM2CODFREQ.out.codfreq              // channel: [ val(meta), [ codfreq ] ]
    mutation_report     = RESISTANCE_REPORT.out.mutation_csv   // channel: [ val(meta), [ mutation_csv ] ]
    resistance_report   = RESISTANCE_REPORT.out.resistance_csv // channel: [ val(meta), [ resistance_csv ] ]

    versions            = ch_versions                          // channel: [ versions.yml ]
}
