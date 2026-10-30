/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: Duplex metrics and consensus reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { clean_chr_names                                                   } from '../../../modules/local/utils'

include { FGUMI_COLLECTDUPLEXSEQMETRICS     as COLLECTSEQMETRICS            } from '../../../modules/local/fgumi/collectduplexseqmetrics/main'
include { FGUMI_COLLECTDUPLEXSEQMETRICS     as COLLECTSEQMETRICSONTARGET    } from '../../../modules/local/fgumi/collectduplexseqmetrics/main'
include { FAMILYSIZEMETRICS                 as FAMILYMETRICS                } from '../../../modules/local/duplex_metrics/familysizes/main'
include { FAMILYSIZEMETRICS                 as FAMILYMETRICSONTARGET        } from '../../../modules/local/duplex_metrics/familysizes/main'
include { COMPUTE_OPTIMAL_SEQUENCING        as OPTIMALSEQUENCING            } from '../../../modules/local/duplex_metrics/compilation/main'
include { UMICOLLISIONS                                                     } from '../../../modules/local/duplex_metrics/umicollisions/main'

workflow DUPLEX_METRICS {

    take:
    grouped_bam
    position_group_info
    interval_list

    main:
    COLLECTSEQMETRICS(grouped_bam, [])

    family_sizes_metrics = COLLECTSEQMETRICS.out.family_sizes
    if (params.split_by_chrom) {
        family_sizes_metrics = clean_chr_names(family_sizes_metrics)
    }

    FAMILYMETRICS(family_sizes_metrics)

    FAMILYMETRICS.out.sample_data.map{it -> it[1]}.collectFile(
        name: "metrics_summary.tsv",
        storeDir: "${params.outdir}/metrics/duplex/familymetrics",
        skip: 1,
        keepHeader: true
    ).set { family_metrics_all_metrics_data }

    FAMILYMETRICS.out.curve_data.map{it -> it[1]}.collectFile(
        name: "curves_summary.tsv",
        storeDir: "${params.outdir}/metrics/duplex/familymetrics",
        skip: 1,
        keepHeader: true
    ).set { family_metrics_all_curve_data }

    COLLECTSEQMETRICSONTARGET(grouped_bam, interval_list)

    family_sizes_metrics_ontarget = COLLECTSEQMETRICSONTARGET.out.family_sizes
    if (params.split_by_chrom) {
        family_sizes_metrics_ontarget = clean_chr_names(family_sizes_metrics_ontarget)
    }

    clean_chr_names(position_group_info)
        .join(clean_chr_names(COLLECTSEQMETRICSONTARGET.out.umi_counts))
        .join(clean_chr_names(COLLECTSEQMETRICSONTARGET.out.duplex_umi_counts))
        .set { umi_collisions_input }

    FAMILYMETRICSONTARGET(family_sizes_metrics_ontarget)

    FAMILYMETRICSONTARGET.out.sample_data.map{it -> it[1]}.collectFile(
        name: "metrics_summary.tsv",
        storeDir: "${params.outdir}/metrics/duplex/familymetricsontarget",
        skip: 1,
        keepHeader: true
    ).set { family_metrics_ontarget_metrics_data }

    FAMILYMETRICSONTARGET.out.curve_data.map{it -> it[1]}.collectFile(
        name: "curves_summary.tsv",
        storeDir: "${params.outdir}/metrics/duplex/familymetricsontarget",
        skip: 1,
        keepHeader: true
    ).set { family_metrics_ontarget_curve_data }

    UMICOLLISIONS(umi_collisions_input)

    UMICOLLISIONS.out.tsv.map{it -> it[1]}.collectFile(
        name: "umi_collisions_summary.tsv",
        storeDir: "${params.outdir}/metrics/duplex/umicollisions",
        skip: 1,
        keepHeader: true
    )

    // OPTIMALSEQUENCING(
    //     family_metrics_all_metrics_data,
    //     family_metrics_all_curve_data,
    //     family_metrics_ontarget_metrics_data,
    //     family_metrics_ontarget_curve_data
    // )

    emit:
    all_metrics         = family_metrics_all_metrics_data
    all_curves          = family_metrics_all_curve_data
    ontarget_metrics    = family_metrics_ontarget_metrics_data
    ontarget_curves     = family_metrics_ontarget_curve_data
}
