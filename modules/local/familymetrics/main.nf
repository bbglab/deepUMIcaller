// TODO
// change the execution of this container in such a way that the script.py file is in the bin folder of the directory and then this can run from here
//  see sample check subworkflow

process FAMILYSIZEMETRICS {
    tag "$meta.id"
    label 'process_medium'
    
    // TODO
    // update this in the nfcore format once the container is available in biocontainers and galaxy singularity
    conda "anaconda::seaborn=0.12.2"
    container "biocontainers/seaborn:0.12.2_cv1"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
    //         'https://depot.galaxyproject.org/singularity/seaborn' : 
    //         'biocontainers/seaborn' }"


    input:
    tuple val(meta), path(groupby_metrics), path(duplex_metrics)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf

    // path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python

    import pandas as pd
    import seaborn as sns
    import numpy as np
    import matplotlib.pyplot as plt


    def compute_duplicates(sample,
                            file,
                            count_var = "count"):

        groupby_umi_families = pd.read_table(file, sep = "\\t", header = 0)
        groupby_umi_families['reads'] = groupby_umi_families[count_var] * groupby_umi_families['family_size']

        unique_molecules = groupby_umi_families[count_var].sum()
        total_molecules = groupby_umi_families['reads'].sum()
        duplicated_molecules = total_molecules - unique_molecules

        return duplicated_molecules / total_molecules    



    def stats_fam_size2plot(sample, groupby_metrics_file, duplex_metrics_file, limx):
        
        # compute duplicate rate from groupby stats
        prop_duplicates = compute_duplicates(sam, groupby_metrics_file)
        percent_duplicates = prop_duplicates * 100
        
        # compute family size distributions from duplex stats data
        data_duplex_families = pd.read_table(f"{duplex_metrics_file}")
        
        data_duplex_families["in_duplex"] = np.logical_and(data_duplex_families["ab_size"] > 0,
                                                            data_duplex_families["ba_size"] > 0)
        data_duplex_families_small = data_duplex_families[["ab_size", "ba_size", "count", "in_duplex"]]
        
        family_size1 = data_duplex_families_small[["ab_size", "count", "in_duplex"]]
        family_size2 = data_duplex_families_small[["ba_size", "count", "in_duplex"]]
        family_size1.columns = ["family_size", "count", "in_duplex"]
        family_size2.columns = ["family_size", "count", "in_duplex"]
        
        data_scss = pd.concat((family_size1, family_size2))
        data_scss = data_scss[data_scss["family_size"] > 0]
        data_scss["count_reads"] = data_scss["family_size"] * data_scss["count"]
        
        data_scss_grouped = data_scss.groupby(["family_size", "in_duplex"]).sum().reset_index()
        data_scss_grouped["family_size"] = data_scss_grouped["family_size"].astype(int)
        data_scss_grouped["fraction"] = data_scss_grouped["count"] / data_scss_grouped["count"].sum()
        data_scss_grouped["fraction_reads"] = data_scss_grouped["count_reads"] / data_scss_grouped["count_reads"].sum()

        total_duplex = data_duplex_families_small["count"][data_duplex_families_small["in_duplex"]].sum()
        total_scss = data_scss["count"].sum()
        total_reads = data_scss_grouped["count_reads"].sum()
        total_non_duplex = data_duplex_families_small["count"][~data_duplex_families_small["in_duplex"]].sum()
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (14, 6))
        sns.lineplot(data = data_scss_grouped,
                    x = "family_size",
                    y = "fraction_reads",
                    hue = "in_duplex",
                    marker = 'o',
                    palette= {True : "k", False : "r"},
                    ax = ax1
                    )

        ax1.set_xlim(limx)
        ax1.set_xlabel("Family size")
        ax1.set_ylabel("Fraction of reads")
        ax1.legend(title = "in duplex\\nfamilies")

        ax2.text(0.1, 0.8, f"Duplicates:             {percent_duplicates:.1f}%\\nNon-duplex SSCs:   {total_non_duplex/(total_non_duplex + total_duplex)*100:.1f}%")
        ax2.text(0.1, 0.5, f"Raw reads:   {total_reads:,}\\nSSCS:            {total_scss:,}\\nNon-duplex:  {total_non_duplex:,}\\nDuplex:         {total_duplex:,}" )
        ax2.text(0.1, 0.2, f"Raw/DCS:       {total_reads/total_duplex:.3f}\\nRaw/SSCS:      {total_reads/total_scss:.3f}\\nSSCS/Duplex:  {total_scss/total_duplex:.3f}")
        ax2.axis('off')

        fig.suptitle(sample)    

        plt.show()

        return fig


    sam = "${prefix}"
    groupby_metrics_file = f"${groupby_metrics}"
    duplex_metrics_file = f"{sam}.duplex_seq_metrics.duplex_family_sizes.txt"
    output_file = f"{sam}.family_sizes_plot_n_stats.pdf"
    x_axis_limits = (0,50)

    figure = stats_fam_size2plot(sam,
                                groupby_metrics_file,
                                duplex_metrics_file,
                                x_axis_limits)

    figure.savefig(output_file, bbox_inches='tight')
    """
}

