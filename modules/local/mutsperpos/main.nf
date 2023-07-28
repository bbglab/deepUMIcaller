process MUTS_PER_POS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::pysam=0.21.0--py38h15b938a_1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.21.0--py38h15b938a_1' :
        'biocontainers/pysam:0.21.0--py38h15b938a_1' }"

    input:
    tuple val(meta), path(bam), path(bam_index), path(vcf)

    output:
    tuple val(meta), path("**.png")     , emit: plots
    path  "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # this has been run here in the past
    /workspace/projects/prominent/analysis/prom10/get_muts_per_position




    vcf is filtered for only low VAF variants in principle, but we could let the python script filter it itself

    python countMutsPerCycle.py \\
    --inFile /workspace/nobackup/prominent/PROMINENT_05/results/deepUMIcaller_160423/sortbamduplexconshigh/B1.sorted.bam
    --inVCF ../data/K_6_1_A_1.noheader.vcf
    -o ../data/K_6_1_A_1
    -l 142
    --filter INCLUDE
    -t 0
    --clonality_limit 0.1
    -c 



    we also need to include functions and classes from the VCFReader file in the same scripts directory



# mutsByReadPos
myPlots = []

colorTrans = {
    "C>A": '#5abdeb',
    "C>G": '#050708',
    "C>T": '#d43c32',
    "T>A": '#cbcacb',
    "T>C": '#aacb72',
    "T>G": '#e7c9c6',
    "C>N": "#750082",
    "T>N": "#4B0082"
}


location = "/workspace/projects/prominent/analysis/prom10/get_muts_per_position/results"
    
    myMutsByCyc_original = pd.read_csv(
        f"{location}/{sample}.onlyVAF010_MutsPerCycle.dat.csv",
        # f"{location}/{sample}.onlyVAF010_original_MutsPerCycle.dat.csv",
        # f"{location}/{sample}_MutsPerCycle.dat.csv",
        sep = ",",
        skiprows = 1,
        names = ["Cycle", "C>T", "C>A", "C>G", "T>A", "T>C", "T>G", "C>N", "T>N", "Count", "Base", "Count_percent"]
    )
        
    myMutsByCyc_original['totMuts'] = myMutsByCyc_original['C>T'] + myMutsByCyc_original['C>G'] + myMutsByCyc_original['C>A'] + \
                        myMutsByCyc_original['T>A'] + myMutsByCyc_original['T>C'] + myMutsByCyc_original['T>G']

    myMutsByCyc = pd.melt(
        myMutsByCyc_original,
        id_vars=["Cycle", "Base", "Count", "Count_percent", "totMuts"],
        value_vars=["C>T", "C>A", "C>G", "T>A", "T>C", "T>G", "C>N", "T>N"],
        var_name="mutType",
        value_name="Number"
    )

    myMutsByCyc[['mutFrom', 'mutTo']] = myMutsByCyc['mutType'].str.split(">", expand=True)
    myMutsByCyc = myMutsByCyc[myMutsByCyc['mutTo'] != "N"]

    maxReads = myMutsByCyc['Count'].max()
    insertSizeCtr = 0
    maxCt = myMutsByCyc[myMutsByCyc['mutTo'] != "N"]['totMuts'].max()

    
    plt.figure(figsize = (13, 5))
    ax = plt.subplot(111)
    
    ax.set_ylim( 0, maxCt + maxCt * 0.05 )
    
    # Initialize a variable to keep track of the bottom positions for each mutation type
    bottoms = pd.Series(0, index=myMutsByCyc['Base'].unique())
    
    for mutType, color in colorTrans.items():
        mut_type_data = myMutsByCyc[myMutsByCyc['mutType'] == mutType].reset_index(drop = True)
        if len(mut_type_data) > 0:
            ax.bar(
                mut_type_data['Base'],
                mut_type_data['Number'],
                bottom=bottoms[mut_type_data['Base']],
                label=mutType,
                color=color
            )
            
        # Update the bottom positions for the next mutation type
        mut_type_data_indexed = mut_type_data.set_index("Base")
        bottoms += mut_type_data_indexed['Number']
        # print(bottoms.shape)

    ax2 = ax.twinx()
    ax2.plot(
        myMutsByCyc_original['Base'],
        myMutsByCyc_original['Count_percent'],
        label=f"Count Percent",
        color="black",
        linestyle='dashed'
    )
    
    ax2.set_ylim(
        min(myMutsByCyc_original['Count_percent']) - 0.05,
        max(myMutsByCyc_original['Count_percent']) + 0.05
    )
    
    
    # Set custom tick positions and labels for the secondary y-axis
    tick_values = [min(myMutsByCyc_original['Count_percent']),
                   max(myMutsByCyc_original['Count_percent'])]
    tick_labels = [f"{val:.2f}" for val in tick_values]
    ax2.set_yticks(tick_values)
    ax2.set_yticklabels(tick_labels)
    
    
    ax.set_ylabel("Count")
    ax2.set_ylabel("Fraction of Total Reads")
    ax.set_title(f"{sample} ({myMutsByCyc_original['totMuts'].sum()} mutated bases)")
    ax.legend()
    plt.show()
    myPlots.append(plt.gcf())
    plt.close()
    

# plt.savefig(f"{args[0]}_summaryMutsByCycle.pdf", dpi=300, bbox_inches='tight')
plt.show()



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.readjusted.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

    // recompute_depth.py \\
    //         --mpileup_file ${pileup_mutations} \\
    //         --vcf_file ${vcf} \\
    //         --output ${prefix}.readjusted.vcf