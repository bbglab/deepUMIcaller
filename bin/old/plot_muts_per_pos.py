#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt


muts_per_cycle_file = sys.argv[1]   # f"{sample}_MutsPerCycle.dat.csv"
sample_name = sys.argv[2]


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



myMutsByCyc_original = pd.read_csv(
    muts_per_cycle_file,
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
ax.set_title(f"{sample_name} ({myMutsByCyc_original['totMuts'].sum()} mutated bases)")
ax.legend()
plt.show()
myPlots.append(plt.gcf())
plt.close()


# plt.savefig(f"{args[0]}_summaryMutsByCycle.pdf", dpi=300, bbox_inches='tight')
plt.show()
