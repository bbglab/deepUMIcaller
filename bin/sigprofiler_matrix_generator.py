#!/usr/bin/env python3

import os, sys, glob, gzip, shutil
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

project_name = sys.argv[1]
genome = sys.argv[2]


file_list = glob.glob(f"./*.vcf*")
input_dir = "./input_mutations"

# create dir
os.mkdir(input_dir)

genome_mapping = {"hg19" : "GRCh37", "hg38": "GRCh38",
                    "GRCh37": "GRCh37", "GRCh38": "GRCh38",
                    "mm9": "mm9", "mm10": "mm10", "mm39": "mm39",
                    "GRCm37": "mm9", "GRCm38": "mm10", "GRCm39": "mm39",
                    }

for f in file_list:
    file = f.split("/")[-1]
    if file.endswith(".gz"):
        with gzip.open(f, 'r') as f_in, open(f'{input_dir}/{file[:-3]}', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        shutil.move(f, f'{input_dir}/{file}')

if genome in ["hg19", "hg38", "GRCh37", "GRCh38", "mm9", "mm10", "mm39", "GRCm37", "GRCm38", "GRCm39"]:
    chosen_genome = genome_mapping[genome]
else:
    print("genome not found, using GRCh38")
    chosen_genome = "GRCh38"

matGen.SigProfilerMatrixGeneratorFunc(project_name,
                                        chosen_genome,
                                        f"{input_dir}",
                                        plot=True, 
                                        exome=False, 
                                        bed_file=None, 
                                        chrom_based=False,
                                        tsb_stat=True,
                                        seqInfo=True,
                                        cushion=100,
                                        )

# from SigProfilerExtractor import sigpro as sig 
# sig.sigProfilerExtractor("vcf", "results_exome", "", reference_genome = chosen_genome, minimum_signatures=1, maximum_signatures=5, nmf_replicates=100) 