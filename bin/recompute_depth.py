#!/usr/bin/env python


import sys
import re
import pandas as pd
import numpy as np
from collections import Counter

def parse_mpu(x):

    # a typical mpileup output line look like (with option: -s -f $refgen)
    # ..+4ACAC.+4ACAC.+2AC => ['.', '.+ACAC', '.+ACAC', '.+AC']
    # .,,-7tttttgtt => ['.', ',', ',-tttttgt', 't']
    # ,$.,,^>. => [',$', '.', ',', ',', '^.']
    # ,.,.*.**.^*. => [',', '.', ',', '.', '*', '.', '*', '*', '.', '^.']
    # ,....,Tn.t => [',', '.', '.', '.', '.', ',', 'T', 'n', '.', 't']
    # A-1N => ['A-N']
    
    reads = x["STATUS"].upper()
    bqlist = x["QNAME"].split(",")

    readlist = []
    i = 0        # input pointer in reads
    j = 0        # output pointer in readlist
    while i < len(reads):
        if reads[i] in "ACGTNacgtn.,*":
            readlist.append(reads[i])
            i += 1
            j += 1
        elif reads[i] in '+-':
            # determine length
            digit = re.findall('[\+-](\d+)[ACGTNacgtn*]+',reads[i:])[0]
            readlist[j-1] += reads[i] + reads[i+1+len(digit):i+1+len(digit)+int(digit)]
            i += 1 + len(digit) + int(digit)
        elif reads[i] == '$':
            readlist[j-1] += '$'
            i += 1
        elif reads[i] == '^':
            # ^Xa, ^Xa$
            # readlist.append(reads[i:i+3])
            readlist.append(reads[i] + reads[i+2])
            i += 3
            j += 1
        else:
            print('*ERROR* mpileup parser: Unknown char {} in {}[{}]'
                    .format(reads[i], reads, i), file=sys.stderr)

    if len(readlist) != len(bqlist):
        print('*ERROR* mpileup parser: length mismatch between BQ string {} '
                'and reads string {} (breakdown={})'
                .format(bqlist, reads, ':'.join(readlist)), file=sys.stderr)

    return readlist, bqlist


def vartype(x,
            letters = ['A', 'T', 'C', 'G'],
            len_SV_lim = 100
            ):
        
    if ">" in (x["REF"] + x["ALT"]) or "<" in (x["REF"] + x["ALT"]):
        return "SV"

    elif len(x["REF"]) > (len_SV_lim+1) or len(x["ALT"]) > (len_SV_lim+1) :
        return "SV"
    
    elif x["REF"] in letters and x["ALT"] in letters:
        return "SNV"
    
    elif len(x["REF"]) == len(x["ALT"]):
        return "MNV"
    
    elif x["REF"] == "-" or ( len(x["REF"]) == 1 and x["ALT"].startswith(x["REF"]) ):
        return "INSERTION"
    
    elif x["ALT"] == "-" or ( len(x["ALT"]) == 1 and x["REF"].startswith(x["ALT"]) ):
        return "DELETION"
    
    return "COMPLEX"


def count_freq_in_row(elements_in_row, searching_elems):
    dp = 0
    for search in searching_elems:
        dp += elements_in_row.get(search, 0)
    return dp


def update_vcf_fields(row, suffix = ''):
    # update the formats
    formats = row["FORMAT"].split(':')
    values = row["SAMPLE"].split(':')
    formats_values = dict( zip(formats, values) )
    formats_values[f"CDP{suffix}"] = f'{int(row["TOT_DP"])}'
    formats_values[f"CAD{suffix}"] = f'{int(row["REF_DP"])},{int(row["ALT_DP"])}'
    formats_values[f"NDP{suffix}"] = f'{int(row["Ns_DP"])}'
    
    if row["NEW_FILTER"] != "":
        if row["FILTER"] == "PASS":
            filt = row["NEW_FILTER"]
        else:
            filt = row["FILTER"] + ";" + row["NEW_FILTER"]
    else:
        filt = row["FILTER"]
    
    updfilt = ";".join( sorted(filt.split(";")) )
    
    return (updfilt, ":".join(formats_values.keys()), ":".join(formats_values.values()) )


def recompute_depth(vcf,
                    mpileup_data,
                    not_supported_filter = "no_pileup_support",
                    not_searched_filter = "not_searched_",
                    suffix_label = ''
                    ):
    """
    vcf file
    mpileup_data
        This should be parsed with parse_mpu(x) function and with the CHROM and POS in the index
    """
    # elements of this list will have 3 sections
    #    total_dp
    #    ref dp
    #    alt dp
    #    Ns dp
    #    filter additions
    
    info_vcf = []
    print("Parsing rows and recomputing depth...", end = "\t")
    for ind, row in vcf.iterrows():
        var_tp = vartype(row[["REF", "ALT"]])

        if var_tp == "SV":
            if (row["CHROM"], row["POS"]) in mpileup_data.index:
                mpileup_row = mpileup_data.loc[(row["CHROM"], row["POS"])]
                elements_in_row = Counter(mpileup_row["SPLIT_bases"])
                count_ns = count_freq_in_row(elements_in_row, ["N"])
                ref_dp = count_freq_in_row(elements_in_row, [",", "."])
                deleted_nucs = count_freq_in_row(elements_in_row, ["*"])
                total_dp = mpileup_row["DEPTH"] - count_ns - deleted_nucs

                info_vcf.append(
                            (total_dp, ref_dp, 0, count_ns, f"{not_searched_filter}{var_tp};{not_supported_filter}")
                        )
                # print(info_vcf[-1])
            else:
                info_vcf.append(
                            (0, 0, 0, 0, f"{not_searched_filter}{var_tp};{not_supported_filter}")
                        )

        elif var_tp == "COMPLEX":
            if (row["CHROM"], row["POS"]) in mpileup_data.index:
                mpileup_row = mpileup_data.loc[(row["CHROM"], row["POS"])]
                elements_in_row = Counter(mpileup_row["SPLIT_bases"])
                count_ns = count_freq_in_row(elements_in_row, ["N"])
                ref_dp = count_freq_in_row(elements_in_row, [",", "."])
                deleted_nucs = count_freq_in_row(elements_in_row, ["*"])
                total_dp = mpileup_row["DEPTH"] - count_ns - deleted_nucs
                
                info_vcf.append(
                            (total_dp, ref_dp, 0, count_ns, f"{not_searched_filter}{var_tp};{not_supported_filter}")
                        )
                # print(info_vcf[-1])
            else:
                info_vcf.append(
                            (0, 0, 0, 0, f"{not_searched_filter}{var_tp};{not_supported_filter}")
                        )


        elif var_tp == "MNV":
            # here there is a possible update to check that the mutated bases appear in the same read

            alt_dp_count = -1
            ref_dp_count = 0
            ns_dp_count = 0
            total_dp_count = 0
            
            pos_count = 0
            not_present = False
            
            for i in range( len(row["ALT"]) ):
                search = row["ALT"][i]
                if (row["CHROM"], row["POS"] + i) in mpileup_data.index:
                    pos_count += 1
                    mpileup_row = mpileup_data.loc[(row["CHROM"], row["POS"] + i)]
                    # splitted_row_read_names = mpileup_row["SPLIT_reads"]

                    elements_in_row = Counter(mpileup_row["SPLIT_bases"])
                    count_ns = count_freq_in_row(elements_in_row, ["N"])
                    deleted_nucs = count_freq_in_row(elements_in_row, ["*"])
                    ns_dp_count += count_ns

                    ref_dp = count_freq_in_row(elements_in_row, [",", "."])
                    ref_dp_count += ref_dp
                    
                    total_dp = mpileup_row["DEPTH"] - count_ns - deleted_nucs
                    total_dp_count += total_dp

                    count_alt = count_freq_in_row(elements_in_row, search)
                    if not_present:
                        # if the variant is not there, keep iterating to update
                        # the total and reference depth
                        # but not the alt depth
                        continue
                        
                    elif count_alt > 0:
                        if alt_dp_count == -1:
                            alt_dp_count = count_alt
                        else:
                            alt_dp_count = min(alt_dp_count, count_alt)
                    else:
                        alt_dp_count = -1
                        not_present = True

            # if the variant has enough support add it to the list otherwise add the depth and the Ns
            if alt_dp_count > 0:
                info_vcf.append(
                        (round(total_dp_count / pos_count),
                         round(ref_dp_count / pos_count),
                         alt_dp_count,
                         round(ns_dp_count / pos_count),
                         f"")
                    )

            elif pos_count == 0:
                info_vcf.append(
                            (0, 0, 0, 0, f"{not_supported_filter}")
                        )

            else:
                info_vcf.append(
                    (round(total_dp_count / pos_count),
                     round(ref_dp_count / pos_count),
                     0,
                     round(ns_dp_count / pos_count),
                     f"{not_supported_filter}")
                )


        elif var_tp == "DELETION":
            # print("DELETION")

            search = f'-{row["REF"][1:]}'
            # print(search)
            alt_dp = 0
            count_ns = 0
            ref_dp = 0
            total_dp = 0

            if (row["CHROM"], row["POS"]) in mpileup_data.index:
                mpileup_row = mpileup_data.loc[(row["CHROM"], row["POS"])]

                elements_in_row = Counter(mpileup_row["SPLIT_bases"])
                alt_dp = count_freq_in_row(elements_in_row, [f".{search}", f",{search}"])
                count_ns = count_freq_in_row(elements_in_row, ["N"])
                ref_dp = count_freq_in_row(elements_in_row, [",", "."])
                deleted_nucs = count_freq_in_row(elements_in_row, ["*"])
                total_dp = mpileup_row["DEPTH"] - count_ns - deleted_nucs                

            if alt_dp > 0:
                info_vcf.append(
                        (total_dp,
                            ref_dp,
                            alt_dp,
                            count_ns,
                            f"")
                    )
            else:
                info_vcf.append(
                        (total_dp,
                            ref_dp,
                            0,
                            count_ns,
                            f"{not_supported_filter}")
                    )

            # print(info_vcf[-1])


        elif var_tp == "INSERTION":
            # print("INSERTION")

            search = f'+{row["ALT"][1:]}'
            # print(search)

            alt_dp = 0
            count_ns = 0
            ref_dp = 0
            total_dp = 0
            
            if (row["CHROM"], row["POS"]) in mpileup_data.index:
                mpileup_row = mpileup_data.loc[(row["CHROM"], row["POS"])]

                elements_in_row = Counter(mpileup_row["SPLIT_bases"])
                alt_dp = count_freq_in_row(elements_in_row, [f".{search}", f",{search}"])
                count_ns = count_freq_in_row(elements_in_row, ["N"])
                ref_dp = count_freq_in_row(elements_in_row, [",", "."])
                deleted_nucs = count_freq_in_row(elements_in_row, ["*"])
                total_dp = mpileup_row["DEPTH"] - count_ns - deleted_nucs

            if alt_dp > 0:
                info_vcf.append(
                        (total_dp,
                            ref_dp,
                            alt_dp,
                            count_ns,
                            f"")
                    )
            else:
                info_vcf.append(
                        (total_dp,
                            ref_dp,
                            0,
                            count_ns,
                            f"{not_supported_filter}")
                    )
            # print(info_vcf[-1])


        # it means it is a SNV
        else:

            alt_dp = 0
            count_ns = 0
            ref_dp = 0
            total_dp = 0
            
            if (row["CHROM"], row["POS"]) in mpileup_data.index:
                mpileup_row = mpileup_data.loc[(row["CHROM"], row["POS"])]

                elements_in_row = Counter(mpileup_row["SPLIT_bases"])
                alt_dp = count_freq_in_row(elements_in_row, [row["ALT"]])
                count_ns = count_freq_in_row(elements_in_row, ["N"])
                ref_dp = count_freq_in_row(elements_in_row, [",", "."])
                deleted_nucs = count_freq_in_row(elements_in_row, ["*"])
                total_dp = mpileup_row["DEPTH"] - count_ns - deleted_nucs

            if alt_dp > 0:
                info_vcf.append(
                        (total_dp,
                            ref_dp,
                            alt_dp,
                            count_ns,
                            f"")
                    )
            else:
                info_vcf.append(
                        (total_dp,
                            ref_dp,
                            0,
                            count_ns,
                            f"{not_supported_filter}")
                    )

            # print(info_vcf[-1])


        # print()
    
    print("Done\nAdding columns to vcf dataframe...", end = "\t")
    vcf[["TOT_DP", "REF_DP", "ALT_DP", "Ns_DP", "NEW_FILTER"]] = info_vcf

    print("Done\nUpdating VCF columns...", end = "")
    vcf[["FILTER", "FORMAT", "SAMPLE"]] = vcf.apply(lambda x: pd.Series(update_vcf_fields(x, suffix = suffix_label)), axis = 1)
    print("\nDONE")
    
    return vcf[["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]]




def main(mpileup_file, vcf_file, output_filename, suffix = ''):
    """
    Your script's main function.
    """

    ###
    # Read the VCF header
    ###
    first_header_lines = []
    header_lines = []
    with open(vcf_file, 'r') as vcf_open_file:
        for line in vcf_open_file:
            if line.startswith('##'):
                if line.strip().strip("#").split("=")[0].islower():
                    first_header_lines.append(line.strip())
                else:
                    header_lines.append(line.strip())
            elif line.startswith('#'):
                single_header = line.strip()
            else:
                break

    if suffix != '':
            not_supported_filter = f"{suffix}_no_pileup_support"
            not_searched_filter = f"{suffix}_not_searched_"
    else:
        not_supported_filter = "no_pileup_support"
        not_searched_filter = "not_searched_"
    
    ###
    # Add your custom header lines
    ###
    # FORMAT/INFO fields
    header_lines.append(f'##FORMAT=<ID=CDP{suffix},Number=1,Type=Integer,Description="Total Depth recomputed using mpileup output. deepUMIcaller.">')
    header_lines.append(f'##FORMAT=<ID=CAD{suffix},Number=2,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed. ONLY ONE ALT allele. Recomputed using mpileup output, this might not sum to CDP. deepUMIcaller.">')
    header_lines.append(f'##FORMAT=<ID=NDP{suffix},Number=1,Type=Integer,Description="Total number of Ns computed using mpileup output. deepUMIcaller.">')

    # FILTER field
    header_lines.append(f'##FILTER=<ID={not_supported_filter},Description="Variant not supported when inspecting the BAM with mpileup. deepUMIcaller.">')
    header_lines.append(f'##FILTER=<ID={not_searched_filter}SV,Description="Structural variant. Recounting of alt depth from mpileup output not done for this type of variants. deepUMIcaller.">')
    header_lines.append(f'##FILTER=<ID={not_searched_filter}COMPLEX,Description="Complex variant. Recounting of alt depth from mpileup output not done for this type of variants. deepUMIcaller.">')

     ###
    # Combine the modified header rows
    ###
    # Convert the list of header lines to a single string
    header_str = '\n'.join(sorted(first_header_lines) + sorted(header_lines))
    header_str = header_str + "\n" + single_header

    ###
    # Read the VCFs per chunk
    ###
    
    # Chunk size for reading
    chunk_size = 1000
    
    ###
    # Read and preprocess the mpileup data
    ###
    mpileup_data_chunks = pd.read_csv(mpileup_file, sep="\t", header=None,
                                dtype={0: str, 1: int, 2: str, 3: int, 4: str, 5: str, 6: str},
                                na_filter=False, chunksize=chunk_size)

    mpileup_data_list = []
    for mpileup_chunk in mpileup_data_chunks:
        mpileup_chunk.columns = ["CHROM", "POS", "REF", "DEPTH", "STATUS", "QUAL", "QNAME"]
        mpileup_chunk[["SPLIT_bases", "SPLIT_reads"]] = mpileup_chunk[["STATUS", "QNAME"]].apply(lambda x: pd.Series(parse_mpu(x)), axis=1)
        mpileup_chunk.drop(["STATUS", "QUAL", "QNAME"], axis=1, inplace=True)
        mpileup_chunk.set_index(["CHROM", "POS"], inplace=True)
        mpileup_data_list.append(mpileup_chunk)

    # Concatenate all chunks into one dataframe
    mpileup_data = pd.concat(mpileup_data_list)


    ###
    # Write the modified VCF header and the updated VCF body into a new file with the appropriate name
    ###
    with open(output_filename, 'w', encoding="utf-8") as new_vcf_file:
        new_vcf_file.write(header_str + '\n')  # Write the modified header

        ###
        # Read and preprocess the VCF file body
        ###
        vcf_chunks = pd.read_csv(vcf_file, sep = '\t', header = None, comment= '#',
                            dtype={0: str, 1: int, 2: str, 3: str, 4: str, 5: int, 6: str, 7: str, 8: str, 9: str},
                            na_filter=False, chunksize=chunk_size)    
    
        for vcf_chunk in vcf_chunks:
            vcf_chunk.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
            updated_vcf_chunk = recompute_depth(vcf_chunk, mpileup_data,
                                        not_supported_filter = not_supported_filter,
                                        not_searched_filter = not_searched_filter,
                                        suffix_label= suffix)
            updated_vcf_chunk.to_csv(new_vcf_file, sep='\t', index=False, header=False, mode='a')

   
    # Print a success message or return a result if needed
    print(f"{output_filename} VCF file created successfully.")

if __name__ == '__main__':
    mpileup_file = sys.argv[1]
    vcf_file = sys.argv[2]
    output_filename = sys.argv[3]

    # main(mpileup_file, vcf_file, output_filename) # 2024-07-27
    # Ferriol commented this line since it seems to be repetitive and a remaining of a previous version of the code
    # that should have been removed after commit 70ed9c2a29530ae2961423e7930fb5cda4049bc3

    if len(sys.argv) >= 5:
        suffix_cmd = sys.argv[4]
    else:
        suffix_cmd = ''
    main(mpileup_file, vcf_file, output_filename, suffix_cmd)


# @click.command()
# @click.option('--mpileup_file', '-m', required=True, type=click.Path(exists=True),
#                 help='Path to the input mpileup file.')
# @click.option('--vcf_file', '-v', required=True, type=click.Path(exists=True),
#                 help='Path to the output VCF file.')
# @click.option('--output_filename', '-o', required=True, type=click.Path(),
#                 help='Path to the output file.')



###
# Draft code to start recounting some complex variants...
###

# print("Deletion with non coincident start")

#             ## search for the deleted part
#             search = f'-{row["REF"][:-1]}'
#             updated_pos = row["POS"] - 1

#             mpileup_row = mpileup_data.loc[(row["CHROM"], updated_pos)]
#             splitted_row = mpileup_row["SPLIT_bases"]#.values[0]
#             splitted_row_read_names = mpileup_row["SPLIT_reads"]#.values[0]
#             elements_in_row = Counter(splitted_row)

#             ref_dp = 0
#             if "," in elements_in_row:
#                 ref_dp += elements_in_row[","]

#             if "." in elements_in_row:
#                 ref_dp += elements_in_row["."]

#             alt_dp_part1 = 0
#             if f".{search}" in elements_in_row:
#                 alt_dp_part1 += elements_in_row[f".{search}"]
#             if f",{search}" in elements_in_row:
#                 alt_dp_part1 += elements_in_row[f",{search}"]



#             ## search for the change of nucleotide part
#             search = f'{row["ALT"]}'
#             updated_pos = row["POS"] + len(row["REF"]) - 1
#             mpileup_row = mpileup_data.loc[(row["CHROM"], updated_pos)]
#             splitted_row = mpileup_row["SPLIT_bases"]#.values[0]
#             splitted_row_read_names = mpileup_row["SPLIT_reads"]#.values[0]
#             elements_in_row = Counter(splitted_row)


#             ref_dp2 = 0
#             if "," in elements_in_row:
#                 ref_dp2 += elements_in_row[","]

#             if "." in elements_in_row:
#                 ref_dp2 += elements_in_row["."]

#             alt_dp_part2 = 0
#             if f"{search}" in elements_in_row:
#                 alt_dp_part2 += elements_in_row[f"{search}"]            


#             if alt_dp_part1 == alt_dp_part2:
#                 alt_dp = alt_dp_part1
#                 ref_dp = (ref_dp + ref_dp2) // 2

#             else:
#                 print("variant not fully reconstructed")
#                 alt_dp = 0


#             list_depth_diff.append(total_dp - row["DEPTH"])
#             total_dp = ref_dp + alt_dp

#             if alt_dp > 0:
#                 info_vcf.append((ref_dp, alt_dp, ""))

#                 print_all = False

#                 if alt_dp == row["ALT_DEPTH"]:
#                     print("agree on alt_depth")
#                     dp_diff = total_dp - row["DEPTH"]
#                     print(f"Mpileup - VarDict: \tDepth: {dp_diff}")

#                 else:
#                     list_disagreeing.append( (row["POS"], row["ALT_DEPTH"], alt_dp ) )
#                     print("samtools mpileup\t", total_dp, ref_dp, alt_dp, alt_dp / total_dp )
#                     print("VarDict dp. values\t", row["DEPTH"], row["DEPTH"] - row["ALT_DEPTH"], row["ALT_DEPTH"], row["ALT_DEPTH"] / row["DEPTH"])
#                     print(row)


#             else:
#                 info_vcf.append((ref_dp, 0, "not_in_bam"))
#                 print("variant does not exist")
#                 print("samtools mpileup\t", total_dp, ref_dp, alt_dp, alt_dp / total_dp )
#                 print(row)
