import sys
import re

def load_variants(vcf_file):
    variants = set()
    with open(vcf_file, "r") as f:
        for line in f:
            if not line.strip():  # skip empty lines
                continue
            if line.startswith("#"):  # skip header
                continue

            cols = re.split(r"\s+", line.strip())
            if len(cols) < 5:  # malformed line, skip
                continue

            chrom = cols[0]
            pos = cols[1]
            ref = cols[3]
            alt = cols[4]
            variants.add((chrom, pos, ref, alt))
    return variants

def main(vcf1, vcf2):
    v1_variants = load_variants(vcf1)
    v2_variants = load_variants(vcf2)

    if not v1_variants:
        print("No variants in vcf1")
        return

    overlap = v1_variants.intersection(v2_variants)
    union = v1_variants.union(v2_variants) #Jaccard index
    percent = int((len(overlap) / len(union)) * 100)

    print(percent)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <vcf1> <vcf2>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
