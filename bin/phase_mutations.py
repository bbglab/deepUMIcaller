#!/usr/bin/env python3
"""Phase reads into two mutation chains from per-position mutated/reference read lists.

Input TSV is expected to contain at least these columns:
	- CHROM
	- POS
	- REF
	- ALT
	- MUTATED_READS
	- REFERENCE_READS

`MUTATED_READS` and `REFERENCE_READS` must be Python-like list strings, e.g.
"['read1', 'read2']" or "[]".

The script builds read-level constraints for each variant:
	- reads in the same allele group (MUTATED or REFERENCE) -> same chain
	- reads in opposite allele groups -> opposite chains

Then it resolves all constraints with a union-find structure with parity, producing
chain labels (chain_0 / chain_1) inside each connected component.
"""

from __future__ import annotations

import argparse
import ast
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple, Set

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

class Component:
	"""Segment with coverage"""

	def __init__(self, number: int, chromosome: str, starting_pos: int, chain_0_name: str, chain_0: Set[str], chain_1: Set[str], variant_id: str) -> None:
		self.number: int = number
		self.starting_pos: int = starting_pos
		self.chromosome: str = chromosome

		self.size: int = 1
		self.variant_id: List[str] = [variant_id]
		
		self.chain_0: Set[str] = chain_0
		self.chain_1: Set[str] = chain_1

		self.size_chain_0: List[int] = [len(self.chain_0)]
		self.size_chain_1: List[int] = [len(self.chain_1)]

		self.chain_0_status: List[str] = [chain_0_name]
		self.chain_1_status: List[str] = ["reference" if chain_0_name == "mutated" else "mutated"]

		self.outliers: List[str] = []

		
	def increase_chains_0_n_1(self, variant_id: str, status_0 :str, status_1 :str, reads_0: Set[str], reads_1: Set[str]) -> None:
		self.variant_id.append(variant_id)
		self.size += 1

		self.chain_0.update(reads_0)
		self.chain_1.update(reads_1)
		
		self.size_chain_0.append(len(self.chain_0))
		self.size_chain_1.append(len(self.chain_1))
		
		self.chain_0_status.append(status_0)
		self.chain_1_status.append(status_1)

	def add_missing_mutation(self, variant_id: str, chain: str, reads: Set[str]) -> None:
		self.outliers.append({"variant": variant_id, "chain": chain, "type": "expected", "reads": sorted(reads)})


	def add_extra_mutation(self, variant_id: str, chain: str, reads: Set[str]) -> None:
		self.outliers.append({"variant": variant_id, "chain": chain, "type": "unexpected", "reads": sorted(reads)})


def chrom_sort_key(chrom: str) -> Tuple[int, int | str]:
	"""Natural-ish chromosome ordering: chr1..chr22, chrX/Y/M, then others."""
	text = str(chrom)
	lower = text.lower()

	if lower.startswith("chr") and lower[3:].isdigit():
		return (0, int(lower[3:]))

	specials = {"chrx": 1000, "chry": 1001, "chrm": 1002, "chrmt": 1002}
	if lower in specials:
		return (1, specials[lower])

	return (2, lower)


def parse_read_list(raw_value: object) -> Set[str]:
	"""Parse list-like string from MUTATED_READS / REFERENCE_READS columns."""
	if pd.isna(raw_value):
		return set()

	text = str(raw_value).strip()
	if text in {"", "[]", "nan", "NaN", "None"}:
		return set()

	try:
		parsed = ast.literal_eval(text)
	except (SyntaxError, ValueError) as exc:
		raise ValueError(f"Could not parse read list: {text}") from exc

	if not isinstance(parsed, list):
		raise ValueError(f"Expected list-like value, got: {type(parsed).__name__}: {text}")

	return {str(x) for x in parsed if str(x).strip()}

def build_variant_constraints(df: pd.DataFrame ) -> pd.DataFrame:
	
	df["MUTATED_READS"] = df["MUTATED_READS"].apply(parse_read_list)
	df["REFERENCE_READS"] = df["REFERENCE_READS"].apply(parse_read_list)

	df = df[df.apply(lambda row: len(row["MUTATED_READS"]) > 0 or len(row["REFERENCE_READS"]) > 0, axis=1)]
	return df

def compute_continuity_assignments(variants_df: pd.DataFrame) -> Tuple[pd.DataFrame, List[dict]]:
	"""
	Assign continuity components by comparing adjacent mutation groups.
	For each variant, we compare the mutated and reference read sets to the previous component's chains.
	We assign the current variant's mutated group to the chain with which it has more overlap, and the reference group to the other chain.
	This way we can track continuity of mutation chains across variants
	"""

	ordered = variants_df.copy()
	ordered["_chrom_key"] = ordered["CHROM"].map(chrom_sort_key)
	ordered["_pos_num"] = pd.to_numeric(ordered["POS"], errors="coerce")
	ordered = ordered.sort_values(["_chrom_key", "_pos_num", "variant_id"], kind="stable").reset_index(drop=True)

	continuity_rows = []

	component_index = 0
	current_component = None

	for _, row in ordered.iterrows():
		variant_id = row["variant_id"]
		chrom = str(row["CHROM"])
		pos = row["POS"]
		mut_reads = row["MUTATED_READS"]
		ref_reads = row["REFERENCE_READS"]

		if current_component is None:
			current_component = Component(component_index, chrom, pos, "mutated", mut_reads, ref_reads, variant_id)
			component_index += 1

		elif chrom != current_component.chromosome:
			continuity_rows.append( current_component.__dict__)
			current_component = Component(component_index, chrom, pos, "mutated", mut_reads, ref_reads, variant_id)
			component_index += 1

		else:
			overlap_1_vs_mut = current_component.chain_1 & mut_reads
			overlap_1_vs_ref = current_component.chain_1 & ref_reads
			overlap_0_vs_mut = current_component.chain_0 & mut_reads
			overlap_0_vs_ref = current_component.chain_0 & ref_reads

			overlap_1_mut_0_ref = len(overlap_1_vs_mut) + len(overlap_0_vs_ref)
			overlap_1_ref_0_mut = len(overlap_1_vs_ref) + len(overlap_0_vs_mut)

			prev_0_status = "mutated" if len(overlap_0_vs_mut) > len(overlap_0_vs_ref) else "reference"
			prev_1_status = "mutated" if len(overlap_1_vs_mut) > len(overlap_1_vs_ref) else "reference"

			prev_0_status = "unknown" if len(overlap_0_vs_mut) == len(overlap_0_vs_ref) else prev_0_status
			prev_1_status = "unknown" if len(overlap_1_vs_mut) == len(overlap_1_vs_ref) else prev_1_status

			if not any(len(x) > 0 for x in [overlap_1_vs_mut, overlap_1_vs_ref, overlap_0_vs_mut, overlap_0_vs_ref]):
				print(f"New")
				# store current component before starting new one, to preserve continuity info
				continuity_rows.append( current_component.__dict__)

				current_component = Component(component_index, chrom, pos, "mutated", mut_reads, ref_reads, variant_id)
				component_index += 1

			elif prev_1_status == prev_0_status and prev_1_status == "mutated":

				# homozygous position where mutated group tracks both chains of previous component
				print(len(overlap_1_vs_mut), len(overlap_1_vs_ref),
		 				len(overlap_0_vs_mut), len(overlap_0_vs_ref))
				print(f"Continuity: {variant_id} mutated group tracks both chains of previous component")
				current_component.increase_chains_0_n_1(variant_id, "mutated", "mutated", set(), set())
				if overlap_1_vs_ref:
					print(f"  - {len(overlap_1_vs_ref)} reads in previous chain 1 are not mutated when they should be: {sorted(overlap_1_vs_ref)}")
					current_component.add_missing_mutation(variant_id, "chain_1", overlap_1_vs_ref)
					# TODO: add it as an outlier and document in which variant is this conflicting read observed
					
				if overlap_0_vs_ref:
					print(f"  - {len(overlap_0_vs_ref)} reads in previous chain 0 are not mutated when they should be: {sorted(overlap_0_vs_ref)}")
					current_component.add_missing_mutation(variant_id, "chain_0", overlap_0_vs_ref)

			elif prev_1_status == prev_0_status and prev_1_status == "reference":

				# homozygous position where mutated group tracks both chains of previous component
				print(len(overlap_1_vs_mut), len(overlap_1_vs_ref),
		 				len(overlap_0_vs_mut), len(overlap_0_vs_ref))
				print(f"Continuity: {variant_id} reference group tracks both chains of previous component")
				current_component.increase_chains_0_n_1(variant_id, "reference", "reference", set(), set())
				if overlap_1_vs_mut:
					print(f"  - {len(overlap_1_vs_mut)} reads in previous chain 1 are mutated when they should not be: {sorted(overlap_1_vs_mut)}")
					# TODO: add it as an outlier and document in which variant is this conflicting read observed
					current_component.add_extra_mutation(variant_id, "chain_1", overlap_1_vs_mut)
				if overlap_0_vs_mut:
					print(f"  - {len(overlap_0_vs_mut)} reads in previous chain 0 are mutated when they should not be: {sorted(overlap_0_vs_mut)}")
					current_component.add_extra_mutation(variant_id, "chain_0", overlap_0_vs_mut)



			elif overlap_1_mut_0_ref > overlap_1_ref_0_mut and prev_1_status == "mutated" and prev_0_status == "reference":
				# heterozygous position where mutated group tracks previous chain 1, and reference group tracks previous chain 0
				print(len(overlap_1_vs_mut), len(overlap_1_vs_ref),
		 			len(overlap_0_vs_mut), len(overlap_0_vs_ref))
				print(f"Continuity: {variant_id} mutated group tracks previous chain 1 (overlap {overlap_1_mut_0_ref} vs {overlap_1_ref_0_mut})")
				current_component.increase_chains_0_n_1(variant_id, "reference", "mutated", ref_reads, mut_reads)

				if overlap_1_vs_ref:
					print(f"  - {len(overlap_1_vs_ref)} reads in previous chain 1 overlap current reference group: {sorted(overlap_1_vs_ref)}")
					current_component.add_missing_mutation(variant_id, "chain_1", overlap_1_vs_ref)
				if overlap_0_vs_mut:
					print(f"  - {len(overlap_0_vs_mut)} reads in previous chain 0 overlap current mutated group: {sorted(overlap_0_vs_mut)}")
					current_component.add_extra_mutation(variant_id, "chain_0", overlap_0_vs_mut)


			elif overlap_1_mut_0_ref < overlap_1_ref_0_mut and prev_1_status == "reference" and prev_0_status == "mutated":
				# heterozygous position where mutated group tracks previous chain 0, and reference group tracks previous chain 1
				print(len(overlap_1_vs_mut), len(overlap_1_vs_ref),
		 			len(overlap_0_vs_mut), len(overlap_0_vs_ref))
				print(f"Continuity: {variant_id} mutated group tracks previous chain 0 (overlap {overlap_1_mut_0_ref} vs {overlap_1_ref_0_mut})")
				current_component.increase_chains_0_n_1(variant_id, "mutated", "reference", mut_reads, ref_reads)
				if overlap_1_vs_mut:
					print(f"  - {len(overlap_1_vs_mut)} reads in previous chain 1 overlap current mutated group: {sorted(overlap_1_vs_mut)}")
					current_component.add_extra_mutation(variant_id, "chain_1", overlap_1_vs_mut)
				if overlap_0_vs_ref:
					print(f"  - {len(overlap_0_vs_ref)} reads in previous chain 0 overlap current reference group: {sorted(overlap_0_vs_ref)}")
					current_component.add_missing_mutation(variant_id, "chain_0", overlap_0_vs_ref)

			else:
				print("Confusing point")
				print(len(overlap_1_vs_mut), len(overlap_1_vs_ref),
		 			len(overlap_0_vs_mut), len(overlap_0_vs_ref))
				print("********************************")
		print()
	
	continuity_rows.append( current_component.__dict__)

	return pd.DataFrame(continuity_rows)


def summarize_major_chain_by_chromosome(variant_summary_df: pd.DataFrame, output_pdf: Path) -> pd.DataFrame:
	"""Summarize major-chain mutated-read proportion per chromosome."""
	# use this code below to explore the cases for possible copy number changes
	
	data = variant_summary_df.copy()
	het_positions = data[(data["proportion_chain_0"] > 0.1) & (data["proportion_chain_0"] < 0.9)].reset_index(drop=True)
	het_positions = het_positions[(het_positions["total_size"] > 100)
				& (het_positions['REF'].str.len() == 1) & (het_positions['ALT'].str.len() == 1)
				].reset_index(drop=True)
	
	with PdfPages(output_pdf) as pdf:
		
		for chromo in het_positions['CHROM'].unique():
			chrom_positions = het_positions[het_positions['CHROM'] == chromo]
			starting_intervals = chrom_positions.drop_duplicates(subset=["number"])
			plt.scatter(chrom_positions.index,
						chrom_positions["minor_proportion"])
			plt.scatter(chrom_positions.index,
						chrom_positions["major_proportion"])
			plt.plot(chrom_positions.index,
					(chrom_positions["total_size"] / chrom_positions["total_size"].max()) * 0.9,
					)
			for indexx in starting_intervals.index.to_list():
				plt.axvline(x = indexx,
							color='red', alpha=0.3, linestyle='--')

			plt.ylim(0, 1)
			plt.title(chromo)
			pdf.savefig()
			plt.close()


def plot_major_chain_proportions(variant_summary_df: pd.DataFrame, output_png: Path) -> None:
	"""Plot overall and per-chromosome distributions of major-chain proportions."""
	plot_df = variant_summary_df[variant_summary_df["n_mut_reads"] > 0].copy()
	if plot_df.empty:
		return

	fig, axes = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)

	axes[0].hist(
		plot_df["mut_prop_major_chain"].dropna(),
		bins=20,
		alpha=0.85,
		color="#2f6f4f",
		edgecolor="black",
	)
	axes[0].set_title("Major-chain mutated-read proportion (all sites)")
	axes[0].set_xlabel("Proportion")
	axes[0].set_ylabel("Number of sites")
	axes[0].set_xlim(0, 1)

	chroms = sorted(plot_df["CHROM"].dropna().unique(), key=str)
	data_by_chrom = [plot_df.loc[plot_df["CHROM"] == c, "mut_prop_major_chain"].dropna().values for c in chroms]
	axes[1].boxplot(data_by_chrom, tick_labels=chroms, showfliers=False)
	axes[1].set_title("Major-chain proportion by chromosome")
	axes[1].set_xlabel("Chromosome")
	axes[1].set_ylabel("Proportion")
	axes[1].set_ylim(0, 1)
	axes[1].tick_params(axis="x", rotation=90)

	output_png.parent.mkdir(parents=True, exist_ok=True)
	fig.savefig(output_png, dpi=160)
	plt.close(fig)

def compute_continuity_metrics(df):
	df["total_size"] = df["size_chain_0"] + df["size_chain_1"]
	df["proportion_chain_0"] = df["size_chain_0"] / df["total_size"]
	df["proportion_chain_1"] = df["size_chain_1"] / df["total_size"]
	major_chain_proportions = df.groupby("number")["proportion_chain_0"].mean()
	major_chain_0 = major_chain_proportions[major_chain_proportions > 0.5].index
	df["major_chain"] = df["number"].isin(major_chain_0).map({True: "chain_0", False: "chain_1"})
	df["major_proportion"] = df["proportion_chain_1"].copy()
	df.loc[df["major_chain"] == 'chain_0', "major_proportion"] = df["proportion_chain_0"]
	df["minor_proportion"] = df["proportion_chain_1"].copy()
	df.loc[df["major_chain"] == 'chain_1', "minor_proportion"] = df["proportion_chain_0"]

	return df



def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description=(
			"Group reads into two mutation chains using per-position mutated/reference read evidence."
		)
	)
	parser.add_argument(
		"-i",
		"--input",
		required=True,
		help="Input TSV with CHROM, POS, REF, ALT, MUTATED_READS, REFERENCE_READS columns.",
	)
	parser.add_argument(
		"-o",
		"--output-prefix",
		default=None,
		help="Prefix for output files. Defaults to input path without extension.",
	)
	parser.add_argument(
		"--plot",
		default=None,
		help="Output PNG path for major-chain proportion plots. Defaults to <output_prefix>.chain_proportions.png",
	)
	return parser.parse_args()


def main() -> None:
	args = parse_args()
	input_path = Path(args.input)
	output_prefix = Path(args.output_prefix) if args.output_prefix else input_path.with_suffix("")

	df = pd.read_csv(input_path, sep="\t")
	required_cols = {"CHROM", "POS", "REF", "ALT", "MUTATED_READS", "REFERENCE_READS"}
	missing = sorted(required_cols - set(df.columns))
	if missing:
		raise ValueError(f"Missing required columns: {', '.join(missing)}")

	df["variant_id"] = df.apply(lambda row: f"{row['CHROM']}:{row['POS']}:{row['REF']}:{row['ALT']}", axis=1)

	detected_variants_df = build_variant_constraints(df)
	# read_chains_df = assign_read_chains(uf)
	# variant_summary_df = summarize_variants(variants_df, read_chains_df)
	continuity_df = compute_continuity_assignments(detected_variants_df)
	continuity_df_no_outliers = continuity_df.drop(columns=["outliers"])
	outliers = continuity_df[continuity_df["outliers"].apply(len) > 0]["outliers"]
	outliers_df = pd.DataFrame(outliers.to_list())
	outliers_df = pd.DataFrame(outliers_df[0].to_dict().values())
	continuity_df_expl = continuity_df_no_outliers.explode(["variant_id", "size_chain_1", "size_chain_0",
		"chain_0_status", "chain_1_status"
		]).reset_index(drop=True)

	continuity_df_with_metrics = compute_continuity_metrics(continuity_df_expl)
	
	variant_summary_df = df.merge(continuity_df_with_metrics, on="variant_id", how="left")
	try:
		summarize_major_chain_by_chromosome(variant_summary_df, output_prefix.parent / f"{output_prefix.name}.chain_proportions_by_chromosome.pdf")
	except Exception as e:
		print(f"Error occurred while summarizing major chain by chromosome: {e}")

	chains_out = output_prefix.parent / f"{output_prefix.name}.read_chains.tsv"
	variants_out = output_prefix.parent / f"{output_prefix.name}.variant_chain_support.tsv"
	outliers_out = output_prefix.parent / f"{output_prefix.name}.outliers.tsv"

	chains_out.parent.mkdir(parents=True, exist_ok=True)

	continuity_df.to_csv(chains_out, sep="\t", index=False)
	variant_summary_df.to_csv(variants_out, sep="\t", index=False)
	outliers_df.to_csv(outliers_out, sep="\t", index=False)

	print(f"Input variants: {len(df)}")



if __name__ == "__main__":
	main()