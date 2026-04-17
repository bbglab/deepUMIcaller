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
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import pandas as pd


class ParityUnionFind:
	"""Union-find with parity (xor) to represent same/opposite relations."""

	def __init__(self) -> None:
		self.parent: Dict[str, str] = {}
		self.rank: Dict[str, int] = {}
		self.parity_to_parent: Dict[str, int] = {}

	def add(self, x: str) -> None:
		if x not in self.parent:
			self.parent[x] = x
			self.rank[x] = 0
			self.parity_to_parent[x] = 0

	def find(self, x: str) -> Tuple[str, int]:
		"""Return (root, parity from x to root)."""
		if self.parent[x] == x:
			return x, 0

		px = self.parent[x]
		root, p_to_root = self.find(px)
		self.parity_to_parent[x] ^= p_to_root
		self.parent[x] = root
		return self.parent[x], self.parity_to_parent[x]

	def union(self, a: str, b: str, relation: int) -> bool:
		"""Union a and b with relation: 0=same, 1=opposite.

		Returns False when a contradiction is found.
		"""
		self.add(a)
		self.add(b)

		ra, pa = self.find(a)
		rb, pb = self.find(b)

		if ra == rb:
			return (pa ^ pb) == relation

		# Attach lower-rank tree under higher-rank tree.
		if self.rank[ra] < self.rank[rb]:
			ra, rb = rb, ra
			pa, pb = pb, pa

		self.parent[rb] = ra
		# We need: parity(a, root_a) ^ parity(b, root_b) ^ parity(root_b, root_a) == relation
		self.parity_to_parent[rb] = pa ^ pb ^ relation

		if self.rank[ra] == self.rank[rb]:
			self.rank[ra] += 1

		return True


def parse_read_list(raw_value: object) -> List[str]:
	"""Parse list-like string from MUTATED_READS / REFERENCE_READS columns."""
	if pd.isna(raw_value):
		return []

	text = str(raw_value).strip()
	if text in {"", "[]", "nan", "NaN", "None"}:
		return []

	try:
		parsed = ast.literal_eval(text)
	except (SyntaxError, ValueError) as exc:
		raise ValueError(f"Could not parse read list: {text}") from exc

	if not isinstance(parsed, list):
		raise ValueError(f"Expected list-like value, got: {type(parsed).__name__}: {text}")

	return [str(x) for x in parsed if str(x).strip()]


def add_same_constraints(uf: ParityUnionFind, reads: Iterable[str], contradictions: List[dict], variant_id: str) -> None:
	reads = list(reads)
	if len(reads) < 2:
		return

	pivot = reads[0]
	for read in reads[1:]:
		if not uf.union(pivot, read, 0):
			contradictions.append(
				{
					"variant_id": variant_id,
					"read_a": pivot,
					"read_b": read,
					"expected_relation": "same",
				}
			)


def build_constraints(df: pd.DataFrame, uf: ParityUnionFind) -> Tuple[pd.DataFrame, List[dict]]:
	rows = []
	contradictions: List[dict] = []

	for _, row in df.iterrows():
		variant_id = f"{row['CHROM']}:{row['POS']}:{row['REF']}:{row['ALT']}"
		mut_reads = parse_read_list(row["MUTATED_READS"])
		ref_reads = parse_read_list(row["REFERENCE_READS"])

		for r in mut_reads:
			uf.add(r)
		for r in ref_reads:
			uf.add(r)

		add_same_constraints(uf, mut_reads, contradictions, variant_id)
		add_same_constraints(uf, ref_reads, contradictions, variant_id)

		# Opposite relation between mutated and reference reads.
		for rm in mut_reads:
			for rr in ref_reads:
				if not uf.union(rm, rr, 1):
					contradictions.append(
						{
							"variant_id": variant_id,
							"read_a": rm,
							"read_b": rr,
							"expected_relation": "opposite",
						}
					)

		rows.append(
			{
				"variant_id": variant_id,
				"CHROM": row["CHROM"],
				"POS": row["POS"],
				"REF": row["REF"],
				"ALT": row["ALT"],
				"mutated_reads": mut_reads,
				"reference_reads": ref_reads,
			}
		)

	return pd.DataFrame(rows), contradictions


def assign_read_chains(uf: ParityUnionFind) -> pd.DataFrame:
	components: Dict[str, List[Tuple[str, int]]] = defaultdict(list)
	for read_id in uf.parent:
		root, parity = uf.find(read_id)
		components[root].append((read_id, parity))

	records = []
	for i, (root, read_infos) in enumerate(sorted(components.items(), key=lambda x: (-len(x[1]), x[0])), start=1):
		component_id = f"component_{i}"
		for read_id, parity in sorted(read_infos):
			records.append(
				{
					"read_id": read_id,
					"component_root": root,
					"component_id": component_id,
					"chain": f"chain_{parity}",
					"chain_index": parity,
				}
			)

	return pd.DataFrame(records)


def summarize_variants(variants_df: pd.DataFrame, read_chains_df: pd.DataFrame) -> pd.DataFrame:
	read2chain = dict(zip(read_chains_df["read_id"], read_chains_df["chain_index"]))

	summary_rows = []
	for _, row in variants_df.iterrows():
		mut_counts = [0, 0]
		ref_counts = [0, 0]

		for read_id in row["mutated_reads"]:
			if read_id in read2chain:
				mut_counts[read2chain[read_id]] += 1
		for read_id in row["reference_reads"]:
			if read_id in read2chain:
				ref_counts[read2chain[read_id]] += 1

		n_mut = len(row["mutated_reads"])
		n_ref = len(row["reference_reads"])
		mut_prop_chain_0 = (mut_counts[0] / n_mut) if n_mut > 0 else float("nan")
		mut_prop_chain_1 = (mut_counts[1] / n_mut) if n_mut > 0 else float("nan")
		ref_prop_chain_0 = (ref_counts[0] / n_ref) if n_ref > 0 else float("nan")
		ref_prop_chain_1 = (ref_counts[1] / n_ref) if n_ref > 0 else float("nan")

		# chain_0/chain_1 labels are arbitrary per connected component; major/minor is orientation-free.
		mut_major = max(mut_counts)
		mut_minor = min(mut_counts)
		mut_prop_major = (mut_major / n_mut) if n_mut > 0 else float("nan")
		mut_prop_minor = (mut_minor / n_mut) if n_mut > 0 else float("nan")

		mutated_support = max(mut_counts)
		if mutated_support == 0 or mut_counts[0] == mut_counts[1]:
			mutated_chain = "undetermined"
		else:
			mutated_chain = f"chain_{0 if mut_counts[0] > mut_counts[1] else 1}"

		summary_rows.append(
			{
				"variant_id": row["variant_id"],
				"CHROM": row["CHROM"],
				"POS": row["POS"],
				"REF": row["REF"],
				"ALT": row["ALT"],
				"n_mut_reads": n_mut,
				"n_ref_reads": n_ref,
				"mut_chain_0": mut_counts[0],
				"mut_chain_1": mut_counts[1],
				"ref_chain_0": ref_counts[0],
				"ref_chain_1": ref_counts[1],
				"mut_prop_chain_0": mut_prop_chain_0,
				"mut_prop_chain_1": mut_prop_chain_1,
				"ref_prop_chain_0": ref_prop_chain_0,
				"ref_prop_chain_1": ref_prop_chain_1,
				"mut_prop_major_chain": mut_prop_major,
				"mut_prop_minor_chain": mut_prop_minor,
				"assigned_mutated_chain": mutated_chain,
			}
		)

	return pd.DataFrame(summary_rows)


def summarize_major_chain_by_chromosome(variant_summary_df: pd.DataFrame) -> pd.DataFrame:
	"""Summarize major-chain mutated-read proportion per chromosome."""
	plot_df = variant_summary_df[variant_summary_df["n_mut_reads"] > 0].copy()
	if plot_df.empty:
		return pd.DataFrame(
			columns=[
				"CHROM",
				"n_sites",
				"mean_mut_prop_major_chain",
				"median_mut_prop_major_chain",
				"std_mut_prop_major_chain",
			]
		)

	summary = (
		plot_df.groupby("CHROM", as_index=False)["mut_prop_major_chain"]
		.agg(["count", "mean", "median", "std"])
		.reset_index()
		.rename(
			columns={
				"count": "n_sites",
				"mean": "mean_mut_prop_major_chain",
				"median": "median_mut_prop_major_chain",
				"std": "std_mut_prop_major_chain",
			}
		)
	)

	return summary.sort_values(["n_sites", "CHROM"], ascending=[False, True])


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

	uf = ParityUnionFind()
	variants_df, contradictions = build_constraints(df, uf)
	read_chains_df = assign_read_chains(uf)
	variant_summary_df = summarize_variants(variants_df, read_chains_df)
	per_chrom_summary_df = summarize_major_chain_by_chromosome(variant_summary_df)

	reads_out = output_prefix.parent / f"{output_prefix.name}.read_chains.tsv"
	variants_out = output_prefix.parent / f"{output_prefix.name}.variant_chain_support.tsv"
	per_chrom_out = output_prefix.parent / f"{output_prefix.name}.major_chain_per_chromosome.tsv"
	contradictions_out = output_prefix.parent / f"{output_prefix.name}.phasing_contradictions.tsv"
	plot_out = Path(args.plot) if args.plot else (output_prefix.parent / f"{output_prefix.name}.chain_proportions.png")

	reads_out.parent.mkdir(parents=True, exist_ok=True)

	read_chains_df.to_csv(reads_out, sep="\t", index=False)
	variant_summary_df.to_csv(variants_out, sep="\t", index=False)
	per_chrom_summary_df.to_csv(per_chrom_out, sep="\t", index=False)
	pd.DataFrame(contradictions).to_csv(contradictions_out, sep="\t", index=False)
	plot_major_chain_proportions(variant_summary_df, plot_out)

	print(f"Input variants: {len(df)}")
	print(f"Phased reads: {len(read_chains_df)}")
	print(f"Connected components: {read_chains_df['component_id'].nunique() if not read_chains_df.empty else 0}")
	print(f"Contradictions: {len(contradictions)}")
	print(f"Wrote: {reads_out}")
	print(f"Wrote: {variants_out}")
	print(f"Wrote: {per_chrom_out}")
	print(f"Wrote: {contradictions_out}")
	if (variant_summary_df["n_mut_reads"] > 0).any():
		print(f"Wrote: {plot_out}")


if __name__ == "__main__":
	main()