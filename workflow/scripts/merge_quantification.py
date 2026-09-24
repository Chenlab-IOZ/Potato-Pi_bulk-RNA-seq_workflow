#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge per-sample StringTie TPM and featureCounts outputs."
    )
    parser.add_argument("--samples", nargs="+", required=True)
    parser.add_argument("--tpm-files", nargs="+", required=True)
    parser.add_argument("--count-files", nargs="+", required=True)
    parser.add_argument("--output-tpm", required=True)
    parser.add_argument("--output-count", required=True)
    return parser.parse_args()


def read_stringtie_tpm(path):
    values = {}
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Gene ID", "TPM"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"Invalid StringTie gene abundance table: {path}")
        for row in reader:
            gene_id = row["Gene ID"]
            if gene_id in values:
                raise ValueError(f"Duplicate gene ID in {path}: {gene_id}")
            values[gene_id] = row["TPM"]
    return values


def read_featurecounts(path):
    with Path(path).open(newline="") as handle:
        lines = (line for line in handle if not line.startswith("#"))
        reader = csv.DictReader(lines, delimiter="\t")
        if reader.fieldnames is None or "Geneid" not in reader.fieldnames:
            raise ValueError(f"Invalid featureCounts table: {path}")
        count_column = reader.fieldnames[-1]
        values = {}
        for row in reader:
            gene_id = row["Geneid"]
            if gene_id in values:
                raise ValueError(f"Duplicate gene ID in {path}: {gene_id}")
            values[gene_id] = row[count_column]
    return values


def write_matrix(path, samples, per_sample_values):
    genes = sorted(set().union(*(values.keys() for values in per_sample_values)))
    with Path(path).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *samples])
        for gene_id in genes:
            writer.writerow(
                [gene_id, *(values.get(gene_id, "0") for values in per_sample_values)]
            )


def main():
    args = parse_args()
    expected = len(args.samples)
    if len(args.tpm_files) != expected or len(args.count_files) != expected:
        raise ValueError(
            "The numbers of samples, TPM files, and count files must match"
        )

    tpm_values = [read_stringtie_tpm(path) for path in args.tpm_files]
    count_values = [read_featurecounts(path) for path in args.count_files]
    write_matrix(args.output_tpm, args.samples, tpm_values)
    write_matrix(args.output_count, args.samples, count_values)


if __name__ == "__main__":
    main()
