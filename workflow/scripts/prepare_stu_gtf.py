#!/usr/bin/env python3
import argparse
from pathlib import Path


GENE_FEATURES = {"gene", "ncRNA_gene"}


def parse_attributes(text):
    attributes = {}
    for field in text.rstrip(";").split(";"):
        if not field:
            continue
        key, separator, value = field.partition("=")
        if separator:
            attributes[key] = value
    return attributes


def strip_prefix(value, prefix):
    return value[len(prefix) :] if value.startswith(prefix) else value


def gtf_attributes(gene_id, transcript_id):
    return f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'


def convert_gff3(input_path, output_path):
    genes = {}
    transcripts = {}
    transcript_count = 0
    exon_count = 0

    with Path(input_path).open() as source, Path(output_path).open("w") as target:
        for line_number, line in enumerate(source, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Invalid GFF3 line {line_number}: expected 9 fields")

            seqname, source_name, feature, start, end, score, strand, phase, raw = fields
            attributes = parse_attributes(raw)
            seqname = seqname if seqname.startswith("Stu_") else f"Stu_{seqname}"

            if feature in GENE_FEATURES:
                source_id = attributes.get("ID")
                gene_id = attributes.get("gene_id")
                if source_id and gene_id:
                    genes[source_id] = gene_id
                continue

            if all(key in attributes for key in ("ID", "Parent", "transcript_id")):
                source_id = attributes.get("ID")
                parent = attributes.get("Parent")
                transcript_id = attributes.get("transcript_id")
                if not source_id or not parent or not transcript_id:
                    raise ValueError(
                        f"Transcript lacks ID, Parent, or transcript_id at line {line_number}"
                    )
                gene_id = genes.get(parent, strip_prefix(parent, "gene:"))
                transcripts[source_id] = (gene_id, transcript_id)
                output = [
                    seqname,
                    source_name,
                    "transcript",
                    start,
                    end,
                    score,
                    strand,
                    ".",
                    gtf_attributes(gene_id, transcript_id),
                ]
                target.write("\t".join(output) + "\n")
                transcript_count += 1
                continue

            if feature != "exon":
                continue

            parents = attributes.get("Parent", "").split(",")
            for parent in filter(None, parents):
                if parent not in transcripts:
                    raise ValueError(
                        f"Exon parent {parent!r} is unknown at line {line_number}"
                    )
                gene_id, transcript_id = transcripts[parent]
                output = [
                    seqname,
                    source_name,
                    "exon",
                    start,
                    end,
                    score,
                    strand,
                    ".",
                    gtf_attributes(gene_id, transcript_id),
                ]
                target.write("\t".join(output) + "\n")
                exon_count += 1

    if transcript_count == 0 or exon_count == 0:
        raise ValueError("Converted GTF has no transcripts or exons")


def main():
    parser = argparse.ArgumentParser(
        description="Convert the SolTub 3.0 GFF3 annotation to Stu-prefixed GTF."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    convert_gff3(args.input, args.output)


if __name__ == "__main__":
    main()
