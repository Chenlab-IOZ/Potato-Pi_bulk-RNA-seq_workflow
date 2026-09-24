import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "merge_quantification.py"


class MergeQuantificationTest(unittest.TestCase):
    def test_merges_tpm_and_readcount_in_requested_sample_order(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            workdir = Path(tmpdir)
            tpm_a = workdir / "sample_a.gene_abund.tsv"
            tpm_b = workdir / "sample_b.gene_abund.tsv"
            count_a = workdir / "sample_a.featurecounts.tsv"
            count_b = workdir / "sample_b.featurecounts.tsv"
            output_tpm = workdir / "gene_tpm.tsv"
            output_count = workdir / "gene_readcount.tsv"

            tpm_a.write_text(
                "Gene ID\tGene Name\tReference\tStrand\tStart\tEnd\tCoverage\tFPKM\tTPM\n"
                "gene_2\tgene_2\tStu_2\t+\t20\t30\t1\t2\t3.5\n"
                "gene_1\tgene_1\tStu_1\t+\t10\t15\t1\t2\t4.5\n"
            )
            tpm_b.write_text(
                "Gene ID\tGene Name\tReference\tStrand\tStart\tEnd\tCoverage\tFPKM\tTPM\n"
                "gene_1\tgene_1\tStu_1\t+\t10\t15\t1\t2\t8\n"
            )
            count_a.write_text(
                "# Program:featureCounts\n"
                "Geneid\tChr\tStart\tEnd\tStrand\tLength\ta.bam\n"
                "gene_1\tStu_1\t10\t15\t+\t6\t11\n"
            )
            count_b.write_text(
                "Geneid\tChr\tStart\tEnd\tStrand\tLength\tb.bam\n"
                "gene_2\tStu_2\t20\t30\t+\t11\t7\n"
            )

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--samples",
                    "sample_b",
                    "sample_a",
                    "--tpm-files",
                    str(tpm_b),
                    str(tpm_a),
                    "--count-files",
                    str(count_b),
                    str(count_a),
                    "--output-tpm",
                    str(output_tpm),
                    "--output-count",
                    str(output_count),
                ],
                check=True,
            )

            with output_tpm.open(newline="") as handle:
                tpm_rows = list(csv.reader(handle, delimiter="\t"))
            with output_count.open(newline="") as handle:
                count_rows = list(csv.reader(handle, delimiter="\t"))

            self.assertEqual(tpm_rows[0], ["gene_id", "sample_b", "sample_a"])
            self.assertEqual(tpm_rows[1:], [["gene_1", "8", "4.5"], ["gene_2", "0", "3.5"]])
            self.assertEqual(count_rows[0], ["gene_id", "sample_b", "sample_a"])
            self.assertEqual(count_rows[1:], [["gene_1", "0", "11"], ["gene_2", "7", "0"]])


if __name__ == "__main__":
    unittest.main()
