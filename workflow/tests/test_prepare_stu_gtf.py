import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "prepare_stu_gtf.py"


class PrepareStuGtfTest(unittest.TestCase):
    def test_converts_transcripts_and_exons_with_stu_contigs_and_gene_ids(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            workdir = Path(tmpdir)
            input_path = workdir / "genes.gff3"
            output_path = workdir / "genes.gtf"
            input_path.write_text(
                "##gff-version 3\n"
                "1\tpgsc\tgene\t10\t30\t.\t+\t.\t"
                "ID=gene:gene_1;gene_id=gene_1\n"
                "1\tpgsc\tmRNA\t10\t30\t.\t+\t.\t"
                "ID=transcript:tx_1;Parent=gene:gene_1;transcript_id=tx_1\n"
                "1\tpgsc\texon\t10\t15\t.\t+\t.\t"
                "Parent=transcript:tx_1\n"
                "1\tpgsc\tncRNA_gene\t40\t50\t.\t-\t.\t"
                "ID=gene:ncrna_1;gene_id=ncrna_1\n"
                "1\tpgsc\tlnc_RNA\t40\t50\t.\t-\t.\t"
                "ID=transcript:nctx_1;Parent=gene:ncrna_1;transcript_id=nctx_1\n"
                "1\tpgsc\texon\t40\t50\t.\t-\t.\t"
                "Parent=transcript:nctx_1\n"
            )

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--input",
                    str(input_path),
                    "--output",
                    str(output_path),
                ],
                check=True,
            )

            self.assertEqual(
                output_path.read_text().splitlines(),
                [
                    'Stu_1\tpgsc\ttranscript\t10\t30\t.\t+\t.\tgene_id "gene_1"; transcript_id "tx_1";',
                    'Stu_1\tpgsc\texon\t10\t15\t.\t+\t.\tgene_id "gene_1"; transcript_id "tx_1";',
                    'Stu_1\tpgsc\ttranscript\t40\t50\t.\t-\t.\tgene_id "ncrna_1"; transcript_id "nctx_1";',
                    'Stu_1\tpgsc\texon\t40\t50\t.\t-\t.\tgene_id "ncrna_1"; transcript_id "nctx_1";',
                ],
            )


if __name__ == "__main__":
    unittest.main()
