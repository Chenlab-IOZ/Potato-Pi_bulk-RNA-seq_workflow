"""Exercise method selection through an actual Snakemake DAG."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

WORKFLOW = Path(__file__).resolve().parents[1]
SNAKEMAKE = shutil.which("snakemake")


@unittest.skipUnless(SNAKEMAKE, "snakemake is required")
class SampleRelationshipSelectionTest(unittest.TestCase):
    def dry_run(self, method):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "03.expression/matrices").mkdir(parents=True)
            (root / "03.expression/matrices/gene_readcount.tsv").touch()
            (root / "samples.tsv").touch()
            options = dict(min_total_count=10, top_variable_genes=500,
                           correlation_method="pearson", clustering_method="average",
                           memory_mb=8192)
            if method is not None:
                options["method"] = method
            snakefile = root / "Snakefile"
            snakefile.write_text(
                "from snakemake.exceptions import WorkflowError\n"
                f"WORK_DIR = {str(root)!r}\n"
                f"SAMPLE_FILE = {str(root / 'samples.tsv')!r}\n"
                f"SCRIPT_PATH = {str(WORKFLOW / 'scripts')!r}\n"
                f"config = {{'sample_relationship': {options!r}}}\n"
                f"include: {str(WORKFLOW / 'rules/sample_relationship.smk')!r}\n"
            )
            return subprocess.run(
                [SNAKEMAKE, "-s", str(snakefile), "sample_relationship_all", "-n", "--cores", "2"],
                cwd=root, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            )

    def test_method_selection(self):
        for method, selected in [("deseq2", ["deseq2"]), ("edgeR", ["edgeR"]),
                                 ("both", ["deseq2", "edgeR"]), (None, ["edgeR"])]:
            with self.subTest(method=method):
                result = self.dry_run(method)
                self.assertEqual(result.returncode, 0, result.stdout)
                for candidate in ["deseq2", "edgeR"]:
                    self.assertEqual(f"wildcards: method={candidate}" in result.stdout,
                                     candidate in selected, result.stdout)

    def test_invalid_method(self):
        result = self.dry_run("invalid")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("sample_relationship.method must be", result.stdout)


if __name__ == "__main__":
    unittest.main()
