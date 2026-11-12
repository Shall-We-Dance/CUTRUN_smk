"""DAG-only compatibility/validation checks; needs Snakemake and PyYAML."""
from copy import deepcopy
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import yaml

ROOT = Path(__file__).resolve().parents[1]
SNAKEMAKE = shutil.which("snakemake")


@unittest.skipUnless(SNAKEMAKE, "Requires Snakemake >=9.3.3 on PATH")
class WorkflowConfigTests(unittest.TestCase):
    def setUp(self):
        with (ROOT / "tests/config.tiny_replicates.yaml").open() as handle:
            self.config = yaml.safe_load(handle)

    def dry_run(self, config, expected_error=None):
        with tempfile.TemporaryDirectory(prefix="cutrun_dag_test_") as directory:
            path = Path(directory) / "config.yaml"
            config = deepcopy(config)
            config["output"]["dir"] = str(Path(directory) / "output")
            path.write_text(yaml.safe_dump(config))
            with tempfile.TemporaryFile(mode="w+t") as log:
                result = subprocess.run(
                    [SNAKEMAKE, "--snakefile", "workflow/Snakefile", "--configfile", str(path),
                     "--dry-run", "--cores", "2"], cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, timeout=60)
                log.seek(0)
                output = log.read()
            if expected_error is None:
                self.assertEqual(result.returncode, 0, output[-12000:])
            else:
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(expected_error, output)

    def test_invalid_settings_fail_before_execution(self):
        cases = []
        cfg = deepcopy(self.config)
        cfg["consensus"] = {"min_replicates": 3}
        cases.append(("support_exceeds_replicates", cfg, "min_replicates"))
        cfg = deepcopy(self.config)
        cfg["idr"] = {"enabled": True}
        cfg["macs3"]["broad"] = True
        cases.append(("broad_idr", cfg, "IDR requires"))
        cfg = deepcopy(self.config)
        cfg["idr"] = {"enabled": True}
        del cfg["samples"]["Tiny_ChIP"]["replicates"]["rep2"]
        cases.append(("single_idr", cfg, "IDR requires"))
        cfg = deepcopy(self.config)
        cfg["idr"] = {"enabled": True}
        cfg["macs3"]["enabled"] = False
        cases.append(("idr_without_macs", cfg, "IDR requires"))
        cfg = deepcopy(self.config)
        cfg["macs3"]["keep_dup"] = 1
        cases.append(("second_duplicate_filter", cfg, "keep_dup must be all"))
        cfg = deepcopy(self.config)
        cfg["fragments"] = {"se_extension": 0}
        cases.append(("zero_extension", cfg, "se_extension"))
        cfg = deepcopy(self.config)
        cfg["consensus"] = {"narrow": {"match": "any_overlap"}}
        cases.append(("ambiguous_match_rule", cfg, "any_overlap"))
        cfg = deepcopy(self.config)
        cfg["read_type"] = "SE"
        cfg["macs3"]["nomodel"] = False
        for rep in cfg["samples"]["Tiny_ChIP"]["replicates"].values():
            del rep["R2"]
        cases.append(("se_model_incompatible", cfg, "SE shared fragments require"))
        for name, cfg, error in cases:
            with self.subTest(case=name):
                self.dry_run(cfg, error)

    def test_optional_and_control_branches_build(self):
        # IDR disabled by default: no new caller is required.
        self.dry_run(self.config)
        cfg = deepcopy(self.config)
        cfg["macs3"]["enabled"] = False
        cfg["cutrun_qc"]["enabled"] = False
        self.dry_run(cfg)
        # Pooled IP with an explicitly named, single-replicate control.
        cfg = deepcopy(self.config)
        reads = deepcopy(cfg["samples"]["Tiny_ChIP"]["replicates"]["rep1"])
        cfg["samples"]["Tiny_INPUT"] = {"call_peaks": False, "replicates": {"input_a": reads}}
        cfg["samples"]["Tiny_ChIP"].update(
            control="Tiny_INPUT", control_replicates={"rep1": "input_a", "rep2": "pooled"})
        cfg["idr"] = {"enabled": True}
        self.dry_run(cfg)


if __name__ == "__main__":
    unittest.main()
