#!/usr/bin/env python3
"""三大组学管线 CLI 组装单测：mock subprocess，校验参数字符串与占位符流转。"""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from gibh_agent.core.workflows.epigenomics_workflow import EpigenomicsWorkflow
from gibh_agent.core.workflows.genomics_workflow import GenomicsWorkflow
from gibh_agent.core.workflows.proteomics_workflow import ProteomicsWorkflow
from gibh_agent.tools import omics_genomics_runner as g_runner
from gibh_agent.tools import omics_epigenomics_runner as e_runner
from gibh_agent.tools.omics_proteomics_pipeline_tools import (
    build_proteomics_database_search_cli,
)


class TestOmicsCliAssembly(unittest.TestCase):
    def test_bwa_mem_command_contains_penalties_and_threads(self) -> None:
        cmd = g_runner.build_bwa_mem_command(
            "/usr/bin/bwa",
            "/ref/hg38.fa",
            "/data/sample.fq.gz",
            threads=16,
            mismatch_penalty=5,
            gap_open_penalty=7,
        )
        self.assertEqual(cmd[0], "/usr/bin/bwa")
        self.assertIn("mem", cmd)
        self.assertIn("-t", cmd)
        self.assertIn("16", cmd)
        self.assertIn("-B", cmd)
        self.assertIn("5", cmd)
        self.assertIn("-O", cmd)
        idx_o = cmd.index("-O")
        self.assertEqual(cmd[idx_o + 1], "7,7")

    def test_gatk_hard_filtering_command_shape(self) -> None:
        cmd = g_runner.build_gatk_hard_filtering_command(
            "/opt/gatk/gatk",
            "/ref/hg38.fa",
            "/tmp/raw.vcf.gz",
            "/tmp/filtered.vcf.gz",
        )
        self.assertEqual(cmd[0], "/opt/gatk/gatk")
        self.assertIn("VariantFiltration", cmd)
        self.assertIn("--filter-expression", cmd)
        self.assertIn("QD < 2.0 || FS > 60.0 || MQ < 40.0", cmd)
        self.assertIn("--filter-name", cmd)
        self.assertIn("HardFiltered", cmd)
        self.assertIn("/tmp/raw.vcf.gz", cmd)
        self.assertIn("/tmp/filtered.vcf.gz", cmd)

    def test_gatk_haplotypecaller_contains_stand_call_conf(self) -> None:
        cmd = g_runner.build_gatk_haplotypecaller_command(
            "/opt/gatk/gatk",
            "/ref/hg38.fa",
            "/data/aln.bam",
            "/tmp/out.vcf.gz",
            stand_call_conf=35.0,
            min_mapping_quality=25,
        )
        self.assertIn("HaplotypeCaller", cmd)
        self.assertIn("--standard-min-confidence-threshold-for-calling", cmd)
        self.assertIn("35.0", cmd)
        self.assertIn("--minimum-mapping-quality", cmd)
        self.assertIn("25", cmd)

    @patch.object(g_runner.subprocess, "run")
    @patch.object(g_runner, "resolve_cli_exe")
    @patch.object(g_runner, "resolve_genomics_reference_fasta", return_value="/ref/hg38.fa")
    def test_try_germline_invokes_gatk_with_conf(
        self,
        _mock_ref_g: MagicMock,
        mock_resolve_cli: MagicMock,
        mock_run: MagicMock,
    ) -> None:
        mock_resolve_cli.side_effect = lambda name: "/bin/gatk" if name == "gatk" else None
        mock_run.return_value = MagicMock(returncode=0)
        with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as tf:
            bam = tf.name
        try:
            out = g_runner._try_germline(
                bam,
                "hg38",
                stand_call_conf=42.0,
                min_base_quality=18,
                min_mapping_quality=22,
            )
            self.assertIsNotNone(out)
            mock_run.assert_called()
            args, kwargs = mock_run.call_args
            cmd_list = args[0]
            joined = " ".join(cmd_list)
            self.assertIn("--standard-min-confidence-threshold-for-calling", joined)
            self.assertIn("42.0", joined)
            self.assertIn("/ref/hg38.fa", joined)
            self.assertIn(bam, joined)
        finally:
            Path(bam).unlink(missing_ok=True)

    def test_macs2_callpeak_command(self) -> None:
        cmd = e_runner.build_macs2_callpeak_command(
            "/usr/bin/macs2",
            "/data/trt.bam",
            "/tmp/macs_prefix",
            qvalue_threshold=0.01,
            broad_peak=True,
        )
        self.assertIn("callpeak", cmd)
        self.assertIn("-q", cmd)
        self.assertIn("0.01", cmd)
        self.assertIn("--broad", cmd)

    def test_bowtie2_command_threads_mp(self) -> None:
        cmd = e_runner.build_bowtie2_single_end_command(
            "/usr/bin/bowtie2",
            "/idx/hg38",
            "/data/atac.fq.gz",
            threads=12,
            mismatch_penalty=6,
        )
        self.assertIn("--threads", cmd)
        self.assertIn("12", cmd)
        self.assertIn("--mp", cmd)
        self.assertIn("6,6", cmd)

    def test_proteomics_search_cli(self) -> None:
        cmd = build_proteomics_database_search_cli(
            "/data/run.mzML",
            fragment_tol_da=0.02,
            missed_cleavages=3,
            peptide_fdr=0.005,
        )
        joined = " ".join(cmd)
        self.assertIn("missed-cleavages", joined)
        self.assertIn("3", joined)
        self.assertIn("peptide-fdr", joined)
        self.assertIn("0.005", joined)
        self.assertIn("0.02", joined)

    def test_workflow_templates_expose_algorithm_keys(self) -> None:
        gw = GenomicsWorkflow()
        meta_align = gw.get_step_metadata("step_genomics_align")
        params_a = meta_align["default_params"]
        self.assertIn("threads", params_a)
        self.assertIn("stand_call_conf", gw.get_step_metadata("step_genomics_germline")["default_params"])

        ew = EpigenomicsWorkflow()
        peak = ew.get_step_metadata("step_epi_peak")["default_params"]
        self.assertIn("qvalue_threshold", peak)
        self.assertIn("broad_peak", peak)

        pw = ProteomicsWorkflow()
        db = pw.get_step_metadata("step_prot_db_search")["default_params"]
        self.assertIn("missed_cleavages", db)
        self.assertIn("peptide_fdr", db)


if __name__ == "__main__":
    unittest.main()
