"""STED-EC 细胞类型探针与自动标注兜底单元测试。"""
import os
import sys
import tempfile
import unittest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)


class TestStedEcAutoCellType(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import anndata  # noqa: F401
            import numpy as np  # noqa: F401
            import scanpy as sc  # noqa: F401
        except ImportError:
            raise unittest.SkipTest("scanpy/anndata 未安装，跳过 STED-EC 自动标注测试")

    def _make_synthetic_h5ad(
        self,
        *,
        with_cell_type: bool,
        cell_type_col: str = "cell_type",
        extra_obs=None,
    ):
        import numpy as np
        import scanpy as sc

        n_obs, n_vars = 120, 80
        rng = np.random.default_rng(42)
        X = rng.poisson(3, size=(n_obs, n_vars)).astype(np.float32)
        genes = [f"GENE{i}" for i in range(n_vars)]
        for g in ("CD3D", "CD3E", "MS4A1", "CD79A", "LYZ", "CD14"):
            if g in genes:
                j = genes.index(g)
                X[:40, j] += 8
                X[40:80, j] += 2
                X[80:, j] += 1

        adata = sc.AnnData(X=X)
        adata.var_names = genes
        adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
        adata.obs["day"] = [0.0] * 40 + [1.0] * 40 + [2.0] * 40
        if with_cell_type:
            adata.obs[cell_type_col] = ["TypeA"] * 40 + ["TypeB"] * 40 + ["TypeC"] * 40
        if extra_obs:
            for k, v in extra_obs.items():
                adata.obs[k] = v
        return adata

    def test_probe_existing_cell_type_unchanged(self):
        from gibh_agent.tools.sted_ec_tools import _probe_cell_type_key, sted_ec_data_validation
        import scanpy as sc

        with tempfile.TemporaryDirectory() as td:
            path = os.path.join(td, "with_ct.h5ad")
            adata = self._make_synthetic_h5ad(with_cell_type=True)
            adata.write(path)
            original_labels = adata.obs["cell_type"].tolist()

            self.assertEqual(_probe_cell_type_key(adata), "cell_type")
            res = sted_ec_data_validation(h5ad_path=path, time_key="day", cell_type_key="cell_type")
            self.assertEqual(res.get("status"), "success")
            self.assertEqual(res.get("cell_type_key"), "cell_type")
            self.assertEqual(os.path.normpath(res["h5ad_path"]), os.path.normpath(path))
            resolution = res.get("cell_type_resolution") or {}
            self.assertTrue(resolution.get("skipped_auto_annotation"))
            self.assertNotIn("auto_annotated", resolution)

            out = sc.read_h5ad(path)
            self.assertEqual(out.obs["cell_type"].tolist(), original_labels)

    def test_probe_alias_skips_auto_annotation_and_preserves_file(self):
        from gibh_agent.tools.sted_ec_tools import sted_ec_data_validation
        import scanpy as sc

        with tempfile.TemporaryDirectory() as td:
            path = os.path.join(td, "alias.h5ad")
            adata = self._make_synthetic_h5ad(with_cell_type=True, cell_type_col="celltype")
            adata.write(path)

            res = sted_ec_data_validation(h5ad_path=path, time_key="day", cell_type_key="cell_type")
            self.assertEqual(res.get("status"), "success")
            self.assertEqual(res.get("cell_type_key"), "celltype")
            self.assertEqual(os.path.normpath(res["h5ad_path"]), os.path.normpath(path))
            self.assertTrue((res.get("cell_type_resolution") or {}).get("skipped_auto_annotation"))

            out = sc.read_h5ad(path)
            self.assertIn("celltype", out.obs.columns)
            self.assertNotIn("cell_type", out.obs.columns)

    def test_missing_cell_type_triggers_auto_annotation(self):
        from gibh_agent.tools.sted_ec_tools import sted_ec_data_validation
        import scanpy as sc

        with tempfile.TemporaryDirectory() as td:
            path = os.path.join(td, "no_ct.h5ad")
            adata = self._make_synthetic_h5ad(with_cell_type=False)
            adata.write(path)

            res = sted_ec_data_validation(h5ad_path=path, time_key="day", cell_type_key="cell_type")
            self.assertEqual(res.get("status"), "success")
            self.assertEqual(res.get("cell_type_key"), "cell_type")
            resolution = res.get("cell_type_resolution") or {}
            self.assertTrue(resolution.get("auto_annotated"))
            self.assertNotEqual(os.path.normpath(res["h5ad_path"]), os.path.normpath(path))

            out = sc.read_h5ad(res["h5ad_path"])
            self.assertIn("cell_type", out.obs.columns)
            self.assertGreater(out.obs["cell_type"].nunique(), 0)

    def test_leiden_only_does_not_skip_auto_annotation(self):
        from gibh_agent.tools.sted_ec_tools import _probe_cell_type_key, sted_ec_data_validation

        with tempfile.TemporaryDirectory() as td:
            path = os.path.join(td, "leiden_only.h5ad")
            adata = self._make_synthetic_h5ad(
                with_cell_type=False,
                extra_obs={"leiden": ["0"] * 60 + ["1"] * 60},
            )
            adata.write(path)

            self.assertIsNone(_probe_cell_type_key(adata))
            res = sted_ec_data_validation(h5ad_path=path, time_key="day", cell_type_key="cell_type")
            self.assertEqual(res.get("status"), "success")
            self.assertEqual(res.get("cell_type_key"), "cell_type")
            self.assertTrue((res.get("cell_type_resolution") or {}).get("auto_annotated"))


if __name__ == "__main__":
    unittest.main()
