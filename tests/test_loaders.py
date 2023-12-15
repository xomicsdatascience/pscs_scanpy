import unittest
from pscs_scanpy.loaders import AnnDataLoadingNode
import anndata as ad
import pandas as pd
import os
from os.path import join

test_data_path = join(os.path.dirname(__file__), "test_data")


class TestAnnData(unittest.TestCase):
    def test_init(self):
        node = AnnDataLoadingNode()
        return

    def test_init_path(self):
        node = AnnDataLoadingNode(path="/path/to/file.h5ad")
        self.assertEqual(node.parameters["path"], "/path/to/file.h5ad")
        return

    def test_run(self):
        node = AnnDataLoadingNode(path=join(test_data_path, "anndata.h5ad"))
        node.run()
        self.assertIsNotNone(node.result[0])
        return

    def test_run_result(self):
        node = AnnDataLoadingNode(path=join(test_data_path, "anndata.h5ad"))
        node.run()
        expected_data = pd.DataFrame(data=[[1,2,3],[2,3,1],[3,1,2]],
                                     index=["first", "second", "third"],
                                     columns=["col0", "col1", "col2"])
        expected_data.index.name = "index"
        expected_data = ad.AnnData(expected_data)
        loaded_data = node.result
        self.assertTrue((expected_data.X == loaded_data.X).all())
        self.assertTrue((expected_data.obs_names == loaded_data.obs_names).all())
        self.assertTrue((expected_data.var_names == loaded_data.var_names).all())
        return

