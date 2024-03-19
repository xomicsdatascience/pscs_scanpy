import anndata
import pandas as pd
import os
from pscs_api import InputNode


class AnnDataLoadingNode(InputNode):
    important_parameters = []

    def __init__(self,
                 path: str = None):
        """
        Loads an AnnData file (.h5ad) into memory.
        Parameters
        ----------
        path : str
            Path to the file to load.
        """
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        # Load h5ad
        data = anndata.read_h5ad(self.parameters["path"])
        self._terminate(result=data)
        return
