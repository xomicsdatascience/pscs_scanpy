from typing import Collection, Optional, Literal, Union, Sequence
import math
import scanpy as sc
from scanpy import preprocessing as pp
from pscs_api import PipelineNode
from pscs_api.base import InteractionList, Interaction, istr

class CalculateQCMetrics(PipelineNode):
    important_parameters = ["expr_type", "var_type"]
    requirements = InteractionList(var=[istr("qc_vars")],
                                   layers=[istr("layer")])
    effects = InteractionList(obs=["total_"+istr("var_type")+"_by_"+istr("expr_type"),
                                   "total_"+istr("expr_type"),
                                   f"pct_{istr('expr_type')}_in_top_{istr('percent_top')}_{istr('var_type')}",
                                   f"total_{istr('expr_type')}_{istr('qc_vars')}",
                                   f"pct_{istr('expr_type')}_{istr('qc_vars')}"],
                              var=[f"total_{istr('expr_type')}",
                                   f"n_genes_by_{istr('expr_type')}",
                                   f"mean_{istr('expr_type')}",
                                   f"n_cells_by_{istr('expr_type')}",
                                   f"pct_dropout_by_{istr('expr_type')}"])

    def __init__(self,
                 expr_type: str = "counts",
                 var_type: str = "genes",
                 qc_vars: Collection[str] = (),
                 percent_top: Optional[Collection[int]] = (50, 100, 200, 500),
                 layer: Optional[str] = None,
                 use_raw: bool = False,
                 log1p: bool = True):
        # Assign all arguments to self.parameters
        super().__init__()
        self.store_vars_as_parameters(**vars(), inplace=True)
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.calculate_qc_metrics(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class FilterCells(PipelineNode):
    important_parameters = ["min_counts", "max_counts", "min_genes", "max_genes"]
    effects = InteractionList(Interaction(obs=["n_genes"]),
                              Interaction(obs=["n_counts"]))
    effects = effects * effects

    def __init__(self,
                 min_counts: Optional[int] = None,
                 min_genes: Optional[int] = None,
                 max_counts: Optional[int] = None,
                 max_genes: Optional[int] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars(), inplace=True)
        return

    def run(self):
        ann_data = self._previous[0].result
        # filter_cells only accepts one arg at a time; go through each one
        other_params = {}
        functional_params = {"min_counts", "min_genes", "max_counts", "max_genes"}
        for par_key, par_val in self.parameters.items():
            if par_key not in functional_params:
                other_params[par_key] = par_val
        for fpar in functional_params:
            fdict = {fpar: self.parameters[fpar]}
            if fdict[fpar] is not None:
                _ = sc.pp.filter_cells(ann_data, **fdict, **other_params)
        self._terminate(ann_data)
        return


class FilterGenes(PipelineNode):
    important_parameters = ["min_counts", "max_counts", "min_cells", "max_cells"]
    effects = InteractionList(Interaction(var=["n_counts"]),
                              Interaction(var=["n_cells"]))
    effects = effects * effects

    def __init__(self,
                 min_counts: Optional[int] = None,
                 min_cells: Optional[int] = None,
                 max_counts: Optional[int] = None,
                 max_cells: Optional[int] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars(), inplace=True)
        return

    def run(self):
        ann_data = self._previous[0].result
        # filter_genes only accepts one arg at a time; go through each one
        other_params = {}
        functional_params = {"min_counts", "min_cells", "max_counts", "max_cells"}
        for par_key, par_val in self.parameters.items():
            if par_key not in functional_params:
                other_params[par_key] = par_val
        for fpar in functional_params:
            fdict = {fpar: self.parameters[fpar]}
            if fdict[fpar] is not None:
                _ = sc.pp.filter_genes(ann_data, **fdict, **other_params)  # inplace == True is in other_params
        self._terminate(ann_data)
        return


class HighlyVariableGenes(PipelineNode):
    important_parameters = ["n_top_genes"]
    requirements = InteractionList(layers=[istr("layer")],
                                   obs=[istr("batch_key")])
    effects = InteractionList(var=["highly_variable", "means", "dispersions", "dispersions_norm", "variances", "variances_norm",
                                   "highly_variable_rank", "highly_variable_nbatches", "highly_variable_intersection"])

    def __init__(self,
                 layer: Optional[str] = None,
                 n_top_genes: Optional[int] = None,
                 min_mean: Optional[float] = 0.0125,
                 max_mea: Optional[float] = 3,
                 min_disp: Optional[float] = 0.5,
                 max_disp: Optional[float] = str(math.inf),
                 span: Optional[float] = 0.3,
                 n_bins: int = 20,
                 flavor: Literal['seurat', 'cell_ranger', 'seurat_v3'] = 'seurat',
                 subset: bool = False,
                 batch_key: Optional[str] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars(), inplace=True)
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.highly_variable_genes(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class Log1p(PipelineNode):
    important_parameters = ["base"]
    requirements = InteractionList(layers=[istr("layer")],
                                   obsm=[istr("obsm")])

    def __init__(self,
                 base: Optional[float] = None,
                 chunked: Optional[bool] = None,
                 chunk_size: Optional[int] = None,
                 layer: Optional[str] = None,
                 obsm: Optional[str] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.log1p(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class PCA(PipelineNode):
    important_parameters = ["n_comps", "zero_center"]
    requirements = InteractionList(layers=[istr("layer")],
                                   var=[istr("mask_var")])
    effects = InteractionList(obsm=["X_pca"],
                              varm=["PCs"],
                              uns=["pca"])

    def __init__(self,
                 n_comps: Optional[int] = None,
                 zero_center: Optional[bool] = True,
                 svd_solver: Literal["arpack", "randomized", "auto", "lobpcg"] = "arpack",
                 random_state: Union[None, int] = 0,
                 use_highly_variable: Optional[bool] = None,
                 dtype: str = 'float32',
                 chunked: bool = False,
                 chunk_size: Optional[int] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.pca(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class NormalizeTotal(PipelineNode):
    important_parameters = ["target_sum", "exclude_highly_expressed"]
    requirements = InteractionList(layers=[istr("layer")])
    effects = InteractionList(obs=[istr("key_added")])

    def __init__(self,
                 target_sum: Optional[float] = None,
                 exclude_highly_expressed: bool = False,
                 max_fraction: float = 0.05,
                 key_added: Optional[str] = None,
                 layer: Optional[str] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars(), inplace=True)
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.normalize_total(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class RegressOut(PipelineNode):
    important_parameters = ["keys"]
    requirements = InteractionList(obs=[istr("keys")],
                                   layers=[istr("layer")])
    effects = InteractionList(layers=[istr("layer")])

    def __init__(self,
                 keys: Union[str, Sequence[str]] = None,
                 n_jobs: Optional[int] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.regress_out(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class Scale(PipelineNode):
    important_parameters = ["zero_center"]
    requirements = InteractionList(layers=[istr("layer")],
                                   obsm=[istr("obsm")],
                                   obs=[istr("mask_obs")])
    effects = InteractionList(layers=[istr("layer")],
                              var=["mean", "std", "var"])

    def __init__(self,
                 zero_center: bool = True,
                 max_value: Optional[float] = None,
                 layer: Optional[str] = None,
                 obsm: Optional[str] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.scale(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class Subsample(PipelineNode):
    important_parameters = ["fraction", "n_obs"]

    def __init__(self,
                 fraction: Optional[float] = None,
                 n_obs: Optional[int] = None,
                 random_state: Union[None, int] = 0):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.subsample(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class DownsampleCounts(PipelineNode):
    important_parameters = ["counts_per_cell", "total_counts"]

    def __init__(self,
                 counts_per_cell: Union[int, Collection[int], None] = None,
                 total_counts: Optional[int] = None,
                 random_state: Union[None, int] = 0,
                 replace: bool = False):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.downsample_counts(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class Combat(PipelineNode):
    important_parameters = ["covariates"]
    requirements = InteractionList(obs=[istr("key"), istr("covariates")])

    def __init__(self,
                 key: str = "batch",
                 covariates: Optional[Collection[str]] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.combat(ann_data, **self.parameters, inplace=True)
        self._terminate(ann_data)
        return


class Neighbors(PipelineNode):
    important_parameters = ["n_neighbors", "n_pcs"]
    requirements = InteractionList(obsm=[istr("use_rep")])
    effects = (InteractionList(uns=["neighbors"],
                              obsp=["distances", "connectivities"]) *
               Interaction(uns=[istr("key_added")],
                           obsp=[istr("key_added")+"_distances", istr("key_added")+"_connectivities"]))

    def __init__(self,
                 n_neighbors: int = 15,
                 n_pcs: Optional[int] = None,
                 use_rep: Optional[str] = None,
                 knn: bool = True,
                 random_state: Union[None, int] = 0,
                 method: Optional[Literal["umap", "gauss", "rapids"]] = "umap",
                 metric: Literal['cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan', 'braycurtis', 'canberra',
                                 'chebyshev', 'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis',
                                 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener',
                                 'sokalsneath', 'sqeuclidean', 'yule'] = "euclidean",
                 key_added: Optional[str] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        pp.neighbors(ann_data, **self.parameters)
        self._terminate(ann_data)
        return
