import os
from pathlib import Path
import scanpy as sc
from scanpy import plotting as pl
from typing import Optional, Sequence, Literal, Collection, Union, Tuple
from pscs_api import OutputNode


class Scatter(OutputNode):
    important_parameters = ["x", "y", "save"]

    def __init__(self,
                 x: Optional[str] = None,
                 y: Optional[str] = None,
                 color: Union[str, Collection[str], None] = None,
                 use_raw: Optional[bool] = None,
                 layers: Union[str, Collection[str], None] = None,
                 basis: Optional[Literal["pca", "tsne", "umap", "diffmap", "draw_graph_fr"]] = None,
                 sort_order: bool = True,
                 groups: Union[str, Collection[str], None] = None,
                 components: Union[str, Collection[str], None] = None,
                 projection: Literal["2d", "3d"] = "2d",
                 size: Union[int, float, None] = None,
                 color_map: str = None,
                 title: Optional[str] = None,
                 save: Union[str, bool, None] = "scatter.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.scatter(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class HeatMap(OutputNode):
    important_parameters = ["var_names", "groupby", "save"]
    def __init__(self,
                 var_names: Union[str, Collection[str]],
                 groupby: Union[str, Collection[str]],
                 use_raw: Optional[bool] = None,
                 log: bool = False,
                 num_categories: int = 7,
                 figsize: Optional[Tuple[float, float]] = None,
                 dendrogram: Union[bool, str] = False,
                 gene_symbols: Optional[str] = None,
                 var_group_positions: Optional[Collection[Tuple[int,int]]] = None,
                 var_group_labels: Optional[Collection[str]] = None,
                 var_group_rotation: Optional[float] = None,
                 layer: Optional[str] = None,
                 standard_scale: Optional[Literal["var", "obs"]] = None,
                 swap_axes: bool = False,
                 show_gene_labels: Optional[bool] = None,
                 save: Union[str, bool, None] = "heatmap.png",
                 vmin: Optional[float] = None,
                 vmax: Optional[float] = None,
                 vcenter: Optional[float] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.heatmap(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class DotPlot(OutputNode):
    important_parameters = ["var_names", "groupy", "save"]
    def __init__(self,
                 var_names: Union[str, Collection],
                 groupby: Union[str, Sequence[str]],
                 use_raw: Optional[bool] = None,
                 log: bool = False,
                 num_categories: int = 7,
                 figsize: Optional[Tuple[float, float]] = None,
                 dendrogram: Union[bool, str] = False,
                 gene_symbols: Optional[str] = None,
                 var_group_positions: Optional[Collection[Tuple[int, int]]] = None,
                 var_group_labels: Optional[Collection[str]] = None,
                 var_group_rotation: Optional[float] = None,
                 layer: Optional[str] = None,
                 title: Optional[str] = None,
                 colorbar_title: Optional[str] = "Mean expression in group",
                 cmap: str = "Reds",
                 standard_scale: Optional[Literal["var", "group"]] = None,
                 swap_axes: Optional[bool] = False,
                 size_title: Optional[str] = "Fraction of cells in group (%)",
                 expression_cutoff: float = 0.0,
                 mean_only_expressed: bool = False,
                 dot_max: Optional[float] = None,
                 dot_min: Optional[float] = None,
                 smallest_dot: Optional[float] = 0.0,
                 save: Union[str, bool, None] = "dotplot.png",
                 vmin: Optional[float] = None,
                 vmax: Optional[float] = None,
                 vcenter: Optional[float] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.dotplot(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class TracksPlot(OutputNode):
    important_parameters = ["var_names", "groupy", "save"]
    def __init__(self,
                 var_names: Union[str, Collection[str]],
                 groupby: Union[str, Collection[str]],
                 use_raw: Optional[bool] = None,
                 log: bool = False,
                 num_categories: int = 7,
                 figsize: Optional[Tuple[float,float]] = None,
                 dendrogram: Union[bool, str] = False,
                 gene_symbols: Optional[str] = None,
                 var_group_positions: Optional[Collection[Tuple[int, int]]] = None,
                 var_group_labels: Optional[Collection[str]] = None,
                 var_group_rotation: Optional[float] = None,
                 layer: Optional[str] = None,
                 save: Union[str, bool, None] = "tracksplot.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.dotplot(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class Violin(OutputNode):
    important_parameters = ["keys", "groupy", "save"]
    def __init__(self,
                 keys: Union[str, Collection[str]],
                 groupby: Optional[str] = None,
                 log: bool = False,
                 use_raw: Optional[bool] = None,
                 stripplot: bool = True,
                 jitter: Union[float, bool] = True,
                 size: int = 1,
                 layer: Optional[str] = None,
                 scale: Literal["area", "count", "width"] = "width",
                 order: Optional[Collection[str]] = None,
                 multi_panel: Optional[bool] = None,
                 xlabel: str = "",
                 ylabel: Union[str, Collection[str], None] = None,
                 rotation: Optional[float] = None,
                 save: Union[str, bool, None] = "violin.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.violin(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class StackedViolin(OutputNode):
    important_parameters = ["var_names", "groupy", "save"]
    def __init__(self,
                 var_names: Union[str, Collection[str]],
                 groupby: Union[str, Sequence[str]],
                 use_raw: Optional[bool] = None,
                 log: bool = False,
                 num_categories: int = 7,
                 figsize: Optional[Tuple[float, float]] = None,
                 dendrogram: Union[bool, str] = False,
                 gene_symbols: Optional[str] = None,
                 var_group_positions: Optional[Collection[Tuple[int,int]]] = None,
                 var_group_labels: Optional[Collection[str]] = None,
                 var_group_rotation: Optional[float] = None,
                 layer: Optional[str] = None,
                 title: Optional[str] = None,
                 colorbar_title: Optional[str] = "Median expression in group",
                 cmap: Optional[str] = "Blues",
                 standard_scale: Optional[Literal["var", "obs"]] = None,
                 swap_axes: bool = False,
                 stripplot: bool = False,
                 jitter: Union[float, bool] = False,
                 size: int = 1,
                 order: Optional[Collection[str]] = None,
                 scale: Literal["area", "count", "width"] = "width",
                 yticklabels: Optional[bool] = False,
                 row_palette: Optional[str] = None,
                 save: Union[bool] = "stackedviolin.png",
                 vmin: Optional[float] = None,
                 vmax: Optional[float] = None,
                 vcenter: Optional[float] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.stacked_violin(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class MatrixPlot(OutputNode):
    important_parameters = ["var_names", "groupy", "save"]

    def __init__(self,
                 var_names: Union[str, Collection[str]],
                 groupby: Union[str, Collection[str]],
                 use_raw: Optional[bool] = None,
                 log: bool = False,
                 num_categories: int = 7,
                 figsize: Optional[Tuple[float, float]] = None,
                 dendrogram: Union[bool, str] = False,
                 gene_symbols: Optional[str] = None,
                 var_group_positions: Optional[Sequence[Tuple[int,int]]] = None,
                 var_group_labels: Optional[Sequence[str]] = None,
                 var_group_rotation: Optional[float] = None,
                 layer: Optional[str] = None,
                 title: Optional[str] = None,
                 colorbar_title: Optional[str] = "Mean expression in group",
                 cmap: Optional[str] = "viridis",
                 standard_scale: Optional[Literal["var", "group"]] = None,
                 swap_axes: bool = False,
                 save: Union[str, bool, None] = "matrixplot.png",
                 vmin: Optional[float] = None,
                 vmax: Optional[float] = None,
                 vcenter: Optional[float] = None):
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.matrixplot(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class ClusterMap(OutputNode):
    important_parameters = ["obs_keys", "save"]

    def __init__(self,
                 obs_keys: Optional[str] = None,
                 use_raw: Optional[bool] = None,
                 save: Union[str, bool, None] = "clustermap.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.clustermap(ann_data)
        self._terminate(ann_data)
        return


class Ranking(OutputNode):
    important_parameters = None

    def __init__(self,
                 attr: Literal["var", "obs", "uns", "varm", "obsm"],
                 keys: Union[str, Collection[str]]):
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.ranking(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class Dendrogram(OutputNode):
    important_parameters = ["groupby", "save"]

    def __init__(self,
                 groupby: str,
                 dendrogram_key: Optional[str] = None,
                 orientation: Literal["top", "bottom", "left", "right"] = "top",
                 remove_labels: bool = False,
                 save: Union[str, bool, None] = "dendrogram.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.dendrogram(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class HighestExprGenes(OutputNode):
    important_parameters = ["n_top", "save"]

    def __init__(self,
                 n_top: int = 30,
                 save: Union[str, bool, None] = "highestexprgenes.png",
                 gene_symbols: Optional[str] = None,
                 log: bool = False):
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.highest_expr_genes(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class FilterGenesDispersion(OutputNode):
    important_parameters = ["save"]

    def __init__(self,
                 log: bool = False,
                 save: Union[str, bool, None] = "filtergenesdisp.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.filter_genes_dispersion(ann_data)
        self._terminate(ann_data)
        return


class HighlyVariableGenes(OutputNode):
    important_parameters = ["save"]

    def __init__(self,
                 log: bool = False,
                 save: Union[str, bool, None] = "highlyvariablegenes.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.highly_variable_genes(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class PCAPlot(OutputNode):
    important_parameters = ["color", "save"]

    def __init__(self,
                 color: Union[str, Collection[str], None] = None,
                 gene_symbols: Optional[str] = None,
                 use_raw: Optional[bool] = None,
                 layer: Optional[str] = None,
                 annotate_var_explained: bool = False,
                 sort_order: bool = True,
                 groups: Optional[str] = None,
                 dimensions: Union[Tuple[int, int], Collection[Tuple[int, int]], None] = None,
                 components: Union[str, Collection[str]] = None,
                 projection: Literal["2d", "3d"] = "2d",
                 legend_loc: str = "right margin",
                 legend_fontsize: Union[int, float, Literal["xx-small", "x-small", "small", "medium", "large", "x-large", "xx-large"], None] = None,
                 legend_fontweight: Union[int, Literal["light", "normal", "medium", "semibold", "bold", "heavy", "black"]] = "bold",
                 legend_fontoutline: Optional[int] = None,
                 colorbar_loc: Optional[str] = "right",
                 size: Union[float, Collection[float], None] = None,
                 color_map: Union[str, None] = None,
                 palette: Union[str, Collection[str], None] = None,
                 na_color: str = "lightgray",
                 na_in_legend: bool = True,
                 frameon: Optional[bool] = None,
                 title: Union[str, Collection[str], None] = None,
                 vmin: Union[str, float] = None,
                 vmax: Union[str, float] = None,
                 vcenter: Union[str, float] = None,
                 add_outline: Optional[bool] = False,
                 outline_color: Tuple[str, str] = ("black", "white"),
                 outline_width: Tuple[float, float] = (0.3, 0.05),
                 ncols: int = 4,
                 wspace: Optional[float] = None,
                 hspace: float = 0.25,
                 save: Union[str, bool, None] = "pcaplot.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.pca(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class PCALoadings(OutputNode):
    important_parameters = ["components", "save"]

    def __init__(self,
                 components: Union[str, Sequence[int], None] = None,
                 include_lowest: bool = True,
                 n_points: Optional[bool] = None,
                 save: Union[str, bool, None] = "pcaloadings.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.pca_loadings(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class PCAVarianceRatio(OutputNode):
    important_parameters = ["n_pcs", "save"]

    def __init__(self,
                 n_pcs: int = 30,
                 log: bool = False,
                 save: Union[str, bool, None] = "pcavarianceratio.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.pca_variance_ratio(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class PCAOverview(OutputNode):
    important_parameters = ["save"]

    def __init__(self,
                 color: Union[str, Collection[str], None] = None,
                 use_raw: Optional[bool] = None,
                 layers: Union[str, Collection[str], None] = None,
                 basis: Optional[Literal["pca", "tsne", "umap", "diffmap", "draw_graph_fr"]] = None,
                 sort_order: bool = True,
                 groups: Union[str, Collection[str], None] = None,
                 components: Union[str, Collection[str], None] = None,
                 projection: Literal["2d", "3d"] = "2d",
                 size: Union[int, float, None] = None,
                 color_map: str = None,
                 title: Optional[str] = None,
                 legend_loc: str = "right margin",
                 legend_fontsize: Union[int, float, Literal[
                     "xx-small", "x-small", "small", "medium", "large", "x-large", "xx-large"], None] = None,
                 legend_fontweight: Union[
                     int, Literal["light", "normal", "medium", "semibold", "bold", "heavy", "black"]] = "bold",
                 legend_fontoutline: Optional[int] = None,
                 colorbar_loc: Optional[str] = "right",
                 palette: Union[str, Collection[str], None] = None,
                 na_color: str = "lightgray",
                 na_in_legend: bool = True,
                 frameon: Optional[bool] = None,
                 vmin: Union[str, float] = None,
                 vmax: Union[str, float] = None,
                 vcenter: Union[str, float] = None,
                 add_outline: Optional[bool] = False,
                 outline_color: Tuple[str, str] = ("black", "white"),
                 outline_width: Tuple[float, float] = (0.3, 0.05),
                 ncols: int = 4,
                 wspace: Optional[float] = None,
                 hspace: float = 0.25,
                 save: Union[str, bool, None] = "pcaoverview.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.pca_overview(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class TSNE(OutputNode):
    important_parameters = ["save"]

    def __init__(self,
                 color: Union[str, Sequence[str]] = None,
                 gene_symbols: Optional[str] = None,
                 use_raw: Optional[bool] = None,
                 layer: Optional[str] = None,
                 edges: bool = False,
                 edges_width: float = 0.1,
                 edges_color: Union[str, Sequence[float], Sequence[str]] = "grey",
                 neighbors_key: Optional[str] = None,
                 arrows: bool = False,
                 sort_order: bool = True,
                 groups: Optional[str] = None,
                 dimensions: Union[Tuple[int, int], Sequence[Tuple[int, int]], None] = None,
                 components: Union[str, Sequence[str]] = None,
                 projection: Literal["2d", "3d"] = "2d",
                 legend_loc: str = "right margin",
                 legend_fontsize: Union[int, float, Literal["xx-small", "x-small", "small", "medium", "large", "x-large", "xx-large"], None] = None,
                 legend_fontweight: Union[int, Literal['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']]  = "bold",
                 legend_fontoutline: Optional[int] = None,
                 colorbar_loc: Optional[str] = "right",
                 size: Union[float, Sequence[float], None] = None,
                 color_map: Union[str, None] = None,
                 palette: Union[str, Sequence[str], None] = None,
                 na_color: str = "lightgray",
                 na_in_legend: bool = True,
                 frameon: Optional[bool] = None,
                 title: Union[str, Sequence[str], None] = None,
                 vmin: Union[str, float, None] = None,
                 vmax: Union[str, float, None] = None,
                 vcenter: Union[str, float, None] = None,
                 add_outline: Optional[bool] = False,
                 outline_color: Tuple[str, str] = ("black", "white"),
                 outline_width: Tuple[float, float] = (0.3, 0.05),
                 ncols: int = 4,
                 wspace: Optional[float] = None,
                 hspace: float = 0.25,
                 save: Union[bool, str, None] = "tsne.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.tsne(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class UMAP(OutputNode):
    important_parameters = ["save"]

    def __init__(self,
                 color: Union[str, Sequence[str], None] = None,
                 gene_symbols: Optional[str] = None,
                 use_raw: Optional[bool] = None,
                 layer: Optional[str] = None,
                 edges: bool = False,
                 edges_width: float = 0.1,
                 edges_color: Union[str, Sequence[float], Sequence[str]] = "grey",
                 neighbors_key: Optional[str] = None,
                 sort_order: bool = True,
                 groups: Optional[str] = None,
                 dimensions: Union[Tuple[int, int], Sequence[Tuple[int, int]], None] = None,
                 projection: Literal["2d", "3d"] = "2d",
                 legend_loc: str = "right margin",
                 legend_fontsize: Union[int, float, Literal["xx-small", "x-small", "small", "medium", "large", "x-large", "xx-large"]] = None,
                 legend_fontweight: Union[int, Literal["light", "normal", "medium", "semibold", "bold", "heavy", "black"]] = "black",
                 legend_fontoutline: Optional[int] = None,
                 colorbar_loc: Optional[str] = "right",
                 size: Union[float, Sequence[float], None] = None,
                 color_map: str = None,
                 palette: Union[str, Sequence[str], None] = None,
                 na_color: str = "lightgray",
                 na_in_legend: bool = True,
                 frameon: Optional[bool] = None,
                 title: Union[str, Sequence[str], None] = None,
                 vmin: Union[str, float] = None,
                 vmax: Union[str, float] = None,
                 vcenter: Union[str,float] = None,
                 add_outline: Optional[bool] = False,
                 outline_color: Tuple[str, str] = ("black", "white"),
                 outline_width: Tuple[float, float] = (0.3, 0.05),
                 ncols: int = 4,
                 wspace: Optional[float] = None,
                 hspace: float = 0.25,
                 save: Union[bool, str, None] = "umap.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.umap(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class DiffMap(OutputNode):
    important_parameters = ["save"]

    def __init__(self,
                 color: Union[str, Sequence[str], None] = None,
                 gene_symbols: Optional[str] = None,
                 use_raw: Optional[bool] = None,
                 layer: Optional[str] = None,
                 edges: bool = False,
                 edges_width: float = 0.1,
                 edges_color: Union[str, Sequence[float], Sequence[str]] = "grey",
                 neighbors_key: Optional[str] = None,
                 sort_order: bool = True,
                 groups: Optional[str] = None,
                 dimensions: Union[Tuple[int, int], Sequence[Tuple[int, int]], None] = None,
                 projection: Literal["2d", "3d"] = "2d",
                 legend_loc: str = "right margin",
                 legend_fontsize: Union[int, float, Literal[
                     "xx-small", "x-small", "small", "medium", "large", "x-large", "xx-large"]] = None,
                 legend_fontweight: Union[
                     int, Literal["light", "normal", "medium", "semibold", "bold", "heavy", "black"]] = "black",
                 legend_fontoutline: Optional[int] = None,
                 colorbar_log: Optional[str] = "right",
                 size: Union[float, Sequence[float], None] = None,
                 color_map: str = None,
                 palette: Union[str, Sequence[str], None] = None,
                 na_color: str = "lightgray",
                 na_in_legend: bool = True,
                 frameon: Optional[bool] = None,
                 title: Union[str, Sequence[str], None] = None,
                 vmin: Union[str, float] = None,
                 vmax: Union[str, float] = None,
                 vcenter: Union[str, float] = None,
                 add_outline: Optional[bool] = False,
                 outline_color: Tuple[str, str] = ("black", "white"),
                 outline_width: Tuple[float, float] = (0.3, 0.05),
                 ncols: int = 4,
                 wspace: Optional[float] = None,
                 hspace: float = 0.25,
                 save: Union[bool, str, None] = "diffmap.png"
                 ):
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.diffmap(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class DrawGraph(OutputNode):
    important_parameters = ["color", "save"]

    def __init__(self,
                 color: Union[str, Sequence[str], None] = None,
                 gene_symbols: Optional[str] = None,
                 use_raw: Optional[bool] = None,
                 layer: Optional[str] = None,
                 layout: Optional[Literal["fa", "fr", "rt", "rt_circular", "drl", "eq_tree", Ellipsis]] = None,
                 edges: bool = False,
                 edges_width: float = 0.1,
                 edges_color: Union[str, Sequence[float], Sequence[str]] = "grey",
                 neighbors_key: Optional[str] = None,
                 arrows: bool = False,
                 sort_order: bool = True,
                 groups: Union[str, Collection[str], None] = None,
                 dimensions: Union[Tuple[int, int], Collection[Tuple[int, int]], None] = None,
                 components: Union[str, Collection[str], None] = None,
                 projection: Literal["2d", "3d"] = "2d",
                 legend_loc: str = "right margin",
                 legend_fontsize: Union[int, float, Literal[
                     "xx-small", "x-small", "small", "medium", "large", "x-large", "xx-large"], None] = None,
                 legend_fontweight: Union[
                     int, Literal["light", "normal", "medium", "semibold", "bold", "heavy", "black"]] = "bold",
                 legend_fontoutline: Optional[int] = None,
                 colorbar_loc: Optional[str] = "right",
                 size: Union[int, float, None] = None,
                 color_map: str = None,
                 palette: Union[str, Collection[str], None] = None,
                 na_color: str = "lightgray",
                 na_in_legend: bool = True,
                 frameon: Optional[bool] = None,
                 title: Union[str, Collection[str], None] = None,
                 vmin: Union[str, float] = None,
                 vmax: Union[str, float] = None,
                 vcenter: Union[str, float] = None,
                 add_outline: Optional[bool] = False,
                 outline_color: Tuple[str, str] = ("black", "white"),
                 outline_width: Tuple[float, float] = (0.3, 0.05),
                 ncols: int = 4,
                 wspace: Optional[float] = None,
                 hspace: float = 0.25,
                 save: Union[str, bool, None] = "drawgraph.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.draw_graph(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class Embedding(OutputNode):
    important_parameters = ["basis", "color", "save"]

    def __init__(self,
                 basis: str,
                 color: Union[str, Sequence[str], None] = None,
                 gene_symbols: Optional[str] = None,
                 use_raw: Optional[bool] = None,
                 edges: bool = False,
                 edges_width: float = 0.1,
                 edges_color: Union[str, Sequence[float], Sequence[str]] = "grey",
                 neighbors_key: Optional[str] = None,
                 arrows: bool = False,
                 sort_order: bool = True,
                 groups: Union[str, Collection[str], None] = None,
                 dimensions: Union[Tuple[int, int], Collection[Tuple[int, int]], None] = None,
                 components: Union[str, Collection[str], None] = None,
                 projection: Literal["2d", "3d"] = "2d",
                 legend_loc: str = "right margin",
                 legend_fontsize: Union[int, float, Literal[
                     "xx-small", "x-small", "small", "medium", "large", "x-large", "xx-large"], None] = None,
                 legend_fontweight: Union[
                     int, Literal["light", "normal", "medium", "semibold", "bold", "heavy", "black"]] = "bold",
                 legend_fontoutline: Optional[int] = None,
                 colorbar_loc: Optional[str] = "right",
                 size: Union[int, float, None] = None,
                 color_map: str = None,
                 palette: Union[str, Collection[str], None] = None,
                 na_color: str = "lightgray",
                 na_in_legend: bool = True,
                 frameon: Optional[bool] = None,
                 title: Union[str, Collection[str], None] = None,
                 vmin: Union[str, float] = None,
                 vmax: Union[str, float] = None,
                 vcenter: Union[str, float] = None,
                 add_outline: Optional[bool] = False,
                 outline_color: Tuple[str, str] = ("black", "white"),
                 outline_width: Tuple[float, float] = (0.3, 0.05),
                 ncols: int = 4,
                 wspace: Optional[float] = None,
                 hspace: float = 0.25,
                 save: Union[str, bool, None] = "embedding.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.embedding(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class EmbeddingDensity(OutputNode):
    important_parameters = ["basis"]

    def __init__(self,
                 basis: str = "umap",
                 key: Optional[str] = None,
                 groupby: Optional[str] = None,
                 group: Union[str, Collection[str], None] = "all",
                 colormap: str = "YlOrRd",
                 bg_dotsize: Optional[int] = 80,
                 fg_dotsize: Optional[int] = 180,
                 vmin: Union[str, float] = None,
                 vmax: Union[str, float] = None,
                 vcenter: Union[str, float] = None,
                 ncols: int = 4,
                 wspace: Optional[float] = None,
                 hspace: float = 0.25,
                 save: Union[bool, str] = "embedding_density.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.embedding_density(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class DPTTimeseries(OutputNode):
    important_parameters = ["save"]

    def __init__(self,
                 color_map: str = None,
                 as_heatmap: bool = True,
                 save: Union[bool, str] = "dpt_timeseries.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.dpt_timeseries(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class PAGA(OutputNode):
    important_parameters = ["threshold", "save"]

    def __init__(self,
                 threshold: Optional[float] = None,
                 color: Union[str] = None,
                 labels: Union[str, Sequence[str], None] = None,
                 pos: Optional[str] = None,
                 layout: Optional[Literal["fa", "fr", "rt", "rt_circular", "drl", "eq_tree", Ellipsis]] = None,
                 random_state: Optional[int] = 0,
                 root: Union[int, str, Sequence[int], None] = 0,
                 transitions: Optional[str] = None,
                 solid_edges: str = "connectivities",
                 dashed_edges: Optional[str] = None,
                 single_component: bool = False,
                 fontsize: Optional[int] = None,
                 fontoutline: Optional[int] = None,
                 node_size_power: float = 0.5,
                 edge_with_scale: float = 1.0,
                 min_edge_width: Optional[float] = None,
                 max_edge_width: Optional[float] = None,
                 arrowsize: int = 30,
                 normalize_to_color: bool = False,
                 cmap: Optional[str] = None,
                 add_pos: bool = True,
                 title: Optional[str] = None,
                 save: Union[bool, str] = "paga.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.paga(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class RankGenesGroups(OutputNode):
    important_parameters = ["n_genes", "sharey"]

    def __init__(self,
                 groups: Union[str, Sequence[str], None] = None,
                 gene_symbols: Optional[str] = None,
                 n_genes: int = 20,
                 fontsize: int = 8,
                 ncols: int = 4,
                 sharey: bool = True,
                 save: Union[bool, str] = "rank_genes_groups.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.rank_genes_groups(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class RankGenesGroupsViolin(OutputNode):
    important_parameters = ["save", "groups", "n_genes", "split"]

    def __init__(self,
                 groups: Optional[Sequence[str]] = None,
                 n_genes: int = 20,
                 gene_names: Optional[Collection[str]] = None,
                 gene_symbols: Optional[str] = None,
                 use_raw: Optional[bool] = None,
                 split: bool = True,
                 scale: str = "width",
                 strip: bool = True,
                 jitter: Union[int, float, bool] = True,
                 size: int = 1,
                 save: Union[bool, str] = "rank_genes_groups_violin.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.rank_genes_groups_violin(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class RankGenesGroupsStackedViolin(OutputNode):
    important_parameters = ["save", "n_genes"]

    def __init__(self,
                 groups: Union[str, Sequence[str]] = None,
                 n_genes: Optional[int] = None,
                 gene_symbols: Optional[str] = None,
                 groupby: Optional[str] = None,
                 min_logfoldchange: Optional[float] = None,
                 key: Optional[str] = None,
                 save: Union[bool, str] = "rank_genes_groups_stacked_violin.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.rank_genes_groups_stacked_violin(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class RankGenesGroupsHeatmap(OutputNode):
    important_parameters = ["save", "n_genes"]

    def __init__(self,
                 groups: Union[str, Sequence[str], None] = None,
                 n_genes: Optional[int] = None,
                 gene_symbols: Optional[str] = None,
                 groupby: Optional[str] = None,
                 min_logfoldchange: Optional[float] = None,
                 key: Optional[str] = None,
                 save: Union[bool, str] = "rank_genes_groups_heatmap.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.rank_genes_groups_heatmap(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class RankGenesGroupsDotPlot(OutputNode):
    important_parameters = ["save", "n_genes"]

    def __init__(self,
                 groups: Union[str, Sequence[str], None] = None,
                 n_genes: Optional[int] = None,
                 gene_symbols: Optional[str] = None,
                 groupby: Optional[str] = None,
                 min_logfoldchange: Optional[float] = None,
                 key: Optional[str] = None,
                 values_to_plot: Optional[Literal["scores", "logfoldchanges", "pvals", "pvals_adj", "log10_pvals", "log10_pvals_adj"]] = None,
                 var_names: Union[Sequence[str], None] = None,
                 save: Union[bool, str] = "rank_genes_groups_dotplot.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars)
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.rank_genes_groups_dotplot(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class RankGenesGroupsMatrixplot(OutputNode):
    important_parameters = ["save", "n_genes"]

    def __init__(self,
                 groups: Union[str, Sequence[str], None] = None,
                 n_genes: Optional[int] = None,
                 gene_symbols: Optional[str] = None,
                 groupby: Optional[str] = None,
                 min_logfoldchange: Optional[float] = None,
                 key: Optional[str] = None,
                 values_to_plot: Optional[Literal["scores", "logfoldchanges", "pvals", "pvals_adj", "log10_pvals", "log10_pvals_adj"]] = None,
                 var_names: Union[Sequence[str], None] = None,
                 save: Union[bool, str] = "rank_genes_groups_matrixplot.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.rank_genes_groups_matrixplot(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class RankGenesGroupsTracksplot(OutputNode):
    important_parameters = ["save", "n_genes"]

    def __init__(self,
                 groups: Union[str, Sequence[str], None] = None,
                 n_genes: Optional[int] = None,
                 gene_symbols: Optional[str] = None,
                 groupby: Optional[str] = None,
                 min_logfoldchange: Optional[float] = None,
                 key: Optional[str] = None,
                 values_to_plot: Optional[Literal[
                     "scores", "logfoldchanges", "pvals", "pvals_adj", "log10_pvals", "log10_pvals_adj"]] = None,
                 var_names: Union[Sequence[str], None] = None,
                 save: Union[bool, str] = "rank_genes_groups_matrixplot.png"):
        super().__init__()
        self.store_vars_as_parameters(**vars())
        return

    def run(self):
        ann_data = self._previous[0].result
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.rank_genes_groups_tracksplot(ann_data, **self.parameters)
        self._terminate(ann_data)
        return
