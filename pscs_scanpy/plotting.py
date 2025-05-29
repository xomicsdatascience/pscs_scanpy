import os
from pathlib import Path
import scanpy as sc
from scanpy import plotting as pl
from typing import Optional, Sequence, Literal, Collection, Union, Tuple
from pscs_api import Interaction, PlottingNode
from pscs_api.base import InteractionList
from pscs_api.base import istr
from matplotlib import pyplot as plt
from typing import Collection


class Scatter(PlottingNode):
    important_parameters = ["x", "y", "save"]
    requirements = (InteractionList(Interaction(obs=[istr("x")]),
                                   Interaction(var=[istr("x")]),
                                   Interaction(var_names=[istr("x")])) *
                    InteractionList(Interaction(obs=[istr("y")]),
                                    Interaction(var=[istr("y")]),
                                    Interaction(var_names=[istr("y")])) *
                    InteractionList(Interaction(layers=[istr("layers")])) *
                    InteractionList(Interaction(obs=[istr("groups")])))
    doc_url = "https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.scatter.html"
    function = sc.pl.scatter
    __doc__ = PlottingNode.set_doc(function, doc_url)

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
        self.function = sc.pl.scatter
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        self.save_path = self.parameters["save"]
        return

    def run(self):
        if isinstance(self.parameters["color"], str) or self.parameters["color"] is None:
            super().run()
        elif isinstance(self.parameters["color"], Collection):
            # Scanpy will create multiple panels; instead, split to keep the figures editable
            all_colors = self.parameters["color"].copy()
            orig_save = self.parameters["save"]
            file_ext_split = os.path.splitext(self.parameters["save"])
            for idx, c in enumerate(all_colors):
                self.parameters["color"] = c
                self.parameters["save"] = file_ext_split[0] + "_" + str(idx) + file_ext_split[1]
                super().run()
            self.parameters["color"] = all_colors
            self.parameters["save"] = orig_save
        return


class HeatMap(PlottingNode):
    important_parameters = ["var_names", "groupby", "save"]
    requirements = InteractionList(var_names=[istr("var_names"), istr("var_group_positions")],
                                   obs=[istr("groupby")],
                                   var=[istr("gene_symbols")],
                                   layers=[istr("layer")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.heatmap.html"
    function = sc.pl.heatmap
    __doc__ = PlottingNode.set_doc(function, doc_url)

    def __init__(self,
                 var_names: Union[str, Collection[str]],
                 groupby: Union[str, Collection[str]] = None,
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
        self.function = sc.pl.heatmap
        if groupby is None:
            groupby = []
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        self.save_path = self.parameters["save"]
        return
    

class DotPlot(PlottingNode):
    important_parameters = ["var_names", "groupby", "save"]
    requirements = InteractionList(var_names=[istr("var_names")],
                                    obs=[istr("groupby")],
                                    var=[istr("gene_symbols")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.dotplot.html"
    function = sc.pl.dotplot
    __doc__ = PlottingNode.set_doc(function, doc_url)

    def __init__(self,
                 var_names: Union[str, Collection[str]],
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
        self.function = sc.pl.dotplot
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return


class TracksPlot(PlottingNode):
    important_parameters = ["var_names", "groupby", "save"]
    requirements = InteractionList(var_names=[istr("var_names")],
                                    obs=[istr("groupby")],
                                    var=[istr("gene_symbols")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.tracksplot.html"
    function = sc.pl.tracksplot
    __doc__ = PlottingNode.set_doc(function, doc_url)

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
        self.function = sc.pl.tracksplot
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return


class Violin(PlottingNode):
    important_parameters = ["keys", "groupby", "save"]
    requirements = (InteractionList(Interaction(var_names=[istr("keys")]),
                                    Interaction(obs=[istr("keys")])) *
                    InteractionList(obs=[istr("groupby")],
                                    layers=[istr("layer")]))
    doc_url = "https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.violin.html"
    function = sc.pl.violin
    __doc__ = PlottingNode.set_doc(function, doc_url)

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
        self.function = sc.pl.violin
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)


class StackedViolin(PlottingNode):
    important_parameters = ["var_names", "groupby", "save"]
    requirements = InteractionList(var_names=[istr("var_names")],
                                   obs=[istr("groupby")],
                                   var=[istr("gene_symbols")],
                                   layers=[istr("layer")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.tracksplot.html"
    function = sc.pl.tracksplot
    __doc__ = PlottingNode.set_doc(function, doc_url)

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
                 save: Union[str, bool] = "stackedviolin.png",
                 vmin: Optional[float] = None,
                 vmax: Optional[float] = None,
                 vcenter: Optional[float] = None):
        self.function = sc.pl.stacked_violin
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)


class MatrixPlot(PlottingNode):
    important_parameters = ["var_names", "groupby", "save"]
    requirements = InteractionList(var_names=[istr("var_names")],
                                   obs=[istr("groupby")],
                                   var=[istr("gene_symbols")],
                                   layers=[istr("layer")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.matrixplot.html"
    function = sc.pl.matrixplot
    __doc__ = PlottingNode.set_doc(function, doc_url)

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
        self.function = sc.pl.matrixplot
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return


class ClusterMap(PlottingNode):
    important_parameters = ["obs_keys", "save"]
    requirements = InteractionList(obs=[istr("obs_keys")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.clustermap.html"
    function = sc.pl.clustermap
    __doc__ = PlottingNode.set_doc(function, doc_url)

    def __init__(self,
                 obs_keys: Optional[str] = None,
                 use_raw: Optional[bool] = None,
                 save: Union[str, bool, None] = "clustermap.png"):
        self.function = sc.pl.clustermap
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return


class Ranking(PlottingNode):
    important_parameters = []
    requirements = InteractionList(Interaction(obs=[istr("keys")]),
                                   Interaction(var=[istr("keys")]),
                                   Interaction(uns=[istr("keys")]),
                                   Interaction(varm=[istr("keys")]),
                                   Interaction(obsm=[istr("keys")]))
    doc_url = "https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.ranking.html"
    function = sc.pl.ranking
    __doc__ = PlottingNode.set_doc(function, doc_url)

    def __init__(self,
                 attr: Literal["var", "obs", "uns", "varm", "obsm"],
                 keys: Union[str, Collection[str]]):
        self.function = sc.pl.ranking
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return


class Dendrogram(PlottingNode):
    important_parameters = ["groupby", "save"]
    requirements = InteractionList(obs=[istr("groupby")],
                                   uns=[istr("dendrogram_key")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.dendrogram.html"
    function = sc.pl.dendrogram
    __doc__ = PlottingNode.set_doc(function, doc_url)

    def __init__(self,
                 groupby: str,
                 dendrogram_key: Optional[str] = None,
                 orientation: Literal["top", "bottom", "left", "right"] = "top",
                 remove_labels: bool = False,
                 save: Union[str, bool, None] = "dendrogram.png"):
        self.function = sc.pl.dendrogram
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return


class HighestExprGenes(PlottingNode):
    important_parameters = ["n_top", "save"]
    requirements = InteractionList(var=[istr("gene_symbols")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.highest_expr_genes.html"
    function = sc.pl.highest_expr_genes
    __doc__ = PlottingNode.set_doc(function, doc_url)

    def __init__(self,
                 n_top: int = 30,
                 save: Union[str, bool, None] = "highestexprgenes.png",
                 gene_symbols: Optional[str] = None,
                 log: bool = False):
        self.function = sc.pl.highest_expr_genes
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return


class FilterGenesDispersion(PlottingNode):
    important_parameters = ["save"]
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.filter_genes_dispersion.html"
    function = sc.pl.filter_genes_dispersion
    __doc__ = PlottingNode.set_doc(function, doc_url)

    def __init__(self,
                 log: bool = False,
                 save: Union[str, bool, None] = "filtergenesdisp.png"):
        self.function = sc.pl.filter_genes_dispersion
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return


class HighlyVariableGenes(PlottingNode):
    important_parameters = ["save"]
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.highly_variable_genes.html"
    function = sc.pl.highly_variable_genes
    __doc__ = PlottingNode.set_doc(function, doc_url)
    def __init__(self,
                 log: bool = False,
                 save: Union[str, bool, None] = "highlyvariablegenes.png"):
        self.function = sc.pl.highest_expr_genes
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return


class PCA(PlottingNode):
    important_parameters = ["color", "save"]
    requirements = (InteractionList(Interaction(obs=[istr("color")]),
                                   Interaction(var_names=[istr("color")])) *
                    InteractionList(var=[istr("gene_symbols")],
                                    layers=[istr("layer")],
                                    obs=[istr("groups")]))
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.pca.html"
    function = sc.pl.pca
    __doc__ = PlottingNode.set_doc(function, doc_url)

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
        self.function = sc.pl.pca
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        if isinstance(self.parameters["color"], str) or self.parameters["color"] is None:
            super().run()
        elif isinstance(self.parameters["color"], Collection):
            # Scanpy will create multiple panels; instead, split to keep the figures editable
            all_colors = self.parameters["color"].copy()
            orig_save = self.parameters["save"]
            file_ext_split = os.path.splitext(self.parameters["save"])
            for idx, c in enumerate(all_colors):
                self.parameters["color"] = c
                self.parameters["save"] = file_ext_split[0] + "_" + str(idx) + file_ext_split[1]
                super().run()
            self.parameters["color"] = all_colors
            self.parameters["save"] = orig_save
        return


class PCALoadings(PlottingNode):
    important_parameters = ["components", "save"]
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.pca_loadings.html"
    function = sc.pl.pca_loadings
    __doc__ = PlottingNode.set_doc(function, doc_url)

    def __init__(self,
                 components: Union[str, Sequence[int], None] = None,
                 include_lowest: bool = True,
                 n_points: Optional[bool] = None,
                 save: Union[str, bool, None] = "pcaloadings.png"):
        self.function = sc.pl.pca_loadings
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return


class PCAVarianceRatio(PlottingNode):
    important_parameters = ["n_pcs", "save"]
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.pca_variance_ratio.html"
    function = sc.pl.pca_variance_ratio
    __doc__ = PlottingNode.set_doc(function, doc_url)

    def __init__(self,
                 n_pcs: int = 30,
                 log: bool = False,
                 save: Union[str, bool, None] = "pcavarianceratio.png"):
        self.function = sc.pl.pca_variance_ratio
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        ann_data = self.input_data[0]
        sc._settings.ScanpyConfig.figdir = Path(os.path.dirname(self.parameters["save"]))
        self.parameters["save"] = os.path.basename(self.parameters["save"])
        pl.pca_variance_ratio(ann_data, **self.parameters)
        self._terminate(ann_data)
        return


class PCAOverview(PlottingNode):
    important_parameters = ["save"]
    requirements = InteractionList(obs=[istr("color"), istr("groups")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.pca_overview.html"
    function = sc.pl.pca_overview
    __doc__ = PlottingNode.set_doc(function, doc_url)

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
        self.function = sc.pl.pca_overview
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        if isinstance(self.parameters["color"], str) or self.parameters["color"] is None:
            super().run()
        elif isinstance(self.parameters["color"], Collection):
            # Scanpy will create multiple panels; instead, split to keep the figures editable
            all_colors = self.parameters["color"].copy()
            orig_save = self.parameters["save"]
            file_ext_split = os.path.splitext(self.parameters["save"])
            for idx, c in enumerate(all_colors):
                self.parameters["color"] = c
                self.parameters["save"] = file_ext_split[0] + "_" + str(idx) + file_ext_split[1]
                super().run()
            self.parameters["color"] = all_colors
            self.parameters["save"] = orig_save
        return


class TSNE(PlottingNode):
    important_parameters = ["save"]
    requirements = InteractionList(obs=[istr("color")],
                                   var=[istr("gene_symbols")],
                                   layers=[istr("layer")],
                                   uns=[istr("neighbors_key")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.tsne.html"
    function = sc.pl.tsne
    __doc__ = PlottingNode.set_doc(function, doc_url)

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
        self.function = sc.pl.tsne
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        if isinstance(self.parameters["color"], str) or self.parameters["color"] is None:
            super().run()
        elif isinstance(self.parameters["color"], Collection):
            # Scanpy will create multiple panels; instead, split to keep the figures editable
            all_colors = self.parameters["color"].copy()
            orig_save = self.parameters["save"]
            file_ext_split = os.path.splitext(self.parameters["save"])
            for idx, c in enumerate(all_colors):
                self.parameters["color"] = c
                self.parameters["save"] = file_ext_split[0] + "_" + str(idx) + file_ext_split[1]
                super().run()
            self.parameters["color"] = all_colors
            self.parameters["save"] = orig_save
        return


class UMAP(PlottingNode):
    important_parameters = ["color", "save"]
    requirements = (InteractionList(Interaction(obs=[istr("color")]),
                                   Interaction(var_names=[istr("color")]),
                                    Interaction(var=[istr("color")])) *
                    InteractionList(obs=[istr("groups")],
                                    obsm=["X_umap"],
                                    var=[istr("gene_symbols")],
                                    layers=[istr("layer")],
                                    uns=[istr("neighbors_key"), "umap"]) *
                    InteractionList(Interaction(obsp=["connectivities"])))
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.umap.html"
    function = sc.pl.umap
    __doc__ = PlottingNode.set_doc(function, doc_url)

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
        self.function = sc.pl.umap
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        # self.save_path = self.parameters["save"]
        return

    def run(self):
        if isinstance(self.parameters["color"], str) or self.parameters["color"] is None:
            super().run()
        elif isinstance(self.parameters["color"], Collection):
            # Scanpy will create multiple panels; instead, split to keep the figures editable
            all_colors = self.parameters["color"].copy()
            orig_save = self.parameters["save"]
            file_ext_split = os.path.splitext(self.parameters["save"])
            for idx, c in enumerate(all_colors):
                self.parameters["color"] = c
                self.parameters["save"] = file_ext_split[0] + "_" + str(idx) + file_ext_split[1]
                super().run()
            self.parameters["color"] = all_colors
            self.parameters["save"] = orig_save
        return


class DiffMap(PlottingNode):
    important_parameters = ["save"]
    requirements = InteractionList(obs=[istr("color"), istr("groups")],
                                   var=[istr("gene_symbols")],
                                   layers=[istr("layer")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.diffmap.html"
    function = sc.pl.diffmap
    __doc__ = PlottingNode.set_doc(function, doc_url)

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
        self.function = sc.pl.diffmap
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        if isinstance(self.parameters["color"], str) or self.parameters["color"] is None:
            super().run()
        elif isinstance(self.parameters["color"], Collection):
            # Scanpy will create multiple panels; instead, split to keep the figures editable
            all_colors = self.parameters["color"].copy()
            orig_save = self.parameters["save"]
            file_ext_split = os.path.splitext(self.parameters["save"])
            for idx, c in enumerate(all_colors):
                self.parameters["color"] = c
                self.parameters["save"] = file_ext_split[0] + "_" + str(idx) + file_ext_split[1]
                super().run()
            self.parameters["color"] = all_colors
            self.parameters["save"] = orig_save
        return


class DrawGraph(PlottingNode):
    important_parameters = ["color", "save"]
    requirements = InteractionList(obs=[istr("color")],
                                   var=[istr("gene_symbols")],
                                   layers=[istr("layer")],
                                   uns=[istr("neighbors_key")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.draw_graph.html"
    function = sc.pl.draw_graph
    __doc__ = PlottingNode.set_doc(function, doc_url)

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
        self.function = sc.pl.draw_graph
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        if isinstance(self.parameters["color"], str) or self.parameters["color"] is None:
            super().run()
        elif isinstance(self.parameters["color"], Collection):
            # Scanpy will create multiple panels; instead, split to keep the figures editable
            all_colors = self.parameters["color"].copy()
            orig_save = self.parameters["save"]
            file_ext_split = os.path.splitext(self.parameters["save"])
            for idx, c in enumerate(all_colors):
                self.parameters["color"] = c
                self.parameters["save"] = file_ext_split[0] + "_" + str(idx) + file_ext_split[1]
                super().run()
            self.parameters["color"] = all_colors
            self.parameters["save"] = orig_save
        return


class Embedding(PlottingNode):
    important_parameters = ["basis", "color", "save"]
    requirements = InteractionList(obsm=[istr("basis")],
                                   obs=[istr("color"), istr("groups")],
                                   var=[istr("gene_symbols")],
                                   layers=[istr("layer")],
                                   uns=[istr("neighbors_key")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.embedding.html"
    function = sc.pl.embedding
    __doc__ = PlottingNode.set_doc(function, doc_url)

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
        self.function = sc.pl.embedding
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        if isinstance(self.parameters["color"], str) or self.parameters["color"] is None:
            super().run()
        elif isinstance(self.parameters["color"], Collection):
            # Scanpy will create multiple panels; instead, split to keep the figures editable
            all_colors = self.parameters["color"].copy()
            orig_save = self.parameters["save"]
            file_ext_split = os.path.splitext(self.parameters["save"])
            for idx, c in enumerate(all_colors):
                self.parameters["color"] = c
                self.parameters["save"] = file_ext_split[0] + "_" + str(idx) + file_ext_split[1]
                super().run()
            self.parameters["color"] = all_colors
            self.parameters["save"] = orig_save
        return


class EmbeddingDensity(PlottingNode):
    important_parameters = ["basis"]
    requirements = InteractionList(obsm=["X_" + istr("basis")],
                                   obs=[istr("key"), istr("groupby"), istr("group")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.embedding_density.html"
    function = sc.pl.embedding_density
    __doc__ = PlottingNode.set_doc(function, doc_url)

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
        self.function = sc.pl.embedding_density
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return


class DPTTimeseries(PlottingNode):
    important_parameters = ["save"]
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.dpt_timeseries.html"
    function = sc.pl.dpt_timeseries
    __doc__ = PlottingNode.set_doc(function, doc_url)

    def __init__(self,
                 color_map: str = None,
                 as_heatmap: bool = True,
                 save: Union[bool, str] = "dpt_timeseries.png"):
        self.function = sc.pl.dpt_timeseries
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return


class PAGA(PlottingNode):
    important_parameters = ["threshold", "save"]
    requirements = InteractionList(Interaction(obs=[istr("color")]),
                                   Interaction(var_names=[istr("color")]))
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.paga.html"
    function = sc.pl.paga
    __doc__ = PlottingNode.set_doc(function, doc_url)

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
                 edge_width_scale: float = 1.0,
                 min_edge_width: Optional[float] = None,
                 max_edge_width: Optional[float] = None,
                 arrowsize: int = 30,
                 normalize_to_color: bool = False,
                 cmap: Optional[str] = None,
                 add_pos: bool = True,
                 title: Optional[str] = None,
                 save: Union[bool, str] = "paga.png"):
        self.function = sc.pl.paga
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return

    def run(self):
        if isinstance(self.parameters["color"], str) or self.parameters["color"] is None:
            super().run()
        elif isinstance(self.parameters["color"], Collection):
            # Scanpy will create multiple panels; instead, split to keep the figures editable
            all_colors = self.parameters["color"].copy()
            orig_save = self.parameters["save"]
            file_ext_split = os.path.splitext(self.parameters["save"])
            for idx, c in enumerate(all_colors):
                self.parameters["color"] = c
                self.parameters["save"] = file_ext_split[0] + "_" + str(idx) + file_ext_split[1]
                super().run()
            self.parameters["color"] = all_colors
            self.parameters["save"] = orig_save
        return


class RankGenesGroups(PlottingNode):
    important_parameters = ["n_genes", "sharey"]
    requirements = InteractionList(obs=[istr("groups")],
                                   var=[istr("gene_symbols")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.rank_genes_groups.html"
    function = sc.pl.rank_genes_groups
    __doc__ = PlottingNode.set_doc(function, doc_url)

    def __init__(self,
                 groups: Union[str, Sequence[str], None] = None,
                 gene_symbols: Optional[str] = None,
                 n_genes: int = 20,
                 fontsize: int = 8,
                 ncols: int = 4,
                 sharey: bool = True,
                 save: Union[bool, str] = "rank_genes_groups.png"):
        self.function = sc.pl.rank_genes_groups
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return


class RankGenesGroupsViolin(PlottingNode):
    important_parameters = ["save", "groups", "n_genes", "split"]
    requirements = InteractionList(obs=[istr("groups")],
                                   var_names=[istr("gene_names")],
                                   var=[istr("gene_symbols")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.rank_genes_groups_violin.html"
    function = sc.pl.rank_genes_groups_violin
    __doc__ = PlottingNode.set_doc(function, doc_url)

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
        self.function = sc.pl.rank_genes_groups_violin
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return


class RankGenesGroupsStackedViolin(PlottingNode):
    important_parameters = ["save", "n_genes"]
    requirements = InteractionList(obs=[istr("groups"), istr("groupby")],
                                   var=[istr("gene_symbols")],
                                   uns=[istr("key")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.rank_genes_groups_stacked_violin.html"
    function = sc.pl.rank_genes_groups_stacked_violin
    __doc__ = PlottingNode.set_doc(function, doc_url)

    def __init__(self,
                 groups: Union[str, Sequence[str]] = None,
                 n_genes: Optional[int] = None,
                 gene_symbols: Optional[str] = None,
                 groupby: Optional[str] = None,
                 min_logfoldchange: Optional[float] = None,
                 key: Optional[str] = None,
                 save: Union[bool, str] = "rank_genes_groups_stacked_violin.png"):
        self.function = sc.pl.rank_genes_groups_stacked_violin
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return


class RankGenesGroupsHeatmap(PlottingNode):
    important_parameters = ["save", "n_genes"]
    requirements = InteractionList(obs=[istr("groups"), istr("groupby")],
                                   var=[istr("gene_symbols")],
                                   uns=[istr("key")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.rank_genes_groups_heatmap.html"
    function = sc.pl.rank_genes_groups_heatmap
    __doc__ = PlottingNode.set_doc(function, doc_url)

    def __init__(self,
                 groups: Union[str, Sequence[str], None] = None,
                 n_genes: Optional[int] = None,
                 gene_symbols: Optional[str] = None,
                 groupby: Optional[str] = None,
                 min_logfoldchange: Optional[float] = None,
                 key: Optional[str] = None,
                 save: Union[bool, str] = "rank_genes_groups_heatmap.png"):
        self.function = sc.pl.rank_genes_groups_heatmap
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return


class RankGenesGroupsDotPlot(PlottingNode):
    important_parameters = ["save", "n_genes"]
    requirements = InteractionList(obs=[istr("groups"), istr("groupby")],
                                   var=[istr("gene_symbols")],
                                   uns=[istr("key")],
                                   var_names=[istr("var_names")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.rank_genes_groups_dotplot.html"
    function = sc.pl.rank_genes_groups_dotplot
    __doc__ = PlottingNode.set_doc(function, doc_url)

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
        self.function = sc.pl.rank_genes_groups_dotplot
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return


class RankGenesGroupsMatrixplot(PlottingNode):
    important_parameters = ["save", "n_genes"]
    requirements = InteractionList(obs=[istr("groups"), istr("groupby")],
                                   var=[istr("gene_symbols")],
                                   uns=[istr("key")],
                                   var_names=[istr("var_names")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.rank_genes_groups_matrixplot.html"
    function = sc.pl.rank_genes_groups_matrixplot
    __doc__ = PlottingNode.set_doc(function, doc_url)

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
        self.function = sc.pl.rank_genes_groups_matrixplot
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return


class RankGenesGroupsTracksplot(PlottingNode):
    important_parameters = ["save", "n_genes"]
    requirements = InteractionList(obs=[istr("groups"), istr("groupby")],
                                   var=[istr("gene_symbols")],
                                   uns=[istr("key")])
    doc_url = "https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pl.rank_genes_groups_tracksplot.html"
    function = sc.pl.rank_genes_groups_tracksplot
    __doc__ = PlottingNode.set_doc(function, doc_url)

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
        self.function = sc.pl.rank_genes_groups_tracksplot
        super().__init__()
        self.store_vars_as_parameters(**vars(), show=False)
        return
