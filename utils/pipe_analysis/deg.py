import scanpy as sc
from utils.data_process.read_cell_type import read_cell_type
import pandas as pd
import os
from config import Config

class DiffereceGene:
    def __init__(self,
                 data_path,
                 cell_type_fn="cell_type.tsv"):
        self.data_path = data_path
        self.cell_type_fn = cell_type_fn
        self.deg_dict = {}
        self.deg_df_dict = {}

    def load_data(self):
        adata = sc.read_10x_mtx(self.data_path, var_names='gene_symbols', cache=False)
        return adata

    def assign_cell_type(self, adata):
        cell_type_fn = os.path.join(self.data_path, self.cell_type_fn)
        cell_types = read_cell_type(cell_type_fn)
        adata.obs['cell_type'] = cell_types
        return adata

    def get_diff_genes(self, adata):
        pvals = Config.deg_threshold_dict["p_value"]
        pvals_adj = Config.deg_threshold_dict["pval_adj"]
        log2fc = Config.deg_threshold_dict["log2fc"]
        sc.tl.rank_genes_groups(adata, "cell_type", method="t-test")
        # get all marker genes for all clusters
        result = adata.uns["rank_genes_groups"]
        groups = result["names"].dtype.names
        marker_genes = pd.DataFrame(
            {group + "__" + key: result[key][group]
             for group in groups
             for key in ["names", "pvals", "pvals_adj", "logfoldchanges"]
             }
        )
        unique_cells = set(adata.obs['cell_type'])
        for cell in unique_cells:
            deg = marker_genes[marker_genes[cell + "__pvals"] < pvals_adj]
            deg = deg[deg[cell + "__pvals"] < pvals]
            deg = deg[deg[cell + "__logfoldchanges"] > log2fc]
            self.deg_dict[cell] = deg[cell + "__names"].tolist()
            # filter all other cell type columns, and rename the columns: trim the "cell__" prefix
            deg = deg.filter(like=cell + "__")
            # deg.columns = [col.replace(cell + "__", "") for col in deg.columns]
            self.deg_df_dict[cell] = deg

    def preprocess(self, adata):
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.tl.rank_genes_groups(adata, "cell_type", method="t-test")

        # get the differentially expressed genes
        sc.tl.filter_rank_genes_groups(adata, groupby="cell_type", key_added="rank_genes_groups_filtered")
        return adata

    def save_deg_df_dict(self, column_names=None):
        deg_dir = Config.deg_fp
        for cell_type, deg_df in self.deg_df_dict.items():
            deg_fn = os.path.join(deg_dir, f"{cell_type}_deg.tsv")
            if column_names is not None:
                deg_df = deg_df[column_names]
            deg_df.to_csv(deg_fn, sep='\t', index=False)

        print("保存成功")


'''
self.deg_dict的结构如下：
cell_type genes
FCGR3A+_Mono CD14,FCGR3A,FCER1G,CD163,CD1C,CD1A,CD1E,CD68
Platelet CD177,PPBP,PF4,ITGA2B,ITGB3,ITGA2,ITGB3,ITGA2B,ITGB3
B CD79A,CD79B,IGHM,IGHD,IGHG1,IGHG2A,IGHG2B,IGHG2C,IGHG3
...


'''

if __name__ == "__main__":
    data_path = r'/mnt/sda/liuyq/ar_data/sc_data/pbmc'
    cell_type_fn = r'/mnt/sda/liuyq/ar_data/sc_data/pbmc/cell_type.tsv'

    diff_gene = DiffereceGene(data_path, cell_type_fn)
    adata = diff_gene.load_data()
    adata = diff_gene.assign_cell_type(adata)
    adata = diff_gene.preprocess(adata)
    diff_gene.get_diff_genes(adata)

    diff_gene.save_deg_df_dict()

    # # 创建空 DataFrame
    # results_df = pd.DataFrame(columns=['cell_type', 'genes'])
    #
    # # 打印每种细胞类型的差异表达基因
    # for cell_type, genes in diff_gene.deg_dict.items():
    #     print(f"Differentially expressed genes for {cell_type}: {genes}")
    #
    # # 将结果添加到 DataFrame
    # rows = []
    # for cell_type, genes in diff_gene.deg_dict.items():
    #     rows.append({'cell_type': cell_type, 'genes': ','.join(genes)})
    #
    # results_df = pd.concat([results_df, pd.DataFrame(rows)], ignore_index=True)
    #
    # results_df.to_csv('diff_genes_results_pbmc.csv', index=False)
    print("保存成功")
