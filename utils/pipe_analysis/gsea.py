'''
基于gseapy的GSEA分析 获得候选基因的富集通路
'''

import multiprocessing
import os, sys
import numpy as np
import pandas as pd
import networkx as nx
import gseapy as gp

from utils.data_process.read_reactome import Gmt
# gene set 文件由用户提供，需要用户输入的文件包括：gene_list， gene_sets，background， outdir

'''
一个富集分析的结果示例：
enr2.results.head()
Gene_set	Term	Overlap	P-value	Adjusted P-value	Odds Ratio	Combined Score	Genes
0	genes.gmt	BvA_UpIN_A	8/139	0.457390	0.568432	1.161982	0.908925	PCSK6;MAP3K5;MBOAT2;MSRB2;IQGAP2;HAL;PADI2;IL1R1
1	genes.gmt	BvA_UpIN_B	12/130	0.026744	0.187208	2.160059	7.822534	FAM65B;MBNL3;GPX8;DYSF;KCTD12;HEBP1;SUOX;ARHGD...
2	genes.gmt	CvA_UpIN_A	1/12	0.481190	0.568432	2.266479	1.657913	MBOAT2
3	genes.gmt	DvA_UpIN_A	16/284	0.426669	0.568432	1.127395	0.960255	PCSK6;FXYD6;IFNGR2;MAP3K5;MBOAT2;VNN1;IQGAP2;H...
4	genes.gmt	DvA_UpIN_D	13/236	0.487227	0.568432	1.084567	0.779830	GNB4;FAM198B;FAM65B;TXNDC5;GLIPR2;MBNL3;GPX8;D...

'''

class GSEA:
    def __init__(self, gene_sets_fp, outdir, background=None):
        self.gene_sets_fp = gene_sets_fp
        self.background = self.build_background(background)
        self.background_genes = None
        self.outdir = outdir
        self.gmt = self.build_gmt(gene_sets_fp)

    #静态方法，读取文件
    @staticmethod
    def read_file(file_path):
        with open(file_path, 'r') as f:
            lines = f.readlines()
            lines = [line.strip() for line in lines]
        return lines

    #构建背景基因
    def build_background(self,background):
        if background is None:
            _background = 30000
        elif isinstance(background, int):
            _background = background
        else:
            self.background_genes = self.read_file(background)
            _background = len(self.background_genes)
        return _background

    #构建gene_list
    def build_gene_list(self, gene_list):
        gene_list = self.read_file(gene_list)
        return gene_list

    def build_gmt(self, gene_sets_fp):
        gmt = Gmt(gene_sets_fp)
        gmt.read_data()
        #将套嵌字典展开为：键为pathway名，值为基因名的列表
        gmt = {key: value['genes'] for key, value in gmt.gmt.items()}
        return gmt

    def enrich(self, gene_list, save_fp=None):
        if isinstance(gene_list, str):
            gene_list = self.build_gene_list(gene_list)
        background = self.background
        gmt = self.gmt
        outdir = self.outdir if save_fp is None else save_fp
        enrichr = gp.enrichr(gene_list=gene_list,
                             gene_sets=gmt,
                             background=background,
                             outdir=outdir)
        return enrichr.results

    def filter_results(self, results, threshold_dict):
        results = results.copy()
        for metric in threshold_dict.keys():
            results = self.filter_by_metric(results, metric, threshold_dict[metric])
        return results

    #如果metric 是term，那么threshold应该是一个通路的名字的列表
    def filter_by_metric(self, results, metric, threshold):
        if metric == "P-value" or metric == "Adjusted P-value":
            results = results[results[metric] <= threshold]
        elif metric == "Term":
            results = results[results[metric] in threshold]
        return results

    def save_results(self, results, save_fp):
        if isinstance(results, pd.DataFrame) and save_fp is not None:
            results.to_csv(save_fp, sep='\t', index=False)
        else:
            raise ValueError("results should be a pandas DataFrame")

    #返回满足要求的通路的名字
    def get_terms(self, results):
        return results["Term"].tolist()



if __name__ == '__main__':
    base_dir = r"/home/liuyq/data/ar_data/interaction_db/reactome"
    gmt_fp = os.path.join(base_dir, "ReactomePathways.gmt")
    gene_list_fp = os.path.join(base_dir, "gene_list.txt")

    base_out_dir = r"/home/liuyq/data/ar_data/test"
    # out_dir = os.path.join(base_dir, "gsea_results.tsv")
    gsea_obj = GSEA(gmt_fp, background=30000, outdir=base_out_dir)
    results = gsea_obj.enrich(gene_list_fp, save_fp=base_out_dir)
    threshold_dict = {"Adjusted P-value": 0.01}
    results_filtered = gsea_obj.filter_results(results, threshold_dict)
    print(results.head())







