import pandas as pd
import numpy as np

from utils.data_process.read_ppi import PPI
from utils.pipe_analysis.subnet_detection import Graph
from utils.data_process.read_reactome import Gmt, InteractionData, RegulationNetwork
from utils.data_process.read_regnetworks import TFData, RegulationNetwork, InteractionData

#创建一个类，将上述的三个网络合并起来，创建一个新的网络
'''
这个构建好的网的结构应该是一个字典，其中键是基因名，值是与之有互作的基因名的集合

'''


class MergeNetwork:
    def __init__(self, ppi_fp, gmt_fp, interaction_fp, tf_fp):
        self.ppi_fp = ppi_fp
        self.gmt_fp = gmt_fp
        self.interaction_fp = interaction_fp
        self.tf_fp = tf_fp
        self.gene_dict = {}

    # 读取并创建reactome网络
    def read_reactome(self):
        reactome = RegulationNetwork(self.gmt_fp, self.interaction_fp)
        reactome.prepare()
        return reactome.get_data()

    # 读取并创建regnetworks网络
    def read_regnetworks(self):
        regnetworks = RegulationNetwork(self.tf_fp, self.interaction_fp)
        regnetworks.create_network()
        return regnetworks.data

    # 读取并创建ppi网络
    def read_ppi(self):
        ppi = PPI(self.ppi_fp, self.ppi_fp)
        ppi.read_mapper()
        ppi.read_ppi()
        ppi.map_id_to_gene()
        return ppi.interactions

    # 合并三个网络
    def merge_network(self):
        ppi = self.read_ppi()
        reactome = self.read_reactome()
        regnetworks = self.read_regnetworks()
        gene_dict = {}
        for index, row in ppi.iterrows():
            if row['gene1'] not in gene_dict:
                gene_dict[row['gene1']] = set()
            if row['gene2'] not in gene_dict:
                gene_dict[row['gene2']] = set()
            gene_dict[row['gene1']].add(row['gene2'])
            gene_dict[row['gene2']].add(row['gene1'])
        for index, row in reactome.iterrows():
            if row['Gene1'] not in gene_dict:
                gene_dict[row['Gene1']] = set()
            if row['Gene2'] not in gene_dict:
                gene_dict[row['Gene2']] = set()
            gene_dict[row['Gene1']].add(row['Gene2'])
            gene_dict[row['Gene2']].add(row['Gene1'])
        for index, row in regnetworks.iterrows():
            if row['TF'] not in gene_dict:
                gene_dict[row['TF']] = set()
            if row['Target'] not in gene_dict:
                gene_dict[row['Target']] = set()
            gene_dict[row['TF']].add(row['Target'])
            gene_dict[row['Target']].add(row['TF'])
        self.gene_dict = gene_dict
        return gene_dict

    # gene_pair_list的数据结构为：[['gene1', 'gene2'], ['gene3', 'gene4'], ['gene5', 'gene6']]
    def validate_ar(self, gene_pair_list):
        ppi = self.read_ppi()
        reactome = self.read_reactome()
        regnetworks = self.read_regnetworks()
        result_valid = {"ppi": None, "reactome": None, "regnetworks": None}
        #设置双重索引
        ppi.set_index(['gene1', 'gene2'], inplace=True)
        reactome.set_index(['Gene1', 'Gene2'], inplace=True)
        regnetworks.set_index(['TF', 'Target'], inplace=True)
        #将gene_pair_list中的每一对作为索引，从ppi, reactome, regnetworks中提取对应的数据，如果提取到了，将其存入result_valid中
        for gene_pair in gene_pair_list:
            antecedent = gene_pair[0]
            consequent = gene_pair[1]
            if (antecedent, consequent) in ppi.index:
                if result_valid["ppi"] is None:
                    result_valid["ppi"] = ppi.loc[(antecedent, consequent)]
                else:
                    result_valid["ppi"] = pd.concat([result_valid["ppi"], ppi.loc[(antecedent, consequent)]])
            if (antecedent, consequent) in reactome.index:
                if result_valid["reactome"] is None:
                    result_valid["reactome"] = reactome.loc[(antecedent, consequent)]
                else:
                    result_valid["reactome"] = pd.concat([result_valid["reactome"], reactome.loc[(antecedent, consequent)]])
            if (antecedent, consequent) in regnetworks.index:
                if result_valid["regnetworks"] is None:
                    result_valid["regnetworks"] = regnetworks.loc[(antecedent, consequent)]
                else:
                    result_valid["regnetworks"] = pd.concat([result_valid["regnetworks"], regnetworks.loc[(antecedent, consequent)]])
        #用一个字典保存统计信息，包括在每个网络中找到的数量，百分比，以及在所有网络中找到的数量和百分比，以及各种未找到的数量和百分比
        result = {}
        for key, value in result_valid.items():
            if value is None:
                result[key] = 0
            else:
                result[key] = len(value)
        result["total"] = sum(result.values())
        result["not_found"] = len(gene_pair_list) - result["total"]
        result["not_found_percentage"] = result["not_found"] / len(gene_pair_list)
        result["found"] = len(gene_pair_list) - result["not_found"]
        result["found_percentage"] = 1 - result["not_found_percentage"]
        for key in result.keys():
            result[key + "_percentage"] = result[key] / result["total"]
        return result, result_valid



















