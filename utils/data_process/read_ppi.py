'''
读取PPI数据集
该数据集是一个无向图，每一行表示一个边，两个节点之间用空格隔开，

互作网络文件
V:\DATA\ontology_db\gene_networks\snap_ppi_decagon\PP-Decagon_ppi.csv


id->基因名转换文件
https://snap.stanford.edu/biodata/datasets/10022/10022-G-SynMiner.html
V:\DATA\ontology_db\gene_networks\snap_ppi_decagon\G-SynMiner_miner-geneHUGO.tsv
'''


import pandas as pd
import os

class PPI:
    def __init__(self, mapping_fp, ppi_fp):
        self.ppi_fp = ppi_fp
        self.mapping_fp = mapping_fp
        self.interactions = None
        self.mapper = None

    #读取转换文件，存为数据框
    def read_mapper(self, file_path=None):
        if file_path is None:
            file_path = self.mapping_fp
        self.mapper = pd.read_csv(file_path, sep='\t', header=0)
        return self.mapper
    #读取互作网络文件，存为数据框
    def read_ppi(self, file_path=None):
        if file_path is None:
            file_path = self.ppi_fp
        self.interactions = pd.read_csv(file_path, sep=',', header=None)
        return self.interactions
    #将id转换为基因名,在ppi文件中的两列中，将id转换为基因名，mapping文件中，基因名为“symbol”，id为“entrez_id”
    def map_id_to_gene(self):
        self.interactions.columns = ['id1', 'id2']
        self.interactions = pd.merge(self.interactions, self.mapper, left_on='id1', right_on='entrez_id', how='left')
        self.interactions = pd.merge(self.interactions, self.mapper, left_on='id2', right_on='entrez_id', how='left')
        self.interactions = self.interactions[['symbol_x', 'symbol_y']]
        self.interactions.columns = ['gene1', 'gene2']
        return self.interactions
'''
该文件没有方向性，所以不需要考虑方向性，只需要考虑两个节点之间是否有边。
样例数据如下：
    gene1 gene2
0  A1BG   A1BG
1  A1BG   A1BG
2  A1BG   A1BG

'''

if __name__ == '__main__':
    mapping_fp = r"/home/liuyq/data/ar_data/interaction_db/ppi/G-SynMiner_miner-geneHUGO.tsv"
    ppi_fp = r"/home/liuyq/data/ar_data/interaction_db/ppi/PP-Decagon_ppi.csv"
    ppi_obj = PPI(mapping_fp, ppi_fp)
    ppi_obj.read_mapper()
    ppi_obj.read_ppi()
    ppi_obj.map_id_to_gene()
    print(ppi_obj.interactions)
    print(ppi_obj.mapper)
























