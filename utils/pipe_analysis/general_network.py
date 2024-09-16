# 导入相关的包
import numpy as np
import pandas as pd
import networkx as nx
from community import community_louvain

# 以类的方式，将上述的三个函数封装起来
class Graph:
    def __init__(self, edge_list_fp):
        self.edge_list_fp = edge_list_fp
        self.G = self.read_edge_list()
    # 读取无向图
    def read_edge_list(self):
        data = pd.read_csv(self.edge_list_fp, sep="\t", header=None)
        G = nx.Graph()
        for i in range(data.shape[0]):
            G.add_edge(data.iloc[i, 0], data.iloc[i, 1], weight=data.iloc[i, 2])
        return G
    # 读取有向图
    def read_edge_list_directed(self):
        data = pd.read_csv(self.edge_list_fp, sep="\t", header=None)
        G = nx.DiGraph()
        for i in range(data.shape[0]):
            G.add_edge(data.iloc[i, 0], data.iloc[i, 1], type=data.iloc[i, 2])
        return G
    # 读取有向加权图
    def read_edge_list_directed_weighted(self):
        data = pd.read_csv(self.edge_list_fp, sep="\t", header=None)
        G = nx.DiGraph()
        for i in range(data.shape[0]):
            G.add_edge(data.iloc[i, 0], data.iloc[i, 1], weight=data.iloc[i, 2])
        return G
    # 读取无向加权图
    def read_edge_list_weighted(self):
        data = pd.read_csv(self.edge_list_fp, sep="\t", header=None)
        G = nx.Graph()
        for i in range(data.shape[0]):
            G.add_edge(data.iloc[i, 0], data.iloc[i, 1], weight=data.iloc[i, 2])
        return G

    # 划分社区, 输出社区划分结果，返回一个字典，键为节点，值为社区编号
    def community_detection(self):
        communities = community_louvain.best_partition(self.G)
        return communities




