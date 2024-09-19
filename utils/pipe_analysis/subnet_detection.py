# 导入相关的包
import numpy as np
import pandas as pd
import networkx as nx
from community import community_louvain



# 以类的方式，将上述的三个函数封装起来
class Graph:
    def __init__(self, edge_list_fp):
        if isinstance(edge_list_fp, str):
            self.edge_list_fp = edge_list_fp
            self.G = self.read_edge_list()
        else:
            self.edge_list_fp = edge_list_fp
            self.G = self.read_edge_list_from_df()
    # 读取无向图,从数据框
    def read_edge_list_from_df(self):
        data = self.edge_list_fp
        G = nx.Graph()
        for i in range(data.shape[0]):
            G.add_edge(data.iloc[i, 0], data.iloc[i, 1])
        return G

    # 读取无向图
    def read_edge_list(self, sep="\t"):
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

    # build subgraph from communities
    def build_subgraph(self, communities):
        subgraphs = {}
        for node, community_id in communities.items():
            if community_id not in subgraphs:
                subgraphs[community_id] = nx.Graph()
            subgraphs[community_id].add_node(node)
        for edge in self.G.edges():
            node1, node2 = edge
            community_id1 = communities[node1]
            community_id2 = communities[node2]
            if community_id1 == community_id2:
                subgraphs[community_id1].add_edge(node1, node2)
        return subgraphs

    # save graph to edge list type
    def save_graph(self, G, out_fp):
        nx.write_edgelist(G, out_fp)
        print(f"save graph to {out_fp}")

    def summary(self, G=None, out_fp=None):
        # 输出图的节点数和边数, 平均度等信息
        if G is None:
            G = self.G
        info = {}
        info["node_num"] = len(G.nodes)
        info["edge_num"] = len(G.edges)
        info["average_degree"] = np.mean(list(dict(G.degree()).values()))
        print(info)
        if out_fp is not None:
            with open(out_fp, "w") as f:
                f.write(f"node_num: {info['node_num']}\n")
                f.write(f"edge_num: {info['edge_num']}\n")
                f.write(f"average_degree: {info['average_degree']}\n")
        return info













