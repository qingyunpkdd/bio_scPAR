from ..pipe_analysis.deg import DiffereceGene
from ..data_process.merge_network import MergeNetwork




class GeneSelection:
    def __init__(self, data):
        self.data = data
        self.gene_dict = {}
        self.data_seed_genes = None
        self.data_bio_network_genes = None
        self.data_candidate_genes = None

    def set_seed_genes(self, genes=None):
        # 如果参数genes为None，则通过差异基因分析获取基因
        if genes is None:
            raise ValueError('Please provide seed genes')
        else:
            self.data_seed_genes = genes
        return self.data_seed_genes

    def set_bio_network_genes(self, genes=None):
        if genes is None:
            merge = MergeNetwork()
            merge.merge_network()
            self.data_bio_network_genes = merge.gene_dict
        else:
            self.data_bio_network_genes = genes
        return self.data_bio_network_genes

    #以self.data_seed_genes作为self.data_bio_network_genes的键，获取self.data_bio_network_genes中的值，并将所有的键和值存储在一个新的集合中。
    def get_candidate_genes(self):
        candidate_genes = set()
        for gene in self.data_seed_genes:
            if gene in self.data_bio_network_genes:
                candidate_genes.add(gene)
                candidate_genes.update(self.data_bio_network_genes[gene])
        self.data_candidate_genes = candidate_genes
        return self.data_candidate_genes






