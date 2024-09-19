#读取基因互作网络,文件的结构如下，第一行是表头，后面的每一行是一个基因对的互作关系
'''

Gene1	Gene2	Annotation	Direction	Score
A1CF	APOBEC1	catalyzed by; complex; input	<-	1.00
A1CF	CELF2	predicted	-	0.90


'''
import pandas as pd
import os

# InteractionData类，读取基因互作网络文件, 互作信息存储在interactions属性中，该属性为一个dataframe
class InteractionData:
    def __init__(self, fp):
        self.fp = fp
        self.interactions = None

    def read_data(self, file_path=None):
        if file_path is None:
            file_path = self.fp
        self.interactions = pd.read_csv(file_path, sep='\t')
        return self.interactions


# 读取gmt文件，文件的结构如下所示, 文件没有表头，第一列表示通路名，第二列是通路的描述，后面的列是通路中的基因
'''
2-LTR circle formation	R-HSA-164843	BANF1	HMGA1	LIG4	PSIP1	XRCC4	XRCC5	XRCC6	gag	gag-pol	rev	vif	vpr	vpu

'''
#读取gmt文件，结果在self.gmt属性中，该属性套嵌的字典，键为通路名，如：{"pathway_name":{description: "pathway_description", genes: ["gene1", "gene2", ...]}}
class Gmt:
    def __init__(self, fp, sepr="\t"):
        self.fp = fp
        self.gmt = None
        self.sepr = sepr

    def read_data(self, file_path=None):
        if file_path is None:
            file_path = self.fp
        with open(file_path, 'r') as f:
            lines = f.readlines()
            self.gmt = {}
            for line in lines:
                line = line.strip().split(self.sepr)
                self.gmt[line[0]] = {"description": line[1], "genes": line[2:]}
        return self.gmt


# 创建一个gmt对象，同时创建一个InteractionData对象，遍历gmt中的每一条通路，遍历InteractionData对象的interactions属性，如果通路中的基因在interactions中的Gene1和Gene2同时出现，将其添加到一个新的dataframe中，结构和interactions相同，所有的dataframe拼接在一起，形成一个新的dataframe，并添加一列“pathway”，将其值设置为通路名。
class ReactomeNetwork:

    def __init__(self, gmt_fp, interaction_fp):
        self.gmt_fp = gmt_fp
        self.interaction_fp = interaction_fp
        self.gmt_data = Gmt(self.gmt_fp)
        self.interaction_data = InteractionData(self.interaction_fp)
        self.data = None
        self.data_by_pathway = None
        self.prepare()


    def prepare(self):
        self.gmt_data.read_data()
        self.interaction_data.read_data()

    def get_data(self):
        self.data = pd.DataFrame()
        for pathway, genes in self.gmt_data.gmt.items():
            genes = genes['genes']
            data = self.interaction_data.interactions[
                (self.interaction_data.interactions['Gene1'].isin(genes)) & (self.interaction_data.interactions['Gene2'].isin(genes))]
            data['pathway'] = pathway
            self.data = pd.concat([self.data, data])

        return self.data

    # 根据pathway的值，将dataframe分成多个子集，返回的结果为一个字典，键为pathway的值，值为一个dataframe
    def split_data(self):
        if self.data is None:
            self.get_data()
        data = self.data
        self.data_by_pathway = {}
        for pathway in data['pathway'].unique():
            self.data_by_pathway[pathway] = data[data['pathway'] == pathway]
        return self.data_by_pathway

    # 获取所有基因的名字, 基因的名字从interactions属性中的Gene1和Gene2中获取
    def get_all_genes_by_interaction(self):
        all_genes = set(self.interaction_data.interactions['Gene1'].tolist() + self.interaction_data.interactions['Gene2'].tolist())
        return all_genes
    # 获取所有基因的名字，基因的名字从gmt文件中获取
    def get_all_genes_by_gmt(self):
        all_genes = set()
        for genes in self.gmt_data.gmt.values():
            all_genes.update(genes['genes'])
        return all_genes

    # 从self.data_by_pathway中对每一个通路获取基因的名字，返回一个字典，键为通路名，值为该通路中的基因名
    def get_genes_by_pathway(self):
        genes_by_pathway = {}
        for pathway, data in self.data_by_pathway.items():
            genes_by_pathway[pathway] = set(data['Gene1'].tolist() + data['Gene2'].tolist())
        return genes_by_pathway

'''
输出文件的结果如下：
    Gene1 Gene2   Annotation  Direction   Score  Pathway
0   A1CF  APOBEC1 catalyzed by; complex; input    <-  1.00    2-LTR circle formation
1   A1CF  CELF2  predicted    -   0.90    MHC class II antigen presentation
2   A1CF  CELF2  predicted    -   0.90    Interferon alpha/beta signaling

这个文件共有1936838行，6列，分别是Gene1, Gene2, Annotation, Direction, Score, Pathway
'''


# 第二种方法，首先构建基因互作的邻接矩阵，将有互作的基因对的值设置为1，没有互作的基因对的值设置为0，
if __name__ == '__main__':
    base_dir = r"/home/liuyq/data/ar_data/interaction_db/reactome"
    gmt_fp = os.path.join(base_dir, "ReactomePathways.gmt")
    interaction_fp = os.path.join(base_dir, "FIsInGene_070323_with_annotations.txt")

    regnetwork = ReactomeNetwork(gmt_fp, interaction_fp)
    data = regnetwork.get_data()
    #print(data)
    data_by_pathway = regnetwork.split_data()
    #print(data_by_pathway)
    all_genes_by_interaction = regnetwork.get_all_genes_by_interaction()
    #print(all_genes_by_interaction)
    all_genes_by_gmt = regnetwork.get_all_genes_by_gmt()
    #print(all_genes_by_gmt)
    genes_by_pathway = regnetwork.get_genes_by_pathway()
    #print(genes_by_pathway)


































