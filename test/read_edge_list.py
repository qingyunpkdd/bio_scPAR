'''
edge_list.tsv文件的格式如下：
PAX2	PPP2R4
hsa-miR-488	RAG2
ATF5	GNB5
hsa-miR-186	DCP2
MAX	ADAR

注意：这里没有表头，第一列是源节点，第二列是目标节点，两列之间用制表符分隔。

'''

def read_edge_list(file_path):
    with open(file_path, 'r') as f:
        edge_list = f.readlines()
        edge_list = [line.strip().split('\t') for line in edge_list]
    return edge_list

# 测试代码
if __name__ == '__main__':
    file_path = r'/home/liuyq/data/ar_data/edge_list_symbole.tsv'
    edge_list = read_edge_list(file_path)
    print(edge_list)
