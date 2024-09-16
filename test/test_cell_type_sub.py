
from utils.data_process.read_cell_type import read_cell_type
from utils.data_process.load_data import LoadMatrixDataReal, TransformDataReal
from utils.algorithms.matrix_association_rule import MatrixRule

# 加载数据
data_path = r'/home/liuyq/data/ar_data/filtered_gene_bc_matrices/hg19'
data_obj = LoadMatrixDataReal(data_path)
_data = data_obj.load_data()
_data = data_obj.filter_genes(_data)

# 读取cell_type.tsv文件
cell_type_path = r'/home/liuyq/data/ar_data/filtered_gene_bc_matrices/hg19/cell_type.tsv'
cell_type = read_cell_type(cell_type_path)
adata_subs = data_obj.split_by_cell_type(_data, cell_type)

#
results = {}


for cell_type, adata in adata_subs.items():
    adata = data_obj.filter_genes(adata)
    transform_data_obj = TransformDataReal(adata)
    data, genes_info, cells_info = transform_data_obj.run()
    mar_obj = MatrixRule(data, genes_info, cells_info)
    metrics = mar_obj.all_metrics_to_dataframe
    # metrics.loc[1,["metric",2]] 可以提取到第三个基因到第二个基因的所有指标
    # 如下所示：
    '''
    1,support,0.108187
    1,confidence,1.000000
    1,lift,9.243243
    1,leverage,0.096483
    1,conviction,inf
    '''

    print("*"*10)
    print(cell_type)
    print(metrics)
    results[cell_type] = metrics

#写一个函数，指定目录， 以result的键为文件名，将result的值写入文件
import os
from utils.data_process.ar_metrics_process import SaveArMetrics, LoadArMetrics, FilterArMetrics
from utils.data_process.read_gene_list import read_gene_list

results_filter = FilterArMetrics(results, {"support":0.1, "confidence":0.5})
res1 = results_filter.filter_result()
res2 = results_filter.merge_cell_type(res1)
gene_list_fp = r"/home/liuyq/data/ar_data/interaction_db/reactome/gene_list.txt"
gene_list = read_gene_list(gene_list_fp)
#测试extract_by_genes方法
res3 = results_filter.extract_by_genes(gene_list, res1)
gene_pair = [("A4GNT", "MUC1"),
             ("A4GNT", "MUC12"),
             ("A4GNT", "MUC13"),
             ("A4GNT", "MUC15"),
             ("AAG1", "CASP7")]

pair_res = results_filter.extract_by_gene_pairs(gene_pair)



out_path =  r'/home/liuyq/data/ar_data/filtered_gene_bc_matrices/result'
save_obj = SaveArMetrics(results, out_path)
save_obj.write_result()
print("finish writing result...")





def write_result(results, path):
    if not os.path.exists(path):
        os.makedirs(path)
    for key, value in results.items():
        value.to_csv(os.path.join(path, key+".csv"))
out_path =  r'/home/liuyq/data/ar_data/filtered_gene_bc_matrices/result'
write_result(results, out_path)
print("finish testing matrix association rule...")








