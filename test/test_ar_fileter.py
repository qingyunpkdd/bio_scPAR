import os
from utils.data_process.ar_metrics_process import SaveArMetrics, LoadArMetrics, FilterArMetrics
from utils.data_process.read_gene_list import read_gene_list

armetrics = LoadArMetrics(r'/home/liuyq/data/ar_data/filtered_gene_bc_matrices/result/sub_test')
results = armetrics.load_result()

results_filter = FilterArMetrics(results, {"support":0.4, "confidence":0.5})
res1 = results_filter.filter_result()
res2 = results_filter.merge_cell_type(res1)
gene_list_fp = r"/home/liuyq/data/ar_data/interaction_db/reactome/gene_list.txt"
gene_list = read_gene_list(gene_list_fp)
#测试extract_by_genes方法
res3 = results_filter.extract_by_genes(gene_list, res1)
gene_pair = [("A4GNT", "MUC1"),
             ("A4GNT", "MUC12"),
             ("A4GNT", "MUC13"),
             ("ISG15", "ENO1"),
             ("RPS8", "ISG15")]

pair_res = results_filter.extract_by_gene_pairs(gene_pair, res1)



out_path =  r'/home/liuyq/data/ar_data/filtered_gene_bc_matrices/result'
save_obj = SaveArMetrics(results, out_path)
save_obj.write_result()
print("finish writing result...")





# def write_result(results, path):
#     if not os.path.exists(path):
#         os.makedirs(path)
#     for key, value in results.items():
#         value.to_csv(os.path.join(path, key+".csv"))
# out_path =  r'/home/liuyq/data/ar_data/filtered_gene_bc_matrices/result'
# write_result(results, out_path)
# print("finish testing matrix association rule...")
