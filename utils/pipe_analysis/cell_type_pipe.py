'''
读取特定文件夹下的数据，该注塑机包括表达数据和细胞类型注释数据
'''
from utils.data_process.read_cell_type import read_cell_type
from utils.data_process.load_data import LoadMatrixDataReal, TransformDataReal
from utils.algorithms.matrix_association_rule import MatrixRule
import os
from utils.data_process.ar_metrics_process import SaveArMetrics, LoadArMetrics, FilterArMetrics
from utils.data_process.read_gene_list import read_gene_list



class CellTypePipe:
    def __init__(self, data_path,
                 ar_result_path,
                 cell_type_fn="cell_type.tsv",
                 threshold_dict={"support":0.1, "confidence":0.5}):
        self.data_path = data_path
        self.ar_result_path = ar_result_path
        self.cell_type_path = os.path.join(data_path, cell_type_fn)
        self.threshold_dict = threshold_dict
        self.min_support = threshold_dict["support"]

    def read_data(self):
        # 加载数据
        data_obj = LoadMatrixDataReal(self.data_path)
        _data = data_obj.load_data()
        _data = data_obj.filter_genes(_data, threshold=self.min_support)

        # 读取cell_type.tsv文件
        cell_type = read_cell_type(self.cell_type_path)
        adata_subs = data_obj.split_by_cell_type(_data, cell_type)
        return adata_subs, data_obj

    def run_ar(self, adata_subs, data_obj, ):
        results = {}
        for cell_type, adata in adata_subs.items():
            adata = data_obj.filter_genes(adata)
            transform_data_obj = TransformDataReal(adata)
            data, genes_info, cells_info = transform_data_obj.run()
            mar_obj = MatrixRule(data, genes_info, cells_info)
            metrics = mar_obj.all_metrics_to_dataframe
            results[cell_type] = metrics
        return results

    def filter_result(self, results):
        results_filter = FilterArMetrics(results, threshold_dict=self.threshold_dict)
        res1 = results_filter.filter_result()
        return res1

    def extract_by_genes(self, results, gene_list):
        results_filter = FilterArMetrics(results, threshold_dict=self.threshold_dict)
        res1 = results_filter.extract_by_genes(gene_list, results)
        return res1

    def extract_by_gene_pairs(self, results, gene_pair):
        results_filter = FilterArMetrics(results, threshold_dict=self.threshold_dict)
        res1 = results_filter.extract_by_gene_pairs(gene_pair, results)
        return res1

    def save_result(self, results):
        save_obj = SaveArMetrics(results, self.ar_result_path)
        save_obj.write_result()
        print("finish writing result...")

    def load_result(self):
        armetrics = LoadArMetrics(self.ar_result_path)
        results = armetrics.load_result()
        return results












































