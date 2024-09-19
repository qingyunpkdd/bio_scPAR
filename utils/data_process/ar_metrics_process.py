'''

将MatrixRule的all_metrics_to_dataframe的结果保存为csv文件

从用户的角度，是想要一个符合要求的基因数据框，用户可以根据自己的需要，对数据进行筛选。

'''


#保存结果
import os
import numpy as np
import pandas as pd
from copy import deepcopy
class SaveArMetrics:
    def __init__(self, results, out_path):
        self.results = results
        self.out_path = out_path

    def write_result(self):
        if not os.path.exists(self.out_path):
            os.makedirs(self.out_path)
        for key, value in self.results.items():
            value.to_csv(os.path.join(self.out_path, key+".csv"))
        print("finish writing result...")

#加载结果
class LoadArMetrics:
    def __init__(self, in_path):
        self.in_path = in_path

    def load_result(self):
        results = {}
        for file in os.listdir(self.in_path):
            if file.endswith(".csv"):
                key = file.split(".")[0]
                value = pd.read_csv(os.path.join(self.in_path, file), index_col=0)
                results[key] = value
        return results


#根据MatrixRule的all_metrics_to_dataframe的support, confidence, lift, leverage, conviction的阈值，过滤结果，参数threshold为一个字典，key为指标名，value为阈值


class FilterArMetrics:
    '''
     #threshold_dict的结构为样例为
    {"support":(0.1, 0.5), "confidence":(0.5, 1.0), "lift":(0.5, 1.0), "leverage":(0.5, 1.0), "conviction":(0.5, 1.0)}
    如果不确定请设置为 None


    #遍历指标，然后，根据阈值产生布尔值，然后根据布尔值过滤结果，将所有值都是False的行删除，将所有值都是False的列删除
    #这里的results的每一个值的结构如下：

         index,metric,gene1,gene2,gene3
         gene1,support,0.108187,0.108187,0.108187
         gene2,support,0.108187,0.108187,0.108187
         gene3,support,0.108187,0.108187,0.108187
         gene1,confidence,1.000000,1.000000,1.000000
         gene2,confidence,1.000000,1.000000,1.000000
         gene3,confidence,1.000000,1.000000,1.000000
         ...
    '''
    def __init__(self, results, threshold_dict):
        #这里results的结构是列->行的关联关系
        self.results = results
        self.threshold_dict = threshold_dict

    #指标名包括：support, confidence, lift, leverage, conviction
    #当AR的结果不是一个字典时，调用这个函数
    def wrapper_single(self, result, name="all"):
        self.results = {name: result}


    # results的结构是{"细胞类型1":df1, "细胞类型2":df2, ...}， 本函数将所有的df合并，并添加一个新的列"cell_type"，将其值设置为细胞类型
    def merge_cell_type(self, results):
        _results = deepcopy(self.results)
        for key, value in _results.items():
            #在第一列的位置插入一个新的列，列名为"cell_type"，列的值为key
            value.insert(loc=0, column="cell_type", value=key)
        result = pd.concat(_results.values())
        return result
    #这里的过滤是在矩阵的基础上进行的过滤。
    def filter_result(self):
        results = deepcopy(self.results)
        for key, value in results.items():
            for metric in self.threshold_dict.keys():
                value = self.filter_by_metric(value, metric, self.threshold_dict[metric])
            results[key] = value
        return results

    #提取指标的值，然后根据阈值过滤结果，首先将"metric"和原来的索引结合作为双重索引，然后将指标对应的行，所有列的值中小于阈值的元素设置为NaN，然后删除所有值都是NaN的行，记录所有值都是NaN的列，然后删除所有值都是NaN的列，最后将双重索引的第一层索引恢复为原来的索引
    def filter_by_metric(self, df, metric, threshold):
        #将"metric"和原来的索引结合作为双重索引
        df = df.set_index(["metric"], append=True)
        #将指标对应的行，所有列的值中小于阈值的元素设置为NaN
        if threshold[0] is not None:
            df.loc[(slice(None), metric), :] = df.loc[(slice(None), metric), :].apply(
                lambda x: x.where(x > threshold[0], np.nan))
        if threshold[1] is not None:
            df.loc[(slice(None), metric), :] = df.loc[(slice(None), metric), :].apply(
                lambda x: x.where(x < threshold[1], np.nan))
        #删除所有值都是NaN的行
        df = df.dropna(axis=0, how="all")
        #记录指标对应的行中所有列的值为NaN的列的索引记录下来
        drop_columns = df.loc[(slice(None), metric), :].isnull().all()
        drop_columns = drop_columns[drop_columns].index
        #删除所有值都是NaN的列
        df = df.drop(drop_columns, axis=1)
        #将双重索引的第一层索引恢复为原来的索引
        df = df.reset_index(level="metric")

        return df

    #根据提供的A->B的基因名，提取所有指标值, gene_type为antecedent或consequent，如果gene_type为antecedent，则从self.results_T中提取，否则从self.results中提取，其中参数genes为一个列表
    def extract_by_genes(self, genes, results, gene_type="antecedent"):
        assert gene_type in ["antecedent", "consequent"] #"gene_type must be antecedent or consequent"
        results = deepcopy(results)
        for key, value in results.items():
            results[key] = self.extract_by_genes_in_df(value, genes, gene_type)
        return results

    #提取指定的基因对应的所有指标值
    def extract_by_genes_in_df(self, df, genes, gene_type):
        genes = set(genes).intersection(set(df.columns))
        columns = ["metric",  *genes]
        if gene_type == "antecedent":
            #提取指定的基因对应的所有指标值
            df = df.loc[:, columns]
        elif gene_type == "consequent":
            #提取指定的基因对应的所有指标值
            df = df.loc[genes, :]
        else:
            raise ValueError("gene_type must be antecedent or consequent")
        return df

    #根据提供的A->B的基因名，提取所有指标值，其中参数gene_pairs为一个列表，列表中的元素为一个元组，元组中的第一个元素为antecedent基因，第二个元素为consequent基因
    def extract_by_gene_pairs(self, gene_pairs, results):
        _results = {}
        results = deepcopy(results)
        for key, value in results.items():
            _results[key] = self.extract_by_gene_pairs_in_df(value, gene_pairs)
        return _results
        # for key, value in self.results_T.items():
        #
        #     results[key] = self.extract_by_gene_pairs_in_df(value, gene_pairs)
        # return results

    #提取指定的基因对应的所有指标值
    def extract_by_gene_pairs_in_df(self, df, gene_pairs):
        #提取指定的基因对应的所有指标值, 首先遍历列表，然后元组中的第一个元素为antecedent基因，以行索引的方式提取，然后元组中的第二个元素为consequent基因，以列索引的方式提取, 最后的结果为一个dataframe
        results = None
        for gene_pair in gene_pairs:
            antecedent = gene_pair[0]
            consequent = gene_pair[1]
            if antecedent not in df.index or consequent not in df.columns:
                continue
            value = df.loc[consequent, ["metric", antecedent]]
            if results is None:
                results = value
            else:
                results = pd.concat([results, value])
        return results

    #将数据框形式的结果转换为基因对的形式，具体为
    '''
    gene1 gene2 support confidence lift leverage conviction
    IGF1  IGF2  0.1    0.2        0.3  0.4     0.5
    '''
    '''
    原数据框的形式为：
       metric gene1 gene2 gene3
    gene1 support 0.1   0.2   0.3
    gene2 support 0.1   0.2   0.3
    
    '''

    def transform_to_pairs(self, results):
        results = deepcopy(results)
        for key, value in results.items():
            results[key] = self.transform_to_pairs_in_df(value)
        return results

    def transform_to_pairs_in_df(self, results):
        results = deepcopy(results)
        for key, value in results.items():
            results[key] = self._transform_to_pairs_in_df(value)
        return results

    def _transform_to_pairs_in_df(self, df):
        if "metric" in df.columns:
            df.set_index(["metric"], append=True, inplace=True)
        #df.set_index(["metric"], append=True, inplace=True)
        df2 = df.stack()
        df2.index.names = ["antecedent", "metric", "consequent"]
        df3 = df2.unstack(level="metric")
        df3.reset_index(inplace=True)
        #删除antecedent 和 antecedent 重复的行
        df3 = df3[df3["antecedent"] != df3["consequent"]]
        return df3

    #从transform_to_pairs_in_df的结果中，删除不满足threshold_dict_tem 的对，threshold_dict_tem的结构为{"metric": threshold}
    def filter_pairs(self, results, threshold_dict_tem):
        results = deepcopy(results)
        for key, value in results.items():
            results[key] = self.filter_pairs_in_df(value, threshold_dict_tem)
        return results

    def filter_pairs_in_df(self, df, threshold_dict_tem):
        df = df.set_index(["antecedent", "consequent"])
        for metric, threshold in threshold_dict_tem.items():
            if threshold[0] is not None:
                df = df[df[metric] > threshold[0]]
            if threshold[1] is not None:
                df = df[df[metric] < threshold[1]]
        df = df.reset_index()
        return df

    #从transform_to_pairs_in_df的结果中提取满足要求的对，并存为两列的数据框，如果out_fp为None，则返回结果，否则将结果保存为csv文件
    def export_edge_list(self, results, out_dir=None):
        results = deepcopy(results)
        for key, value in results.items():
            out_fp = os.path.join(out_dir, key + ".csv")
            results[key] = self.export_edge_list_in_df(value, out_fp)
        return results

    def export_edge_list_in_df(self, df, out_fp=None):
        df = df.loc[:, ["antecedent", "consequent"]]
        if out_fp is not None:
            df.to_csv(out_fp, index=False)
        else:
            return df




    #构建新的all_metrics_to_dataframe的结果，将原数据框转置，然后将metric行作为列插入，并将原metric行删除
    # def build_new_result(self):
    #     for key, value in self.results.items():
    #         #提取并删除metric列
    #         metric = value.loc[:, "metric"]
    #         value = value.drop("metric", axis=1)
    #         #转置
    #         value = value.T
    #         #在第一列插入metric列
    #         value.insert(0, "metric", metric)
    #         self.results_T[key] = value
    #     return self.results_T


