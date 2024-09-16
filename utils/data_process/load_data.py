# import scanpy package
import os

import scanpy as sc
import tensorflow as tf
import scipy.sparse as sp
import numpy as np
import gc
import sys
import anndata

'''
cell_type.tsv 是一个细胞类型注释文件，只有一列，和barcodes.tsv中的barcode一一对应，每个barcode对应一个细胞类型。

'''


def read_cell_type(file_path):
    with open(file_path, 'r') as f:
        cell_type = f.readlines()
        # 去掉每行末尾的换行符
        cell_type = [line.strip() for line in cell_type]
    return cell_type


class LoadMatrixDataReal:
    def __init__(self, data_path, cache_path="/home/liuyq/data/ar_data"):
        self.data_path = data_path
        cache_dir = cache_path
        os.environ['SCANPY_TEMPDIR'] = cache_dir

    def load_data(self):
        adata = sc.read_10x_mtx(self.data_path, var_names='gene_symbols', cache=False)
        # 设置一个表达阈值，例如1
        threshold = 1

        # 获取adata.X的行和列索引
        rows, cols = adata.X.nonzero()

        # 获取非零元素的值
        values = adata.X.data

        # 应用阈值，将值转换为0和1
        binary_values = np.where(values >= threshold, 1, 0)

        # 创建一个新的稀疏矩阵，其形状与原始矩阵相同
        binary_matrix = sp.csr_matrix((binary_values, (rows, cols)), shape=adata.X.shape)

        # 创建一个新的AnnData对象或更新现有的AnnData对象
        adata_binary = anndata.AnnData(binary_matrix, obs=adata.obs, var=adata.var)

        return adata_binary

    # 过滤部分基因，只保留在细胞中表达占比超过0.1的基因
    def filter_genes(self, adata, threshold=0.1):
        # 计算每个基因在每个细胞中的表达占比
        gene_expression_ratio = adata.X.sum(axis=0) / adata.X.shape[0]
        # 过滤基因
        adata = adata[:, gene_expression_ratio > threshold]
        return adata

    # 根据提供的基因名列表，过滤adata中的基因
    def filter_genes_by_name(self, adata, gene_names):
        adata = adata[:, gene_names]
        return adata

    # cell_type 为和adata.obs中的cell_type列对应的cell type名称，根据cell type将adata分成多个子集
    def split_by_cell_type(self, adata, cell_type):
        adata.obs['cell_type'] = cell_type
        adata_subs = {}
        for cell_type in adata.obs['cell_type'].unique():
            adata_sub = adata[adata.obs['cell_type'] == cell_type]
            adata_subs[cell_type] = adata_sub
        return adata_subs


# class TransformData, transform data to tensorflow tensor
class TransformDataReal:
    def __init__(self, adata):
        self.adata = adata

    def transform_data(self):
        # 将scipy稀疏矩阵转换为TensorFlow的SparseTensor
        sys.getsizeof(self.adata)
        sparse_tensor = tf.sparse.from_dense(self.adata.X.todense())
        # sparse_tensor = tf.constant(self.adata.X.toarray())
        # self.adata = None
        print("稀疏矩阵大小：", tf.size(sparse_tensor))
        data = tf.cast(sparse_tensor, tf.float32)
        print("float32大小：", tf.size(data))
        # 将data转置
        data = tf.sparse.transpose(data)
        # data = tf.cast(data>0, tf.float32)
        del sparse_tensor
        gc.collect()
        # 将data转置
        print("load data success")
        # print(data)
        # # 如果需要，将SparseTensor转换为密集Tensor
        dense_tensor = tf.sparse.to_dense(data)
        # dense_tensor = sparse_tensor
        # data = dense_tensor
        # #data = tf.cast(dense_tensor, tf.float32)
        # #将data 矩阵转化为0,1矩阵
        # #data = tf.where(data < 0.5, tf.zeros_like(data), tf.ones_like(data))
        # #将data转置
        # data = tf.cast(data, tf.bool)
        # gc.collect()
        # data = tf.transpose(data)

        return dense_tensor

    def get_data(self):
        return self.data

    def get_genes_info(self):
        return self.adata.var_names

    def get_cells_info(self):
        return self.adata.obs_names

    def run(self):
        self.tf_data = self.transform_data()
        return self.tf_data, self.adata.var_names, self.adata.obs_names


# 测试加载数据
if __name__ == '__main__':
    data_path = r'/mnt/sda/liuyq/ar_data/sc_data/liver'
    data_obj = LoadMatrixDataReal(data_path)
    _data = data_obj.load_data()

    # 读取cell_type.tsv文件
    cell_type_path = r'/mnt/sda/liuyq/ar_data/sc_data/liver/cell_type.tsv'
    cell_type = read_cell_type(cell_type_path)
    adata_subs = data_obj.split_by_cell_type(_data, cell_type)

    for cell_type, adata in adata_subs.items():
        adata = data_obj.filter_genes(adata)
        print(adata)

        # 测试TransformData类
        transform_data_obj = TransformDataReal(adata)
        data, genes_info, cells_info = transform_data_obj.run()
        print(data)
        print(genes_info)
        print(cells_info)
