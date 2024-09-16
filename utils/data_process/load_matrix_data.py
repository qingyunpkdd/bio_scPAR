'''
加载样例数据，如一个0,1矩阵，从tsv文件当中，下面是一个数据的例子：
        c1	c2	c3	c4	c5
    g1	1	1	0	1	1
    g2	1	0	1	1	1
    g3	1	1	0	1	1
    g4	1	1	1	0	1
    g5	1	1	1	0	0
'''

import pandas as pd

class LoadMatrixData:

    def __init__(self, file_path):
        self.file_path = file_path

    def load_data(self):
        _data = pd.read_csv(self.file_path, sep='\t', index_col=0)
        return _data

# 将pandas dataframe转化为tensorflow tensor, 列名赋给自身的cells_info属性, 行名赋给自身的genes_info属性, 将cells_info属性和genes_info属性转换为字符串类型
import tensorflow as tf


class TransformData:
    def __init__(self, data):
        self.data = data
        self.genes_info = data.index.astype(str).to_list()
        self.cells_info = data.columns.astype(str).to_list()
        self.tf_data = None

    def transform_data(self):
        data = tf.convert_to_tensor(self.data.values)
        data = tf.cast(data, tf.float32)
        return data

    def get_data(self):
        return self.data

    def get_genes_info(self):
        return self.genes_info

    def get_cells_info(self):
        return self.cells_info

    # 运行transform_data、get_data、get_genes_info、get_cells_info方法，并将结果赋给self.data, self.genes_info, self.cells_info
    def run(self):
        self.tf_data = self.transform_data()
        return self.tf_data, self.genes_info, self.cells_info





# 生成main函数，从目录F:\project\teacher_shi_task\data\example_data1.tsv中加载数据

if __name__ == '__main__':

    data_path = r'/home/liuyq/data/ar_data/example_data1.tsv'
    data_obj = LoadMatrixData(data_path)
    data = data_obj.load_data()
    print(data)

    #测试TransformData类
    transform_data_obj = TransformData(data)
    data, genes_info, cells_info = transform_data_obj.run()
    print(data)
    print(genes_info)
    print(cells_info)




