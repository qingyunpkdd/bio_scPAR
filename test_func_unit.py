'''
通过对测试代码的打包和封装，实现测试过程的可复用性
'''

# 将测试测试AssociationRule的方法和MatrixAssociationRule类的所有代码封装到一个类中
from utils.algorithms.association_rule import AssociationRule
from utils.algorithms.matrix_association_rule import MatrixRule
from utils.data_process.load_matrix_data import LoadMatrixData, TransformData
from utils.data_process.load_data import LoadMatrixDataReal, TransformDataReal

class TestAssociationRule:
    def __init__(self, data_path):
        data_obj = LoadMatrixData(data_path)
        _data = data_obj.load_data()
        transform_data_obj = TransformData(self._data)
        self.data, self.genes_info, self.cells_info = transform_data_obj.run()
        print("start to test...")
        #self.ar_obj = AssociationRule(self.data, self.genes_info, self.cells_info)

    #打印字典里面的键值对结果
    def print_dict(self, dict):
        for key, value in dict.items():
            print(key, ":", "\n",value)
        print("\n"*2)

    def test_association_rule(self, i, j):
        # # 测试AssociationRule类的方法
        ar_obj = AssociationRule(self.data, self.genes_info, self.cells_info)
        # # 测试support_single方法
        # support_single = ar_obj.support_single(i)
        # # 测试support_multi方法
        # support_multi = ar_obj.support_multi(i, j)
        # # 测试confidence方法
        # confidence = ar_obj.confidence(i, j)
        # # 测试lift方法
        # lift = ar_obj.lift(i, j)
        # # 测试leverage方法
        # leverage = ar_obj.leverage(i, j)
        # # 测试conviction方法
        # conviction = ar_obj.conviction(i, j)
        # 测试all_metrics方法
        metrics = ar_obj.all_metrics_to_python(i, j)
        print("*"*10)
        print("start to test association rule...")
        print("--"*20)
        self.print_dict(metrics)
        print("--"*20)
        print("finish testing association rule...")
        print("\n"*6)



    def test_matrix_rule(self):
        # 测试MatrixAssociationRule类的方法

        mar_obj = MatrixRule(self.data, self.genes_info, self.cells_info)

        # 将tensorflow tensor转化为numpy array

        def tensor_to_numpy(tensor):
            return tensor.numpy()

        # support_all = mar_obj.support_all
        # support_all = tensor_to_numpy(support_all)
        #
        # support_pair_all = mar_obj.support_pair_all
        # support_pair_all = tensor_to_numpy(support_pair_all)
        #
        # confidence_pair_all = mar_obj.confidence_pair_all
        # confidence_pair_all = tensor_to_numpy(confidence_pair_all)
        #
        # lift_pair_all = mar_obj.lift_pair_all
        # lift_pair_all = tensor_to_numpy(lift_pair_all)
        #
        # leverage_pair_all = mar_obj.leverage_pair_all
        # leverage_pair_all = tensor_to_numpy(leverage_pair_all)
        #
        # conviction_pair_all = mar_obj.conviction_pair_all
        # conviction_pair_all = tensor_to_numpy(conviction_pair_all)

        metrics =mar_obj.all_metrics_to_python
        print("*"*10)
        print("start to test matrix association rule...")
        print("--"*20)
        self.print_dict(metrics)
        print("--"*20)
        print("finish testing matrix association rule...")

#基于real数据的测试
class TestAssociationRuleReal(TestAssociationRule):
    def __init__(self, data_path):
        data_obj = LoadMatrixDataReal(data_path)
        _data = data_obj.load_data()
        _data = data_obj.filter_genes(_data)
        transform_data_obj = TransformDataReal(_data)
        self.data, self.genes_info, self.cells_info = transform_data_obj.run()
        print("start to test...")
        #self.ar_obj = AssociationRule(self.data, self.genes_info, self.cells_info)



# 测试代码
if __name__ == '__main__':
    # 测试example_data1.tsv数据
    data_path = r'/home/liuyq/data/ar_data/example_data1.tsv'
    test_obj = TestAssociationRule(data_path)
    test_obj.test_association_rule(0, 1)
    test_obj.test_matrix_rule()

    # 测试example_data2.tsv数据
    data_path = r'/home/liuyq/data/ar_data/example_data2.tsv'
    test_obj = TestAssociationRule(data_path)
    test_obj.test_association_rule(0, 1)
    test_obj.test_matrix_rule()

    # 测试真实数据
    data_path = r'/home/liuyq/data/ar_data/filtered_gene_bc_matrices/hg19'
    test_obj = TestAssociationRuleReal(data_path)
    test_obj.test_association_rule(0, 1)
    test_obj.test_matrix_rule()
    print("finish testing...")



