'''
通过加载样例数据来测试association_rules.py中的类和函数

'''


# 导入LoadMtrixData类
from utils.data_process.load_matrix_data import LoadMatrixData, TransformData
from utils.algorithms.association_rule import AssociationRule

data_path = r'/home/liuyq/data/ar_data/example_data1.tsv'
data_obj = LoadMatrixData(data_path)
data = data_obj.load_data()

# 将data转化为tensorflow tensor
transform_data_obj = TransformData(data)
data, genes_info, cells_info = transform_data_obj.run()
ar_obj = AssociationRule(data, genes_info, cells_info)


# 测试AssociationRule类的方法

#测试ar_obj的所有方法
#测试support_single方法
support = ar_obj.support_single(0)

#测试support_multi方法
support_mulpy = ar_obj.support_multi(0, 1)
#测试confidence方法
confidence = ar_obj.confidence(0, 1)
#测试lift方法
lift = ar_obj.lift(0, 1)
#测试leverage方法
leverage = ar_obj.leverage(0, 1)
#测试conviction方法
conviction = ar_obj.conviction(0, 1)
#测试all_metrics方法
metrics = ar_obj.all_metrics_to_python(0, 1)
print(metrics)


# 下面的代码为测试MatrixAssociationRule类的方法
from utils.algorithms.matrix_association_rule import MatrixRule

#测试MatrixAssociationRule类的方法

mar_obj = MatrixRule(data, genes_info, cells_info)

# 将tensorflow tensor转化为numpy array

def tensor_to_numpy(tensor):
    return tensor.numpy()


support_all = mar_obj.support_all
support_all = tensor_to_numpy(support_all)

support_pair_all = mar_obj.support_pair_all
support_pair_all = tensor_to_numpy(support_pair_all)

confidence_pair_all = mar_obj.confidence_pair_all
confidence_pair_all = tensor_to_numpy(confidence_pair_all)

lift_pair_all = mar_obj.lift_pair_all
lift_pair_all = tensor_to_numpy(lift_pair_all)

leverage_pair_all = mar_obj.leverage_pair_all
leverage_pair_all = tensor_to_numpy(leverage_pair_all)

conviction_pair_all = mar_obj.conviction_pair_all
conviction_pair_all = tensor_to_numpy(conviction_pair_all)

print("123")





# 将测试测试AssociationRule的方法和MatrixAssociationRule类的所有代码封装到一个类中
from utils.algorithms.association_rule import AssociationRule
from utils.algorithms.matrix_association_rule import MatrixRule
from utils.data_process.load_matrix_data import LoadMatrixData, TransformData

class TestAssociationRule:
    def __init__(self, data_path):
        self.data_path = data_path
        self.data_obj = LoadMatrixData(data_path)
        self.data = self.data_obj.load_data()
        self.transform_data_obj = TransformData(self.data)
        self.data, self.genes_info, self.cells_info = self.transform_data_obj.run()
        self.ar_obj = AssociationRule(self.data, self.genes_info, self.cells_info)

    def test_all_methods(self):
        # 测试AssociationRule类的方法
        # 测试support_single方法
        support = self.ar_obj.support_single(0)

        # 测试support_multi方法
        support_mulpy = self.ar_obj.support_multi(0, 1)
        # 测试confidence方法
        confidence = self.ar_obj.confidence(0, 1)
        # 测试lift方法
        lift = self.ar_obj.lift(0, 1)
        # 测试leverage方法
        leverage = self.ar_obj.leverage(0, 1)
        # 测试conviction方法
        conviction = self.ar_obj.conviction(0, 1)
        # 测试all_metrics方法
        metrics = self.ar_obj.all_metrics_to_python(0, 1)
        print(metrics)

    def test_matrix_rule(self):
        # 测试MatrixAssociationRule类的方法

        mar_obj = MatrixRule(self.data, self.genes_info, self.cells_info)

        # 将tensorflow tensor转化为numpy array

        def tensor_to_numpy(tensor):
            return tensor.numpy()

        support_all = mar_obj.support_all
        support_all = tensor_to_numpy(support_all)

        support_pair_all = mar_obj.support_pair_all
        support_pair_all = tensor_to_numpy(support_pair_all)

        confidence_pair_all = mar_obj.confidence_pair_all
        confidence_pair_all = tensor_to_numpy(confidence_pair_all)

        lift_pair_all = mar_obj.lift_pair_all
        lift_pair_all = tensor_to_numpy(lift_pair_all)

        leverage_pair_all = mar_obj.leverage_pair_all
        leverage_pair_all = tensor_to_numpy(leverage_pair_all)

        conviction_pair_all = mar_obj.conviction_pair_all
        conviction_pair_all = tensor_to_numpy(conviction_pair_all)

        print("123")







