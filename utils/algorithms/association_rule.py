'''
基于TensorFlow框架，实现一个计算关联规则的算法
'''

# 导入TensorFlow
import tensorflow as tf


# # import 10x single-cell RNA-seq data
# data = tf.random.normal((10, 1000))
# data = tf.cast(data, tf.float32)
#
# #将data 矩阵转化为0,1矩阵
# data[data>0] = 1
# data[data<0] = 0


# # 统计data的列数C, 对矩阵的每一行求和，并除以C
# data = tf.math.divide_no_nan(tf.math.reduce_sum(data, axis=1), 1000)


class AbstractAssociationRule:
    def __init__(self, data, genes_info, cells_info):
        self.data = data
        self.genes_info = genes_info
        self.cells_info = cells_info


# 定义一个计算每一个关联规则算法的类

class AssociationRule(AbstractAssociationRule):
    '''
    data 为输入的表达矩阵，注意这里的元素是count而不是归一化，或者标准化后的结果。
    每一行代表一个基因，每一列代表一个样本。
    genes_info 为一个字典，键为基因名，值为对应的行号；
    cell_info 为一个字典，键为样本名，值为对应的列号。
    '''

    # 单个基因的支持度
    def support_single(self, i):
        # 获取self.data 的列数
        C = self.data.shape[1]
        # 计算第i行的元素个数
        count = tf.math.count_nonzero(self.data[i, :])
        # 计算第i列的元素个数占所有元素个数的比例
        support = count / C
        #将support 转换为float32类型
        support = tf.cast(support, tf.float32)
        return support

    # 多个基因的支持度
    def support_multi(self, i, j):
        # 获取self.data 的列数
        C = self.data.shape[1]

        # 提取self.data 的第i行和第j行
        row_i = self.data[i, :]
        # 将row_i 转置为 1行 n列的张量
        row_i = tf.expand_dims(row_i, axis=0)

        row_j = self.data[j, :]
        # 将row_j 转置为 n行 1列的张量
        row_j_transpose = tf.expand_dims(row_j, axis=1)

        # 计算两个矩阵的矩阵乘积，即为 i 和 j 的支持度
        support_ij = tf.matmul(row_i, row_j_transpose) / C
        return support_ij

    # 计算基因i 和基因j的置信度
    def confidence(self, i, j):
        # 获取self.data 的列数
        C = self.data.shape[1]

        # 首先计算 i 和 j 的支持度
        support_ij = self.support_multi(i, j)
        # 计算 i 的支持度
        support_i = self.support_single(i)
        confidence = support_ij / support_i
        return confidence

    # 计算基因i 和基因j的提升度
    def lift(self, i, j):
        # 获取self.data 的列数
        C = self.data.shape[1]

        # 首先计算 i 和 j 的置信度
        confidence = self.confidence(i, j)

        # 计算 j 的支持度
        support_j = self.support_single(j)

        # 计算提升度
        lift = confidence / support_j
        return lift

    # 计算 基因i 和基因j 的 leverage
    def leverage(self, i, j):
        # 获取self.data 的列数
        C = self.data.shape[1]
        # 首先计算 i 和 j 的支持度
        support_ij = self.support_multi(i, j)
        # 计算 i 的支持度
        support_i = self.support_single(i)
        # 计算 j 的支持度
        support_j = self.support_single(j)
        leverage = support_ij - support_i * support_j
        return leverage

    # 计算 基因i 和基因j 的 conviction
    def conviction(self, i, j):
        # 获取self.data 的列数
        C = self.data.shape[1]
        # 首先计算 i 和 j 的置信度
        confidence = self.confidence(i, j)
        # 计算 j 的支持度
        support_j = self.support_single(j)
        conviction = (1 - support_j) / (1 - confidence)
        return conviction

    # 同时计算多个基因的所有指标
    def all_metrics(self, i, j):
        support_ij = self.support_multi(i, j)
        support_i = self.support_single(i)
        support_j = self.support_single(j)
        confidence = self.confidence(i, j)
        lift = self.lift(i, j)
        leverage = self.leverage(i, j)
        conviction = self.conviction(i, j)
        # 返回一个字典，包含所有指标
        return {
            'support_ij': support_ij,
            'support_i': support_i,
            'support_j': support_j,
            'confidence': confidence,
            'lift_ij': lift,
            'leverage_ij': leverage,
            'conviction_ij': conviction
        }
    #将all_metrics里面所有的Tensor的结果转换为python的数据类型
    def all_metrics_to_python(self, i, j):
        metrics = self.all_metrics(i, j)
        metrics = {k: v.numpy() for k, v in metrics.items()}
        return metrics



class AssociationRuleTrans(AbstractAssociationRule):
    '''
    data 为输入的表达矩阵，注意这里的元素是count而不是归一化，或者标准化后的结果。
    每一行代表一个样本，每一列代表一个基因。
    genes_info 为一个字典，键为基因名，值为对应的行号；
    cell_info 为一个字典，键为样本名，值为对应的列号。
    '''

    # 单个基因的支持度
    def support_single(self, i):
        # 获取self.data 的列数
        C = self.data.shape[0]
        # 计算第i行的元素个数
        count = tf.math.count_nonzero(self.data[:, i])
        # 计算第i列的元素个数占所有元素个数的比例
        support = count / C
        return support

    # 多个基因的支持度
    def support_multi(self, i, j):
        # 获取self.data 的列数
        C = self.data.shape[0]

        # 提取self.data 的第i列和第j列
        col_i = self.data[:, i]
        col_j = self.data[:, j]
        # col_i 为 n行 1列的张量，col_j 为 n行 1列的张量， col
        # 计算两个矩阵的矩阵乘积，即为 i 和 j 的支持度
        support_ij = tf.matmul(col_i, col_j, transpose_a=True) / C
        return support_ij

    # 计算基因i 和基因j的置信度
    def confidence(self, i, j):
        # 获取self.data 的列数
        C = self.data.shape[0]

        # 首先计算 i 和 j 的支持度
        support_ij = self.support_multi(i, j)
        # 计算 i 的支持度
        support_i = self.support_single(i)
        confidence = support_ij / support_i
        return confidence

    # 计算基因i 和基因j的提升度
    def lift(self, i, j):
        # 获取self.data 的列数
        C = self.data.shape[0]

        # 首先计算 i 和 j 的置信度
        confidence = self.confidence(i, j)

        # 计算 j 的支持度
        support_j = self.support_single(j)

        # 计算提升度
        lift = confidence / support_j
        return lift

    # 计算 基因i 和基因j 的 leverage
    def leverage(self, i, j):
        # 获取self.data 的列数
        C = self.data.shape[0]
        # 首先计算 i 和 j 的支持度
        support_ij = self.support_multi(i, j)
        # 计算 i 的支持度
        support_i = self.support_single(i)
        # 计算 j 的支持度
        support_j = self.support_single(j)
        leverage = support_ij - support_i * support_j
        return leverage

    # 计算 基因i 和基因j 的 conviction
    def conviction(self, i, j):
        # 获取self.data 的列数
        C = self.data.shape[0]
        # 首先计算 i 和 j 的置信度
        confidence = self.confidence(i, j)
        # 计算 j 的支持度
        support_j = self.support_single(j)
        conviction = (1 - support_j) / (1 - confidence)
        return conviction

    # 同时计算多个基因的所有指标
    def all_metrics(self, i, j):
        support_ij = self.support_multi(i, j)
        support_i = self.support_single(i)
        support_j = self.support_single(j)
        confidence = self.confidence(i, j)
        lift = self.lift(i, j)
        leverage = self.leverage(i, j)
        conviction = self.conviction(i, j)
        # 返回一个字典，包含所有指标
        return {
            'support_ij': support_ij,
            'support_i': support_i,
            'support_j': support_j,
            'confidence': confidence,
            'lift_ij': lift,
            'leverage_ij': leverage,
            'conviction_ij': conviction
        }
