# 从a 引入 Ab
from .association_rule import AbstractAssociationRule
import tensorflow as tf
from ..util_class.utilclass import lazyproperty
import pandas as pd


# 创建一个基于矩阵运算，一次能返回所有结果的关联规则算法
class MatrixRule(AbstractAssociationRule):
    '''
    该类实现的功能是一次性的将所有可能的基因对，通过矩阵运算，一次性的求出
    这里需要注意的是，该算法需要不断地筛选从而保证计算量降低到可以接受的范围内；
    需要注意的时返回的矩阵 i -> j 的值，在矩阵中是第i列，第j行的值，也就是从列到行的值
    '''

    # 计算所有基因的支持度
    @lazyproperty
    def support_all(self):
        # 对self.data 每一行求和
        sum = tf.reduce_sum(self.data, axis=1)
        support = tf.divide(sum, self.data.shape[1])
        return support

    # 计算所有基因对的支持度
    @lazyproperty
    def support_pair_all(self):
        # self.data 矩阵乘 self.data的转置
        support_m = tf.matmul(self.data, self.data, transpose_b=True)
        #support = self.support_all()
        support_m = tf.divide(support_m, self.data.shape[1])
        return support_m


    # 计算所有基因对的置信度, 这里假设support_m 的第i行基因和第j列的支持度，为j -> i 的支持度, 因为广播机制的原理是右端对齐。
    @lazyproperty
    def confidence_pair_all(self):
        # 所有基因的支持度，除以所有基因对的支持度
        support = self.support_all
        support_m = self.support_pair_all
        confidence = tf.divide(support_m, support)
        return confidence

    # 计算所有基因对的提升度
    @lazyproperty
    def lift_pair_all(self):
        # 所有基因对的置信度， 除以所有基因的支持度, 注意这里的广播机制，通过首先创建了一个维度为1的张量 shape=(n_gene, 1)，然后通过广播机制，将其扩展到了所有基因对的数量 confidence shape=(n_gene, n_gene)
        support = self.support_all
        support = tf.expand_dims(support, axis=1)
        confidence = self.confidence_pair_all
        lift = tf.divide(confidence, support)
        return lift

    # 计算所有基因对的leverage levarage(A→C)=support(A→C)−support(A)×support(C),range: [−1,1]
    @lazyproperty
    def leverage_pair_all(self):
        #one = tf.constant(1, dtype=tf.float32)
        support = self.support_all
        support_A = tf.expand_dims(support, axis=1)
        support_B = tf.transpose(support_A)

        support_m = self.support_pair_all
        leverage = support_m - tf.matmul(support_A, support_B)
        return leverage

    # 计算所有基因对的conviction conviction(A→C)=(1−support(C))/(1−confidence(A→C)),range: [0,∞]
    @lazyproperty
    def conviction_pair_all(self):
        one = tf.constant(1, dtype=tf.float32)
        support = self.support_all
        support = tf.expand_dims(support, axis=1)
        confidence = self.confidence_pair_all
        conviction = tf.divide(one - support, one - confidence)
        return conviction

    @lazyproperty
    def all_metrics(self):
        return {
            'support': self.support_pair_all,
            'confidence': self.confidence_pair_all,
            'lift': self.lift_pair_all,
            'leverage': self.leverage_pair_all,
            'conviction': self.conviction_pair_all
        }

    @lazyproperty
    def all_metrics_to_python(self):
        return {
            'support': self.support_pair_all.numpy(),
            'confidence': self.confidence_pair_all.numpy(),
            'lift': self.lift_pair_all.numpy(),
            'leverage': self.leverage_pair_all.numpy(),
            'conviction': self.conviction_pair_all.numpy()
        }

    #将所有的指标转化为dataframe，注意这里的每一个指标都是一个矩阵，思路是将每一个矩阵先转换为dataframe,添加一列作为指标名，然后将dataframe合并
    # 返回的结果为列名和索引名为基因名，因为是多个指标的矩阵，按照行的方式合并，所以行索引存在重复，这样能够根据需要提取对应的值，如，需要提取第i个基因的所有指标值，可以用
    #all_metrics_df.loc[i, :], 这样就能提取第i个基因的所有指标值

    
    @lazyproperty
    def all_metrics_to_dataframe(self, with_name=True):
        all_metrics = self.all_metrics_to_python
        all_metrics_df = None
        for key, value in all_metrics.items():

            #如果with_name 为True, 则将行索引设置为基因名，列名也设置为基因名
            if with_name:
                value_df = pd.DataFrame(value, index=self.genes_info, columns=self.genes_info)
            else:
                value_df = pd.DataFrame(value)
            # 添加一列作为指标名, 并放在第一列
            value_df.insert(0, 'metric', key)

            if all_metrics_df is None:
                all_metrics_df = value_df
            else:
                # 将value_df 添加到all_metrics_df的下面
                all_metrics_df = pd.concat([all_metrics_df, value_df], axis=0, ignore_index=False)

        return all_metrics_df


