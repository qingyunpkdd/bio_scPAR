from .matrix_association_rule import MatrixRule
from .association_rule import AssociationRule


'''
在实际部署中，显存的消耗是一个很大的问题，因此需要对矩阵关联规则进行优化，这里的优化是将矩阵关联规则的计算过程进行分块，这样可以减少显存的消耗，提高计算效率
'''



class MatrixRuleBlock(MatrixRule):
