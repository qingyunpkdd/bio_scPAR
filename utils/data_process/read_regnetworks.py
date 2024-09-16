





#首先读取第一个文件,文件的结构如下所示：
'''
#TF ID Target ID Up or Down or Unknown
AHR 196 CDKN1B 1027 -->
APLNR 187 PIK3C3 5289 -->
APLNR 187 PIK3R4 30849 -->
AR 367 KLK3 354 -->

'''
#其中第一行为注释行，从第二行开始为数据行，每一行有四列，分别为TF ID, Target ID, Up or Down or Unknown
import pandas as pd
class InteractionData:
    def __init__(self, fp, sepr=' '):
        self.fp = fp
        #存储每一行的数据,其中键为TF， 值为Target
        self.interactions = None
        self.sepr = sepr

    def __str__(self):
        return self.tf + " " + self.target + " " + self.up_down_unknown

    # 以pandas 数据框方式读取数据，数据分5列，分别为TF ID, Target ID, Up or Down or Unknown。第一行为注释行
    def read_data(self, file_path=None):
        if file_path is None:
            file_path = self.fp

        data = pd.read_csv(file_path, sep=self.sepr, skiprows=1, header=None)
        data.columns = ['TF', 'ID1', 'Target', "ID2", "Up_or_Down_or_Unknown"]
        self.interactions = data
        return self.interactions





# 读取所有的TF name 文件,文件的结构如下所示：
'''
AC008770.3
AC023509.3
AC092835.1
AC138696.1
ADNP
ADNP2
AEBP1
AEBP2
AHCTF1
'''

# 读取所有的Target name 文件

class TFData:
    def __init__(self, fp):
        self.fp = fp
        self.tf = []

    def read_data(self, file_path=None):
        if file_path is None:
            file_path = self.fp

        with open(file_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                self.tf.append(line.strip())

        return self.tf


#新建一个dataframe，创建TargetData对象和InteractionData对象，遍历InteractionData对象的interaction属性创建一个新的dataframe，添加一列“type”，并根据其值是否在TargetData对象的targets属性中，将其值设置为“TF”或“Gene”

import pandas as pd
import os

class RegulationNetwork:
    def __init__(self, tf_fp, interaction_fp):
        self.tf_fp = tf_fp
        self.interaction_fp = interaction_fp
        self.tf_data = TFData(self.tf_fp)
        self.interaction_data = InteractionData(self.interaction_fp)
        self.data = None
        self.data_by_type = None

    def create_network(self,return_df=True):
        tf = self.tf_data.read_data(self.tf_fp)
        interactions = self.interaction_data.read_data()
        data = interactions[['TF', 'Target']]
        data['Type'] = 'Gene'
        data.loc[data['TF'].isin(tf), 'Type'] = 'TF'
        self.data = data
        if return_df:
            return data

    def save_network(self, save_fp):
        if self.data is None:
            data = self.create_network()
        else:
            data = self.data
        data.to_csv(save_fp, index=False, sep='\t')

    # 根据type的值，将dataframe分成多个子集,返回的结果为一个字典，键为type的值，值为一个dataframe
    def split_network(self):
        if self.data is None:
            self.create_network()
        data = self.data
        self.data_by_type = {}
        for _type in data['Type'].unique():
            self.data_by_type[_type] = data[data['Type'] == _type]
        return self.data_by_type

    # 从self.interaction_data中获取所有的TF和Target的名字
    # def get_all_genes(self):
    #     all_genes = set(self.interaction_data.interactions.keys())
    #     all_genes.update(set(self.interaction_data.interactions.values()))
    #     return all_genes
    def get_all_genes(self):
        all_genes = set(self.data['TF'].tolist() + self.data['Target'].tolist())
        return all_genes

'''
下面是输出的结果文件的样式。
TF	Target	Type
AHR	CDKN1B	TF
APLNR	PIK3C3	Gene
APLNR	PIK3R4	Gene
AR	KLK3	TF
ARNT	ALDOA	TF
ARNT	ANGPT1	TF
ARNT	ANGPT2	TF

'''

if __name__ == '__main__':
    base_dir = r"/home/liuyq/data/ar_data/interaction_db/regnetworks"
    interaction_fp = r'new_kegg.human.reg.direction.txt'
    tf_fp = r'TF_names_v_1.01.txt'
    save_fp = r'dataframe.tsv'

    interaction_fp = os.path.join(base_dir, interaction_fp)
    tf_fp = os.path.join(base_dir, tf_fp)
    save_fp = os.path.join(base_dir, save_fp)

    network = RegulationNetwork(tf_fp, interaction_fp)
    network.create_network()
    network.save_network(save_fp)
    sp_data = network.split_network()
    network.get_all_genes()
    print("finished!")





