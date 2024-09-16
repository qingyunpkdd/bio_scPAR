'''
cell_type.tsv 是一个细胞类型注释文件，只有一列，和barcodes.tsv中的barcode一一对应，每个barcode对应一个细胞类型。

'''
def read_cell_type(file_path):
    with open(file_path, 'r') as f:
        cell_type = f.readlines()
        #去掉每行末尾的换行符
        cell_type = [line.strip() for line in cell_type]
    return cell_type

