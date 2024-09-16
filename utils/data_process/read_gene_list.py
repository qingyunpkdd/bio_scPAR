def read_gene_list(file_path):
    with open(file_path, 'r') as f:
        gene_list = f.readlines()
        gene_list = [line.strip() for line in gene_list]
    return gene_list