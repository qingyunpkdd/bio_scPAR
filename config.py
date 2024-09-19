import os
import re



class Config:
    # 输入数据文件路径
    #############################################################################################
    # 数据文件路径，所有的数据文件都放在这个目录下
    base_data_dir = r"/home/liuyq/data/ar_data"
    #生物网络数据文件路径
    #######################
    bionetwork_med_dir ="interaction_db"
    #reactome数据文件路径
    reactome_dir = "reactome"
    gmt_fp = os.path.join(base_data_dir, bionetwork_med_dir, reactome_dir, "ReactomePathways.gmt")
    reactome_fp = os.path.join(base_data_dir, bionetwork_med_dir, reactome_dir, "FIsInGene_070323_with_annotations.txt")

    #ppi数据文件路径
    ppi_dir = "ppi"
    ppi_fp = os.path.join(base_data_dir, bionetwork_med_dir, ppi_dir, "PP-Decagon_ppi.csv")
    ppi_mapping_fp = os.path.join(base_data_dir, bionetwork_med_dir, ppi_dir, "G-SynMiner_miner-geneHUGO.tsv")

    #regnetwork数据文件路径
    regnetwork_dir = "regnetworks"
    tf_fp = os.path.join(base_data_dir, bionetwork_med_dir, regnetwork_dir, "TF_names_v_1.01.txt")
    regnetwork_fp = os.path.join(base_data_dir, bionetwork_med_dir, regnetwork_dir, "new_kegg.human.reg.direction.txt")

    #数据文件路径
    #######################
    data_med_dir = "sc_data"
    #数据文件路径
    data_dir_name = "pbmc"
    data_fp = os.path.join(base_data_dir, data_med_dir, data_dir_name)
    cell_type_fp = os.path.join(base_data_dir, data_med_dir, data_dir_name, "cell_type.tsv")

    #富集分析文件路径
    #######################
    gsea_med_dir = "gsea"
    go_gmt_fp = os.path.join(base_data_dir, gsea_med_dir, "go.gmt")
    reactome_gmt_fp = os.path.join(base_data_dir, gsea_med_dir, "reactome.gmt")

    #############################################################################################

    # Intermediate result storage path
    #############################################################################################
    # seed genes to cadiates genes mapping file
    intermediate_dir = "intermediate"
    intermediate_dir_fp = os.path.join(base_data_dir, intermediate_dir)
    seed2cand_fp = os.path.join(base_data_dir, intermediate_dir, "seed2cand.obj")
    ppi_interaction_df_fp = os.path.join(base_data_dir, intermediate_dir, "ppi_interaction_df.df")
    regnetwork_interaction_df_fp = os.path.join(base_data_dir, intermediate_dir, "regnetwork_interaction_df.df")
    reactome_interaction_df_fp = os.path.join(base_data_dir, intermediate_dir, "reactome_interaction_df.df")

    # graph Cache
    gene_dict_each_cell_fp = os.path.join(base_data_dir, intermediate_dir, "gene_dict_each_cell.obj")
    graph_dict_each_cell_fp = os.path.join(base_data_dir, intermediate_dir, "graph_dict_each_cell.obj")

    cache_community_detection_dir = os.path.join(base_data_dir, intermediate_dir, "community_detection")

    cahe_result_edge_style_fp = os.path.join(base_data_dir, intermediate_dir, "result_edge_style.obj")
    cahe_result_edge_style_filted_fp = os.path.join(base_data_dir, intermediate_dir, "result_edge_style_filted.obj")

    #############################################################################################

    # 输出以及中间文件路径
    #############################################################################################
    # 输出文件都放在这个目录下
    base_output_dir = r"/home/liuyq/data/ar_result"
    bionetwork_dir = "bionetwork"

    tensor_dir = "ar_tensor"
    ar_dir = "ar"
    deg_dir = "deg"
    gsea_dir = "gsea"
    graph_dir = "graph"

    tensor_fp = os.path.join(base_output_dir, tensor_dir)
    ar_fp = os.path.join(base_output_dir, ar_dir)
    deg_fp = os.path.join(base_output_dir, deg_dir)
    gsea_fp = os.path.join(base_output_dir, gsea_dir)
    graph_fp = os.path.join(base_output_dir, graph_dir)

    #############################################################################################



    # 阈值设置
    #############################################################################################
    min_support = 0.1
    threshold_dict = {"support":[0.05, None], "confidence":[0.5, None], "lift":[1, None], "leverage":[None, None], "conviction":[1, None]}
    gsea_threshold_dict = {"Adjusted P-value": 0.01}
    gesea_meta = {
    "go_gmt_fp":go_gmt_fp,
    "reactome_gmt_fp":reactome_gmt_fp,
    "background": 30000,
    "outdir": gsea_fp,
    }



    # 阈值设置--差异表达基因
    #############################################################################################
    deg_threshold_dict = {"log2fc":0.25, "p_value":0.05, "pval_adj":0.05}

    #############################################################################################





















