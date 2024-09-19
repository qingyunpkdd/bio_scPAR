import os

from config import Config
from utils.pipe_analysis.cell_type_pipe import CellTypePipe

thread_dict = Config.threshold_dict
cell_type_pipe = CellTypePipe(data_path=Config.data_fp,
                              ar_result_path=Config.ar_fp,
                              cell_type_fp=Config.cell_type_fp,
                              threshold_dict=Config.threshold_dict)


#import data
adata_subs, data_obj = cell_type_pipe.read_data(min_support=Config.min_support)

#select genes
bio_net_fp_dict = {"ppi_fp":Config.ppi_fp,
                  "ppi_mapping_fp":Config.ppi_mapping_fp,
                  "gmt_fp":Config.gmt_fp,
                  "reactome_fp":Config.tf_fp,
                  "tf_fp":Config.tf_fp,
                  "regnetwork_fp":Config.regnetwork_fp}
# filter genes by network
adata_subs_gene_selected = cell_type_pipe.subset_by_genes(adata_subs, bio_net_fp_dict, deg_sel=True, network_sel=True)

#run ar
results = cell_type_pipe.run_ar(adata_subs = adata_subs_gene_selected)
#filter results
results_filted = cell_type_pipe.filter_ar_results(results)

result_edge_style = cell_type_pipe.transform_to_edge_style(results_filted)
result_edge_style_filted = cell_type_pipe.filter_edge_style(result_edge_style)

#export edge style results to file
output_path = Config.ar_fp
unfilterd_output_path = os.path.join(output_path, "edge_style","unfilterd")
try:
    os.makedirs(unfilterd_output_path)
except FileExistsError:
    print("FileExistsError")
except Exception as e:
    print(F"Make dir error: {e}")

cell_type_pipe.export_edge_style(result_edge_style, unfilterd_output_path)


filterd_output_path = os.path.join(output_path, "edge_style","filterd")
try:
    os.makedirs(filterd_output_path)
except FileExistsError:
    print("FileExistsError")
except Exception as e:
    print(F"Make dir error: {e}")

cell_type_pipe.export_edge_style(result_edge_style_filted, filterd_output_path)


#get all result_edge_style_filted genes and run gsea
gene_dict = cell_type_pipe.get_all_genes(result_edge_style_filted)


#community_detection
graph_meta_dict = {
    "graph_outdir": Config.graph_fp,
    "gsea_outdir": Config.gsea_fp,
    "gsea_meta": Config.gesea_meta,
    "gsea_threshold_dict":Config.gsea_threshold_dict
}

cell_type_pipe.graph_and_gsea_steps(results = result_edge_style_filted,
                                    meta_dict = graph_meta_dict)


