'''
读取特定文件夹下的数据，该注塑机包括表达数据和细胞类型注释数据
'''
import config
from utils.data_process.read_cell_type import read_cell_type
from utils.data_process.load_data import LoadMatrixDataReal
from utils.data_process.load_data import TransformDataReal
from utils.algorithms.matrix_association_rule import MatrixRule
import os
from utils.data_process.ar_metrics_process import SaveArMetrics, LoadArMetrics, FilterArMetrics
from utils.data_process.read_gene_list import read_gene_list
from utils.pipe_analysis.subnet_detection import Graph
from utils.pipe_analysis.gsea import GSEA
from utils.data_process.merge_network import MergeNetwork
from utils.data_process.gene_selection import GeneSelection
from utils.pipe_analysis.deg import DiffereceGene
from config import Config
from utils.data_process.ar_metrics_process import SaveArMetrics, LoadArMetrics
from utils.data_process.pickle_unpicle import pickle_data, unpickle_data
class CellTypePipe:
    def __init__(self, data_path,
                 ar_result_path,
                 cell_type_fp="cell_type.tsv",
                 threshold_dict=None):
        if threshold_dict is None:
            threshold_dict = {"support": [0.1, None], "confidence": [0.5, None], "lift": [0.5, None],
                              "leverage": [0.5, None], "conviction": [0.5, None]}
        self.data_path = data_path
        self.ar_result_path = ar_result_path
        self.cell_type_path = cell_type_fp
        self.threshold_dict = threshold_dict
        self.min_support = threshold_dict["support"][0]

    def reset_threshold(self, threshold_dict):
        self.threshold_dict = threshold_dict
        self.min_support = threshold_dict["support"][0]

    def read_data(self, min_support=None):
        # 加载数据
        data_obj = LoadMatrixDataReal(self.data_path)
        _data = data_obj.load_data()
        _data = data_obj.filter_genes(_data, threshold=min_support)

        # 读取cell_type.tsv文件
        cell_type = read_cell_type(self.cell_type_path)
        adata_subs = data_obj.split_by_cell_type(_data, cell_type)
        return adata_subs, data_obj

    def get_degs(self, ):
        diff_gene = DiffereceGene(data_path=self.data_path, cell_type_fn=self.cell_type_path)
        adata = diff_gene.load_data()
        adata = diff_gene.assign_cell_type(adata)
        adata = diff_gene.preprocess(adata)
        diff_gene.get_diff_genes(adata)
        degs = diff_gene.deg_dict
        diff_gene.save_deg_df_dict()
        return degs

    def bio_network(self, bionet_fp_dict=None):
        ppi_fp = bionet_fp_dict["ppi_fp"]
        ppi_mapping_fp = bionet_fp_dict["ppi_mapping_fp"]
        gmt_fp = bionet_fp_dict["gmt_fp"]
        interaction_fp = bionet_fp_dict["reactome_fp"]
        tf_fp = bionet_fp_dict["tf_fp"]
        regnetwork_fp = bionet_fp_dict["regnetwork_fp"]
        merge_network = MergeNetwork(
            ppi_fp=ppi_fp,
            ppi_mapping_fp=ppi_mapping_fp,
            gmt_fp=gmt_fp,
            reactome_fp=interaction_fp,
            tf_fp=tf_fp,
            regnetwork_fp=regnetwork_fp
        )
        return merge_network

    def subset_by_genes(self, adata_subs, bionet_fp_dict, deg_sel=True, network_sel=True):
        for cell_type, adata in adata_subs.items():
            gene_selector = GeneSelection(adata)
            if deg_sel:
                # get seed genes from degs
                degs = self.get_degs()
                seed_genes = degs[cell_type]
                gene_selector.set_seed_genes(seed_genes)
                target_genes = seed_genes
            if network_sel:
                if not deg_sel:
                    raise ValueError("Please set deg_sel=True, before setting network_sel=True")
                # get seed genes from bio_network
                merge_network = self.bio_network(bionet_fp_dict)
                bio_network_genes_dict = merge_network.merge_network()
                gene_selector.set_bio_network_genes(bio_network_genes_dict)
                candidate_genes = gene_selector.get_candidate_genes()
                target_genes = candidate_genes
            target_genes = adata.var_names.intersection(target_genes)
            adata = adata[:, target_genes]
            adata_subs[cell_type] = adata
        return adata_subs

    def run_ar(self, adata_subs):
        results = {}
        for cell_type, adata in adata_subs.items():
            transform_data_obj = TransformDataReal(adata)
            data, genes_info, cells_info = transform_data_obj.run()
            mar_obj = MatrixRule(data, genes_info, cells_info)
            metrics = mar_obj.all_metrics_to_dataframe
            results[cell_type] = metrics
        return results

    def filter_result(self, results):
        results_filter = FilterArMetrics(results, threshold_dict=self.threshold_dict)
        res1 = results_filter.filter_result()
        return res1

    def transform_to_edge_style(self, results):
        results_filter = FilterArMetrics(results, threshold_dict=self.threshold_dict)
        new_df = results_filter.transform_to_pairs_in_df(results)
        return new_df

    def filter_transform_to_edge_style(self, results):
        results_filter = FilterArMetrics(results, threshold_dict=self.threshold_dict)
        new_df = results_filter.filter_pairs(results=results, threshold_dict_tem=self.threshold_dict)
        return new_df

    # 这里的输入为filter_transform_to_edge_style 或 transform_to_edge_style
    def export_edge_style(self, results, out_dir=None):
        results_filter = FilterArMetrics(results, threshold_dict=self.threshold_dict)
        new_df = results_filter.export_edge_list(results, out_dir=out_dir)
        return new_df

    def extract_by_genes(self, results, gene_list):
        results_filter = FilterArMetrics(results, threshold_dict=self.threshold_dict)
        res1 = results_filter.extract_by_genes(gene_list, results)
        return res1

    def extract_by_gene_pairs(self, results, gene_pair):
        results_filter = FilterArMetrics(results, threshold_dict=self.threshold_dict)
        res1 = results_filter.extract_by_gene_pairs(gene_pair, results)
        return res1

    def save_result(self, results):
        save_obj = SaveArMetrics(results, self.ar_result_path)
        save_obj.write_result()
        print("finish writing result...")

    def load_result(self):
        armetrics = LoadArMetrics(self.ar_result_path)
        results = armetrics.load_result()
        return results

    # 这里的输入为边列表的样式，
    def get_all_genes(self, results):
        gene_dict = {}
        for cell_type, df in results.items():
            gene_list = []
            gene_list.extend(df["antecedent"].tolist())
            gene_list.extend(df["consequent"].tolist())
            gene_list = list(set(gene_list))
            gene_dict[cell_type] = gene_list
        return gene_dict

    # 这里的输入为边列表的样式，export_edge_style的输出
    def get_all_graph(self, results):
        graph_dict = {}
        for cell_type, df in results.items():
            graph = Graph(df)
            graph_dict[cell_type] = graph
        return graph_dict

    def community_detection(self, G):
        communities = G.community_detection()
        return communities

    # communities content: {node: community_id}, community_to_trans_dict convert to : {community_id: [node1, node2, ...]}
    def community_to_trans_dict(self, communities):
        community_trans_dict = {}
        for node, community_id in communities.items():
            if community_id not in community_trans_dict:
                community_trans_dict[community_id] = {node}
            else:
                community_trans_dict[community_id].add(node)
        return community_trans_dict

    def gesa(self,
             gene_list,
             gsea_threshold_dict,
             gsea_meta,
             gmt_type = "reactome_gmt_fp",
             save_fp = None,
            ):
        gmt_fp = gsea_meta[gmt_type]
        background = gsea_meta["background"]
        out_dir = gsea_meta["outdir"]
        gsea_obj = GSEA(gmt_fp, background=background, outdir=out_dir)
        results = gsea_obj.enrich(gene_list, save_fp=None)
        results_filterd = gsea_obj.filter_results(results, gsea_threshold_dict)
        gsea_obj.save_results(results_filterd, save_fp)
        return results_filterd

    def graph_and_gsea_steps(self, results, meta_dict):
        if not os.path.exists(Config.gene_dict_each_cell_fp):
            gene_dict = self.get_all_genes(results)
            pickle_data(gene_dict, Config.gene_dict_each_cell_fp)
        else:
            gene_dict = unpickle_data(Config.gene_dict_each_cell_fp)
        if not os.path.exists(Config.graph_dict_each_cell_fp):
            graph_dict = self.get_all_graph(results)
            pickle_data(graph_dict, Config.graph_dict_each_cell_fp)
        else:
            graph_dict = unpickle_data(Config.graph_dict_each_cell_fp)
        if not os.path.exists(Config.cache_community_detection_dir):
            os.makedirs(Config.cache_community_detection_dir)
        else:
            print("dir exists")
        for cell_type, gene_list in gene_dict.items():
            G = graph_dict[cell_type]
            communities_fp = os.path.join(Config.cache_community_detection_dir, cell_type + "_communities.obj")
            if os.path.exists(communities_fp):
                communities = unpickle_data(communities_fp)
            else:
                communities = self.community_detection(G)
                pickle_data(communities, communities_fp)

            communities_trans_fp = os.path.join(Config.cache_community_detection_dir, cell_type + "_communities_trans.obj")
            if os.path.exists(communities_trans_fp):
                community_trans_dict = unpickle_data(communities_trans_fp)
            else:
                community_trans_dict = self.community_to_trans_dict(communities)
                pickle_data(community_trans_dict, communities_trans_fp)
            # save graph to edge list type file
            save_graph_fp = os.path.join(meta_dict["graph_outdir"], cell_type + ".tsv")
            if not os.path.exists(save_graph_fp):
                G.save_graph(G=G.G, out_fp=save_graph_fp)
            save_graph_summary_fp = os.path.join(meta_dict["graph_outdir"], cell_type + "_summary.tsv")
            if not os.path.exists(save_graph_summary_fp):
                G.summary(G=G.G, out_fp=save_graph_summary_fp)
            save_fp = os.path.join(meta_dict["gsea_outdir"], cell_type + "reactome_enrich.tsv")

            gsea_results = self.gesa(gene_list=gene_list,
                                     gsea_threshold_dict=meta_dict["gsea_threshold_dict"],
                                     gsea_meta=meta_dict["gsea_meta"],
                                     gmt_type="reactome_gmt_fp",
                                     save_fp=save_fp)
            save_fp = os.path.join(meta_dict["gsea_outdir"], cell_type + "go_enrich.tsv")
            gsea_results = self.gesa(gene_list=gene_list,
                                        gsea_threshold_dict=meta_dict["gsea_threshold_dict"],
                                        gsea_meta=meta_dict["gsea_meta"],
                                        gmt_type="go_gmt_fp",
                                        save_fp=save_fp)
            # run gsea for each community
            for community_id, nodes in community_trans_dict.items():
                community_gene_list = list(set(gene_list).intersection(nodes))
                # reactome enrichment
                save_fp = os.path.join(meta_dict["gsea_outdir"],
                                       cell_type + "_gsea_reactome_community_id_" + str(community_id) + ".tsv")
                gsea_results = self.gesa(gene_list=community_gene_list,
                                         gsea_threshold_dict=meta_dict["gsea_threshold_dict"],
                                         gsea_meta=meta_dict["gsea_meta"],
                                         gmt_type="reactome_gmt_fp",
                                         save_fp=save_fp)
                # go enrichment
                save_fp = os.path.join(meta_dict["gsea_outdir"],
                                       cell_type + "_gsea_go_community_id_" + str(community_id) + ".tsv")
                gsea_results = self.gesa(gene_list=community_gene_list,
                                         gsea_threshold_dict=meta_dict["gsea_threshold_dict"],
                                         gsea_meta=meta_dict["gsea_meta"],
                                         gmt_type="go_gmt_fp",
                                         save_fp=save_fp)

            # save community to edge list type file
            subgraph_dict = G.build_subgraph(communities)
            for community_id, subgraph in subgraph_dict.items():
                save_sub_grahp_fp = os.path.join(meta_dict["graph_outdir"],
                                                 cell_type + "_community_id_" + str(community_id) + ".tsv")
                G.save_graph(G=subgraph, out_fp=save_sub_grahp_fp)

                save_sub_grahp_summary_fp = os.path.join(meta_dict["graph_outdir"],
                                                         cell_type + "_community_id_" + str(community_id) + "_summary.tsv")
                G.summary(G=subgraph, out_fp=save_sub_grahp_summary_fp)

        print(f"cell_type: {cell_type} finished...")

if __name__ == '__main__':
    import os
    from config import Config
    from utils.pipe_analysis.cell_type_pipe import CellTypePipe

    thread_dict = Config.threshold_dict
    cell_type_pipe = CellTypePipe(data_path=Config.data_fp,
                                  ar_result_path=Config.ar_fp,
                                  cell_type_fp=Config.cell_type_fp,
                                  threshold_dict=Config.threshold_dict)

    if os.path.exists(Config.cahe_result_edge_style_fp) and os.path.exists(Config.cahe_result_edge_style_filted_fp):
        print("load from cache...")
        result_edge_style = LoadArMetrics(Config.cahe_result_edge_style_fp).load_result()
        result_edge_style_filted = LoadArMetrics(Config.cahe_result_edge_style_filted_fp).load_result()
    else:

        # import data
        adata_subs, data_obj = cell_type_pipe.read_data(min_support=Config.min_support)

        # select genes
        bio_net_fp_dict = {"ppi_fp": Config.ppi_fp,
                           "ppi_mapping_fp": Config.ppi_mapping_fp,
                           "gmt_fp": Config.gmt_fp,
                           "reactome_fp": Config.reactome_fp,
                           "tf_fp": Config.tf_fp,
                           "regnetwork_fp": Config.regnetwork_fp}
        adata_subs_gene_selected = cell_type_pipe.subset_by_genes(adata_subs, bio_net_fp_dict, deg_sel=True,
                                                                  network_sel=True)

        # run ar
        results = cell_type_pipe.run_ar(adata_subs=adata_subs_gene_selected)

        #filter results
        #results_filted = cell_type_pipe.filter_result(results)

        result_edge_style = cell_type_pipe.transform_to_edge_style(results)
        res_fp = os.path.join(config.Config.ar_fp, "result_edge_style")
        save_obj = SaveArMetrics(results=result_edge_style, out_path=res_fp)
        save_obj.write_result()
        print("finish writing result_edge_style")

        result_edge_style_filted = cell_type_pipe.filter_transform_to_edge_style(result_edge_style)
        res_filtered_fp = os.path.join(config.Config.ar_fp, "result_edge_style_filted")
        save_obj = SaveArMetrics(results=result_edge_style_filted, out_path=res_filtered_fp)
        save_obj.write_result()
        print("finish writing result_edge_style_filted")

        #export edge style results to file
        output_path = Config.ar_fp
        unfilterd_output_path = os.path.join(output_path, "edge_style","unfilterd")
        try:
            os.makedirs(unfilterd_output_path)
        except FileExistsError:
            print("FileExists, edgle style output will be overwrited...")
        except Exception as e:
            print(F"Make dir error: {e}")

        cell_type_pipe.export_edge_style(result_edge_style, unfilterd_output_path)


        filterd_output_path = os.path.join(output_path, "edge_style","filterd")
        try:
            os.makedirs(filterd_output_path)
        except FileExistsError:
            print("FileExists, filtered edgle style output will be overwrited...")
        except Exception as e:
            print(F"Make dir error: {e}")

        cell_type_pipe.export_edge_style(result_edge_style_filted, filterd_output_path)

        #cache result
        cache_result_fp = os.path.join(Config.cahe_result_edge_style_fp)
        save_obj = SaveArMetrics(results=result_edge_style, out_path=cache_result_fp)
        save_obj.write_result()
        print("result_edge_style has been cached...")

        cache_result_filted_fp = os.path.join(Config.cahe_result_edge_style_filted_fp)
        save_filterd_obj = SaveArMetrics(results=result_edge_style_filted, out_path=cache_result_filted_fp)
        save_filterd_obj.write_result()
        print("result_edge_style_filted has been cached...")


    #get all result_edge_style_filted genes and run gsea
    #gene_dict = cell_type_pipe.get_all_genes(result_edge_style_filted)


    #community_detection
    graph_meta_dict = {
        "graph_outdir": Config.graph_fp,
        "gsea_outdir": Config.gsea_fp,
        "gsea_meta": Config.gesea_meta,
        "gsea_threshold_dict":Config.gsea_threshold_dict
    }

    cell_type_pipe.graph_and_gsea_steps(results = result_edge_style_filted,
                                        meta_dict = graph_meta_dict)
    print("all task finish...")



























