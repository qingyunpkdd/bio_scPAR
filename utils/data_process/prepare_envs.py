import os
from config import Config
from utils.data_process.pickle_unpicle import pickle_data, unpickle_data
import shutil

class DirManager:
    def __init__(self):
        self.cache_dirs = Config.intermediate_dir_fp
        self.results_dirs = {
            "ar": Config.ar_fp,
            "gsea": Config.gsea_fp,
            "graph": Config.graph_fp,
            "tensor": Config.tensor_fp,
            "deg": Config.deg_fp
        }

    def mk_cache_dir(self):
        if not os.path.exists(Config.cache_dir):
            os.makedirs(Config.cache_dir)
    #remve self.cache_dirs and create new cache dirs
    def clean_cache_dir(self):
        if os.path.exists(self.cache_dirs):
            print("remove cache dirs")
            shutil.rmtree(self.cache_dirs)
        print("create cache dirs")
        os.makedirs(self.cache_dirs)

    def mk_results_dir(self):
        for key, value in self.results_dirs.items():
            if not os.path.exists(value):
                os.makedirs(value)

    def clean_results_dir(self):
        for key, value in self.results_dirs.items():
            if os.path.exists(value):
                print(f"remove {key} dir")
                shutil.rmtree(value)
        print("create results dirs")
        self.mk_results_dir()


if __name__ == '__main__':
    dir_manager = DirManager()
    dir_manager.clean_cache_dir()
    dir_manager.mk_results_dir()
    print("done")





















