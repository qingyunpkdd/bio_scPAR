import json
import pickle
import pandas as pd
import os
def pickle_data(data, path):
    if not os.path.exists(path):
        if isinstance(data, pd.DataFrame):
            save_pd(data, path)
        else:
            pickle_python_object(data, path)
    else:
        print("file:", path, "already exists")


def pickle_python_object(obj, path):
    with open(path, 'wb') as f:
        pickle.dump(obj, f)
    print("finish pickle object to", path)

def save_pd(df, path):
    if path.endswith(".df"):
        df.to_pickle(path)
    else:
        raise ValueError("file path should end with .df")



def unpickle_data(path):
    if os.path.exists(path):
        print("start load data from:", path)
        if path.endswith(".df"):
            return pd.read_pickle(path)
        else:
            return unpickle_python_object(path)
    else:
        print("file:", path, "does not exists")


def unpickle_python_object(path):
    with open(path, 'rb') as f:
        obj = pickle.load(f)
    return obj


if __name__ == '__main__':
    data = {"a": 1, "b": 2, "c": 3}
    path = "/home/liuyq/data/ar_data/test.obj"
    pickle_data(data=data, path=path)
    test = unpickle_data(path)
    print(test)

    data_pd = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    data_pd.index = ["a", "b", "c"]
    path = "/home/liuyq/data/ar_data/test.df"
    pickle_data(data=data_pd, path=path)
    test = unpickle_data(path)
    print(test)
