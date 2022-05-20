from MMPBSA_mods import API as MMPBSA_API
import pickle
import os


def pickle_mmgbsa_data(mmgbsa_info_file, output_dir=".", output_filename="MMPBSA_data.pickle"):
    wd = os.getcwd()
    print()
    print(mmgbsa_info_file)
    print()
    os.chdir(os.path.dirname(mmgbsa_info_file))
    data = MMPBSA_API.load_mmpbsa_info("_MMPBSA_info")
    os.chdir(wd)
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, output_filename), 'wb') as f:
        pickle.dump(data, f)


def load_pickle_data(pcike_file="MMPBSA_data.pickle"):
    with open(pcike_file, 'rb') as f:
        data = pickle.load(f)
    return data
