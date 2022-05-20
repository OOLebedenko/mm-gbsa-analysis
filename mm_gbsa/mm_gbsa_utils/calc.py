from mm_gbsa.mm_gbsa_utils.dump import load_pickle_data
from tqdm import tqdm
import pandas as pd
import numpy as np
import os


def calc_ddg_decomposition(pickle_mmgbsa_data,
                           first_frame=1,
                           output_directory=".",
                           output_filename="decomposition.csv",
                           energy_threshold=1):
    decomp_data = load_pickle_data(pcike_file=pickle_mmgbsa_data)

    ### read total energy
    receptor_residue_pairs = np.array([*decomp_data['decomp']["gb"]['receptor']["TDC"]])
    ligand_residue_pairs = np.array([*decomp_data['decomp']["gb"]['ligand']["TDC"]])
    complex_residue_pairs = np.array([*decomp_data['decomp']["gb"]['complex']["TDC"]])


    ### get interacting residues lcb+rbd
    mask_intra_residues = np.isin(complex_residue_pairs, receptor_residue_pairs) | np.isin(
        complex_residue_pairs, ligand_residue_pairs)
    complex_residue_pairs = complex_residue_pairs[~(mask_intra_residues)]

    out_decomp_df = pd.DataFrame()
    for residue_pair in tqdm(complex_residue_pairs):
        rid_first, rid_second = residue_pair.split("-")

        ### decomp_matrix is pseido-symmetry. sum of the per-residue contributions will equal the total binding free energy
        ### Gohlke, H., C. Kiel, and D. A. Case. 2003 https://pubmed.ncbi.nlm.nih.gov/12850155/
        total_energy = decomp_data['decomp']["gb"]['complex']["TDC"][f"{rid_first}-{rid_second}"]['tot'].copy() + \
                       decomp_data['decomp']["gb"]['complex']["TDC"][f"{rid_second}-{rid_first}"]['tot'].copy()

        total_energy_avg, total_energy_std = total_energy[first_frame:].mean(), total_energy[first_frame:].std()

        ### retaining only those pairs, whose trajectory-average energy exceeded 1 kcal/mol
        if abs(total_energy_avg) > energy_threshold:
            temp = {"rId_1-rId_2": residue_pair,
                    "total_energy_avg": total_energy_avg, "total_energy_std": total_energy_std}
            out_decomp_df = pd.concat([out_decomp_df, pd.DataFrame(temp, index=[0])])

    os.makedirs(output_directory, exist_ok=True)
    out_decomp_df.to_csv(os.path.join(output_directory, output_filename), index=False)


if __name__ == "__main__":
    pickle_mmgbsa_data = "/home/olebedenko/cov2/handling/mm_gbsa/lcb3/igb8_salt150_decomp/MMPBSA_data_wt.pickle"

    calc_ddg_decomposition(pickle_mmgbsa_data)
