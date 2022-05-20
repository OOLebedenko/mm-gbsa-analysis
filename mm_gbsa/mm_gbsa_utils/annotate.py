from pyxmolpp2 import PdbFile, mName, rId
import pandas as pd
import os


def annotate_ddg_decomposition(path_to_decomposition_csv,
                               path_to_pdb_reference,
                               partner_A="partner_A",
                               partner_B="partner_B",
                               output_directory=".",
                               output_filename="decomposition_annotated.csv"
                               ):
    ### read reference
    reference = PdbFile(path_to_pdb_reference).frames()[0]
    molecule_partner_A = reference.molecules.filter(mName == "A")
    molecule_partner_B = reference.molecules.filter(mName == "B")
    rids_partner_A = [residue.id.serial for residue in molecule_partner_A.residues]
    rids_partner_B = [residue.id.serial for residue in molecule_partner_B.residues]

    ### set first rid in case end-to-end numbering (relevant for mmgbsa analysis)
    partner_A_start_rid = rids_partner_A[0]
    partner_B_start_rid = rids_partner_A[-1] + 1

    ### read decopmosition data
    df_decomp = pd.read_csv(path_to_decomposition_csv)
    interacting_residues = df_decomp["rId_1-rId_2"]

    ### get upper triangular matrix
    mask_pairs = [True if int(residue_pairs.split("-")[0]) <= partner_B_start_rid else False for residue_pairs in
                  interacting_residues]
    interacting_residues = interacting_residues[mask_pairs]

    ### create annotated decomposition data
    df_decomp_annotated = pd.DataFrame()

    for residue_pair in interacting_residues:
        rid_first, rid_second = residue_pair.split("-")
        rid_partner_A = rids_partner_A[int(rid_first) - partner_A_start_rid]
        rid_partner_B = rids_partner_B[int(rid_second) - partner_B_start_rid]
        rname_partner_A = molecule_partner_A.residues.filter(rId == rid_partner_A)[0].name
        rname_partner_B = molecule_partner_B.residues.filter(rId == rid_partner_B)[0].name

        temp = {f"rId_{partner_A}": rid_partner_A, f"rName_{partner_A}": rname_partner_A,
                f"rId_{partner_B}": rid_partner_B, f"rName_rbd_{partner_B}": rname_partner_B,
                "total_energy_avg": df_decomp["total_energy_avg"][df_decomp["rId_1-rId_2"] == residue_pair].values[0],
                "total_energy_std": df_decomp["total_energy_std"][df_decomp["rId_1-rId_2"] == residue_pair].values[0]
                }
        df_decomp_annotated = pd.concat([df_decomp_annotated, pd.DataFrame(temp, index=[0])])

    os.makedirs(output_directory, exist_ok=True)
    df_decomp_annotated.to_csv(os.path.join(output_directory, output_filename), index=False)
