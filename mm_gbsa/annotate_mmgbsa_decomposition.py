from mm_gbsa.mm_gbsa_utils.annotate import annotate_ddg_decomposition
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calc DD-CSA ccr func')
    parser.add_argument('--path-to-decomposition-csv', required=True)
    parser.add_argument('--path-to-pdb-reference', required=True)
    parser.add_argument('--output-directory', default=".")
    parser.add_argument('--output-filename', default="decomposition_annotated.csv")
    args = parser.parse_args()

    annotate_ddg_decomposition(path_to_decomposition_csv=args.path_to_decomposition_csv,
                               path_to_pdb_reference=args.path_to_pdb_reference,
                               output_directory=args.output_directory,
                               output_filename=args.output_filename,
                               partner_A="lcb",
                               partner_B="rbd"
                               )
