from mm_gbsa.mm_gbsa_utils.dump import pickle_mmgbsa_data
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calc DD-CSA ccr func')
    parser.add_argument('--mmgbsa-info_file', required=True)
    parser.add_argument('--output-directory', default=".")
    parser.add_argument('--output-filename', default="MMPBSA_data.pickle")
    args = parser.parse_args()

    pickle_mmgbsa_data(mmgbsa_info_file=args.mmgbsa_info_file,
                       output_dir=args.output_directory,
                       output_filename=args.output_filename)
