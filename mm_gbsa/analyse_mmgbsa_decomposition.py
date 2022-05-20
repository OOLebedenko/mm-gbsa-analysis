from mm_gbsa.mm_gbsa_utils.calc import calc_ddg_decomposition
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calc DD-CSA ccr func')
    parser.add_argument('--pickle-mmgbsa-data', required=True)
    parser.add_argument('--first-frame', default=0, type=int)
    parser.add_argument('--energy-threshold', default=1, type=int)
    parser.add_argument('--output-directory', default=".")
    parser.add_argument('--output-filename', default="decomposition.csv")
    args = parser.parse_args()

    calc_ddg_decomposition(pickle_mmgbsa_data=args.pickle_mmgbsa_data,
                           output_directory=args.output_directory,
                           output_filename=args.output_filename,
                           first_frame=args.first_frame,
                           energy_threshold=args.energy_threshold
                           )
