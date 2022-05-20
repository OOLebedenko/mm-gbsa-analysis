include ../../common.mk

all: data/decomposition_annotated.csv

clean:
	rm -rf data

data/decomposition_annotated.csv: ${SCRIPT_DIR}/annotate_mmgbsa_decomposition.py data/decomposition.csv
	python3 "$<" \
	    --path-to-decomposition-csv=data/decomposition.csv \
	    --path-to-pdb-reference="${REFERENCE_PDB_PATH}" \
	    --output-directory="data" \

data/decomposition.csv: ${SCRIPT_DIR}/analyse_mmgbsa_decomposition.py data/MMPBSA_data.pickle
	python3 "$<" \
	    --pickle-mmgbsa-data=data/MMPBSA_data.pickle \
	    --first-frame="${FIRST_FRAME_FOR_MMGBSA}" \
	    --energy-threshold="${ENERGY_THRESHOLD}" \
	    --output-directory="data" \

data/MMPBSA_data.pickle: ${SCRIPT_DIR}/pickle_mmgbsa.py
	python3 "$<" \
	    --mmgbsa-info_file="${PATH_TO_MMGBSA_INFO}" \
	    --output-directory="data" 
