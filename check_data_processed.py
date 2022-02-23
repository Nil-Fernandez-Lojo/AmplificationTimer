from pathlib import Path
import pandas as pd

path_data = Path("../memoire/data")
path_summary_table = path_data / "summary_table_combined_annotations_v4.txt" 
summary_table = pd.read_csv(path_summary_table, sep='\t')
N_samples = summary_table.shape[0]

for i in range(N_samples):
	name = summary_table['samplename'][i]
	path_file = Path('preprocessed_data')/(name+'.json')
	if not path_file.exists():
		print(i,name)