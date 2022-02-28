import pickle
import pymc3 as pm
from pathlib import Path
import matplotlib.pyplot as plt

path_folder = Path("inferred_times")
sample_name = "02e5c36f-5bec-45e2-a048-875653b85ca1"
amplification_index = 0
subclonality_modelled = True 

if subclonality_modelled:
	path_folder = path_folder / "subclonality_modelled"
else:
	path_folder = path_folder / "subclonality_not_modelled"

for model_idx in range(1,5):
	path_file = path_folder / ("Model"+str(model_idx)) / (sample_name + "_amp_"+str(amplification_index)+"_trace.pkl")
	trace = pickle.load(open( path_file, "rb" ))
	pm.plots.plot_posterior(trace)
	path_figure = path_folder / ("Model"+str(model_idx)) / (sample_name + "_amp_"+str(amplification_index)+"_plot_posterior.png")
	plt.savefig(path_figure, bbox_inches='tight')
