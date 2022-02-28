import pymc3 as pm
import pickle
import os
import sys
import arviz as avz
from pathlib import Path

path_data = Path(sys.argv[1])

summarry_inference = dict()
for file in os.listdir(path_data):
	if file.endswith("trace.pkl"):
		print(file)
		trace = pickle.load( open(path_data/file, "rb"))
		try:
			summary = avz.summary(trace)
		except:
			summary = None
		a = file.split('_')
		sample_name = a[0]
		amp_index = int(a[2])
		if sample_name not in summarry_inference.keys():
			summarry_inference[sample_name] = dict()
		summarry_inference[sample_name][amp_index] = summary
	pickle.dump(summarry_inference,open(path_data/"summary.pkl", "wb" ))
