import numpy as np
import arviz as az
from models_dev import Model1, Model2, Model3, Model4
import matplotlib.pyplot as plt
import pickle
import pymc3 as pm

def inference_t(sample,
	amplification_idx,
	model_idx,
	subclonality_modelled,
	n_MCMC_iterations,
	folder_models,
	cores,
	true_t = None,
	true_u = None,
	suffix_name = "",
	save = True,
	filter_APOBEC = True):
	# Model1
	amplification = sample.amplifications[amplification_idx]

	if model_idx == 1:
		model = Model1(amplification,sample.mutation_rate,sample.subclonal_structure,subclonality_modelled,cores,filter_APOBEC)

	elif model_idx == 4:
		model = Model4(amplification,sample.mutation_rate,sample.subclonal_structure,subclonality_modelled,cores)
	
	