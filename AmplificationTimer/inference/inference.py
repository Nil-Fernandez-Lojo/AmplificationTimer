import numpy as np
import arviz as az
from models import Model1, Model2, Model3, Model4
import matplotlib.pyplot as plt
import pickle
import pymc3 as pm
import arviz as avz

def inference_t(amplification,
	model_idx,
	subclonality_modelled,
	filter_APOBEC,
	true_t = None,
	true_u = None,
	save = True):

	if subclonality_modelled:
		path_subclonality_modelled = "subclonality_modelled"
	else:
		path_subclonality_modelled = "subclonality_not_modelled"

	if filter_APOBEC:
		path_filter_APOBEC = "filter_APOBEC"
	else:
		path_filter_APOBEC = "no_filter_APOBEC"

	samplename = amplification.clinical_data['samplename']
	trace_path = amplification.config["path_folder_inferred_times"]/"traces"/path_filter_APOBEC/path_subclonality_modelled/("Model"+str(model_idx))
	figure_path = amplification.config["path_folder_inferred_times"]/"figures"/path_filter_APOBEC/path_subclonality_modelled/("Model"+str(model_idx))

	if model_idx == 1:
		model = Model1(amplification,
					   amplification.mutation_rate,
					   amplification.subclonal_structure,
					   subclonality_modelled,
					   amplification.config["cores"],
					   filter_APOBEC)
		if save:
			t = np.linspace(0,1,100)
			analytical_posterior = model.get_analytical_posterior(t)
			analytical_posterior_path = trace_path / (samplename+'_amp_'+str(amplification.chromosome)+amplification.arm+'_analytical_posterior.pkl')
			with open(analytical_posterior_path, 'wb') as buff:
				pickle.dump(analytical_posterior, buff)
	elif model_idx == 2:
		model = Model2(amplification,
					   amplification.mutation_rate,
					   amplification.subclonal_structure,
					   subclonality_modelled,
					   amplification.config["cores"],
					   filter_APOBEC)

	elif model_idx == 3:
		model = Model3(amplification,
					   amplification.mutation_rate,
					   amplification.subclonal_structure,
					   subclonality_modelled,
					   amplification.config["cores"],
					   filter_APOBEC)

	elif model_idx == 4:
		model = Model4(amplification,
					   amplification.mutation_rate,
					   amplification.subclonal_structure,
					   subclonality_modelled,
					   amplification.config["cores"],
					   filter_APOBEC)
	
	trace = model.get_MCMC_samples_posterior(amplification.config["n_MCMC_iterations"])
	summary = avz.summary(trace)
	
	if save:
		ax = pm.plots.plot_posterior(trace)
		if true_t is not None:
			if model_idx <=2:
				ax.axvline(x=true_t)
			else:
				print(ax.flatten())
				print(ax.flatten()[0])
				ax.flatten()[0].axvline(x=true_t)
		if true_u is not None:
			if not hasattr(true_u, '__iter__'):
				true_u = [true_u]
			for i,u in enumerate(true_u):
				print(i,u)
				print(ax.flatten()[1+i])
				ax.flatten()[1+i].axvline(x=u)
		save_file_prefix = samplename+str(amplification.chromosome)+amplification.arm
		trace_path = trace_path / (save_file_prefix+'_trace.pkl')
		with open(trace_path, 'wb') as buff:
			pickle.dump(trace, buff)
		figure_path = figure_path / (save_file_prefix+"_plot_posterior.png")
		plt.savefig(figure_path, bbox_inches='tight')
	return summary