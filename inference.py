import numpy as np
import arviz as az
from models import Model1, Model2, Model3, Model4
import matplotlib.pyplot as plt
import pickle
import pymc3 as pm
import arviz as avz

def inference_t(sample,
	amplification_idx,
	model_idx,
	subclonality_modelled,
	n_MCMC_iterations,
	folder_models,
	cores,
	filter_APOBEC,
	true_t = None,
	true_u = None,
	suffix_name = "",
	save = True):
	# Model1
	amplification = sample.amplifications[amplification_idx]
	
	if subclonality_modelled:
		trace_path = folder_models/"traces"/"subclonality_modelled"
		figure_path = folder_models/"figures"/"subclonality_modelled"
	else:
		trace_path = folder_models/"traces"/"subclonality_not_modelled"
		figure_path = folder_models/"figures"/"subclonality_not_modelled"

	if model_idx == 1:
		model = Model1(amplification,sample.mutation_rate,sample.subclonal_structure,subclonality_modelled,cores,filter_APOBEC)
		trace_path = trace_path / "Model1"
		figure_path = figure_path/"Model1"
		if save:
			t = np.linspace(0,1,100)
			analytical_posterior = model.get_analytical_posterior(t)
			analytical_posterior_path = trace_path / (sample.name+'_amp_'+str(amplification_idx)+suffix_name+'_analytical_posterior.pkl')
			with open(analytical_posterior_path, 'wb') as buff:
				pickle.dump(analytical_posterior, buff)
	elif model_idx == 2:
		model = Model2(amplification,sample.mutation_rate,sample.subclonal_structure,subclonality_modelled,cores,filter_APOBEC)
		trace_path = trace_path / "Model2" 
		figure_path = figure_path/"Model2"

	elif model_idx == 3:
		model = Model3(amplification,sample.mutation_rate,sample.subclonal_structure,subclonality_modelled,cores,filter_APOBEC)
		trace_path = trace_path / "Model3" 
		figure_path = figure_path/"Model3"

	elif model_idx == 4:
		model = Model4(amplification,sample.mutation_rate,sample.subclonal_structure,subclonality_modelled,cores,filter_APOBEC)
		trace_path = trace_path / "Model4"
		figure_path = figure_path/"Model4"
	
	trace = model.get_MCMC_samples_posterior(n_MCMC_iterations)
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

		trace_path = trace_path / (sample.name+'_amp_'+str(amplification_idx)+suffix_name+'_trace.pkl')
		with open(trace_path, 'wb') as buff:
			pickle.dump(trace, buff)
		figure_path = figure_path / (sample.name+'_amp_'+str(amplification_idx)+suffix_name+"_plot_posterior.png")
		plt.savefig(figure_path, bbox_inches='tight')


	return summary