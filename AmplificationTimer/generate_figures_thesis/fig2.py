import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

model_idx_list = [1,4]
path_data = "fig_2_data/"

list_df = []
for file in os.listdir(path_data):
	if file.endswith('.tsv'):
		list_df.append(pd.read_csv(path_data+file, sep = '\t'))
number_entries = 0
for df in list_df:
	number_entries+=df.shape[0]
sim_results = pd.concat(list_df)

i = 0
rmse = []
for model_idx in model_idx_list:
	df_model = sim_results.loc[sim_results["model_idx"] == model_idx]
	e_mut_list = df_model["expected_mutations"].unique()
	fig, ax = plt.subplots(2, 2, sharey=True,sharex=True)
	ax = ax.flatten()
	for i,e_mut in enumerate(e_mut_list):
		if (e_mut == 0): continue
		df = df_model.loc[df_model["expected_mutations"] == e_mut]
		rmse.append([model_idx,e_mut,np.sqrt(np.mean((df["t"].values-df["inferred_t"].values)**2))])
		ax[i].title.set_text("E[N]= "+str(int(e_mut)))
		ax[i].plot([0,1],[0,1],color='k',linewidth = 0.5)
		ax[i].plot(df['t'],df['inferred_t'],marker='.', color='k', linestyle='None')
		ax[i].grid()
		if i>=2:
			ax[i].set_xlabel("t",fontsize=16)
		if i%2==0:
			h = ax[i].set_ylabel(r'$\hat{t}$',fontsize=16)
			h.set_rotation(0)
		i+=1
	fig.suptitle("Model "+str(model_idx))
rmse = pd.DataFrame(data= rmse, columns = ['model_idx',"expected_mutations","rmse"])
plt.figure(20)
for model_idx in model_idx_list:
	df = rmse.loc[rmse['model_idx'] == model_idx]
	plt.plot(df["expected_mutations"], df["rmse"],label = "Model "+str(model_idx))
	plt.xlabel("Expected number of mutations")
	plt.ylabel("RMSE")
	plt.xscale("log")
	plt.legend()
plt.grid()
plt.show()
