import matplotlib.pyplot as plt
import pandas as pd
import json
import matplotlib.pyplot as plt
from scipy.stats import pearsonr,spearmanr

file = "../snvs_1_copy_vs_major_cn.json"
with open(file) as json_file:
    x = json.load(json_file)

df = pd.DataFrame.from_dict(x)
y = df['n_snvs_1']/(df['len_seg']*(df['major_cn']+df['minor_cn']))
print(df)
print(y)
print(pearsonr(df['major_cn'],y))
print(spearmanr(df['major_cn'],y))
plt.scatter(df['major_cn'],y)
plt.yscale('log')
plt.xlabel('Major CN')
plt.ylabel('Number snvs on fewer than 1 copy / Total copy number')
plt.show()
