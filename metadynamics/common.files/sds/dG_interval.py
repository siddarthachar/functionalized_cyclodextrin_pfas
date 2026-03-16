import numpy as np

a=np.load('kb_dg.npy')
print('dG intervals: ', a[2])
print('dG std intervals: ', a[3])
dg_m = np.mean(a[2])
final_std = np.sqrt(np.mean(a[3]**2 + (a[2] - dg_m)**2))
print(dg_m,' +- ',final_std,' kJ/mol')
