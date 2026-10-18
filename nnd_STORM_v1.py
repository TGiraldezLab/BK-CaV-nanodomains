# El programa calcula la minima distancia de cada particula
# con su nearest-neightbor (mismo canal y entre canales)
# y lo muesta en un histograma

import os
import pandas as pd
from scipy import spatial
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from time import time
from matplotlib.ticker import FuncFormatter

tiempo_inicial = time()

directorio = 'C:\\NND_txt\\'   # Cambiar directorio
mol_list = '20240411_488_BK_647_2B WT_4006_CORREGIDO_DUDA.txt'   # Cambiar nombre archivo
os.chdir(directorio)

df = pd.read_csv(mol_list, sep="	")
df_647 = df.loc[df['Channel Name'] == 647]
df_647 = df_647[['Xwc', 'Ywc']].copy()
df_488 = df.loc[df['Channel Name'] == 488]
df_488 = df_488[['Xwc', 'Ywc']].copy()

ch1 = df_647.values
ch2 = df_488.values

mins_ch1 = []
mins_ch2 = []
mins_ch1_ch2 = []
mins_ch2_ch1 = []

# Nomes ch1
print(' ch1')
for i in range(len(ch1)):
    llista = [ch1[i]]
    distance,index = spatial.cKDTree(llista).query(ch1)
    distance = distance[distance != 0]
    mins_ch1.append(min(distance))
    
# Nomes ch2
print(' ch2')
for i in range(len(ch2)):
    llista = [ch2[i]]
    distance,index = spatial.cKDTree(llista).query(ch2)
    distance = distance[distance != 0]
    mins_ch2.append(min(distance))

# Entre ch1 i ch2
print(' ch1 i ch2')
for i in range(len(ch1)):
    llista = [ch1[i]]
    distance,index = spatial.cKDTree(llista).query(ch2)
    mins_ch1_ch2.append(min(distance))  
    
# Entre ch2 i ch1
print(' ch2 i ch1')
for i in range(len(ch2)):
    llista = [ch2[i]]
    distance,index = spatial.cKDTree(llista).query(ch1)
    mins_ch2_ch1.append(min(distance))
    
#print('')
print('nº de particulas canal 647:', len(mins_ch1))
print('nº de particulas canal 488:', len(mins_ch2))
print('nº de particulas canal 647 con 488:', len(mins_ch1_ch2))
print('nº de particulas canal 488 con 647:', len(mins_ch2_ch1))

print('Valor maximo canal 647:', max(mins_ch1))
print('Valor maximo canal 488:', max(mins_ch2))
print('Valor maximo canal 647 con 488:', max(mins_ch1_ch2))
print('Valor maximo canal 488 con 647:', max(mins_ch2_ch1))

# Histograma
step_bins = 100  # <-- CAMBIAR
xmax_bins = 400  # <-- CAMBIAR
bins = np.linspace(0, xmax_bins, step_bins)
#plt.hist(mins_ch1, bins, weights=np.ones_like(mins_ch1) / len(mins_ch1), alpha=0.5, label = '647')
#plt.hist(mins_ch2, bins, weights=np.ones_like(mins_ch2) / len(mins_ch2), alpha=0.5, label = '488')
plt.hist(mins_ch1_ch2, bins, weights=np.ones_like(mins_ch1_ch2) / len(mins_ch1_ch2), alpha=0.5, label = '647 & 488')
plt.hist(mins_ch2_ch1, bins, weights=np.ones_like(mins_ch2_ch1) / len(mins_ch2_ch1), alpha=0.5, label = '488 & 647')
plt.xlabel('Distance (nm)')
plt.ylabel('Number of counts (%)')
plt.legend()
formatter = FuncFormatter(lambda y, _: str(y * 100) + '%')
plt.gca().yaxis.set_major_formatter(formatter)
plt.show()

# Guardar a csv
np.savetxt("C:\\NND_txt\\mins_647.csv", mins_ch1, delimiter=",", fmt='%d')
np.savetxt("C:\\NND_txt\\mins_488.csv", mins_ch2, delimiter=",", fmt='%d')
np.savetxt("C:\\NND_txt\\mins_647_488.csv", mins_ch1_ch2, delimiter=",", fmt='%d')
np.savetxt("C:\\NND_txt\\mins_488_647.csv", mins_ch2_ch1, delimiter=",", fmt='%d')

# Computo tiempo procesamiento
tiempo_final = time() 
tiempo_ejecucion = tiempo_final - tiempo_inicial
print('')
print('El tiempo de ejecucion fue:',tiempo_ejecucion, 's')