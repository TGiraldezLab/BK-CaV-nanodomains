import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
import itertools
from scipy.spatial import ConvexHull
from matplotlib.ticker import FuncFormatter
from scipy.optimize import curve_fit

'''
Poner en una misma carpeta los archivos de molecule lists .txt y este script.
Se crearan 3 archivos .csv con todas las areas de todos los molecule lists del 647,
el 488 y el 647 con el 488. Los archivos se nombran como 'archivo'_647.csv, 
'archivo'_488.csv y 'archivo'_647-488.csv respectivamente. El nombre del 'archivo'
se pone a continuacion:
'''

#===MODIFICAR========================================================
# Nombre archivos csv
archivo = ''   # El nombre de archivo tiene que ir entre comillas
# Parametros DBSCAN
epsilon = 60
minimum_samples = 10
#Step del bin del histograma
step = 2000
# Transparencia barras histograma
transp = 0.6
# Colores de los canales
color_ch1 = 'red'
color_ch2 = 'green'
color_ch3 = 'yellow'
#====================================================================


llista  = (glob.glob("*.txt"))

ch1_tots = []
ch2_tots = []
ch3_tots = []; ch3_new_tots = []


for file in llista:
    print(file)
    df = pd.read_csv(file, sep="\t")
    df_647 = df.loc[df['Channel Name'] == 647]
    df_647 = df_647[['Xwc', 'Ywc']].copy()
    df_488 = df.loc[df['Channel Name'] == 488]
    df_488 = df_488[['Xwc', 'Ywc']].copy()
    df_647_488 = df[['Xwc', 'Ywc']].copy()

    ch1_data = df_647.values
    ch2_data = df_488.values
    ch1ch2_data = df_647_488.values
    
    ch2_data = np.roll(ch2_data, 0, axis=0)

    # Ch1: 647
    db1 = DBSCAN(eps=epsilon, min_samples=minimum_samples, metric='euclidean').fit(ch1_data)
    db1_labels = db1.labels_
    db1n_clusters_ = len(set(db1_labels)) - (1 if -1 in db1_labels else 0)
    print('Estimated number of clusters in 647: %d' % db1n_clusters_)
    labels_ch1 = db1.labels_
    ch1_counts = np.bincount(labels_ch1[labels_ch1>=0])
    #print('Number of particles per cluster in 647:', ch1_counts)
    #print('')

    # Ch2: 488
    db2 = DBSCAN(eps=epsilon, min_samples=minimum_samples, metric='euclidean').fit(ch2_data)
    db2_labels = db2.labels_
    db2n_clusters_ = len(set(db2_labels)) - (1 if -1 in db2_labels else 0)
    print('Estimated number of clusters in 488: %d' % db2n_clusters_)
    labels_ch2 = db2.labels_
    ch2_counts = np.bincount(labels_ch2[labels_ch2>=0])
    #print('Number of particles per cluster in 488:', ch2_counts)

    # ch1 + ch2
    db3 = DBSCAN(eps=epsilon, min_samples=minimum_samples, metric='euclidean').fit(ch1ch2_data)
    db3_labels = db3.labels_
    db3n_clusters_ = len(set(db3_labels)) - (1 if -1 in db3_labels else 0)
    print('Estimated number of clusters in 647 + 488: %d' % db3n_clusters_)
    labels_ch3 = db3.labels_
    ch3_counts = np.bincount(labels_ch3[labels_ch3>=0])
    #print('Number of particles per cluster in 647 + 488:', ch3_counts)


    # Area de los clusters
    def Area(vertices):
        n = len(vertices) # of corners
        a = 0.0
        for i in range(n):
            j = (i + 1) % n
            a += (vertices[i][0] * vertices[j][1]-vertices[j][0] * vertices[i][1])
        result = abs(a / 2.0)
        return result

    # ch1
    clusters_ch1 = [ch1_data[db1_labels == i] for i in range(db1n_clusters_)]
    areas_ch1 = []
    for i,item in enumerate(clusters_ch1):
        points = np.array(clusters_ch1[i])
        hull = ConvexHull(points)
        A = hull.points
        B = hull.vertices
        ixgrid = np.ix_(B)
        AA = A[ixgrid]
        areas_ch1.append(Area(AA))    
    #print(areas_ch1)

    # ch2
    clusters_ch2 = [ch2_data[db2_labels == i] for i in range(db2n_clusters_)]
    areas_ch2 = []
    for i,item in enumerate(clusters_ch2):
        points = np.array(clusters_ch2[i])
        hull = ConvexHull(points)
        A = hull.points
        B = hull.vertices
        ixgrid = np.ix_(B)
        AA = A[ixgrid]
        areas_ch2.append(Area(AA))    
    #print(areas_ch2)

    # ch1 ch2
    clusters_ch1ch2 = [ch1ch2_data[db3_labels == i] for i in range(db3n_clusters_)]
    areas_ch1ch2 = []
    for i,item in enumerate(clusters_ch1ch2):
        points = np.array(clusters_ch1ch2[i])
        hull = ConvexHull(points)
        A = hull.points
        B = hull.vertices
        ixgrid = np.ix_(B)
        AA = A[ixgrid]
        areas_ch1ch2.append(Area(AA))    
    #print(areas_ch1ch2)


    #-----------------
    xch1 = ch1_counts
    ych1 = areas_ch1
    xch2 = ch2_counts
    ych2 = areas_ch2
    xch3 = ch3_counts
    ych3 = areas_ch1ch2


    #-------------------
    ch1 = []
    for combination in itertools.zip_longest(xch1, ych1):
        ch1.append(combination)

    ch2 = []
    for combination in itertools.zip_longest(xch2, ych2):
        ch2.append(combination)

    ch3 = []
    for combination in itertools.zip_longest(xch3, ych3):
        ch3.append(combination)

    eliminar_ch1ch3 = list(set(ch1).intersection(ch3))
    eliminar_ch2ch3 = list(set(ch2).intersection(ch3))
    #print('Elementos ch1 (647): ', len(ch1))
    #print('Elementos ch2 (488): ', len(ch2))
    #print('Elementos ch1 + ch2 (647 + 488): ', len(ch3))
    #print('Elementos comunes a eliminar ch1-ch3: ', len(eliminar_ch1ch3))
    #print('Elementos comunes a eliminar ch2-ch3: ', len(eliminar_ch2ch3))


    # Eliminar de ch3: ch1 - ch3
    index_elim_ch1ch3 = []

    for item in eliminar_ch1ch3:
        index_elim_ch1ch3.append(ch3.index(item))

    for index in sorted(index_elim_ch1ch3, reverse=True):
        del ch3[index]


    # Eliminar de ch3: ch2 - ch3
    index_elim_ch2ch3 = []

    for item in eliminar_ch2ch3:
        index_elim_ch2ch3.append(ch3.index(item))

    for index in sorted(index_elim_ch2ch3, reverse=True):
        del ch3[index]


    #print('Elementos ch1 + ch2 tras las substracciones: ', len(ch3))



    #------------------
    xch3_new = [x[0] for x in ch3]
    ych3_new = [x[1] for x in ch3]
    print('Final number of 647 + 488 clusters: ', len(ych3_new))
    #----------------
    ch1_tots.extend(ych1)
    ch2_tots.extend(ych2)
    ch3_new_tots.extend(ych3_new)



# Cambiar unidades de nm a nm
ch1_tots[:] = [x / 1 for x in ch1_tots]
ch2_tots[:] = [x / 1 for x in ch2_tots]
ch3_new_tots[:] = [x / 1 for x in ch3_new_tots]

# Guardar a un csv
np.savetxt(archivo+"_647.csv", ch1_tots, delimiter=",")
np.savetxt(archivo+"_488.csv", ch2_tots, delimiter=",")
np.savetxt(archivo+"_647-488.csv", ch3_new_tots, delimiter=",")

# Plots
tipo = 'bar'
tots = ch1_tots+ch2_tots+ch3_new_tots 




#647----------------------------------
Bins=np.arange(min(ch1_tots), max(ch1_tots), step)
bin_heights_ch1, bin_borders_ch1, _ = plt.hist(ch1_tots, fc=color_ch1, alpha=transp, bins=Bins, histtype=tipo, label='V618G', weights=np.ones_like(ch1_tots) / len(tots), ec='black')
#647----------------------------------


#647-488------------------------------
Bins=np.arange(min(ch3_new_tots), max(ch3_new_tots), step)
bin_heights_ch3, bin_borders_ch3, _ = plt.hist(ch3_new_tots, fc=color_ch3, alpha=transp, bins=Bins, histtype=tipo, label='BKα + V618G', weights=np.ones_like(ch3_new_tots) / len(tots), ec='black')
#647-488------------------------------




#488----------------------------------
Bins=np.arange(min(ch2_tots), max(ch2_tots), step)
bin_heights_ch2, bin_borders_ch2, _ = plt.hist(ch2_tots, fc=color_ch2, alpha=transp, bins=Bins, histtype=tipo, label='BKα', weights=np.ones_like(ch2_tots) / len(tots), ec='black')
#488----------------------------------




formatter = FuncFormatter(lambda y, _: str(y * 100))
plt.gca().yaxis.set_major_formatter(formatter)
plt.legend(fontsize=12, frameon=False)
plt.xlabel('Cluster area (nm²)', fontsize=14)
plt.ylabel('Clusters (%)', fontsize=14)
plt.title(archivo)
plt.xlim(0,80000)
plt.ylim(0,0.06)
plt.savefig(archivo+'_histograma.png', bbox_inches='tight', dpi=1200)
plt.savefig(archivo+'_histograma.pdf', bbox_inches='tight')

plt.show()