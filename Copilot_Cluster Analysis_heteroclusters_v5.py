# -*- coding: utf-8 -*-
"""
Created on Fri Nov 28 17:28:08 2025

@author: eldei
"""



import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os

# ==============================
# CONFIGURACIÓN
# ==============================
EPSILON = 60        # Distancia máxima (nm)
MIN_SAMPLES = 10    # Mínimo de puntos por cluster
GENERAR_GRAFICAS = True

# ==============================
# FUNCIONES
# ==============================
def calcular_area_cluster(puntos):
    if len(puntos) < 3:
        return 0.0
    hull = ConvexHull(puntos)
    return hull.volume  # Área en 2D

def aplicar_dbscan(coords):
    if len(coords) == 0:
        return np.array([])  # Sin etiquetas
    clustering = DBSCAN(eps=EPSILON, min_samples=MIN_SAMPLES).fit(coords)
    return clustering.labels_

def procesar_clusters(df, etiquetas):
    if len(etiquetas) == 0:
        return pd.DataFrame()
    df = df.copy()
    # ✅ Ajustar Cluster_ID para que empiece en 1
    etiquetas = np.where(etiquetas == -1, -1, etiquetas + 1)
    df['Cluster_ID'] = etiquetas
    resultados = []
    for cluster_id in sorted(df['Cluster_ID'].unique()):
        if cluster_id == -1:
            continue
        puntos = df[df['Cluster_ID'] == cluster_id][['Xwc', 'Ywc']].values
        area = calcular_area_cluster(puntos)
        resultados.append({
            'Cluster_ID': cluster_id,
            'Area_nm2': area,
            'Num_Particles': len(puntos)
        })
    return pd.DataFrame(resultados)

def procesar_heteroclusters(df647, df488):
    if len(df647) == 0 or len(df488) == 0:
        return pd.DataFrame()
    df_comb = pd.concat([df647, df488], ignore_index=True)
    coords = df_comb[['Xwc', 'Ywc']].values
    etiquetas = aplicar_dbscan(coords)
    if len(etiquetas) == 0:
        return pd.DataFrame()
    # ✅ Ajustar Cluster_ID para que empiece en 1
    etiquetas = np.where(etiquetas == -1, -1, etiquetas + 1)
    df_comb['Cluster_ID'] = etiquetas
    resultados = []
    for cluster_id in sorted(df_comb['Cluster_ID'].unique()):
        if cluster_id == -1:
            continue
        cluster_points = df_comb[df_comb['Cluster_ID'] == cluster_id]
        puntos = cluster_points[['Xwc', 'Ywc']].values
        area = calcular_area_cluster(puntos)
        num_total = len(cluster_points)
        num_647 = len(cluster_points[cluster_points['Channel Name'] == 647])
        num_488 = len(cluster_points[cluster_points['Channel Name'] == 488])
        # ✅ Filtrar heteroclusters verdaderos
        if num_647 > 0 and num_488 > 0:
            resultados.append({
                'Cluster_ID': cluster_id,
                'Area_nm2': area,
                'Num_Particles': num_total,
                'Num_647': num_647,
                'Num_488': num_488
            })
    return pd.DataFrame(resultados)

# ==============================
# PROCESAMIENTO EN LOTE
# ==============================
archivos = glob.glob("*.txt")

# ✅ Listas para acumulación global
global_clusters_647 = []
global_clusters_488 = []
global_heteroclusters = []

for archivo in archivos:
    df = pd.read_csv(archivo, sep="\t")
    df647 = df[df['Channel Name'] == 647]
    df488 = df[df['Channel Name'] == 488]

    clusters_647 = procesar_clusters(df647, aplicar_dbscan(df647[['Xwc', 'Ywc']].values))
    clusters_488 = procesar_clusters(df488, aplicar_dbscan(df488[['Xwc', 'Ywc']].values))
    heteroclusters = procesar_heteroclusters(df647, df488)

    base_name = os.path.splitext(archivo)[0]
    # ✅ Exportar CSV con separador ';' para Excel
    clusters_647.to_csv(f"{base_name}_647.csv", index=False, sep=';')
    clusters_488.to_csv(f"{base_name}_488.csv", index=False, sep=';')
    heteroclusters.to_csv(f"{base_name}_647-488.csv", index=False, sep=';')

    # ✅ Acumular para globales
    if not clusters_647.empty:
        clusters_647['Archivo'] = base_name
        global_clusters_647.append(clusters_647)
    if not clusters_488.empty:
        clusters_488['Archivo'] = base_name
        global_clusters_488.append(clusters_488)
    if not heteroclusters.empty:
        heteroclusters['Archivo'] = base_name
        global_heteroclusters.append(heteroclusters)

    # ✅ Generar gráficas si hay datos
    if GENERAR_GRAFICAS and not clusters_647.empty:
        plt.figure(figsize=(8, 5))
        sns.histplot(clusters_647['Area_nm2'], bins=30, kde=True)
        plt.title(f"Distribución Áreas - Canal 647 ({base_name})")
        plt.savefig(f"{base_name}_hist_647.png")
        plt.close()

    if GENERAR_GRAFICAS and not heteroclusters.empty:
        plt.figure(figsize=(8, 5))
        sns.histplot(heteroclusters['Area_nm2'], bins=30, kde=True)
        plt.title(f"Distribución Áreas - Heteroclusters ({base_name})")
        plt.savefig(f"{base_name}_hist_hetero.png")
        plt.close()

# ✅ Generar ficheros globales
if global_clusters_647:
    pd.concat(global_clusters_647, ignore_index=True).to_csv("global_647.csv", sep=';', index=False)
if global_clusters_488:
    pd.concat(global_clusters_488, ignore_index=True).to_csv("global_488.csv", sep=';', index=False)
if global_heteroclusters:
    pd.concat(global_heteroclusters, ignore_index=True).to_csv("global_647-488.csv", sep=';', index=False)

print("✅ Procesamiento completado para todos los archivos. Ficheros globales generados.")
