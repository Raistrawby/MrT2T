import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from collections import defaultdict

# === PARAMÈTRES ===
FILENAME = 'MAP12AY07_H0_bin50kb_R4_corrected.ginteractions.tsv'
BIN_SIZE = 50000  # taille du bin en bp
MIN_CONTACT = 0   # seuil pour ignorer les bins vides
OUTDIR = 'hic_centromere_output'

# Créer le dossier de sortie
os.makedirs(OUTDIR, exist_ok=True)

# Charger le fichier Hi-C
df = pd.read_csv(FILENAME, sep='\t', header=None,
                 names=['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'count'])

# Rassembler les chromosomes présents
chromosomes = sorted(set(df['chr1'].tolist()) | set(df['chr2'].tolist()))

# Fonction pour créer et analyser la matrice d'un chromosome
def analyze_chromosome(chrom):
    # Garder uniquement les interactions intra-chromosomiques
    d = df[(df['chr1'] == chrom) & (df['chr2'] == chrom)].copy()

    # Obtenir tous les bins
    bins = sorted(set(d['start1'].tolist() + d['start2'].tolist()))
    bin_index = {v: i for i, v in enumerate(bins)}
    size = len(bins)

    # Construire la matrice
    mat = np.zeros((size, size))
    for _, row in d.iterrows():
        i = bin_index[row['start1']]
        j = bin_index[row['start2']]
        mat[i, j] = row['count']
        mat[j, i] = row['count']  # symétrique

    # Sauvegarder la heatmap
    plt.figure(figsize=(6, 5))
    plt.imshow(np.log1p(mat), cmap='hot', origin='lower')
    plt.title(f'{chrom} - Hi-C contact matrix')
    plt.xlabel('50kb bins')
    plt.ylabel('50kb bins')
    plt.colorbar(label='log1p(contact count)')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, f'{chrom}_heatmap.png'))
    plt.close()

    # Calculer la densité de contact
    bin_strength = mat.sum(axis=0) + mat.sum(axis=1)

    # Sauvegarder le profil de densité
    plt.figure(figsize=(7, 4))
    plt.plot(np.arange(size) * BIN_SIZE / 1e6, bin_strength, lw=1.5)
    plt.title(f'{chrom} - contact density per 50kb bin')
    plt.xlabel('Genomic position (Mb)')
    plt.ylabel('Total contact count')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, f'{chrom}_density.png'))
    plt.close()

    # Trouver le pic (position centromérique candidate)
    peak_index = np.argmax(bin_strength)
    peak_start = bins[peak_index]
    peak_end = peak_start + BIN_SIZE
    print(f"[{chrom}] Centromere candidate at: {peak_start // 1000}-{peak_end // 1000} kb")

# === Appliquer sur chaque chromosome ===
for chrom in chromosomes:
    analyze_chromosome(chrom)

print(f"\n📁 All output files saved in: {OUTDIR}/")
