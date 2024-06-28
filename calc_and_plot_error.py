import os
import tqdm
import subprocess
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter
from itertools import product
from matplotlib.colors import LinearSegmentedColormap

def read_reference_sequences(refer_path):
    return [str(item.seq) for idx, item in enumerate(SeqIO.parse(refer_path, 'fasta'))]

def get_motifs(base_alphabet, step):
    return [''.join(_) for _ in product(base_alphabet, repeat=step)]

def align_sequences(query, ref, step, base_alphabet, motifs):
    error_type_array = np.zeros((len(base_alphabet) ** step, len(base_alphabet) ** step))
    with open('temp_align.fasta', 'w') as f:
        f.write(f'>query\n{query}\n>ref\n{ref}\n')
    
    with open(os.devnull, 'w') as devnull:
        subprocess.run('muscle -align temp_align.fasta -output temp_aligned.fasta', shell=True, stdout=devnull, stderr=devnull, check=True)
    
    aligned_query, aligned_ref = '', ''
    for record in SeqIO.parse('temp_aligned.fasta', 'fasta'):
        if record.id == 'query':
            aligned_query = str(record.seq)
        if record.id == 'ref':
            aligned_ref = str(record.seq)
    
    for i in range(len(aligned_ref) - step + 1):
        a, b = aligned_query[i:i+step], aligned_ref[i:i+step]
        if a in motifs and b in motifs:
            error_type_array[motifs.index(a), motifs.index(b)] += 1

    return error_type_array

def calculate_error_rates(refs, motifs, max_sample_num, step, base_alphabet):
    calculating_info = pd.DataFrame(columns=['seq_index', 'copy_index'] + [f'{i[0]}>{i[1]}' for i in product(motifs, repeat=2)])
    
    with tqdm.trange(len(refs), desc='Calculating error rate...') as pbar:
        for idx, ref in enumerate(refs):
            reads = [str(item.seq) for item in SeqIO.parse(f'/data/biolab-nvme-pool1/zhaoxy/guppy_basecalled/smeseqs_2307/grouping_res_allbase/grouping_res_{idx//8}_{idx%8}.fasta', 'fasta')]
            if len(reads) > max_sample_num:
                reads = random.sample(reads, max_sample_num)
            
            for copy_index, query in enumerate(reads):
                error_type_array = align_sequences(query, ref, step, base_alphabet, motifs)
                temp_dict = {'seq_index': idx, 'copy_index': copy_index}
                for row in range(error_type_array.shape[0]):
                    for column in range(error_type_array.shape[1]):
                        temp_dict[f'{motifs[row]}>{motifs[column]}'] = error_type_array[row, column]
                
                calculating_info = pd.concat([calculating_info, pd.DataFrame(temp_dict, index=[0])], ignore_index=True)
            
            pbar.update()
    
    calculating_info.fillna(0, inplace=True)
    return calculating_info

def create_confusion_matrix(calculating_info, base_alphabet, step):
    all_type_error = calculating_info.sum()[2:]
    error_type_array = np.zeros((len(base_alphabet)**step, len(base_alphabet)**step))

    for idx, error_type in enumerate(all_type_error):
        if idx // (len(base_alphabet)**step) != idx % (len(base_alphabet)**step):
            error_type_array[idx // (len(base_alphabet)**step)][idx % (len(base_alphabet)**step)] = error_type

    return error_type_array * 100 / np.sum(error_type_array)

def plot_heatmap(confusion_matrix, motifs, colors, output_file):
    df_cm = pd.DataFrame(confusion_matrix.T, index=motifs, columns=motifs)
    np.fill_diagonal(df_cm.values, np.nan)

    custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=100)

    plt.figure(figsize=(5, 5))
    ax = sns.heatmap(df_cm, annot=True, fmt='.2f', cmap=custom_cmap, annot_kws={"size": 12}, cbar_kws={'label': 'Substitution (%)'})

    cbar = ax.collections[0].colorbar
    cbar.ax.yaxis.label.set_size(12)
    cbar.ax.tick_params(labelsize=12)

    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('Substitution', fontsize=12)
    plt.ylabel('Reference', fontsize=12)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()

def main():
    refer_path = r'/data/nas-shared/zhaoxy/CompositeHedges/all_seqs20230512/zxy_sme_seqs.fasta'
    refs = read_reference_sequences(refer_path)

    MAX_SAMPLE_NUM = 20
    STEP = 1
    base_alphabet = 'ACGT'
    motifs = get_motifs(base_alphabet, STEP)

    calculating_info = calculate_error_rates(refs, motifs, MAX_SAMPLE_NUM, STEP, base_alphabet)
    # calculating_info.to_csv('biased_error_type.csv', index=False)

    confusion_matrix = create_confusion_matrix(calculating_info, base_alphabet, STEP)
    
    orange = ["#f8f8f8", "#ED6C00"] #SUSTech Orange
    azure = ["#f8f8f8", "#2BB7B3"] #SUSTech Azure
    darkgreen = ["#f8f8f8", "#003F43"] #SUSTech Dark green
    plot_heatmap(confusion_matrix, motifs, darkgreen, './error_type_heatmap.svg')

if __name__ == '__main__':
    main()
