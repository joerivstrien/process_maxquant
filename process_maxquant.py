import pandas as pd
import numpy as np

import argparse

EXACT_MATCHES = [
    # 'Majority protein IDs',
    'Number of proteins',
    'Razor + unique peptides',
    'Fasta headers',
    'PEP',
    'Intensity',
    'Mol. weight [kDa]',
    'Sequence length',
    'MS/MS Count'
    'iBAQ',
]

CONTAINS = [
    'iBAQ ',
]

PROTEIN_FILTERS = [
    'REV',
    'CON',
]

def process_maxquant(filename):
    """"""
    data = pd.read_csv(filename, sep='\t', index_col = 'Majority protein IDs')
    print(data.shape)

    # select columns to keep
    columns = list(data.columns)
    columns = select_columns(columns)
    data=data[columns]
    print(data.shape)

    # filter unwanted proteins
    proteins = select_proteins(data.index)
    data = data.loc[proteins]

    print(data.shape)

    data = data.dropna()
    print(data.shape)

    # fetch id's

    data['identifier'] = data['Fasta headers'].apply(parse_identifier)

    print(data.head())

    return data

def parse_identifier(entry):
    try:
        result = entry.split('|')[1]
    except:
        print(f'didnt work!: {fasta_entry}')
        return np.nan
    return result    

def select_columns(columns):
    """selects columns to keep"""

    keep_cols = []

    for match in EXACT_MATCHES:
        keep = [col for col in columns if col == match]
        keep_cols += keep

    for match in CONTAINS:
        keep = [col for col in columns if col.startswith(match)] 
        keep_cols += keep
    
    return keep_cols


def select_proteins(index):
    """selects proteins to keep"""
    keep_proteins = index
    for protfilt in PROTEIN_FILTERS:
        keep_proteins = [prot for prot in keep_proteins if protfilt not in prot]

    return keep_proteins

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str)
    
    args = parser.parse_args()

    processed = process_maxquant(args.filename)
    
    # processed.to_csv('processed_out.tsv', sep='\t')