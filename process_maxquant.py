import pandas as pd
import numpy as np

import json

import argparse

# columns whose name matches these strings
# exactly are kept
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

# columns containing these strings are kept
CONTAINS = [
    'iBAQ ',
]

# proteins whose name contains one of these
#  strings it is filtered out
PROTEIN_FILTERS = [
    'REV',
    'CON',
]
def load_json(json_filepath):
    """
    input:
    json_filepath = string
    output:
    json_object = dict{}
    parameters:
    json_file = IO text object
    """
    with open(json_filepath) as json_file:
        json_object = json.load(json_file)
    
    return json_object

def process_maxquant(filename):
    """
    main function that processes maxquant protein-groups file

    Args:
        filename: path/filename of proteingroups file
    
    Returns:
        processed dataframe
    """
 
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
    """
    parse protein identifier from the fasta header
    """
    try:
        result = entry.split('|')[1]
    except:
        print(f'didnt work!: {entry}')
        return np.nan
    return result    

def select_columns(columns):
    """
    selects columns to keep from list of colnames
    
    Args:
        columns: list, containing all column names

    Returns:
        list, containing columns to retain
    """

    keep_cols = []

    for match in EXACT_MATCHES:
        keep = [col for col in columns if col == match]
        keep_cols += keep

    for match in CONTAINS:
        keep = [col for col in columns if col.startswith(match)] 
        keep_cols += keep
    
    return keep_cols


def select_proteins(index):
    """
    selects proteins to keep from index
    
    Args:
        dataframe index, list-like    
    
    Returns:
        list of protein identifiers to keep
    """
    keep_proteins = index
    for protfilt in PROTEIN_FILTERS:
        keep_proteins = [prot for prot in keep_proteins if protfilt not in prot]

    return keep_proteins

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="Group proteins file name", type=str, default="20190115_HEKwt_and_MICS1ko_proteinGroups.txt")
    parser.add_argument("settings_filename", help="Settings file name", type=str, default="maxquant_settings.json")
    development_arguments = ["20190115_HEKwt_and_MICS1ko_proteinGroups.txt", "maxquant_settings.json"]
    args = parser.parse_args(development_arguments)

    #get settings
    settings_dict = load_json(args.settings_filename)

    processed = process_maxquant(args.filename)
    
    # processed.to_csv('processed_out.tsv', sep='\t')
