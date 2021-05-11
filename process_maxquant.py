import sys
import pandas as pd
import numpy as np

import json

import argparse

### columns whose name matches these strings
### exactly are kept
##EXACT_MATCHES = [
##    # 'Majority protein IDs',
##    'Number of proteins',
##    'Razor + unique peptides',
##    'Fasta headers',
##    'PEP',
##    'Intensity',
##    'Mol. weight [kDa]',
##    'Sequence length',
##    'MS/MS Count'
##    'iBAQ',
##]
##
### columns containing these strings are kept
##CONTAINS = [
##    'iBAQ ',
##]
##
### proteins whose name contains one of these
###  strings it is filtered out
##PROTEIN_FILTERS = [
##    'REV',
##    'CON',
##]
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
        try:
            json_object = json.load(json_file)
        except FileNotFoundError as file_not_found_error:
            print(file_not_found_error)
        except Exception as error:
            print("A unhandled error has occurred while reading {json_file_path}, please see error belows:")
            print(error)
    return json_object
def read_in_protein_groups_file(filename):
    """
    input:
    string, path of proteingroups file
    output:
    protein_groups_dataframe = pd.DataFrame
    """
    print(f"Read in \'{filename}\'")
    try:
        protein_groups_dataframe = pd.read_csv(filename, sep='\t', index_col = 'Majority protein IDs')
        number_of_rows = protein_groups_dataframe.shape[0]
        number_of_columns = protein_groups_dataframe.shape[1]
    except FileNotFoundError as file_not_found_error:
        print(file_not_found_error)
    except IOError as io_error:
        print(io_error)
    except Exception as error:
        print(f"Something went wrong while reading in \'{filename}\', please see the error below:")
        print(error)
    print(f"Succesfully read in \'{filename}\' and created a dataframe. The dataframe has {number_of_rows} rows and {number_of_columns} columns")
    print("-"*40)
    return protein_groups_dataframe
    
def process_maxquant(protein_groups_dataframe, settings_dict):
    """
    main function that processes maxquant protein-groups file
    input:
    protein_groups_dataframe = pd.DataFrame()
    settings_dict = dict{setting : value}, dict containing parameter settings for this search. 
    output:
    processed_dataframe = pd.DataFrame()
    """
    
    # select columns to keep
    print("Start removing unwanted columns from the dataframe:")
    all_dataframe_columns = list(protein_groups_dataframe.columns)
    selected_dataframe_columns = select_columns(all_dataframe_columns, settings_dict["EXACT_MATCHES"], "exact_matches")
    selected_dataframe_columns = select_columns(all_dataframe_columns, settings_dict["CONTAINS"], "contains")
    
    protein_groups_dataframe = protein_groups_dataframe[selected_dataframe_columns]
    number_of_rows = protein_groups_dataframe.shape[0]
    number_of_columns = protein_groups_dataframe.shape[1]
    print(f"Finished removing unwanted columns from the dataframe. The filtered dataframe has {number_of_rows} rows an {number_of_columns} columns")
    print("-"*40)

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

def select_columns(all_column_names, selected_column_names, method):
    """
    selects columns to keep from list of colnames
    input:
    all_column_names = list, list of all column names of the dataframe
    selected_column_names = list, list of column names to select
    method = string, in what way does the match have to match? 
    output:
    keep_columns = list, list containing columns to retain
    """
    keep_columns = []

    for match in selected_column_names:
        if "exact_matches" == method:
            keep_column = [column for column in all_column_names if column == match]
        elif "contains" == method:
            keep_column = [column for column in all_column_names if column.startswith(match)]
        else:
            print(f"In function 'select_columns' the variable name {method} has been passed but isn't a valid method name. Change this to a valid method name")
            sys.exit()
        keep_columns += keep_column
        
    return keep_columns

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
    protein_groups_dataframe = read_in_protein_groups_file(args.filename)

    processed = process_maxquant(protein_groups_dataframe, settings_dict)
    
    # processed.to_csv('processed_out.tsv', sep='\t')
