import sys
import pandas as pd
import numpy as np

import json

import argparse

def get_user_arguments():
    """
    Use argparse to get user arguments
    input:
    None
    output:
    args = argparse object
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="Group proteins file name", type=str, default="20190115_HEKwt_and_MICS1ko_proteinGroups.txt")
    parser.add_argument("settings_filename", help="Settings file name", type=str, default="maxquant_settings.json")
    development_arguments = ["20190115_HEKwt_and_MICS1ko_proteinGroups.txt", "maxquant_settings.json"]

    args = parser.parse_args(development_arguments)
    return args

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
        protein_groups_dataframe = pd.read_csv(filename, sep='\t', index_col = 'Majority protein IDs', low_memory=False)
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

def filter_dataframe_columns(protein_groups_dataframe, settings_dict):
    """
    input:
    protein_groups_dataframe = pd.DataFrame
    settings_dict = dict{parameter: [values]}
    output:
    protein_groups_dataframe = pd.DataFrame
    """
    print("Start removing unwanted columns from the dataframe:")
    all_dataframe_columns = list(protein_groups_dataframe.columns)
    selected_dataframe_columns = []
    for key_word, method in zip(["EXACT_MATCHES", "CONTAINS"], ["exact_matches", "contains"]):
        selected_columns = select_columns(all_dataframe_columns, settings_dict[key_word], method)
        selected_dataframe_columns += selected_columns
    protein_groups_dataframe = protein_groups_dataframe[selected_dataframe_columns]
    
    number_of_rows = protein_groups_dataframe.shape[0]
    number_of_columns = protein_groups_dataframe.shape[1]
    print(f"Finished removing unwanted columns from the dataframe. The filtered dataframe has {number_of_rows} rows and {number_of_columns} columns")
    print("-"*40)
    return protein_groups_dataframe

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

def filter_dataframe_rows(protein_groups_dataframe, settings_dict):
    """
    input:
    protein_groups_dataframe = pd.DataFrame
    settings_dict = dict{parameter: [values]}
    output:
    protein_groups_dataframe = pd.DataFrame
    """
    print("Start filtering out unwanted proteins and drop proteins containing NA values: ")
    non_applying_rows, applying_rows = select_proteins(protein_groups_dataframe.index, settings_dict["PROTEIN_FILTERS"])
    filtered_groups_dataframe = protein_groups_dataframe.loc[applying_rows]
    protein_groups_dataframe = protein_groups_dataframe.loc[non_applying_rows]
    #drop any row containing NA values:
    protein_groups_dataframe = protein_groups_dataframe.dropna(axis="index")
    number_of_rows = protein_groups_dataframe.shape[0]
    number_of_columns = protein_groups_dataframe.shape[1]
    print(f"Finished filtering out unwanted proteins. The filtered dataframe has {number_of_rows} rows and {number_of_columns} columns")
    print("-"*40)
    return protein_groups_dataframe, filtered_groups_dataframe

def select_proteins(row_labels_index, protein_filters):
    """
    selects proteins to keep using index and protein_filter, but keep proteins which do not apply to the filter.
    input:
    row_labels_index = pd.DataFrame().index, the row labels(in this case 'Majority protein IDs') of the dataframe
    protein_filters = list, list of filters 
    output:
    keep_proteins = list, list with rows to keep
    """
    for protein_filter in protein_filters:
        non_applying_proteins = [row for row in row_labels_index if not protein_filter in row]
        applying_proteins = [row for row in row_labels_index if protein_filter in row]
        
    return non_applying_proteins, applying_proteins

def fetch_identifiers(protein_groups_dataframe):
    """
    input:
    protein_groups_dataframe = pd.DataFrame
    output:
    protein_groups_dataframe = pd.DataFrame
    """
    print("Start fetching identifiers: ")
    protein_groups_dataframe['identifier'] = protein_groups_dataframe['Fasta headers'].apply(parse_identifier)
    print("Finished fetching identifiers. Please see the first 5 rows of the dataframe: ")
    print(protein_groups_dataframe.head())
    print("-"*40)
    return protein_groups_dataframe

def parse_identifier(fasta_header):
    """
    parse protein identifier from the fasta header
    input:
    fasta_header = string
    output:
    fasta_identifier = string, a protein identifier encoded into the fasta header
    """
    try:
        fasta_identifier = fasta_header.split('|')[1]
        return fasta_identifier
    except Exception as error:
        print(f'Parsing a fasta header to get the identifier did not work for:\n{entry}')
        return np.nan

if __name__ == "__main__":
    args = get_user_arguments()

    #Read in files:
    settings_dict = load_json(args.settings_filename)
    protein_groups_dataframe = read_in_protein_groups_file(args.filename)

    #process maxquant file:
    protein_groups_dataframe = filter_dataframe_columns(protein_groups_dataframe, settings_dict)
    protein_groups_dataframe, filtered_groups_dataframe = filter_dataframe_rows(protein_groups_dataframe, settings_dict)
    protein_groups_dataframe = fetch_identifiers(protein_groups_dataframe)
    
    



    
    # processed.to_csv('processed_out.tsv', sep='\t')
