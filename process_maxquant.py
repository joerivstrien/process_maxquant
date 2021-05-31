#Import standard python libraries:
import sys
import requests
import time
import re
import json
import argparse
import urllib.parse
import urllib.request
#Import third-part libraries
import pandas as pd
import numpy as np


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

def fetch_uniprot_annotation(identifiers, sleep_time=2, batch_length=100):
    """
    input:
    identifiers = pd.Series
    output:
    protein_data_list = list, list of protein_data_dict
    """
    print(f"Start fetching data from uniprot. After each query {sleep_time} seconds pass before a new protein is queried to uniprot.\nThis safety feature is put in place to prevent being blacklisted.")
    print("There are {n_identifier} identifiers so fetching data from uniprot will take atleast {total_seconds} seconds.".format(n_identifier=str(len(identifiers)), total_seconds=str(int(sleep_time)*len(identifiers))))
    protein_data_dict = {}
    base_url = "https://www.ebi.ac.uk/proteins/api/proteins?"
    request_base_url = "offset=0&size=100&accession="

    function_dict = construct_function_dict()

    #split identifiers into multiple sub arrays of length batch_length:
    identifiers = np.split(identifiers, range(batch_length,len(identifiers), batch_length))
    for n_batch, batch in enumerate(identifiers):
        print(f"Start fetching uniprot data for batch number {n_batch} of the total {len(identifiers)} batches")
        requestURL = base_url+request_base_url+",".join(batch)
        request = requests.get(requestURL, headers={"Accept" : "application/json"})
        if not request.ok:
            print(f"Something went wrong with batch {n_batch} when getting the uniprot data request. Proteins of this batch will be ignored")
            request.raise_for_status()
            for identifier in batch:
                protein_data_dict[identifier] = {"gene_name":np.nan, "protein_name":np.nan, "organism_name":np.nan, "hyperlink":np.nan, "cell_compartment":np.nan, "string_linkout":np.nan}
            continue
        else:
            print(f"Succesfully fetched uniprot data for batch {n_batch}:")
            protein_data_dict = update_protein_data_dict(request.json(), batch, function_dict, protein_data_dict)
        #Let the program sleep for a bit else uniprot is going to be overloaded and I get a problem.
        time.sleep(sleep_time)
    return protein_data_dict
    print("Finished fetching data from uniprot")

def construct_function_dict():
    """
    input:
    None
    output:
    function_dict = dict{function_name : function}
    """
    function_dict = {"gene_name":get_uniprot_gene_name, "protein_name":get_protein_name, "organism_name":get_organism_name, "cell_compartment":get_cell_compartment}
    return function_dict

def update_protein_data_dict(uniprot_output_dict, identifiers, function_dict, protein_data_dict):
    """
    This function is meant to loop over the identifiers and update a dict with values per protein.
    input:
    identifiers = list, list of uniprot identifiers
    function_dict = dict{function_name(string) : function(function)}
    protein_data_dict = dict{protein : {gene_name, protein_name, organism_name, uniprot_hyperlink, cell_comparment}}
    output:
    protein_data_dict = dict{protein : {gene_name, protein_name, organism_name, uniprot_hyperlink, cell_comparment}}
    """
    uniprot_base_url = "https://www.uniprot.org/uniprot/"
    for request_number in range(len(identifiers)):
        identifier = identifiers[request_number]
        uniprot_data_dict = uniprot_output_list[request_number]
        for function_name, function in function_dict.items():
            if not identifier in protein_data_dict:
                protein_data_dict[identifier] = {}
            protein_data_dict[identifier].update({function_name : function(uniprot_data_dict)})
        #Per protein the uniprot_hyperlink needs to be added separately because, the other functions depend on the uniprot_data_dict while this function depends on the identifier. 
        protein_data_dict[identifier].update({"uniprot_hyperlink": uniprot_base_url+identifier})
        print(protein_data_dict[identifier])
        print("-"*40)
    return protein_data_dict
def filter_uniprot_query(uniprot_data_dict, identifier):
    """
    Because the input has so little certainty this function looks like a mess. In order to be error prone each and every variable needs to be checked if it is present. 
    input:
    uniprot_data_dict = [{accession: "", id:"", proteinExistence:"", info:{}, organism:{}, protein:{}, gene:{}, features:{}, dbReferences:{}, keywords:[], references:[], sequence:{}}], this is the best case scenario.
    fields can be missing.
    identifier = string
    output:
    gene_name = string
    protein_name = string
    organism_name = string
    hyperlink = string
    cell_compartment = np.nan
    string_linkout = string
    """
    hyperlink_base_url = "https://www.uniprot.org/uniprot/"
    try:
        gene_name = get_uniprot_gene_name(uniprot_data_dict)
        protein_name = get_protein_name(uniprot_data_dict)
        organism_name = get_organism_name(uniprot_data_dict)
        cell_compartment = get_cell_compartment(uniprot_data_dict)
        second_string_linkout = get_database_reference_element(uniprot_data_dict, "STRING")
        hyperlink = hyperlink_base_url+identifier
        print(string_linkout, second_string_linkout)
    except IndexError as index_error:
        print(f"An index error occured , see below for more information\n{index_error}")
        print(f"Protein {identifier} uniprot output will be ignored")
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    except Exception as error:
        print(f"An error occured while filtering the uniprot query, see below for more information\n{error}")
        print(f"Protein {identifier} uniprot output will be ignored")
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    return gene_name, protein_name, organism_name, hyperlink, cell_compartment, string_linkout

def get_uniprot_gene_name(uniprot_data_dict):
    """
    input:
    uniprot_data_dict = [{accession: "", id:"", proteinExistence:"", info:{}, organism:{}, protein:{}, gene:{}, features:{}, dbReferences:{}, keywords:[], references:[], sequence:{}}], this is the best case scenario.
    output:
    gene_name = string
    """
    if "gene" in uniprot_data_dict.keys():
        gene_name = uniprot_data_dict["gene"][0]["name"]["value"]
    else:
        gene_name = np.nan
    return gene_name

def get_protein_name(uniprot_data_dict):
    """
    input:
    uniprot_data_dict = [{accession: "", id:"", proteinExistence:"", info:{}, organism:{}, protein:{}, gene:{}, features:{}, dbReferences:{}, keywords:[], references:[], sequence:{}}], this is the best case scenario.
    output:
    protein_name = string
    """
    if "protein" in uniprot_data_dict.keys():
        if "recommendedName" in uniprot_data_dict["protein"].keys(): 
            protein_name = uniprot_data_dict["protein"]["recommendedName"]["fullName"]["value"]
        elif "submittedName" in uniprot_data_dict["protein"].keys():
            protein_name = uniprot_data_dict["protein"]["submittedName"][0]["fullName"]["value"]
        elif "alternativeName" in uniprot_data_dict["protein"].keys():
            protein_name = uniprot_data_dict["protein"]["alternativeName"]["fullName"]["value"]
        else:
            print("Found a protein name which is not recommonededName or submittedName but: "+", ".join(uniprot_data_dict["protein"].keys()))
            protein_name = uniprot_data_dict["protein"].values()[0]["fullName"]["value"]
    return protein_name

def get_organism_name(uniprot_data_dict):
    """
    input:
    uniprot_data_dict = [{accession: "", id:"", proteinExistence:"", info:{}, organism:{}, protein:{}, gene:{}, features:{}, dbReferences:{}, keywords:[], references:[], sequence:{}}], this is the best case scenario.
    output:
    organism_name = string 
    """
    if "organism" in uniprot_data_dict.keys():
        organism_name = uniprot_data_dict["organism"]["names"][0]["value"]
    else:
        organism_name = np.nan
    return organism_name

def get_cell_compartment(uniprot_data_dict):
    """
    input:
    uniprot_data_dict = {accession: "", id:"", proteinExistence:"", info:{}, organism:{}, protein:{}, gene:{}, features:{}, dbReferences:{}, keywords:[], references:[], sequence:{}}, this is the best case scenario.
    output:
    cell_compartment=string
    """
    cell_compartment = ""
    if "comments" in uniprot_data_dict.keys():
        for reference in uniprot_data_dict["comments"]:
            if "SUBCELLULAR_LOCATION" == reference["type"]:
                cell_compartment = ""
                for location in reference["locations"]:
                    cell_compartment += str(location["location"]["value"])+";"
            else:
                if cell_compartment == "":
                    cell_compartment = np.nan
    else:
        cell_compartment = np.nan
    return cell_compartment

def get_database_reference_element(uniprot_data_dict, reference_type):
    """
    In the uniprot_data_dict there is a dict with database linkouts. This function enables you to get a specific linkout, mainly used to get a STRING database linkout. 
    input:
    uniprot_data_dict = {accession: "", id:"", proteinExistence:"", info:{}, organism:{}, protein:{}, gene:{}, features:{}, dbReferences:{}, keywords:[], references:[], sequence:{}}, this is the best case scenario.
    reference_type = string
    output:
    string_linkout = string
    """
    element = ""
    if "dbReferences" in uniprot_data_dict.keys():
        for reference in uniprot_data_dict["dbReferences"]:
            if reference_type == reference["type"]:
                element = reference["id"]
                return element
    else:
        return np.nan

def get_string_linkout(identifiers):
    """
    Use the uniprot identifier mapping service to get a linkout to the string database per uniprot identifier.
    input:
    identifiers = list, list of uniprot identifiers
    output:
    string_linkout_dict = dict{identifier : string_linkout}
    """
    base_string_url = "https://string-db.org/network/"
    uniprot_mapping_service_url = 'https://www.uniprot.org/uploadlists/'
    string_linkout_dict = {}

    formatted_identifier_string = process_uniprot_identifier_input(identifiers)
    parameters = {'from': 'ACC+ID', 'to': 'STRING_ID', 'format': 'tab', 'query': formatted_identifier_string}

    data = urllib.parse.urlencode(parameters)
    data = data.encode('utf-8')
    request = urllib.request.Request(uniprot_mapping_service_url, data)
    with urllib.request.urlopen(request) as uniprot_mapping_request:
       response = uniprot_mapping_request.read()
    uniprot_mapped_proteins_dict = process_uniprot_mapping_service_output(response.decode('utf-8'))
    for identifier in identifiers:
        if not identifier in uniprot_mapped_proteins_dict.keys():
            string_linkout_dict[identifier] = np.nan
        else:
            string_linkout_dict[identifier] = base_string_url+uniprot_mapped_proteins_dict[identifier]
    return string_linkout_dict

def process_uniprot_identifier_input(identifiers):
    """
    input:
    identifiers = list, list of uniprot identifiers
    output:
    formatted_identifier_string = string, should be "identifier identifier etc."
    """
    pattern = "\-[0-9]{1}$"
    formatted_identifier_string = ""
    for identifier in identifiers:
        #A regex function is used to identify proteins with for example: "-2" as suffix and removes the suffix.
        if None != re.search(pattern, identifier): identifier = identifier[:-2]
        formatted_identifier_string += identifier+" "
    return formatted_identifier_string
    
def process_uniprot_mapping_service_output(uniprot_mapped_proteins):
    """
    This function enables to process more than one mapped identifier from the uniprot mapping service.
    input:
    uniprot_mapped_proteins = string, should look like: From\tTo\nuniprot_identifier\tstring_identifier\nuniprot_identifier\tstring_identifier\n
    output:
    uniprot_mapped_proteins_dict = dict{uniprot identifier:string identifier}
    """
    uniprot_mapped_proteins_dict = {}
    uniprot_mapped_proteins = uniprot_mapped_proteins.replace("From\tTo\n", "")
    for comparison in uniprot_mapped_proteins.split("\n"):
        if not "" == comparison or None == comparison:
            uniprot_id, string_id = comparison.split("\t")
            uniprot_mapped_proteins_dict[uniprot_id] = string_id
    return uniprot_mapped_proteins_dict

def append_uniprot_data_to_dataframe(protein_groups_dataframe, protein_data_list):
    """
    input:
    protein_groups_dataframe = pd.DataFrame
    protein_data_list = dict{identifier: {"gene_name":"", "protein_name":"", "organism_name":"", "hyperlink":"", "cell_compartment":np.nan, "string_linkout":""}}
    output:
    protein_groups_dataframe = pd.DataFrame
    """
    print("Start adding the uniprot columns to the existing dataframe")
    uniprot_column_names = ["gene_name", "protein_name", "organism_name", "hyperlink", "cell_compartment", "string_linkout"]
    for uniprot_column_name in uniprot_column_names:
        uniprot_column_values = get_uniprot_column_values(protein_groups_dataframe["identifier"], uniprot_column_name, protein_data_list)
        protein_groups_dataframe[uniprot_column_name] = uniprot_column_values
    print("Finished adding the uniprot columns to the existing dataframe")
    print("-"*40)
    return protein_groups_dataframe

def get_uniprot_column_values(identifiers, uniprot_column_name, protein_data_list):
    """
    input:
    identifiers = pd.Series
    uniprot_column_name = string
    protein_data_list = list, list of dicts
    output:
    uniprot_column_values = list
    """
    uniprot_column_values = []
    for identifier in identifiers:
        uniprot_column_value = protein_data_list[0][identifier][uniprot_column_name]
        uniprot_column_values.append(uniprot_column_value)
    return uniprot_column_values
    
if __name__ == "__main__":
    """
    ToDo
    1, add boolean values to the second step in order to control the output. 
    """
    args = get_user_arguments()
    
    #Read in files:
    settings_dict = load_json(args.settings_filename)
    protein_groups_dataframe = read_in_protein_groups_file(args.filename)

    #process group proteins file by filtering columns and rows:
    protein_groups_dataframe = filter_dataframe_columns(protein_groups_dataframe, settings_dict)
    protein_groups_dataframe, filtered_groups_dataframe = filter_dataframe_rows(protein_groups_dataframe, settings_dict)
    protein_groups_dataframe = fetch_identifiers(protein_groups_dataframe)

    #fetch annotation for uniprot identifiers:
    protein_data_list = fetch_uniprot_annotation(protein_groups_dataframe["identifier"])
    with open("example_proteins_group_data.json", 'r') as inputfile:
        protein_data_list = json.load(inputfile)
    protein_groups_dataframe = append_uniprot_data_to_dataframe(protein_groups_dataframe, protein_data_list)
    
    #map proteins to mitocarta:
    
    #apply hierarchical cluster analysis

    #write away dataframe to an excel file:
    



    
    # processed.to_csv('processed_out.tsv', sep='\t')
