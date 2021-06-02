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
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as spd
import openpyxl
import xlsxwriter

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
        protein_groups_dataframe = pd.read_csv(filename, sep='\t', index_col = "Majority protein IDs", low_memory=False)
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
def read_in_excel_file(filename, sheet_name):
    """
    input:
    filename = string
    output:
    dataframe = pd.DataFrame()
    """
    print(f"Read in {filename} for sheet {sheet_name}.")
    dataframe = pd.read_excel(io=filename, sheet_name=sheet_name, index_col=None, na_filter=False)
    print(f"Succesfully read in\'{filename}\' for sheet {sheet_name} and created a dataframe.")
    print("-"*40)
    return dataframe

def filter_dataframe_columns(protein_groups_dataframe, settings_dict):
    """
    input:
    protein_groups_dataframe = pd.DataFrame
    settings_dict = dict{parameter: [values]}
    output:
    protein_groups_dataframe = pd.DataFrame
    non_selected_dataframe = pd.DataFrame
    """
    print("Start removing unwanted columns from the dataframe:")
    all_dataframe_columns = list(protein_groups_dataframe.columns)
    selected_dataframe_columns = []
    non_selected_dataframe_columns = []
    for key_word, method in zip(["EXACT_MATCHES", "CONTAINS"], ["exact_matches", "contains"]):
        selected_columns = select_columns(all_dataframe_columns, settings_dict[key_word], method)
        selected_dataframe_columns += selected_columns
    non_selected_dataframe_columns = list(set(protein_groups_dataframe.columns).difference(selected_dataframe_columns))
    non_selected_dataframe = protein_groups_dataframe[non_selected_dataframe_columns]
    protein_groups_dataframe = protein_groups_dataframe[selected_dataframe_columns]
    
    number_of_rows = protein_groups_dataframe.shape[0]
    number_of_columns = protein_groups_dataframe.shape[1]
    print(f"Finished removing unwanted columns from the dataframe. The filtered dataframe has {number_of_rows} rows and {number_of_columns} columns")
    print("-"*40)
    return protein_groups_dataframe, non_selected_dataframe

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
    non_applying_proteins = list, list with rows not to keep
    applying_proteins = list, list of rows to keep
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
        print(f'Parsing a fasta header to get the identifier did not work for header {fasta_header}')
        return np.nan

def fetch_uniprot_annotation(identifiers, settings_dict):
    """
    input:
    identifiers = pd.Series
    output:
    protein_data_list = list, list of protein_data_dict
    """
    print("Start fetching data from uniprot. After each batch {sleep_time} seconds pass before a new barch is queried to uniprot.\nThis safety feature is put in place to prevent being blacklisted."
          .format(sleep_time=settings_dict["request_idle_time"]))
    protein_data_dict = {}
    function_dict = construct_function_dict(settings_dict)

    #split identifiers into multiple sub arrays of length batch_length:
    identifier_batches = np.split(identifiers, range(settings_dict["batch_amount"],len(identifiers), settings_dict["batch_amount"]))
    for n_batch, identifiers_batch in enumerate(identifier_batches):
        print(f"Start fetching uniprot data for batch number {n_batch} of the total {len(identifier_batches)} batches")
        requestURL = settings_dict["uniprot_base_url"]+settings_dict["uniprot_request_url"]+",".join(identifiers_batch)
        request = requests.get(requestURL, headers={"Accept" : "application/json"})
        if not request.ok:
            print(f"Something went wrong with query {n_batch} while querying the uniprot database. The {len(identifiers)} proteins of this batch will be ignored")
            request.raise_for_status()
            for identifier in identifiers_batch:
                protein_data_dict[identifier] = {"gene_name":np.nan, "protein_name":np.nan, "organism_name":np.nan, "hyperlink":np.nan, "cell_compartment":np.nan, "string_linkout":np.nan}
        else:
            print(f"Succesfully fetched uniprot data for batch {n_batch}:")
            protein_data_dict = update_protein_data_dict(request.json(), identifiers_batch, function_dict, protein_data_dict, settings_dict)
            protein_data_dict = add_uniprot_hyperlink(protein_data_dict, settings_dict, identifiers_batch)
            protein_data_dict = add_string_linkout(protein_data_dict, settings_dict, identifiers_batch)
        #Let the program sleep for a bit else uniprot is going to be overloaded and I get a problem.
        time.sleep(settings_dict["request_idle_time"])
    return protein_data_dict
    print("Finished fetching data from uniprot")

def construct_function_dict(settings_dict):
    """
    Based on the user input apply which functions should be used.
    input:
    settings_dict = dict, dict containing settings
    output:
    function_dict = dict{function_name : function}
    """
    function_dict = {}
    for function_setting, function_value in settings_dict["uniprot_options"].items():
        if function_setting == "get_gene_name" and function_value == True:
            function_dict.update({"gene_name":get_uniprot_gene_name})
        if function_setting == "get_protein_name" and function_value == True:
            function_dict.update({"protein_name":get_protein_name})
        if function_setting == "get_organism_name" and function_value == True:
            function_dict.update({"organism_name":get_organism_name})
        if function_setting == "get_cell_compartment" and function_value == True:
            function_dict.update({"cell_compartment":get_cell_compartment})
    return function_dict

def update_protein_data_dict(uniprot_output_list, identifiers, function_dict, protein_data_dict, settings_dict):
    """
    This function is meant to loop over the identifiers and update a dict with values per protein.
    input:
    identifiers = list, list of uniprot identifiers
    function_dict = dict{function_name(string) : function(function)}
    protein_data_dict = dict{protein : {gene_name, protein_name, organism_name, cell_comparment}}
    output:
    protein_data_dict = dict{protein : {gene_name, protein_name, organism_name, cell_comparment}}
    """
    for request_number in range(len(identifiers)):
        identifier = identifiers[request_number]
        uniprot_data_dict = get_matching_uniprot_query(uniprot_output_list, identifier)
        if not identifier in protein_data_dict.keys():
            protein_data_dict[identifier] = {}
        for function_name, function in function_dict.items():
            if None == uniprot_data_dict:
                protein_data_dict[identifier].update({function_name : np.nan})
            else:
                try:
                    protein_data_dict[identifier].update({function_name : function(uniprot_data_dict)})
                except KeyError as key_error:
                    print(f"A key error occured for protein {identifier}, the field {function_name} will be ignored for this protein. Please see error message below: ")
                    print(key_error)
                    protein_data_dict[identifier].update({function_name : np.nan})
                except IndexError as index_error:
                    print(f"A key index occured for protein {identifier}, the field {function_name} will be ignored for this protein. Please see error message below: ")
                    print(index_error)
                    protein_data_dict[identifier].update({function_name : np.nan})
                except Exception as exception:
                    print(f"A general exception occured for protein {identifier}, the field {function_name} will be ignored for this protein. Please see error message below: ")
                    print(exception)
                    protein_data_dict[identifier].update({function_name : np.nan})
    return protein_data_dict

def get_matching_uniprot_query(uniprot_output_list, identifier):
    """
    input:
    uniprot_output_list = list, list of uniprot entries
    identifier = string
    output:
    uniprot_data_dict
    """
    for uniprot_data_dict in uniprot_output_list:
        if "accession" in uniprot_data_dict:
            if uniprot_data_dict["accession"] == identifier:
                return uniprot_data_dict
    return None

def add_uniprot_hyperlink(protein_data_dict, settings_dict, identifiers):
    """
    input:
    protein_data_dict = dict{identifier : {gene_name, protein_name, organism_name, cell_comparment}}
    settings_dict = dict, dictionary containing settings
    identifiers = list, list of strings
    output:
    protein_data_dict = dict{identifier : {gene_name, protein_name, organism_name, uniprot_hyperlink, cell_comparment}}
    """
    if settings_dict["uniprot_options"]["get_uniprot_hyperlink"] == True:
        for identifier in identifiers:
            protein_data_dict[identifier].update({"uniprot_hyperlink":settings_dict["uniprot_protein_base_url"]+identifier})
    return protein_data_dict

def add_string_linkout(protein_data_dict, settings_dict, identifiers):
    """
    input:
    protein_data_dict = dict{identifier : {gene_name, protein_name, organism_name, cell_comparment, uniprot_hyperlink}}
    settings_dict = dict, dictionary containing settings
    identifiers = list, list of strings
    output:
    protein_data_dict = dict{identifier : {gene_name, protein_name, organism_name, cell_comparment, uniprot_hyperlink, string_linkout}}
    """
    if settings_dict["uniprot_options"]["get_string_linkout"] == True:
        string_linkout_dict = get_string_linkout(identifiers, settings_dict["string_linkout_parameters"])
        for identifier in identifiers:
            protein_data_dict[identifier].update({"string_linkout":string_linkout_dict[identifier]})
            print(protein_data_dict[identifier])
            print("-"*40)
    return protein_data_dict 

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

def get_string_linkout(identifiers, settings_dict):
    """
    Use the uniprot identifier mapping service to get a linkout to the string database per uniprot identifier.
    input:
    identifiers = list, list of uniprot identifiers
    settings_dict = dict, dictionary with the specific parameters for the string linkout. 
    output:
    string_linkout_dict = dict{identifier : string_linkout}
    """
    string_linkout_dict = {}

    formatted_identifier_string = process_uniprot_identifier_input(identifiers, settings_dict["regex_pattern"])
    parameters = {'from': 'ACC+ID', 'to': 'STRING_ID', 'format': 'tab', 'query': formatted_identifier_string}

    data = urllib.parse.urlencode(parameters)
    data = data.encode('utf-8')
    request = urllib.request.Request(settings_dict["uniprot_mapping_service_url"], data)
    with urllib.request.urlopen(request) as uniprot_mapping_request:
       response = uniprot_mapping_request.read()
    uniprot_mapped_proteins_dict = process_uniprot_mapping_service_output(response.decode('utf-8'))
    for identifier in identifiers:
        if not identifier in uniprot_mapped_proteins_dict.keys():
            string_linkout_dict[identifier] = np.nan
        else:
            string_linkout_dict[identifier] = settings_dict["string_base_url"]+uniprot_mapped_proteins_dict[identifier]
    return string_linkout_dict

def process_uniprot_identifier_input(identifiers, regex_pattern):
    """
    input:
    identifiers = list, list of uniprot identifiers
    regex_pattern = string, pattern that is used to detect "-[0-9]" suffixes behind proteins. 
    output:
    formatted_identifier_string = string, should be "identifier identifier etc."
    """
    formatted_identifier_string = ""
    for identifier in identifiers:
        #A regex function is used to identify proteins with for example: "-2" as suffix and removes the suffix.
        if None != re.search(regex_pattern, identifier): identifier = identifier[:-2]
        #each identifier must be spaced besides each other according to the example python3 query from: https://www.uniprot.org/help/api_idmapping, accessed at 31-05-2021
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

def append_uniprot_data_to_dataframe(protein_groups_dataframe, protein_data_dict, uniprot_options_dict):
    """
    input:
    protein_groups_dataframe = pd.DataFrame
    protein_data_dict = dict{identifier: {"gene_name":"", "protein_name":"", "organism_name":"", "hyperlink":"", "cell_compartment":np.nan, "string_linkout":""}}
    uniprot_options_dict = dict, dictionary containing information about which information is gained from the uniprot server
    output:
    protein_groups_dataframe = pd.DataFrame
    """
    print("Start adding the uniprot columns to the existing dataframe")
    uniprot_column_names = get_column_names(uniprot_options_dict)
    for uniprot_column_name in uniprot_column_names:
        uniprot_column_values = get_uniprot_column_values(protein_groups_dataframe["identifier"], uniprot_column_name, protein_data_dict)
        protein_groups_dataframe[uniprot_column_name] = uniprot_column_values
    print("Finished adding the uniprot columns to the existing dataframe")
    print("-"*40)

    return protein_groups_dataframe

def get_column_names(uniprot_options_dict):
    """
    input:
    uniprot_options_dict = dict, dictionary containing information about which information is gained from the uniprot server
    output:
    uniprot_column_names = list, list of column names
    """
    uniprot_column_names = []
    for option_name, option_value in uniprot_options_dict.items():
        if option_name == "get_gene_name" and option_value == True:
            uniprot_column_names.append("gene_name")
        if option_name == "get_protein_name" and option_value == True:
            uniprot_column_names.append("protein_name")
        if option_name == "get_organism_name" and option_value == True:
            uniprot_column_names.append("organism_name")
        if option_name == "get_uniprot_hyperlink" and option_value == True:
            uniprot_column_names.append("uniprot_hyperlink")
        if option_name == "get_cell_compartment" and option_value == True:
            uniprot_column_names.append("cell_compartment")
        if option_name == "get_string_linkout" and option_value == True:
            uniprot_column_names.append("string_linkout")
    return uniprot_column_names

def get_uniprot_column_values(identifiers, uniprot_column_name, protein_data_dict):
    """
    input:
    identifiers = pd.Series
    uniprot_column_name = string
    protein_data_list = dict{identifier: {"gene_name":"", "protein_name":"", "organism_name":"", "uniprot_hyperlink":"", "cell_compartment":np.nan, "string_linkout":""}}
    output:
    uniprot_column_values = list
    """
    uniprot_column_values = []
    for identifier in identifiers:
        uniprot_column_value = protein_data_dict[identifier][uniprot_column_name]
        uniprot_column_values.append(uniprot_column_value)
    return uniprot_column_values
def evaluate_uniprot_settings(uniprot_options):
    """
    Evaluate whether the user didn't accidentally put all settings to 0 while the step is 1. 
    input:
    uniprot_options = dict{"option":boolean}
    output:
    boolean
    """
    #any() returns True whenever one element is True and False if all elements are False.
    return any(uniprot_options.values())

def is_protein_in_mitocarta(protein_groups_dataframe, mitocarta_species_dataframe, species_name, new_column_name):
    """
    input:
    protein_groups_dataframe = pd.DataFrame()
    mitocarta_species_dataframe = pd.DataFrame()
    species_name = string
    new_column_name = string
    output:
    protein_groups_dataframe = pd.DataFrame()
    """
    print(f"Start finding proteins found in the {species_name} mitocarta dataset")
    protein_groups_dataframe["gene_name"].fillna('', inplace=True)
    mitocarta_species_dataframe["Synonyms"].fillna('', inplace=True)
    #In the mitocarta_mouse_data there are several synonyms which are floats, in order to remove these the following line of code is needed:
    mitocarta_species_dataframe["Synonyms"].where(cond=[True if type(x) == str else False for x in mitocarta_species_dataframe["Synonyms"].str.split("|")], other="-", inplace=True)
    
    organism_rows = protein_groups_dataframe["organism_name"] == species_name
    symbol_rows = protein_groups_dataframe["gene_name"].isin([symbol for symbol in protein_groups_dataframe["gene_name"] if mitocarta_species_dataframe["Symbol"].str.contains(symbol).any() == True
                                                              or [symbol in series_list for series_list in mitocarta_species_dataframe["Synonyms"].str.split("|")]])
    
    presency_index = protein_groups_dataframe.index[(organism_rows) & (symbol_rows)]
    protein_groups_dataframe[new_column_name] = [1.0 if index in presency_index else 0.0 for index in range(protein_groups_dataframe.shape[0])]
    
    print(f"Finished finding proteins found in the {species_name} mitocarta dataset")
    print("-"*40)
    return protein_groups_dataframe


def cluster_reorder(sample_specific_dataframe, method = 'average', metric = 'correlation'):
    """
    clusters proteins with hierarchical clustering. Determines optimal order, resulting in minimal distance between adjacent leaves
    Returns: dict with mapping of identifier to position in ordered protein list
    input:
    sample_specific_dataframe = pd.DataFrame()
    method = string
    metric = string
    output:
    order = dict{protein_identifier : ordered_index}
    clustered = 
    """
    condensed_distance_matrix = spd.pdist(np.array(sample_specific_dataframe))
    
    clustered = sch.linkage(condensed_distance_matrix, method = method, metric = metric, optimal_ordering = True)
    dendrogram = sch.dendrogram(clustered,labels = sample_specific_dataframe.index.values, orientation = 'right')
    ordered_index = sample_specific_dataframe.iloc[dendrogram['leaves'],:].index.values
    order = {label:ix for ix,label in enumerate(ordered_index)}
    return order, clustered

def dump_data_to_excel(protein_groups_dataframe, non_selected_dataframe, settings_dict):
    """
    The last part of this script, dump the complexome profiling data into an excel file.
    input:
    protein_groups_dataframe = pd.DataFrame()
    non_selected_dataframe = pd.DataFrame()
    settings_dict = dict, dictionary with parameters for this function
    output:
    None
    """
    print("Start writing away the data to file {file_name}".format(file_name=settings_dict["make_excel_file_step"]["excel_file_name"]))
    writer = pd.ExcelWriter(settings_dict["make_excel_file_step"]["excel_file_name"], engine='xlsxwriter', mode="w")
    workbook = writer.book

    complexome_profiling_dataframe, data_dataframe = split_dataframe(protein_groups_dataframe)
    data_dataframe['uniprot_hyperlink'].transform(lambda x: make_hyperlink(x))
    data_dataframe['string_linkout'].transform(lambda x: make_hyperlink(x))
    complexome_profiling_dataframe = order_complexome_profiling_dataframe(complexome_profiling_dataframe)

    non_selected_dataframe.to_excel(writer, sheet_name = "non selected columns", index=False)
    data_dataframe.to_excel(writer, sheet_name = 'data', index=False)
    complexome_profiling_dataframe.to_excel(writer, sheet_name= "complexome profiles", index=False)
    worksheet = writer.sheets['complexome profiles']
    
    positions = get_sample_positions(complexome_profiling_dataframe.columns.tolist())
    apply_conditional_formating_per_sample(complexome_profiling_dataframe, positions, writer, worksheet, workbook)
        
    writer.save()
    print("Finished writing away the data to file {file_name}".format(file_name=settings_dict["make_excel_file_step"]["excel_file_name"]))

def make_hyperlink(hyperlink):
    """
    input:
    hyperlink = string
    output:
    hyperlink = string
    """
    if "" != hyperlink and pd.isnull(hyperlink) == False:
        identifier = hyperlink.split("/")[-1:][0]
        return f'=HYPERLINK("{hyperlink}", "{identifier}")'
    else:
        return ""

def split_dataframe(protein_groups_dataframe):
    """
    Split the dataframe into columns containing complexome profiling information and all the other data
    input:
    protein_groups_dataframe = pd.DataFrame()
    output:
    complexome_profiling_dataframe = pd.DataFrame()
    data_dataframe = pd.DataFrame()
    """
    applying_column_names = []
    non_applying_column_names = []
    complexome_profiling_sample_filter = "iBAQ"
    clustered_sample_filter = "clustered"
    
    for column_name in protein_groups_dataframe.columns:
        if complexome_profiling_sample_filter in column_name:
            applying_column_names.append(column_name)
        elif clustered_sample_filter in column_name:
            applying_column_names.append(column_name)
        else:
            non_applying_column_names.append(column_name)
    complexome_profiling_dataframe = protein_groups_dataframe[applying_column_names]
    data_dataframe = protein_groups_dataframe[non_applying_column_names]
    return complexome_profiling_dataframe, data_dataframe
def order_complexome_profiling_dataframe(complexome_profiling_dataframe):
    """
    Order the complexome profiling dataframe so the pattern emerges: sample iBAQ columns - sample clustered column - etc. - global clustered column
    input:
    complexome_profiling_dataframe = pd.DataFrame()
    output:
    complexome_profiling_dataframe = pd.DataFrame()
    """
    ordered_columns = []
    cluster_column = ""
    cluster_key_word = "clustered"
    global_key_word = "global"

    sample_names = list(set([i[1].split("_")[0] for i in complexome_profiling_dataframe.columns.str.split(" ") if i[0] == "iBAQ" and len(i) > 1]))
    for sample_name in sample_names:
        sample_specific_columns = protein_groups_dataframe.columns[protein_groups_dataframe.columns.to_series().str.contains(sample_name)]
        for column_value in sample_specific_columns:
            if cluster_key_word in column_value:
                cluster_column = column_value
            else:
                ordered_columns.append(column_value)
        if cluster_column != "": ordered_columns.append(cluster_column)
    ordered_columns.extend([x for x in complexome_profiling_dataframe.columns[complexome_profiling_dataframe.columns.str.contains(global_key_word)]])
    complexome_profiling_dataframe = complexome_profiling_dataframe.reindex(columns=ordered_columns)
    return complexome_profiling_dataframe

def get_sample_positions(column_names):
    """
    Per sample, get the start and end column positions
    input:
    column_names = pd.Index(), iterable with the column names of the complexome profiling dataframe
    output:
    positions = list, [[start_position, stop_position], [start_position, stop_position]]
    """
    positions = []
    start_position = 0
    end_position = 0
    for column_name in column_names:
        if "iBAQ" in column_name:
            if start_position == 0:
                start_position = column_names.index(column_name)
            end_position = column_names.index(column_name)
        else:
            if end_position != 0:
                positions.append([start_position, end_position])
                start_position, end_position = 0, 0 
    return positions

def apply_conditional_formating_per_sample(complexome_profiling_dataframe, positions, writer, worksheet, workbook):
    """
    Per sample, apply conditional formatting using the positions
    input:
    complexome_profiling_dataframe = pd.DataFrame()
    positions = list, [[start_position, end_position],[start_position, end_position]]
    writer = pd.ExcelWriter object
    worksheet = writer.sheets object
    workbook = writer.book object
    output:
    complexome_profiling_dataframe = pd.DataFrame()
    """
    for start_position, end_position in positions:
        worksheet.set_column(start_position,end_position, 0.5)
        apply_cond_format(complexome_profiling_dataframe,start_position,end_position,writer,worksheet,workbook)

def apply_cond_format(dataframe,startcol,endcol,writer,worksheet,workbook):
    """
    Apply conditional formatting to relevant columns
    input:
    dataframe = pd.DataFrame()
    startcol = int
    endcol = int
    writer = pd.ExcelWriter
    worksheet = writer.worksheets
    workbook = writer.workbook
    output:
    None
    """
    row_numbers = dataframe.shape[0]
    worksheet.conditional_format(1,startcol,row_numbers, endcol,
                                              {'type'      : '3_color_scale',
                                               'min_color' : "#000000",
                                               'mid_type'  :  'percentile',
                                               'mid_value' : 95,
                                               'mid_color' : "#FFFF00",
                                               'max_color' : "#FF0000"})
    
if __name__ == "__main__":
    """
    ToDo
    each step should be put in a separate file
    the link should be transformed to something that can be clicked upon.
    """
    args = get_user_arguments()
    
    #Read in files:
    settings_dict = load_json(args.settings_filename)
    protein_groups_dataframe = read_in_protein_groups_file(args.filename)

    #process group proteins file by filtering columns and rows:
    if settings_dict["steps_dict"]["filtering_step"] == True:
        protein_groups_dataframe, non_selected_dataframe = filter_dataframe_columns(protein_groups_dataframe, settings_dict["filtering_step"])
        protein_groups_dataframe, filtered_groups_dataframe = filter_dataframe_rows(protein_groups_dataframe, settings_dict["filtering_step"])
        protein_groups_dataframe = fetch_identifiers(protein_groups_dataframe)
        #Reset the index, through this way new columns can be added to the dataframe 
        protein_groups_dataframe.reset_index(inplace=True)
    else:
        print("The proteins will not be filtered")
        print("-"*50)
    #fetch annotation for uniprot identifiers:
    if settings_dict["steps_dict"]["uniprot_step"] == True and evaluate_uniprot_settings(settings_dict["uniprot_step"]["uniprot_options"]) == True:
##        protein_data_dict = fetch_uniprot_annotation(protein_groups_dataframe["identifier"], settings_dict["uniprot_step"])
        with open("example_proteins_group_data.json", 'r') as inputfile:
            protein_data_dict = json.load(inputfile)
        protein_groups_dataframe = append_uniprot_data_to_dataframe(protein_groups_dataframe, protein_data_dict, settings_dict["uniprot_step"]["uniprot_options"])
    else:
        print("Uniprot will not be queried for information due to the step being disabled or none of the fields are set to True.")
        print("-"*50)
    #map proteins to mitocarta:
    if settings_dict["steps_dict"]["mitocarta_step"] == True and settings_dict["steps_dict"]["uniprot_step"] == True and \
       settings_dict["uniprot_step"]["uniprot_options"]["get_gene_name"] == True and settings_dict["uniprot_step"]["uniprot_options"]["get_organism_name"] == True:
        mitocarta_mouse_dataframe = read_in_excel_file(settings_dict["mitocarta_step"]["mitocarta_mouse_ftp_link"], settings_dict["mitocarta_step"]["mouse_sheet_name"])
        mitocarta_human_dataframe = read_in_excel_file(settings_dict["mitocarta_step"]["mitocarta_human_ftp_link"], settings_dict["mitocarta_step"]["human_sheet_name"])
        protein_groups_dataframe = is_protein_in_mitocarta(protein_groups_dataframe, mitocarta_mouse_dataframe, "Mus musculus", "mitocarta_mouse_presency")
        protein_groups_dataframe = is_protein_in_mitocarta(protein_groups_dataframe, mitocarta_human_dataframe, "Homo sapiens", "mitocarta_human_presency")
        
    else:
        if settings_dict["steps_dict"]["mitocarta_step"] == False and settings_dict["steps_dict"]["uniprot_step"] == True:
            print("The final dataframe will not be enriched with information from mitocarta because the mitocarta step is disabled.")
        elif settings_dict["steps_dict"]["uniprot_step"] == True and settings_dict["uniprot_step"]["uniprot_options"]["get_gene_name"] == False or settings_dict["uniprot_step"]["uniprot_options"]["get_organism_name"] == False:
            print("Uniprot has been queried but the gene name or the organism name has not been retrieved. Set the organism and gene_name to 1 in order to get the information from mitocarta")
        elif settings_dict["steps_dict"]["mitocarta_step"] == True and settings_dict["steps_dict"]["uniprot_step"] == False:
            print("The final dataframe will not be enriched with information from mitocarta because the uniprot step is disabled.\nUniprot is needed to get per protein the gene symbol, the key to find values in uniprot.")
        elif settings_dict["steps_dict"]["mitocarta_step"] == False and settings_dict["steps_dict"]["uniprot_step"] == False:
            print("The final dataframe will not be enriched with information from mitocarta because the mitocarta and uniprot steps are disabled")
        print("-"*50)
    
    #apply hierarchical cluster analysis
    #get the different sample names:
    if settings_dict["steps_dict"]["clustering_step"] == True:
        sample_names = list(set([i[1].split("_")[0] for i in protein_groups_dataframe.columns.str.split(" ") if i[0] == "iBAQ" and len(i) > 1]))
        for sample_name in sample_names:
            print(f"Start hierarchical clustering for sample {sample_name}")
            sample_specific_dataframe = pd.DataFrame(protein_groups_dataframe[protein_groups_dataframe.columns[protein_groups_dataframe.columns.to_series().str.contains(sample_name)]], dtype="float64")
            order_mapping, clustered = cluster_reorder(sample_specific_dataframe, settings_dict["clustering_step"]["method"], settings_dict["clustering_step"]["metric"])
            protein_groups_dataframe[f'sample_{sample_name}_clustered'] = pd.Series(order_mapping)
            print(f"Finished hierarchical clustering for sample {sample_name}")
            print("-"*40)
        print("Start hierarchical clustering for all samples") 
        global_order_mapping, global_clustered = cluster_reorder(protein_groups_dataframe[protein_groups_dataframe.columns[protein_groups_dataframe.columns.to_series().str.contains("iBAQ")]])
        protein_groups_dataframe['global_clustered'] = pd.Series(global_order_mapping)
        print("Finished hierarchical clustering for all samples")
        print("-"*40)
    else:
        print("Hierarchical clustering will not be applied for the available clusters because the clustering_step has been disabled")

    #for testing purposes create a csv with the current dataframe:
    #protein_groups_dataframe.to_csv("excel_input.csv", sep=",", index=False) 
   
    #write away dataframe to an excel file:
    if settings_dict["steps_dict"]["make_excel_file_step"] == True:
        #read in an example dataframe:
        protein_groups_dataframe = pd.read_csv("excel_input.csv", sep=',', index_col = None, low_memory=False)
        dump_data_to_excel(protein_groups_dataframe, non_selected_dataframe, settings_dict)
        
    else:
        print("The data will not be written away to an excel file because the make_excel_file_step has been disabled")
    



    
    # processed.to_csv('processed_out.tsv', sep='\t')
