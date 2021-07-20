#Import standard python libraries:
import logging
import os.path
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
    parser = argparse.ArgumentParser(description="Process maxquant files, please use complete file paths because else errors occur",
                                     usage="USAGE: [process_maxquant.exe] --maxquant_file_path [maxquant_file_path] --settings_file_path [settings_file_path]")
    parser.add_argument("--maxquant_file_path", help="The file path to the maxquant file", type=str)
    parser.add_argument("--settings_file_path", help="The path to the settings dict, should be a json file", type=str)
    #Test files which are more usefull to have like this instead of having to type them in over and over again:
    #development_options = [["20190417_Methanop_DDMfull_BLZ2_proteinGroups.txt", "maxquant_settings.json"], ["20200123_Scer_Daniel_NB40_BY4742_DFM6_DFM9_proteinGroups.txt", "maxquant_settings.json"],["--maxquant_file_path", "C:\\Users\\ariel\\Documents\\Bio_informatica_bestanden\\maxquant_project\\pythonProject\\20190115_HEKwt_and_MICS1ko_proteinGroups.txt", "--settings_file_path", "C:\\Users\\ariel\\Documents\\Bio_informatica_bestanden\\maxquant_project\\pythonProject\\maxquant_settings.json"]]
    #development_arguments = development_options[2]

    #args = parser.parse_args(development_arguments)
    args = parser.parse_args()
    #Detect whether none arguments were submitted:
    if not len(sys.argv) > 1:
        parser.print_help()
        sys.exit()
    return args


def check_user_input(gui_object, settings_file_path, maxquant_file_path):
    """
    input:
    gui_object = PyQt5 Qapplication
    settings_file_path = string
    maxquant_file_path = string
    output:
    boolean, True == both file_paths are valid and False == both or one file is not valid.
    """
    expected_settings_file_extension = ".json"
    expected_maxquant_file_extension = ".txt"
    if is_user_input_valid(gui_object, settings_file_path, maxquant_file_path) == False:
        return False
    if is_extension_valid(settings_file_path, expected_settings_file_extension) == False:
        gui_object.report_error(f"The settings file has a different file extension than the expected {expected_settings_file_extension}\nAre you sure it is the correct file?")
        return False
    if is_extension_valid(maxquant_file_path, expected_maxquant_file_extension) == False:
        gui_object.report_error(f"The maxquant file has a different file extension than the expected {expected_maxquant_file_extension}\nAre you sure it is the correct file?")
        return False
    return True


def is_extension_valid(file_path, prefered_extension):
    file_extension = os.path.splitext(file_path)[1]
    if file_extension == prefered_extension:
        return True
    else:
        return False


def is_user_input_valid(gui_object, settings_file_path, maxquant_file_path):
    """
    Check whether the user submitted two file paths
    input:
    gui_object = PyQt5 Qapplication
    settings_file_path = string
    maxquant_file_path = string
    output:
    boolean, True == both file_paths are not empty and False == both or one file path is empty
    """
    if "" == settings_file_path and "" == maxquant_file_path:
        gui_object.report_error("No arguments were submitted. Select the settings file and a maxquant file in order to let the program work")
        return False
    elif "" == settings_file_path:
        gui_object.report_error("The settings file path was not selected.")
        return False
    elif "" == maxquant_file_path:
        gui_object.report_error("The maxquant file path was not selected.")
        return False
    else:
        return True


def load_json(gui_object, json_filepath):
    """
    input:
    gui_object = PyQt5 Qapplication
    json_filepath = string
    output:
    json_object = dict{}
    parameters:
    json_file = IO text object
    """
    gui_object.report_status("Start reading the settings file")
    try:
        with open(json_filepath) as json_file:
            json_object = json.load(json_file)
    except FileNotFoundError as file_not_found_error:
        log_error(gui_object, f"The settings file at {json_filepath} was not found, please try again", file_not_found_error)
        return None, False
    except Exception as error:
        log_error(gui_object, f"A unhandled error has occurred while reading {json_filepath}", error)
        return None, False
    gui_object.report_status("Finished reading in the settings file")
    return json_object, True


def read_in_protein_groups_file(gui_object, filename):
    """
    input:
    gui_object = PyQt5 Qapplication
    string, path of proteingroups file
    output:
    protein_groups_dataframe = pd.DataFrame
    """
    gui_object.report_status(f"Start reading in \'{filename}\'")
    try:
        protein_groups_dataframe = pd.read_csv(filename, sep='\t', index_col="Majority protein IDs", low_memory=False)
    except FileNotFoundError as file_not_found_error:
        log_error(gui_object, f"The maxquant file \'{filename}\' was not found, please try again.", file_not_found_error)
        return None, False
    except IOError as io_error:
        log_error(gui_object, "An IO(input/output) error occurred while reading in the maxquant file.", io_error)
        return None, False
    except Exception as error:
        log_error(gui_object, f"Something went wrong while reading in \'{filename}\', please see the error below:", error)
        return None,False
    logging.info(f"Successfully read in \'{filename}\' and created a dataframe. The dataframe has {len(protein_groups_dataframe)}"
                 f" rows and {len(protein_groups_dataframe.columns)} columns")
    gui_object.report_status("Finished reading in the maxquant file")
    return protein_groups_dataframe, True

def log_error(gui_object, error_message, exception):
    """"
    This function is meant to generalize the code under each exception clause
    input:
    gui_object = PyQt5, Qapplication
    error_message = string, user defined message
    exception = string, the error python detected
    output:
    None
    """
    logging.error(error_message)
    logging.error(exception)
    gui_object.report_error(f"{error_message}\n{exception}")

def read_in_excel_file(gui_object, filename, sheet_name):
    """
    input:
    gui_object = PyQt5, Qapplication
    filename = string
    sheet_name = string
    output:
    dataframe = pd.DataFrame()
    """
    gui_object.report_status(f"Start reading in the excel file {filename} for sheet {sheet_name}.")
    try:
        dataframe = pd.read_excel(io=filename, sheet_name=sheet_name, index_col=None, na_filter=False)
    except FileNotFoundError as file_not_found_error:
        log_error(gui_object, f"The excel file \'{filename}\' has not been found at the specified location", file_not_found_error)
    except Exception as error:
        log_error(gui_object, f"Something went wrong while reading in the excel file \'{filename}\'", error)
    gui_object.report_status(f"Successfully created a dataframe from excel file \'{filename}\' with sheet {sheet_name}.")
    return dataframe

def filter_dataframe_columns(gui_object, protein_groups_dataframe, settings_dict):
    """
    input:
    gui_object = PyQt5, Qapplication
    protein_groups_dataframe = pd.DataFrame
    settings_dict = dict{parameter: [values]}
    output:
    protein_groups_dataframe = pd.DataFrame
    non_selected_dataframe = pd.DataFrame
    """
    logging.info("Start removing unwanted columns from the dataframe.")
    all_dataframe_columns = list(protein_groups_dataframe.columns)
    selected_dataframe_columns = []
    try:
        for key_word, method in zip(["EXACT_MATCHES", "CONTAINS"], ["exact_matches", "contains"]):
            selected_columns = select_columns(all_dataframe_columns, settings_dict[key_word], method)
            selected_dataframe_columns.extend(selected_columns)
        protein_groups_dataframe = protein_groups_dataframe[selected_dataframe_columns]
    except IndexError as index_error:
        log_error(gui_object, f"While filtering the columns the key word {key_word} was not found in the settings", index_error)
    logging.info(f"Finished removing unwanted columns from the dataframe.\n"
                 f"The filtered dataframe has {len(protein_groups_dataframe)} rows and {len(protein_groups_dataframe.columns)} columns")
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
            logging.debug(f"In function 'select_columns' the variable name {method} has been passed but isn't a valid method name. Change this to a valid method name")
            return
        keep_columns += keep_column

    return keep_columns

def filter_dataframe_rows(protein_groups_dataframe, settings_dict):
    """
    input:
    protein_groups_dataframe = pd.DataFrame
    settings_dict = dict{parameter: [values]}
    output:
    protein_groups_dataframe = pd.DataFrame, proteins applying to the condition
    filtered_groups_dataframe = pd.DataFrame, proteins not applying to the condition but should still be retained
    """
    logging.info("Start filtering out unwanted proteins and drop proteins containing NA values: ")
    non_applying_rows, applying_rows = select_proteins(protein_groups_dataframe.index, settings_dict["PROTEIN_FILTERS"])
    filtered_groups_dataframe = protein_groups_dataframe.loc[non_applying_rows]
    protein_groups_dataframe = protein_groups_dataframe.loc[applying_rows]

    logging.info(f"Finished filtering out proteins not applying to the set protein_filters. "
                 f"The filtered dataframe has {len(protein_groups_dataframe)} rows and {len(protein_groups_dataframe.columns)} columns")
    return protein_groups_dataframe, filtered_groups_dataframe

def select_proteins(row_labels_index, protein_filters):
    """
    selects proteins to keep using index(should be Majority_proteins_id) and protein_filter, but keep proteins which do not apply to the filter.
    input:
    row_labels_index = pd.DataFrame().index, the row labels(in this case 'Majority protein IDs') of the dataframe
    protein_filters = list, list of filters 
    output:
    non_applying_proteins = list, list with rows not to keep
    applying_proteins = list, list of rows to keep
    """
    non_applying_proteins = []
    applying_proteins = []
    for row_label in row_labels_index:
        for protein_filter in protein_filters:
            if protein_filter in row_label and not row_label in non_applying_proteins:
                non_applying_proteins.append(row_label)
        if not row_label in non_applying_proteins:
            applying_proteins.append(row_label)

    return non_applying_proteins, applying_proteins

def fetch_identifiers(protein_groups_dataframe):
    """
    input:
    protein_groups_dataframe = pd.DataFrame
    output:
    identifiers = pd.Series
    """
    logging.info("Start fetching identifiers based on the Fasta headers")
    identifiers = protein_groups_dataframe['Fasta headers'].apply(parse_identifier)
    logging.info("Finished fetching identifiers from the fasta headers")
    return identifiers

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
        logging.error(f'Parsing a fasta header to get the identifier did not work for header {fasta_header}\nPython error: {error}')
        return np.nan


def are_identifiers_not_available(identifiers):
    """
    Are all identifiers np.nan or are they strings?
    input:
    identifiers = pd.Series()
    output:
    boolean, True == all identifiers are NA, and False == some identifiers are not NA
    """
    boolean = identifiers.isna().all()
    return boolean


def fetch_uniprot_annotation(gui_object, identifiers, settings_dict):
    """
    input:
    gui_object = PyQt5, Qapplication
    identifiers = pd.Series
    settings_dict = dict["uniprot_step"], dictionary containing information about the uniprot user definable parameters.
    output:
    protein_data_list = list, list of protein_data_dict
    """
    gui_object.report_status("Step 2, fetching data from uniprot. The data is fetched from uniprot in batches, after each batch {sleep_time} seconds pass\n"
                             "before a new batch is queried to uniprot. This feature is implement to prevent being blacklisted".format(sleep_time=settings_dict["request_idle_time"]))
    protein_data_dict = {}
    function_dict = construct_function_dict(settings_dict)

    #split identifiers into multiple sub arrays of length batch_length:
    identifier_batches = np.split(identifiers, range(settings_dict["batch_amount"],len(identifiers), settings_dict["batch_amount"]))
    for n_batch, identifiers_batch in enumerate(identifier_batches):
        gui_object.report_status(f"Start fetching uniprot data for batch number {n_batch + 1} of the total {len(identifier_batches)} batches")
        requestURL = settings_dict["uniprot_base_url"]+settings_dict["uniprot_request_url"]+",".join(identifiers_batch)
        request = requests.get(requestURL, headers={"Accept" : "application/json"})
        if not request.ok:
            logging.error(f"Something went wrong with query {n_batch} while querying the uniprot database. The {len(identifiers)} proteins of this batch will be ignored")
            request.raise_for_status()
            for identifier in identifiers_batch:
                protein_data_dict[identifier] = {"gene_name":np.nan, "protein_name":np.nan, "organism_name":np.nan, "hyperlink":np.nan, "cell_compartment":np.nan, "string_linkout":np.nan}
        else:
            protein_data_dict = update_protein_data_dict(gui_object, request.json(), identifiers_batch, function_dict, protein_data_dict, settings_dict)
            protein_data_dict = add_uniprot_hyperlink(protein_data_dict, settings_dict, identifiers_batch)
            protein_data_dict = add_string_linkout(protein_data_dict, settings_dict, identifiers_batch)
            logging.info(f"Succesfully fetched and saved data from uniprot for batch {n_batch + 1}:")
        time.sleep(settings_dict["request_idle_time"])
    gui_object.report_status("Step 2, fetching data from uniprot, is finished.")
    return protein_data_dict

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

def update_protein_data_dict(gui_object, uniprot_output_list, identifiers, function_dict, protein_data_dict, settings_dict):
    """
    This function is meant to loop over the identifiers and update a dict with values per protein.
    input:
    gui_object = PyQt5, Qapplication
    identifiers = list, list of uniprot identifiers
    function_dict = dict{function_name(string) : function(function)}
    protein_data_dict = dict{protein : {gene_name, protein_name, organism_name, cell_comparment}}
    output:
    protein_data_dict = dict{protein : {gene_name, protein_name, organism_name, cell_comparment}}
    """
    for identifier in identifiers:
        uniprot_data_dict = get_matching_uniprot_query(uniprot_output_list, identifier)
        if not identifier in protein_data_dict.keys():
            protein_data_dict[identifier] = {}
        for function_name, function in function_dict.items():
            if None == uniprot_data_dict:
                protein_data_dict[identifier].update({function_name : np.nan})
            else:
                try:
                    protein_data_dict[identifier].update({function_name : function(uniprot_data_dict, settings_dict)})
                except KeyError as key_error:
                    log_error(gui_object, f"While getting the {function_name} for protein {identifier} a key error occured thus this field will be ignored", key_error)
                    protein_data_dict[identifier].update({function_name : np.nan})
                except IndexError as index_error:
                    log_error(gui_object, f"While getting the {function_name} for protein {identifier} an index error occured thus this field will be ignored", index_error)
                    protein_data_dict[identifier].update({function_name : np.nan})
                except Exception as exception:
                    log_error(gui_object, f"While getting the {function_name} for protein {identifier} an exception error occured thus this field will be ignored", exception)
                    protein_data_dict[identifier].update({function_name : np.nan})
    return protein_data_dict

def get_matching_uniprot_query(uniprot_output_list, identifier):
    """
    Per batch not all proteins are found. This function returns the uniprot data for the given identifier if it is present in the uniprot output.
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
            protein_data_dict[identifier].update({"uniprot_hyperlink":make_hyperlink(settings_dict["uniprot_protein_base_url"]+identifier)})
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
    return protein_data_dict

def get_uniprot_gene_name(uniprot_data_dict, settings_dict):
    """
    input:
    uniprot_data_dict = [{accession: "", id:"", proteinExistence:"", info:{}, organism:{}, protein:{}, gene:{}, features:{}, dbReferences:{}, keywords:[], references:[], sequence:{}}], this is the best case scenario.
    settings_dict = dict, dictionary containing user defined comments
    output:
    gene_name = string
    """
    if "gene" in uniprot_data_dict.keys():
        for uniprot_gene_dict in uniprot_data_dict["gene"]:
            if isinstance(uniprot_gene_dict[list(uniprot_gene_dict.keys())[0]], dict):
                if list(uniprot_gene_dict.keys())[0] in settings_dict["known_gene_names"]:
                    gene_name = uniprot_data_dict["gene"][0][list(uniprot_gene_dict.keys())[0]]["value"]
                else:
                    logging.debug(f'A new gene name has been found: {list(uniprot_gene_dict.keys())[0]}. You should add this gene name in {settings_dict["known_gene_names"]}')
                    gene_name = uniprot_data_dict["gene"][0][list(uniprot_gene_dict.keys())[0]]["value"]
            elif isinstance(uniprot_gene_dict[list(uniprot_gene_dict.keys())[0]], list):
                if list(uniprot_gene_dict.keys())[0] in settings_dict["known_gene_names"]:
                    gene_name = uniprot_data_dict["gene"][0][list(uniprot_gene_dict.keys())[0]][0]["value"]
                else:
                    logging.debug(f'A new gene name has been found: {list(uniprot_gene_dict.keys())[0]}. You should add this gene name in {settings_dict["known_gene_names"]}')
                    gene_name = uniprot_data_dict["gene"][0][list(uniprot_gene_dict.keys())[0]][0]["value"]
            else:
                logging.debug(f"In the extract gene name from uniprot output the type {type(uniprot_gene_dict[list(uniprot_gene_dict.keys())[0]])} has been found while it is not expected")
    else:
        gene_name = np.nan
    return gene_name

def get_protein_name(uniprot_data_dict, settings_dict):
    """
    input:
    uniprot_data_dict = [{accession: "", id:"", proteinExistence:"", info:{}, organism:{}, protein:{}, gene:{}, features:{}, dbReferences:{}, keywords:[], references:[], sequence:{}}], this is the best case scenario.
    settings_dict = dict, dictionary containing user defined comments
    output:
    protein_name = string
    """
    if "protein" in uniprot_data_dict.keys():
        for uniprot_protein_name in uniprot_data_dict["protein"].keys():
            if uniprot_protein_name in settings_dict["known_protein_names"]:
                if isinstance(uniprot_data_dict["protein"][uniprot_protein_name], list):
                    protein_name = uniprot_data_dict["protein"][uniprot_protein_name][0]["fullName"]["value"]
                elif isinstance(uniprot_data_dict["protein"][uniprot_protein_name], dict):
                    protein_name = uniprot_data_dict["protein"][uniprot_protein_name]["fullName"]["value"]
                return protein_name
            else:
                logging.debug(f"A new protein name has been found: {str(uniprot_protein_name)}. You should add this protein name to {settings_dict['known_protein_names']}")
                protein_name = list(uniprot_data_dict["protein"].values()[0])["fullName"]["value"]
    else:
        protein_name = np.nan
    return protein_name

def get_organism_name(uniprot_data_dict, settings_dict):
    """
    input:
    uniprot_data_dict = [{accession: "", id:"", proteinExistence:"", info:{}, organism:{}, protein:{}, gene:{}, features:{}, dbReferences:{}, keywords:[], references:[], sequence:{}}], this is the best case scenario.
    settings_dict = dict, dictionary containing user defined comments
    output:
    organism_name = string 
    """
    if "organism" in uniprot_data_dict.keys():
        organism_name = uniprot_data_dict["organism"]["names"][0]["value"]
    else:
        organism_name = np.nan
    return organism_name

def get_cell_compartment(uniprot_data_dict, settings_dict):
    """
    input:
    uniprot_data_dict = {accession: "", id:"", proteinExistence:"", info:{}, organism:{}, protein:{}, gene:{}, features:{}, dbReferences:{}, keywords:[], references:[], sequence:{}}, this is the best case scenario.
    settings_dict = dict, dictionary containing user defined comments
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
            string_linkout_dict[identifier] = make_hyperlink(settings_dict["string_base_url"]+uniprot_mapped_proteins_dict[identifier])
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
    logging.info("Start adding the uniprot columns to the existing dataframe")

    uniprot_column_names = get_column_names(uniprot_options_dict)
    for uniprot_column_name in uniprot_column_names:
        uniprot_column_values = get_uniprot_column_values(protein_groups_dataframe["identifier"], uniprot_column_name, protein_data_dict)
        protein_groups_dataframe[uniprot_column_name] = uniprot_column_values

    logging.info("Finished adding the uniprot columns to the existing dataframe")
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
    protein_data_dict = dict{identifier: {"gene_name":"", "protein_name":"", "organism_name":"", "uniprot_hyperlink":"", "cell_compartment":np.nan, "string_linkout":""}}
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
    Evaluate whether the user didn't accidentally put all settings to 0 while the uniprot step is set to 1.
    input:
    uniprot_options = dict{"option":boolean}
    output:
    boolean, True whenever at least 1 element is True and False if all elements are False
    """
    return any(uniprot_options.values())

def is_protein_in_mitocarta(gui_object, protein_groups_dataframe, mitocarta_species_dataframe, species_name, new_column_name, symbol_column, additional_symbol_column):
    """
    input:
    gui_object = PyQt5, Qapplication
    protein_groups_dataframe = pd.DataFrame()
    mitocarta_species_dataframe = pd.DataFrame()
    species_name = string
    new_column_name = string
    symbol_column = string, initial column to evaluate whether proteins are also found in the mitocarta data
    additional_symbol_column = string, secondary column which should contain multiple symbols for the same gene.
    output:
    protein_groups_dataframe = pd.DataFrame()
    """
    gui_object.report_status(f"Start elucidating which proteins are present in the {species_name} mitocarta dataset")
    protein_groups_dataframe["gene_name"].fillna('', inplace=True)
    mitocarta_species_dataframe[additional_symbol_column].fillna('', inplace=True)
    #In the mitocarta_mouse_data there are several synonyms which are not string, in order to remove these I wrote the  following line of code:
    mitocarta_species_dataframe[additional_symbol_column].where(cond=[True if type(x) == str else False for x in mitocarta_species_dataframe[additional_symbol_column].str.split("|")],
                                                  other="-", inplace=True)
    try:
        organism_rows = protein_groups_dataframe["organism_name"] == species_name
        symbol_rows = protein_groups_dataframe["gene_name"].isin([symbol for symbol in protein_groups_dataframe["gene_name"]
                                                                  if mitocarta_species_dataframe[symbol_column].str.contains(symbol).any() == True
                                                                  or [symbol in series_list for series_list in mitocarta_species_dataframe[additional_symbol_column].str.split("|")]])

        presency_index = protein_groups_dataframe.index[(organism_rows) & (symbol_rows)]
        protein_groups_dataframe[new_column_name] = [1.0 if index in presency_index else 0.0 for index in range(protein_groups_dataframe.shape[0])]
    except IndexError as index_error:
        log_error(gui_object, f"An index error occurred while evaluating whether proteins are found in {species_name} mitocarta data set. The mitocarta step will be ignored.", index_error)
    except ValueError as value_error:
        log_error(gui_object, f"An value error occurred while evaluating whether proteins are found in {species_name} mitocarta data set. The mitocarta step will be ignored.", value_error)
    except Exception as error:
        log_error(gui_object, f"An exception occurred while evaluating whether proteins are found in {species_name} mitocarta data set. The mitocarta step will be ignored.", error)
    gui_object.report_status(f"Finished elucidating which proteins are found in the {species_name} mitocarta dataset")
    return protein_groups_dataframe


def cluster_reorder(gui_object, sample_specific_dataframe, method = 'average', metric = 'correlation'):
    """
    The complexome profiling data is transformed to a condensed distance matrix and the proteins are clustered using hierarchical clustering. 
    The optimal order is determined resulting in minimal distance between adjacent leaves. 
    input:
    gui_object = PyQt5, Qapplication
    sample_specific_dataframe = pd.DataFrame()
    method = string
    metric = string
    output:
    order = dict{protein_identifier : ordered_index}
    clustered = np.array(), encoded as linkage matrix 
    """
    try:
        condensed_distance_matrix = spd.pdist(np.array(sample_specific_dataframe))

        clustered = sch.linkage(condensed_distance_matrix, method=method, metric=metric, optimal_ordering=True)
        dendrogram = sch.dendrogram(clustered, labels=sample_specific_dataframe.index.values, orientation='right', no_plot=True)
        ordered_index = sample_specific_dataframe.iloc[dendrogram['leaves'],:].index.values
        order = {label: index_x for index_x, label in enumerate(ordered_index)}
        return order, clustered
    except Exception as error:
        log_error(gui_object, "An exception occured while applying clustering on a sample", error)
        return {}, np.empty([0,0], dtype="float64")

def dump_data_to_excel(gui_object, protein_groups_dataframe, non_selected_dataframe, settings_dict, original_column_order):
    """
    The last part of this script, dump the complexome profiling data into an excel file.
    input:
    gui_object = PyQt5, Qapplication
    protein_groups_dataframe = pd.DataFrame()
    non_selected_dataframe = pd.DataFrame(), proteins that have been filtered away in the analysis
    settings_dict = dict, dictionary with parameters for this function
    output:
    None
    """
    gui_object.report_status("Step 5, start writing away the data to the excel file {file_name}".format(file_name=settings_dict["make_excel_file_step"]["excel_file_name"]))
    try:
        writer = pd.ExcelWriter(settings_dict["make_excel_file_step"]["excel_file_name"], engine='xlsxwriter', mode="w")
        workbook = writer.book

        ordered_columns = get_ordered_sample_columns(protein_groups_dataframe)
        protein_groups_dataframe = order_complexome_profiling_dataframe(protein_groups_dataframe, ordered_columns, original_column_order, settings_dict)
        non_selected_dataframe = order_complexome_profiling_dataframe(non_selected_dataframe, [], original_column_order, settings_dict)

        protein_groups_dataframe.to_excel(writer, sheet_name = 'data', index=False)
        non_selected_dataframe.to_excel(writer, sheet_name = 'filtered away proteins', index=False)
        worksheet = writer.sheets['data']

        positions = get_sample_positions(protein_groups_dataframe.columns.tolist())
        apply_conditional_formating_per_sample(protein_groups_dataframe, positions, writer, worksheet, workbook)

        writer.save()
        gui_object.report_status("Finished writing away the data to the excel file {file_name}".format(file_name=settings_dict["make_excel_file_step"]["excel_file_name"]))
    except Exception as error:
        log_error(gui_object, "An error occured while trying to write away the data. The data will be written away as .csv file", error)
        protein_groups_dataframe.to_csv("maxquant_saved_result.csv", sep=",", index=False)


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

def order_complexome_profiling_dataframe(protein_groups_dataframe, ordered_columns, original_column_order, settings_dict):
    """
    Define the order of the output dataframe. Set the first column to be a identifier column and the rest is not interesting. 
    input:
    protein_groups_dataframe = pd.DataFrame()
    ordered_columns = list, list of column names
    original_column_order = pd.Series(), original column order
    settings_dict = dict, dictionary with user defined settings
    output:
    protein_groups_dataframe = pd.DataFrame
    """
    original_columns = set(original_column_order).difference(set(ordered_columns))
    ordered_columns = list(original_columns) + list(ordered_columns)
    for column in protein_groups_dataframe.columns:
        if column in settings_dict["make_excel_file_step"]["identifier_column_names"] and column in ordered_columns:
            ordered_columns.remove(column)
            ordered_columns.insert(0, column)
        elif column not in ordered_columns:
            ordered_columns.append(column)

    protein_groups_dataframe = protein_groups_dataframe.reindex(columns=ordered_columns)
    return protein_groups_dataframe

def get_ordered_sample_columns(complexome_profiling_dataframe):
    """
    Order the complexome profiling dataframe so the pattern emerges: sample iBAQ columns - sample clustered column - etc. - global clustered column
    input:
    complexome_profiling_dataframe = pd.DataFrame()
    output:
    ordered_columns = list, list of column names
    """
    ordered_columns = []
    global_cluster_column = "global_clustered"
    global_total_protein_abundance_column = "global_summed_iBAQ_value"
    sample_names = get_sample_names(complexome_profiling_dataframe)

    for sample_name in sample_names:
        sample_columns = complexome_profiling_dataframe.filter(regex=f"iBAQ {sample_name}_").columns
        ordered_columns.extend(sample_columns)
        ordered_columns.append(f'sample_{sample_name}_clustered')
        ordered_columns.append(f"{sample_name}_summed_iBAQ_value")

    #add global clustering column to the end of the ordered_columns list:
    ordered_columns.append(global_cluster_column)
    ordered_columns.append(global_total_protein_abundance_column)
    return ordered_columns

def get_sample_positions(column_names):
    """
    Per sample, get the start and end column positions
    input:
    column_names = pd.Index(), iterable with the column names of the complexome profiling dataframe
    output:
    positions = list, [[start_position, stop_position], [start_position, stop_position]]
    """
    positions = []
    start_position, end_position = 0, 0
    for column_name in column_names:
        if "iBAQ " in column_name:
            if start_position == 0:
                start_position = column_names.index(column_name)
            end_position = column_names.index(column_name)
        else:
            if end_position != 0:
                positions.append([start_position, end_position])
                start_position, end_position = 0, 0
    if end_position != 0:
        positions.append([start_position, end_position])
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

def filter_dataframe_step(gui_object, protein_groups_dataframe, settings_dict):
    """
    process group proteins file by filtering columns and rows:
    input:
    gui_object = PyQt5, Qapplication
    protein_groups_dataframe = pd.DataFrame()
    settings_dict = dict, dictionary with user defined settings
    output:
    protein_groups_dataframe = pd.DataFrame()
    non_selected_dataframe = pd.DataFrame()
    """
    if settings_dict["steps_dict"]["filtering_step"] == True:
        gui_object.report_status("Step 1, filtering the dataframe, has started")

        protein_groups_dataframe = filter_dataframe_columns(gui_object, protein_groups_dataframe, settings_dict["filtering_step"])
        protein_groups_dataframe, filtered_groups_dataframe = filter_dataframe_rows(protein_groups_dataframe, settings_dict["filtering_step"])
        #Reset the index, through this way new columns can be added to the dataframe
        protein_groups_dataframe.reset_index(inplace=True)
        protein_groups_dataframe.dropna(axis="index", inplace=True)
        identifiers = fetch_identifiers(protein_groups_dataframe)
        protein_groups_dataframe['identifier'] = identifiers

        add_protein_abundance_columns(protein_groups_dataframe)

        gui_object.report_status("Step 1, filtering the dataframe, is finished")
        return protein_groups_dataframe, filtered_groups_dataframe, protein_groups_dataframe.columns
    else:
        gui_object.report_status("Step 1 (filtering the columns and rows of the main dataframe)\nhas been disabled and will be skipped")
        return protein_groups_dataframe, pd.DataFrame(), pd.Series()


def add_protein_abundance_columns(protein_groups_dataframe):
    """
    Per sample for each protein, add the protein abundances of each fraction together.
    Additionally, add the protein abundances of each fraction from all samples together.
    input:
    protein_groups_dataframe = pd.Dataframe
    output:
    protein_groups_dataframe = pd.Dataframe
    """
    sample_names = get_sample_names(protein_groups_dataframe)
    for sample_name in sample_names:
        sample_columns = protein_groups_dataframe.filter(regex=f"iBAQ {sample_name}_").columns
        sample_protein_abundances = protein_groups_dataframe[sample_columns].apply(summing_protein_abundances, axis=1)
        protein_groups_dataframe[f"{sample_name}_summed_iBAQ_value"] = sample_protein_abundances
    protein_abundance_sample_columns = [f"{sample_name}_summed_iBAQ_value" for sample_name in sample_names]
    global_protein_abundances = protein_groups_dataframe[protein_abundance_sample_columns].apply(summing_protein_abundances, axis=1)
    protein_groups_dataframe["global_summed_iBAQ_value"] = global_protein_abundances


def get_sample_names(protein_groups_dataframe):
    """
    Get the available sample names from the main dataframe
    input:
    protein_groups_dataframe = pd.Dataframe()
    output:
    sample_names = list, list of strings
    """
    all_sample_columns = protein_groups_dataframe.columns.str.contains(pat="^iBAQ ", case=True, na=False, regex=True)
    sample_names_with_fraction_ids = protein_groups_dataframe.columns[all_sample_columns].str.replace(pat="^iBAQ ", repl="", case=True, regex=True)
    sample_names = sample_names_with_fraction_ids.str.replace(pat="_[$0-9]{2}", repl="", case=True, regex=True)
    sample_names = set(sample_names.tolist())

    return sample_names


def summing_protein_abundances(protein_abundances):
    """
    For the sample, add the protein abundances per fraction together for each protein.
    input:
    protein_abundances = pd.Series()
    output:
    sample_protein_abundances = pd.Series(), the summed protein abundances for the given sample where the value corresponds to the index of the main dataframe
    """
    return pd.Series.sum(protein_abundances)


def fetch_uniprot_annotation_step(gui_object, protein_groups_dataframe, settings_dict):
    """
    fetch annotation for uniprot identifiers:
    input:
    protein_groups_dataframe = pd.DataFrame()
    settings_dict = dict, dictionary with user defined settings
    output:
    protein_groups_dataframe = pd.DataFrame()
    """
    if settings_dict["steps_dict"]["uniprot_step"] == True and evaluate_uniprot_settings(settings_dict["uniprot_step"]["uniprot_options"]) == True:
        if are_identifiers_not_available(protein_groups_dataframe["identifier"]) == False:
            protein_data_dict = fetch_uniprot_annotation(gui_object, protein_groups_dataframe["identifier"], settings_dict["uniprot_step"])
            protein_groups_dataframe = append_uniprot_data_to_dataframe(protein_groups_dataframe, protein_data_dict, settings_dict["uniprot_step"]["uniprot_options"])
        else:
            gui_object.report_status("Uniprot will not be queried because no uniprot identifiers were found in the \'Fasta headers\' column.")
    else:
        if settings_dict["steps_dict"]["uniprot_step"] == False:
            gui_object.report_status("Uniprot will not be queried for information due to the step being disabled")
        elif evaluate_uniprot_settings(settings_dict["uniprot_step"]["uniprot_options"]) == False:
            gui_object.report_status("Uniprot will not be queried for information because the step is enabled, but all the fields are disabled")
    return protein_groups_dataframe

def is_protein_in_mitocarta_step(gui_object, settings_dict, protein_groups_dataframe):
    """
    Evaluate per protein if it is found in the mouse or human dataset. 
    input:
    settings_dict = dict, dictionary with user defined settings
    protein_groups_dataframe = pd.DataFrame()
    output:
    protein_groups_dataframe = pd.DataFrame()
    """
    if settings_dict["steps_dict"]["mitocarta_step"] == True and settings_dict["steps_dict"]["uniprot_step"] == True and \
       settings_dict["uniprot_step"]["uniprot_options"]["get_gene_name"] == True and settings_dict["uniprot_step"]["uniprot_options"]["get_organism_name"] == True and \
       are_identifiers_not_available(protein_groups_dataframe["identifier"]) == False:
        mitocarta_mouse_dataframe = read_in_excel_file(gui_object, settings_dict["mitocarta_step"]["mitocarta_mouse_ftp_link"], settings_dict["mitocarta_step"]["mouse_sheet_name"])
        mitocarta_human_dataframe = read_in_excel_file(gui_object, settings_dict["mitocarta_step"]["mitocarta_human_ftp_link"], settings_dict["mitocarta_step"]["human_sheet_name"])
        protein_groups_dataframe = is_protein_in_mitocarta(gui_object, protein_groups_dataframe, mitocarta_mouse_dataframe,
                                                           settings_dict["mitocarta_step"]["mitocarta_mouse_organism"], "mitocarta_mouse_presency",
                                                           settings_dict["mitocarta_step"]["mitocarta_symbol_column"], settings_dict["mitocarta_step"]["mitocarta_additional_symbol_column"])
        protein_groups_dataframe = is_protein_in_mitocarta(gui_object, protein_groups_dataframe, mitocarta_human_dataframe,
                                                           settings_dict["mitocarta_step"]["mitocarta_human_organism"], "mitocarta_human_presency",
                                                           settings_dict["mitocarta_step"]["mitocarta_symbol_column"], settings_dict["mitocarta_step"]["mitocarta_additional_symbol_column"])

    else:
        if settings_dict["steps_dict"]["mitocarta_step"] == False and settings_dict["steps_dict"]["uniprot_step"] == True:
            gui_object.report_status("The final dataframe will not be enriched with information from mitocarta because the mitocarta step is disabled.")
        elif settings_dict["steps_dict"]["uniprot_step"] == True and settings_dict["uniprot_step"]["uniprot_options"]["get_gene_name"] == False or settings_dict["uniprot_step"]["uniprot_options"]["get_organism_name"] == False:
            gui_object.report_status("Uniprot has been queried but the gene name or the organism name has not been retrieved. Set the organism and gene_name to 1 in order to get the information from mitocarta")
        elif settings_dict["steps_dict"]["mitocarta_step"] == True and settings_dict["steps_dict"]["uniprot_step"] == False:
            gui_object.report_status("The final dataframe will not be enriched with information from mitocarta because the uniprot step is disabled.\nUniprot is needed to get per protein the gene symbol, the key to find values in uniprot.")
        elif settings_dict["steps_dict"]["mitocarta_step"] == False and settings_dict["steps_dict"]["uniprot_step"] == False:
            gui_object.report_status("The final dataframe will not be enriched with information from mitocarta because the mitocarta and uniprot steps are disabled")
        elif are_identifiers_not_available(protein_groups_dataframe["identifier"]) == True:
            gui_object.report_status("Step 2 was skipped because no uniprot identifiers were found in the Fasta header column.\nSo, step 3 is also skipped because the gene name column is unavailable.")

    return protein_groups_dataframe

def apply_clustering_step(gui_object, settings_dict, protein_groups_dataframe):
    """
    Apply hierarchical cluster analysis per sample and a global one as well 
    input:
    protein_groups_dataframe = pd.DataFrame()
    settings_dict = dict, dictionary with user defined settings
    output:
    protein_groups_dataframe = pd.DataFrame()
    """
    if settings_dict["steps_dict"]["clustering_step"] == True:
        gui_object.report_status("Step 4, cluster the fractions per sample using hierarchical clustering.")
        sample_names = get_sample_names(protein_groups_dataframe)

        for sample_name in sample_names:
            logging.info(f"Start hierarchical clustering for sample {sample_name}")
            sample_specific_dataframe = pd.DataFrame(protein_groups_dataframe[protein_groups_dataframe.columns[protein_groups_dataframe.columns.to_series().str.contains(pat=f"iBAQ {sample_name}", regex=True)]], dtype="float64")
            order_mapping, clustered = cluster_reorder(gui_object, sample_specific_dataframe, settings_dict["clustering_step"]["method"], settings_dict["clustering_step"]["metric"])
            protein_groups_dataframe[f'sample_{sample_name}_clustered'] = pd.Series(order_mapping)
            logging.info(f"Finished hierarchical clustering for sample {sample_name}")
        logging.info("Start hierarchical clustering for all samples")
        global_order_mapping, global_clustered = cluster_reorder(gui_object, protein_groups_dataframe[protein_groups_dataframe.columns[protein_groups_dataframe.columns.to_series().str.contains("iBAQ ")]])
        protein_groups_dataframe['global_clustered'] = pd.Series(global_order_mapping)
        gui_object.report_status("Step 4, finished clustering the fractions per sample using hierarchical clustering.")
    else:
        gui_object.report_status("Step 4, clustering the fractions per sample using hierarchical clustering has been disabled.")

    return protein_groups_dataframe
def dump_to_excel_step(gui_object, protein_groups_dataframe, filtered_groups_dataframe, settings_dict, original_column_order):
    """
    write away dataframe to an excel file:
    input:
    protein_groups_dataframe = pd.DataFrame(), dataframe containing all the selected proteins
    filtered_groups_dataframe = pd.DataFrame(), dataframe containing all the filtered proteins
    settings_dict = dict, dictionary with user defined settings
    output:
    None
    """
    if settings_dict["steps_dict"]["make_excel_file_step"] == True:
        dump_data_to_excel(gui_object, protein_groups_dataframe, filtered_groups_dataframe, settings_dict, original_column_order)
    else:
        logging.info("Step 5, writing away the data to an excel file, will not be executed because the user has disabled the step.")
        return

if __name__ == "__main__":
    args = get_user_arguments()

    settings_dict = load_json(args.settings_file_path)
    protein_groups_dataframe = read_in_protein_groups_file(args.maxquant_file_path)

    protein_groups_dataframe, filtered_groups_dataframe, original_column_order = filter_dataframe_step(protein_groups_dataframe, settings_dict)

    protein_groups_dataframe = fetch_uniprot_annotation_step(protein_groups_dataframe, settings_dict)

    protein_groups_dataframe = is_protein_in_mitocarta_step(settings_dict, protein_groups_dataframe)

    protein_groups_dataframe = apply_clustering_step(settings_dict, protein_groups_dataframe)

    dump_to_excel_step(protein_groups_dataframe, filtered_groups_dataframe, settings_dict, original_column_order)
