<h3>Description</h3>
<p>
This repository serves as a storage place for code and an executable to process the output from maxquant used for complexome profiling experiments. The maxquant file is processed in 5 steps. 
In the first step unwanted columns are removed and proteins not of interest (these proteins contain the words "REV" or "CON" in the 'Majority protein IDs' column) are filtered away. 
In the second step proteins are queried in batches to uniprot and per protein the gene name, protein name, organism name, the cell compartment, a uniprot specific protein hyperlink and a hyperlink to the string database are retrieved if they are available. The user can specify which elements should be retrieved from uniprot. 
In the third step each protein is compared with the proteins from the mitocarta 'database' and in a column it is denoted if the protein is present in the human or mouse mitocarta database.
In the fourth step the different samples are clustered based on hierarchical clustering. Additionally, all samples together (global) are clustered. 
In the fifth and final step the processed output is written away to an excel file with the 2 different sheets where the first sheet contains the filtered away proteins and the second sheet contains the proteins with information as well as the accompanying complexome profiling data(conditional formatting has been applied on these columns).
</p>
<p>
This program uses a settings file in which the user can disable individual steps, change output behaviour and more. However, it should be noted that this program can easily not work when the input parameters are not as expected. 
Additionally, whenever the uniprot step is disabled the mitocarta step will also not be executed, because the input for the mitocarta step comes from the uniprot output. 
</p>
<h3>Usage</h3>
In order to use this script you first need to download an executable and user defined settings file (the help file is also nice to have).
First step, go to the releases web page. 
![Alt text](images/find_releases.png?raw=true "Go to the releases webpage on github")

Second step, download from the newest release the executable file and the 'maxquant_settings.json' file. 
![Alt text](images/download_this_file.png?raw=true "Download the .exe file and maxquant_settings.json")

Once you have downloaded the executable and settings file place them in a folder of your choosing. 
This folder should contain (1) the maxquant file you wish to process, (2) the settings file and (3) the executable. 
Now is a good moment to open the settings file and based on the description of the user defined parameters down below change parameters to change the behaviour of how the maxquant output is processed.  
An example setup: 
![Alt text](images/example_folder_structure.png?raw=true "An example folder structure")

Now the program is ready to run. Start the program by executing the executable which will start a shell script and after 1 second a Graphical User Interface (GUI) will appear. 
In order to run the program it needs the settings file and maxquant file. These two files do not have to be in the same folder as the executable. 
Anyway, you can select the maxquant output file by clicking 'Select maxquant file' which will open a file selection menu where you can select the file. The file path to the maxquant file is put in the textbox above the 'Select maxquant file' button. 
The settings file is selected in a similar way as the maxquant file, but you should select the maxquant file by clicking 'Select the settings file' button. 
![Alt text](images/file_selection_procedure.png?raw=true "An example folder structure")

When both files have been selected the program can be execute by clicking 'Process maxquant'. The GUI will print status messages and error messages whenever something goes wrong. 

<h4>User defined parameters:</h4>
1. steps_dict -> which steps would you like to execute? 1 means that the step is executed and 0 means the step is not executed. These parameters shouldn't be anything else than 1 or 0. 
   1. filtering_step
   2. uniprot_step
   3. mitocarta_step
   4. clustering_step 
   5. make_excel_file_step

2. filtering_step -> parameters for the filtering step
   1. EXACT_MATCHES -> Retain columns which match exactly the value in this list
   2. CONTAINS -> Columns containing elements from this list should be retained
   3. PROTEIN_FILTERS -> Rows/proteins containing elements from this list will not be retained(will be written to a separate excel sheet) 

3. uniprot_step -> parameters for querying uniprot
   1. uniprot_options -> a dictionary similiar to steps_dict which enables the user to define which elements should be retrieved from uniprot. These parameters shouldn't be anything else than 1 or 0
      1. get_gene_name
      2. get_protein_name
      3. get_organism_name
      4. get_uniprot_hyperlink
      5. get_cell_compartment
      6. get_string_linkout
   2. uniprot_base_url -> This url is the basis for querying uniprot
   3. uniprot_request_url -> This url has parameters necessary to query uniprot and is combined with 'uniprot_base_url'.
   4. request_idle_time -> In order to prevent being blacklisted between batches the program waits a few seconds. How many seconds should the program wait? Change this parameters at you own risk 
   5. batch_amount -> Uniprot is queried in batches of size n. How large should n be? Uniprot will not accept batches larger than 100. 
   6. uniprot_protein_base_url -> The hyperlink to each protein of the uniprot database has a generic part + protein name, this url encodes this generic path. 
   7. known_gene_names -> In order to retrieve the gene name from uniprot several keys are known, this list contains the known gene names.
   8. known_protein_names -> similar to known_gene_names but for protein_names.
   9. string_linkout_parameters -> dictionary containing information to query the uniprot mapping service to get string linkouts
      1. string_base_url -> The base url to query string 
      2. uniprot_mapping_service_url -> the base url to query the uniprot mapping service
      3. regex_pattern -> a regex pattern is needed to locate proteins with a suffix like "-2" which needs to be removed from proteins in order to query them to the mapping service of uniprot. 

4. mitocarta_step -> parameters for querying the mitocarta database
   1. mitocarta_human_ftp_link -> A link to the human mitocarta excel file/database
   2. mitocarta_mouse_ftp_link -> A link to the mouse mitocarta excel file/database
   3. human_sheet_name -> the sheet name that should be used for the human excel file
   4. mouse_sheet_name -> the sheet name that should be used for the mouse excel file
   5. mitocarta_symbol_column -> What is the main column in the mitocarta database to examine whether the protein is in the mitocarta database?
   6. mitocarta_additional_symbol_column -> Proteins have synonyms, which column in the mitocarta database contains synonyms for the gene names? 
   7. mitocarta_mouse_organism -> Which organism should be selected in the mouse mitocarata database? It can be something different from "Mus musculus" like 'Bos taurus'. 
   8. mitocarta_human_organism -> Which organism should be selected in the mouse mitocarata database? It can be something different from "Homo sapiens" like 'Xerophyta viscosa'. 

5. clustering_step -> parameters for clustering the complexome profiling samples
   1. method -> method for how the clustering is performed. Possible options are: 'single', 'complete', 'average', 'weighted', 'centroid', 'median' or 'ward'.
   2. metric -> which distance metric should be used. Possible distance metrics are: braycurtis’, ‘canberra’, ‘chebyshev’, ‘cityblock’, ‘correlation’, ‘cosine’, ‘dice’, ‘euclidean’, ‘hamming’, ‘jaccard’, ‘jensenshannon’, ‘kulsinski’, ‘mahalanobis’, ‘matching’, ‘minkowski’, ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘yule’.

6. make_excel_file_step -> parameters for writing away the data into an excel file 
   1. excel_file_name -> the name of the excel file. This file name can also be the path  
   2. identifier_column_names -> 

<h3>Authors</h3>
Ariel Komen and Joeri van Strien
<h3>Requirements</h3>
The requirements can be found in a separate requirements file. 
