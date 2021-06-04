<h3>Description</h3>
This repository serves as a place to store code in order to process a maxquant file used for complexome profiling. The maxquant file is processed in 5 steps. Firstly, proteins not of interest are written to a separate excel sheet and columns containing the words "REV" or "CON" are filtered away. Secondly, each protein is queried to uniprot and per protein the gene name, protein name, organism name, the cell compartment, a uniprot specific protein hyperlink and a hyperlink to the string database are retrieved if they are available. The user can specify which elements should be retrieved from uniprot. Thirdly, per protein it is denoted if the protein is found in the human or mouse mitocarta database. Fourthly, each sample is clustered based on hierarchical clustering. Fifthly, an excel file is made with the 2 different sheets containing the filtered away proteins and proteins with information as well as the accompanying complexome profiling data(conditional formatting has been applied on these columns). The user is able to disable any step they wish not to do, it should be noted that whenever the uniprot step is disabled the mitocarta step doesn't work because the input for the mitocarta step comes from the uniprot output. 
<h3>Usage</h3>
In order to use this script you first need to download an executable and user defined settings file.
First step, go to the releases web page:
![Alt text](find_releases.png?raw=true "Go to the releases webpage on github")

Second step, download from the newest release the executable file and the 'maxquant_settings.json' file. 
![Alt text](download_this_file.png?raw=true "Download the .exe file and maxquant_settings.json")

Once you have downloaded the executable and settings file place them in a folder of your choosing. 
This folder should contain (1) the maxquant file you wish to process, (2) the settings file and (3) the executable. 
Once the executable has run an excel file will be generated in the same folder as the executable(unless you specify the path the excel file should be written to).
Third step, open a powershell window in the folder where the maxquant processing files can be found. 
![Alt text](get_to_powershell.png?raw=true "Open powershell in windows for specified folder")

Fourth step, type '.\'[executable_name.exe] and the parameters for the maxquant file as shown here with the file names:
![Alt text](alternative_way_to_execute_command.png?raw=true "Alternative way to execute the command")

If windows doesn't recognize the files or it is more convenient to keep the files in another folder try to execute the executable with the complete filepaths, like so:
![Alt text](how_to_execute_command.png?raw=true "How to execute the command")

<h4>User defined parameters:</h4>
1. steps_dict -> which steps should be taken? 1 means active 0 means non-active, shouldn't be anything else than 1 or 0. 
   1. filtering_step
   2. uniprot_step
   3. mitocarta_step
   4. clustering_step 
   5. make_excel_file_step

2. filtering_step -> parameters for the filtering step
   1. EXACT_MATCHES -> Retain columns which match exactly the values in this list
   2. CONTAINS -> Columns containing elements from this list should be retained
   3. PROTEIN_FILTERS -> Columns containing elements from this list will not be retained

3. uniprot_step -> parameters for querying uniprot
   1. uniprot_options -> a dictionary similiar to steps_dict which enables the user to define which elements should be retrieved from uniprot.
      1. get_gene_name
      2. get_protein_name
      3. get_organism_name
      4. get_uniprot_hyperlink
      5. get_cell_compartment
      6. get_string_linkout
   2. uniprot_base_url -> This url is the basis for querying uniprot
   3. uniprot_request_url -> This url has parameters necessary to query uniprot and is combined with 'uniprot_base_url'.
   4. request_idle_time -> how many seconds should be waited between querying uniprot? This number should be low but not 0 to prevent to becoming blacklisted.
   5. batch_amount -> how many proteins should be retrieved from uniprot? There is a maximum of 100 else the program doesn't work
   6. uniprot_protein_base_url -> The hyperlink to each protein of the uniprot database has a generic part + protein name, this url encodes this generic path. 
   7. string_linkout_parameters -> dictionary containing information to query the uniprot mapping service to get string linkouts
      1. string_base_url -> The base url to query string + string identifier
      2. uniprot_mapping_service_url -> the base url to query the uniprot mapping service
      3. regex_pattern -> a regex pattern is needed to locate proteins with a suffix like "-2"

4. mitocarta_step -> parameters for querying the mitocarta database
   1. mitocarta_human_ftp_link -> A link to the human mitocarta excel file/database
   2. mitocarta_mouse_ftp_link -> A link to the mouse mitocarta excel file/database
   3. human_sheet_name -> the sheet name that should be used for the human excel file
   4. mouse_sheet_name -> the sheet name that should be used for the mouse excel file

5. clustering_step -> parameters for clustering the complexome profiling steps
   1. method -> method that should be used. Can be: 'single', 'complete', 'average', 'weighted', 'centroid', 'median' or 'ward'.
   2. metric -> distance metric that should be used, can be: braycurtis’, ‘canberra’, ‘chebyshev’, ‘cityblock’, ‘correlation’, ‘cosine’, ‘dice’, ‘euclidean’, ‘hamming’, ‘jaccard’, ‘jensenshannon’, ‘kulsinski’, ‘mahalanobis’, ‘matching’, ‘minkowski’, ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘yule’.

6. make_excel_file_step -> parameters for writing away the data into an excel file 
   1. excel_file_name -> the name the excel of the excel file

<h3>Authors</h3>
Ariel Komen and Joeri van Strien
<h3>Requirements</h3>
The requirements can be found in a separate requirements file. 
