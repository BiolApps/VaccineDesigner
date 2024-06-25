# Vaccine Designer
Vaccine Designer is an R Shiny application the constructs vaccine sequences based on a streamlined multi-epitope design workflow. The Web interface is available at: bioinformatics.med.auth.gr/VaccineDesigner

The source code is built on R and Python scripts. It is recommended to run to Linux system for the dependencies. 
The application has dependencies on the software of:
1) BepiPred 3.0
2) netMHCpan 4.1
3) netMHCIIpan 4.0
4) ToxinPred2
5) Algpred2
6) IEDB MHC Class I
7) IEDB MHC Class II
8) IEDB Population Coverage
9) DeepVacPred Source code
10) NetChop 3.1

The additional tools used in the pipeline, that do not require local installation are:
1) VaxiJen 2.0
2) ProtParam

# 1. Installation of software dependencies
First install the software dependencies that listed above. Some of the software need academic licence to be installed.

# 2. Run and configure the installed software
Ensure that the files can be executed and run correctly through command prompt.

# 3. Make a Local_software and a temp_data folder in the VaccineDesigner folder
Make two new folders inside the root directory. Place the installed software (Bepipred,netMHCpan,netMHCIIpan,ToxinPred,Algpred) in the folder "Local_Tools" (the whole tool folders, not just the executables).
The second folder is named temp_data and is the directory that some of the data will be generated.

# 4. Install Python packages
Install the required packages by executing the the command:

pip install pandas numpy bs4 torch requests_html fair-esm plotly importlib-resources pytz six tenacity typing-extensions zipp joblib 
pip install scikit-learn==1.2.2

# 5. Install R libraries
Install the required R library dependencies with the command:

RUN R -e "install.packages(c('shinydashboard','shiny','DT','shinyFiles','reticulate','dplyr','shinyjs'))"
RUN R -e "install.packages('BiocManager')"
RUN R -e "library(BiocManager);BiocManager::install('Biostrings')"

# 6. Make executable .sh files
Change direactory (cd) to the folder where the files are stored and run:
chmod +x update_envfile_toxin_alg.sh
chmod +x update_path_net.sh
Use sudo in case of "permission denied" errors.

# 7. Run the executable files
Run the following commands:
update_envfile_toxin_alg.sh path_envfile_toxinpred2 path_blastp ./toxinpred2/toxinpred2/Database/data ./toxinpred2/toxinpred2/progs/MERCI_motif_locator.pl ./toxinpred2/toxinpred2/Database/pos_motif.txt

update_envfile_toxin_alg.sh ./algpred2_new/algpred2_new/envfile path_blastp ./algpred2_new/algpred2_new/Database/data ./algpred2_new/algpred2_new/progs/MERCI_motif_locator.pl ./algpred2_new/algpred2_new/Database/pos_ige_motifs.txt

update_path_net.sh path_to_netMHCpan path_to_folder_netMHCpan-4.1
pdate_path_net.sh path_to_netMHCIIpan path_to_folder_netMHCIIpan-4.1

Change the paths according to the structure of your file system.

# 8. Edit server.R
Edit the server.R script, and change the file paths appropriately. Maybe you will need to make new directories for results etc., be sure to be placed inside the root directory.

# 9. Check software paths
Make sure that the software netMHCpan and netMHCIIpan, are moved to the default location for executable files of your file system e.g. /usr/local/bin.




