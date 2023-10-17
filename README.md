# Vaccine Designer
Vaccine Designer is an R Shiny application aiming for the construction of vaccine sequences  based on multi epitope design workflow. The web-tool can be accessible through : bioinformatics.med.auth.gr/VaccineDesigner

The source code of the Project is based on R and Python programming languages. It is recommended to run to Linux system for the dependencies. 
The application has dependencies on the software of:
1) BepiPred 3.0
2) netMHCpan 4.1
3) netMHCIIpan 4.0
4) ToxinPred2
5) Algpred2

# 1. Installation of Software Dependencies.
The first step to run the code is to have installed the software dependencies that listed above. Some of the software need of academic licence to be installed.

# 2. Run and configure the installed Software.
After downloading the software, you need to ensure that this files can be executed and run correctly through Command Prompt.

# 3. Make a Local_software and a temp_data folder in the Vaccine Designer folder
Make two new folders inside the root directory of the Vaccine Designer. Place the installed softwared (Bepipred,netMHCpan,netMHCIIpan,ToxinPred,Algpred) inside the folder Local_Tools. Don't place only the executables but also the whole folders. 
The second folder is named temp_data and is the directory that some of the data will be generated.

# 4. Install Python packages.
The third step is to install the required packages for the Vaccine Designer scripts with the command:

pip install pandas numpy bs4 torch requests_html fair-esm plotly importlib-resources pytz six tenacity typing-extensions zipp joblib 
pip install scikit-learn==1.2.2

# 5. Install R libraries.
After installing python packages, following is to install the required R libraries with the command:

RUN R -e "install.packages(c('shinydashboard','shiny','DT','shinyFiles','reticulate','dplyr','shinyjs'))"
RUN R -e "install.packages('BiocManager')"
RUN R -e "library(BiocManager);BiocManager::install('Biostrings')"

# 6. Make executable the .sh files
You need to run the following commands in your Command Prompt. Move with the "cd" command to the folder that the files are stored and run:
RUN chmod +x update_envfile_toxin_alg.sh
RUN chmod +x update_path_net.sh

# 7. Run the executable files

RUN update_envfile_toxin_alg.sh path_envfile_toxinpred2 path_blastp ./toxinpred2/toxinpred2/Database/data ./toxinpred2/toxinpred2/progs/MERCI_motif_locator.pl ./toxinpred2/toxinpred2/Database/pos_motif.txt

RUN update_envfile_toxin_alg.sh ./algpred2_new/algpred2_new/envfile path_blastp ./algpred2_new/algpred2_new/Database/data ./algpred2_new/algpred2_new/progs/MERCI_motif_locator.pl ./algpred2_new/algpred2_new/Database/pos_ige_motifs.txt

RUN update_path_net.sh path_to_netMHCpan path_to_folder_netMHCpan-4.1
RUN update_path_net.sh path_to_netMHCIIpan path_to_folder_netMHCIIpan-4.1

You need to change the paths based on your system.

# 8. Go to server.R and change the according paths
Go to server.R script, and change the according paths with the place of the scripts or executables. Maybe you will need to make new directories for results etc., be sure to be placed inside the root directory.

# 9. Make sure that the software netMHCpan and netMHCIIpan, are moved to /usr/local/bin folder, to be executed without path.




