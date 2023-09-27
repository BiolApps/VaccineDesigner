# Vaccine Designer
Vaccine Designer is an R Shiny application aiming for the construction of vaccine sequences  based on multi epitope design workflow.

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

# 3. Install Python packages.
The third step is to install the required packages for the Vaccine Designer scripts with the command:

pip install pandas numpy bs4 torch requests_html fair-esm plotly importlib-resources pytz six tenacity typing-extensions zipp joblib 
pip install scikit-learn==1.2.2

# 4. Install R libraries.
After installing python packages, following is to install the required R libraries with the command:

RUN R -e "install.packages(c('shinydashboard','shiny','DT','shinyFiles','reticulate','dplyr','shinyjs'))"
RUN R -e "install.packages('BiocManager')"
RUN R -e "library(BiocManager);BiocManager::install('Biostrings')"

# 5. Make executable the .sh files



