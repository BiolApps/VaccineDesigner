library(shiny)
library(shinydashboard)
library(DT)
library(Biostrings)
library(shinyFiles)
library(reticulate)
library(dplyr)
library(shinyjs)
library(shinydashboardPlus)
#Sidebar
options(rsconnect.max.bundle.files = 150000)
sidebar <- dashboardSidebar(
  hr(),
  useShinyjs(),
  sidebarMenu(id="tabs",
              menuItem("B cell epitopes", tabName="B_cell_epitopes", icon=icon("tag"),
                       menuSubItem("Upload",tabName="B_upload",icon=icon("file")),
                       menuSubItem("Parameters",tabName="B_params",icon=icon("gear")),
                       menuSubItem("Execute",tabName = "B_res",icon=icon("terminal")),
                       menuSubItem("Filter Results",tabName="B_ranked",icon=icon("filter"))
              ),
              menuItem("CTL epitopes", tabName = "CTL_epitopes", icon=icon("tag"),
                       menuSubItem("Upload",tabName="C_upload",icon=icon("file")),
                       menuSubItem("Parameters",tabName="C_params",icon=icon("gears")),
                       menuSubItem("Execute",tabName = "C_res",icon=icon("terminal")),
                       menuSubItem("Filter Results",tabName="C_ranked",icon=icon("filter"))
              ),
              menuItem("HTL epitopes",  icon = icon("tag"),
                       menuSubItem("Upload",tabName="H_upload",icon=icon("file")),
                       menuSubItem("Parameters",tabName="H_params",icon=icon("gears")),
                       menuSubItem("Execute",tabName = "H_res",icon=icon("terminal")),
                       menuSubItem("Filter Results",tabName="H_ranked",icon=icon("filter"))
              ),
              menuItem("Vaccine Designer",tabName="Revacc",icon=icon("shield-virus"),
                       menuSubItem("Upload",tabName="Re_upload",icon=icon("file")),
                       menuSubItem("Evaluation",tabName="Re_eval",icon=icon("gears")),
                       menuSubItem("Final Results",tabName="Re_results",icon=icon("square-poll-vertical"))
              ),
              menuItem("Information", tabName = "readme", icon=icon("mortar-board")),
              menuItem("About", tabName = "about", icon = icon("question"),selected= TRUE)
  ),
  hr()
  
  
)

#Header
header <- dashboardHeader()
header$children[[2]]$children <-  tags$a(href='#',
                                         tags$img(src='logo.png',height='20',width='200'))

#Body
body <- dashboardBody(
  useShinyjs(),
  tags$head(
    tags$style(
      HTML(
        "
        .buttons {
          font-size: 16px;
          padding: 10px 20px;
        }
        "
      )
    )
  ),
  tags$head(
    tags$style(
      HTML("
        /* Custom styles for the waitbar */
        .shiny-notification {
          position: fixed; top: 85%; right: 65%;
          width : 20%;
          height=30%
          background-color: #eee;
          border: 4px solid #ccc;
        }
        
        .progress-bar {
          width: 100%;
          height: 100%;
          background-color:#607d8b;
          transition: width 0.3s ease-in-out;
        }
      ")
    )
  ),

  tabItems(
    #Upload B cell
    tabItem(tabName = "B_upload",
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Upload protein for B-cell epitopes",
                fileInput("B_protein", "Protein sequence (.fasta format)",multiple=FALSE,accept=c('.fasta','.fa')),
                actionButton("b_example_fasta","Load example")),
            br(),
            HTML("<h3> <b> Select one sequence of your choice </b></h3>"),
            box(width=NULL,status=NULL,solidHeader = TRUE, title = "Sequence Table",dataTableOutput('B_Table')),
            br(),br(),h3(textOutput("b_res_check"))
            
            
    ),   
    #Bcell epitope params
    tabItem(tabName = "B_params",
            fluidRow(
              column(width = 12, 
                     tabBox(width = NULL,
                            tabPanel(h4("BepiPred"),
                                     column(12,
                                            selectInput("bepipred_opt",h3("Select Type of Analysis"),choices=c("No 1"= 1,"No 2"=2),selected=1),
                                            p("Based on the BepiPred scores, Analysis No 1 performs the formation B cell epitopes with different length and Analysis No 2 enables the formation of B cell epitopes of standard length"),br(),
                                            numericInput("bepi_thres_num",h3("BepiPred Threshold"),value=0.1512,min=0,max=0.3),
                                        conditionalPanel(condition="input.bepipred_opt == 1",
                                            sliderInput("bepi_thres",NULL,min = 0, max = 0.3, value = 0.1512,step=0.0001),br(),
                                            numericInput("base_neg",h3('Subthreshold Amino Acid Inclusion Count'),value=0,min=0),
                                            p(" This parameter determines the maximum number of amino acids between predicted epitopes that can be included in the output, allowing the merging of separate epitopes into a larger epitopic region (for example if its' value is equal to “1” and there are two amino acids in a row with lower score than the threshold then these amino acids are not being taken into account and the high scored regions can not be joined, but if this parameter was set to '2' then those regions would join together as a larger epitope."),
                                            br(),br(),
                                            hidden(numericInput("second_thres_num",h3("Secondary threshold"),value=0.12,min=0,max=0.3),
                                                   sliderInput("second_thres",h3("Secondary threshold"),min=0,max=0.3,value=0.12,step=0.0001),
                                                   p(id="secondary_thres_message","Utilize the 'Secondary Threshold' to refine BepiPred's epitope predictions. Only amino acids with scores below the primary threshold will be included if their score exceeds this secondary threshold. For example, if you set the secondary threshold at 0.1, BepiPred will consider amino acids with scores below the primary threshold only if they score above 0.1, enhancing epitope accuracy")
                                            ),br(),
                                            numericInput("amino_region",h3("Minimum epitope length"),value=10),
                                            p("Set the 'Minimum Epitope Length' to specify the minimum consecutive amino acids required to define an epitope region."),
                                        ),
                                        conditionalPanel(condition="input.bepipred_opt==2", numericInput("amino_region",h3("Epitope length"),value=10),
                                                         p("Set the 'Epitope Length' to specify the required number of consecutive amino acids to be considered as an epitope region.")
                                        )
                                        ,br(),br(),br(),br(),HTML("<p><u> Reference </u></p>"),br(),p("Joakim Clifford, Magnus Haraldson Høie, Sebastian Deleuran, Bjoern Peters, Morten Nielsen and Paolo Marcatili BepiPred-3.0: Improved B-cell epitope prediction using protein language models doi: https://doi.org/10.1002/pro.4497")
                                     )),
                            tabPanel(h4("VaxiJen"),column(12,
                                                          selectInput("vax_target",h3("VaxiJen Target"),choices=c('bacteria','virus','tumour','parasite','fungal'),selected='bacteria'),
                                                          p("Use VaxiJen Target to specify the nature of the vaccine you are interested in, whether it's for bacteria, tumors, viruses, or other targets, for tailored antigenicity prediction"),
                                                          br(),
                                                          numericInput("vax_thres_num",h3("VaxiJen Threshold"),value=0.4,min=0,max=1),
                                                          sliderInput("vax_thres",NULL,min=0,max=1,value=0.4),
                                                          p("Set the VaxiJen Threshold to adjust the sensitivity of antigen prediction, controlling the minimum antigenicity score required for protein sequences to be considered as potential vaccine candidates.")
                                                          ,br(),br(),br(),br(),HTML("<p><u> Reference </u></p>"),br(),p("Doytchinova, I.A., Flower, D.R. VaxiJen: a server for prediction of protective antigens, tumour antigens and subunit vaccines. BMC Bioinformatics 8, 4 (2007). https://doi.org/10.1186/1471-2105-8-4")
                            )),
                            tabPanel(h4("ToxinPred"),column(12,
                                                            radioButtons("model_tox",h3("Model"),choices=c('1','2'),selected = '2'),
                                                            p("Model 1 stands for Machine Learning Prediction and 2 stands for Hybrid approach (ML + BLAST + MERCI)"),
                                                            br(),
                                                            numericInput('tox_thres_num',h3('Threshold ToxinPred'),value=0.6,min=0,max=1),
                                                            sliderInput('tox_thres',NULL,min=0,max=1,value=0.6),
                                                            p("Adjust the ToxinPred Threshold to fine-tune the sensitivity of toxin prediction. Higher threshold values will yield more stringent predictions, while lower values may include a broader range of potential toxins"),
                                                            
                                     br(),br(),br(),br(),HTML("<p><u> Reference </u></p>"),br(),p("Sharma N, Naorem LD, Jain S, Raghava GPS. ToxinPred2: an improved method for predicting toxicity of proteins. Brief Bioinform. 2022 Sep 20;23(5):bbac174. doi: 10.1093/bib/bbac174. PMID: 35595541.")
                            )),
                            tabPanel(h4("AlgPred"),column(12,
                                                          radioButtons("model_alg",h3("Model"),choices=c('1','2'),selected='2'),br(),
                                                          p("Model 1 stands for Machine Learning Prediction and 2 stands for Hybrid approach (ML + BLAST + MERCI)"),
                                                          br(),
                                                          numericInput('alg_thres_num',h3('Threshold AlgPred'),value=0.3,min=0,max=1),
                                                          sliderInput('alg_thres',NULL,min=0,max=1,value=0.3),
                                                          p("Adjust the AlgPred Threshold to control the sensitivity of allergenicity prediction. Higher threshold values yield more stringent predictions, while lower values may include a wider range of potential allergenic proteins."),
                                                          
                                                          br(),br(),br(),br(),HTML("<p><u> Reference </u></p>"),br(),p("Sharma N, Patiyal S, Dhall A, Pande A, Arora C, Raghava GPS. AlgPred 2.0: an improved method for predicting allergenic proteins and mapping of IgE epitopes. Brief Bioinform. 2021 Jul 20;22(4):bbaa294. doi: 10.1093/bib/bbaa294. PMID: 33201237.")
                                                          ))
                            
                     )))
    ),
    ##B cell execute
    tabItem(tabName='B_res',
            column(6,align='center',
                   checkboxGroupInput("b_softwares",h3("Softwares to include "),choices=list("VaxiJen" = 1,"ToxinPred" = 2,"AlgPred" = 3),selected=c(1,2,3)),
                   br(),br(),br(),br(),actionButton("b_exe","Execute!",class="buttons"),br(),br(),textOutput("bepi_text"),textOutput("vax_text"),textOutput("tox_text"),textOutput("alg_text")),
            column(6,align='center',box(width=NULL,status=NULL,solidHeader = TRUE, title = "Sequence Table",dataTableOutput('sel_B_Table')),hr(),br(),
                   HTML("<p><u>Parameters </u></p>"),
                   div(
                     style = "text-align: left;",
                     HTML("<p><u> BepiPred</u>: </p>"),
                     textOutput("bepi_params"),br(),
                     HTML("<p><u> VaxiJen</u>: </p> "),
                     textOutput("vaxi_params"),br(),
                     HTML("<p><u>ToxinPred</u>: </p>"),
                     textOutput("tox_params"),br(),
                     HTML("<p><u>AlgPred</u>: </p>"),
                     textOutput("alg_params")
                   )),br(),br(),br(),
            column(id='b_res_col1',width=12,align='center',HTML('<h4><u>Final Results</u></h4>'),br(),
                   dataTableOutput('res_df_b'),br(),br(),br(),
                   p("Move to Final Results Filter tab for results filtering"))
            
            
    ),
    tabItem(tabName = "B_ranked",
            fluidRow(
              column(width = 12, 
                     tabBox(width = NULL,
                            tabPanel(h4("BepiPred"),
                                     column(12,
                                            column(12,align='center',HTML('<h4><u>B cell epitope sequences</u></h4>'),br(),
                                                   dataTableOutput('bepi_df'),br(),
                                                   #HTML('<h4><u> BepiPred Interactive figure</u> </h4>'),
                                                   #uiOutput("figureOutput"),
                                                   HTML('<h4><u> BepiPred raw output</u></h4>'),br(),
                                                   dataTableOutput('bepi_csv')
                                            ))),
                            tabPanel(h4("VaxiJen"),
                                     column(12,align='center',HTML('<h4><u>VaxiJen Results</u></h4>'),br(),
                                            dataTableOutput('vaxi_df'),br())
                            ),
                            tabPanel(h4("ToxinPred"),
                                     column(12,align='center',HTML('<h4><u>ToxinPred Results</u></h4>'),br(),
                                            dataTableOutput('toxi_df'),br())
                            ),
                            tabPanel(h4("AlgPred"),
                                     column(12,align='center',HTML('<h4><u>AlgPred Results</u></h4>'),br(),
                                            dataTableOutput('alg_df'),br())
                            ),
                            tabPanel(h4('Filter'),
                                     column(12,align='center',HTML('<h4><u>Final Results</u></h4>'),br(),
                                            dataTableOutput('filter_df_b'),br(), uiOutput("sel_b_filter_ui"),actionButton("b_fin_fil","Filter!",class="buttons"),
                                            uiOutput("sel_b_ranking_ui"),actionButton("b_fin_rank","Rank!",class="buttons"),br(),br(),p("Download the results in CSV format for the Vaccine Designer pipeline"),dataTableOutput('filtered_df_b')))
                     )
              )
            )
    ),
    #Upload CTL epitopes
    tabItem(tabName = "C_upload",
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Upload protein for CTL epitopes",
                fileInput("C_protein", "Protein sequence (.fasta format)",multiple=FALSE,accept=c('.fasta','.fa')),
                actionButton("c_example_fasta","Load example")),
            br(),
            HTML("<h3> <b> Select one sequence of your choice </b></h3>"),
            box(width=NULL,status=NULL,solidHeader = TRUE, title = "Sequence Table",dataTableOutput('C_Table')),
            conditionalPanel(
              condition = "!is.null(input.C_protein)",br(),br(),p(textOutput("c_res_check")))
    ),   
    #CTL epitopes params
    tabItem(tabName = "C_params",
            fluidRow(
              column(width = 12, 
                     tabBox(width = NULL,
                            tabPanel(h4("netMHCpan"),
                                     column(12,column(4,align='center',
                                                      numericInput("c_length",h3("Length of CTL epitopes"),value= 9,min=8,max=12),
                                                      br(),p("Choose the MHC class I alleles of interest for peptide binding prediction by clicking one by one"),
                                                      selectInput("allel_opt",h3("Select species"),choices=c("HLA supertype representative"="HLA","Chimpanzee (Patr)"="Patr","Rhesus macaque (Mamu)"="Mamu","Pig (SLA)"="SLA","Mouse (H-2)"="H-2","Gorilla lowlands (Gogo)"="Gogo","BoLA"="BoLA","DLA"="DLA","EQCA"="Eqca","all"="all"),selected="all"),
                                                      h3("Select Allelles"),dataTableOutput("allelle_ctl_tab")),
                                            column(4,br(),br(),br(),br(),br(),br(),align='center',actionButton("ctl_allelle_btn", "Select!",class="buttons")),
                                            conditionalPanel(condition="!is.null(input.ctl_allelle_btn)",column(4,br(),br(),br(),br(),br(),br(),
                                                                                                                div(
                                                                                                                  style = "display: flex; justify-content: center;",
                                                                                                                  dataTableOutput("ctl_textarea"))))
                                            ),br(),br(),br(),br(),HTML("<p> Reference </p>"),br(),p("Birkir Reynisson, Bruno Alvarez, Sinu Paul, Bjoern Peters, Morten Nielsen, NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data, Nucleic Acids Research, Volume 48, Issue W1, 02 July 2020, Pages W449–W454, https://doi.org/10.1093/nar/gkaa379")

),
                            tabPanel(h4("VaxiJen"),column(12,
                                                          selectInput("c_vax_target",h3("VaxiJen Target"),choices=c('bacteria','virus','tumour','parasite','fungal'),selected='bacteria'),
                                                          p("Use VaxiJen Target to specify the nature of the vaccine you are interested in, whether it's for bacteria, tumors, viruses, or other targets, for tailored antigenicity prediction"),br(),
                                                          br(),
                                                          numericInput("c_vax_thres_num",h3("VaxiJen Threshold"),value=0.4,min=0,max=1),
                                                          sliderInput("c_vax_thres",NULL,min=0,max=1,value=0.4),
                                                          p("Set the VaxiJen Threshold to adjust the sensitivity of antigen prediction, controlling the minimum antigenicity score required for protein sequences to be considered as potential vaccine candidates."),
                                                          
                                                          br(),br(),br(),br(),HTML("<p><u> Reference </u></p>"),br(),p("Doytchinova, I.A., Flower, D.R. VaxiJen: a server for prediction of protective antigens, tumour antigens and subunit vaccines. BMC Bioinformatics 8, 4 (2007). https://doi.org/10.1186/1471-2105-8-4")
                            )),
                            tabPanel(h4("ToxinPred"),column(12,
                                                            radioButtons("c_model_tox",h3("Model"),choices=c('1','2'), selected = '2'),
                                                            p("Model 1 stands for Machine Learning Prediction and 2 stands for Hybrid approach (ML + BLAST + MERCI)"),
                                                            br(),
                                                            numericInput("c_tox_thres_num",h3("Threshold ToxinPred"),value=0.6,min=0,max=1),
                                                            sliderInput("c_tox_thres",NULL,min=0,max=1,value=0.6),
                                                            p("Adjust the ToxinPred Threshold to fine-tune the sensitivity of toxin prediction. Higher threshold values will yield more stringent predictions, while lower values may include a broader range of potential toxins"),
                                                            br(),br(),br(),br(),HTML("<p> Reference </p>"),br(),p("Sharma N, Naorem LD, Jain S, Raghava GPS. ToxinPred2: an improved method for predicting toxicity of proteins. Brief Bioinform. 2022 Sep 20;23(5):bbac174. doi: 10.1093/bib/bbac174. PMID: 35595541.")
                            )),
                            tabPanel(h4("AlgPred"),column(12,
                                                          radioButtons("c_model_alg",h3("Model"),choices=c('1','2'),selected = '2'),
                                                          p("Model 1 stands for Machine Learning Prediction and 2 stands for Hybrid approach (ML + BLAST + MERCI)"),
                                                          br(),
                                                          numericInput("c_alg_thres_num",h3('Threshold AlgPred'),value=0.3,min=0,max=1),
                                                          sliderInput("c_alg_thres",NULL,min=0,max=1,value=0.3),
                                                          p("Adjust the AlgPred Threshold to control the sensitivity of allergenicity prediction. Higher threshold values yield more stringent predictions, while lower values may include a wider range of potential allergenic proteins."),
                                                          br(),br(),br(),br(),HTML("<p> Reference </p>"),br(),p("Sharma N, Patiyal S, Dhall A, Pande A, Arora C, Raghava GPS. AlgPred 2.0: an improved method for predicting allergenic proteins and mapping of IgE epitopes. Brief Bioinform. 2021 Jul 20;22(4):bbaa294. doi: 10.1093/bib/bbaa294. PMID: 33201237.")
                            ))
                          
                     )))
    ),
    ##CTL execute
    tabItem(tabName='C_res',
            column(6,align='center',
                   checkboxGroupInput("c_softwares",h3("Softwares to include "),choices=list("VaxiJen" = 1,"ToxinPred" = 2,"AlgPred" = 3),selected=c(1,2,3)),
                   br(),br(),br(),br(),
                   actionButton("c_exe","Execute!",class="buttons"),br(),br(),textOutput('net_I_text'),textOutput('c_vax_text'),textOutput('c_tox_text'),textOutput('c_alg_text'),textOutput('immuno_text')),
            column(6,align='center',box(width=NULL,status=NULL,solidHeader = TRUE, title = "Sequence Table",dataTableOutput('sel_C_Table')),hr(),br(),
                   HTML("<p><u>Parameters </u><p>"),
                   div(
                     style = "text-align: left;",
                     HTML("<p><u> netMHCpan</u>: <p>"),
                     textOutput("c_net_params"),br(),
                     HTML("<p><u> VaxiJen</u>: <p> "),
                     textOutput("c_vaxi_params"),br(),
                     HTML("<p><u>ToxinPred</u>: <p>"),
                     textOutput("c_tox_params"),br(),
                     HTML("<p><u>AlgPred</u>: <p>"),
                     textOutput("c_alg_params")
                   )),
            br(),br(),br(),
            column(id='c_res_col1',width=12,align='center',HTML('<h4><u>Final Results</u></h4>'),br(),
                   dataTableOutput('res_df_c'),br(),br(),br(),
                   p("Move to Final Results Filter tab for results filtering"))
            
    ),
    tabItem(tabName = "C_ranked",
            fluidRow(
              column(width = 12, 
                     tabBox(width = NULL,
                            tabPanel(h4("netMHCpan"),
                                     column(12,align='center',HTML('<h4><u>CTL epitope sequences</u></h4>'),br(),
                                            dataTableOutput('net_I_df'),br()
                                     )),
                            tabPanel(h4("VaxiJen"),
                                     column(12,align='center',HTML('<h4><u>VaxiJen Results</u></h4>'),br(),
                                            dataTableOutput('c_vax_df'),br())
                            ),
                            tabPanel(h4("ToxinPred"),
                                     column(12,align='center',HTML('<h4><u>ToxinPred Results</u></h4>'),br(),
                                            dataTableOutput('c_toxi_df'),br()
                                     )
                            ),
                            tabPanel(h4("AlgPred"),
                                     column(12,align='center',HTML('<h4><u>AlgPred Results</u></h4>'),br(),
                                            dataTableOutput('c_alg_df'),br()
                                     )
                            ),
                            
                            tabPanel(h4('Filter'),
                                     column(12,align='center',HTML('<h4><u>Final Results</u></h4>'),br(),
                                            dataTableOutput('filter_df_c'),br(), uiOutput("sel_c_filter_ui"),actionButton("c_fin_fil","Filter!",class="buttons"),
                                            uiOutput("sel_c_ranking_ui"),actionButton("c_fin_rank","Rank!",class="buttons"),br(),br(),p("Download the results in CSV format for the Vaccine Designer pipeline"),dataTableOutput('filtered_df_c')
                                     )
                            )
                     )
              )
            )
    ),
    
    #Upload HTL epitopes
    tabItem(tabName = "H_upload",
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Upload protein for HTL epitopes",
                fileInput("H_protein", "Protein sequence (.fasta format)",multiple=FALSE,accept=c('.fasta','.fa')),
                actionButton("h_example_fasta","Load example")),
            br(),
            HTML("<h3> <b> Select one sequence of your choice </b></h3>"),
            box(width=NULL,status=NULL,solidHeader = TRUE, title = "Sequence Table",dataTableOutput('H_Table')),
            conditionalPanel(
              condition = "!is.null(input.H_protein)",br(),br(),p(textOutput("h_res_check")))
    ),   
    #HTL epitopes params
    tabItem(tabName="H_params",
            fluidRow(
              column(width = 12, 
                     tabBox(width = NULL,
                            tabPanel(h4("netMHCIIpan"),column(12,
                                                              column(4,align='center',
                                                                     numericInput("h_length",h3("Length of HTL epitopes"),value= 15),
                                                                     br(),p("Choose the MHC class II alleles of interest for peptide binding prediction by clicking one by one"),
                                                                     selectInput("allel_opt_h",h3("Select species"),choices=c("DRB1"="DRB1","DRB3"="DRB3","DRB4"="DRB4","DRB5"="DRB5","DP"="DP","DQ"="DQ","all"="all"),selected="all"),
                                                                     h3("Select Allelles"),dataTableOutput("allelle_htl_tab")),
                                                              column(4,align='center',br(),br(),br(),br(),br(),br(),actionButton("htl_allelle_btn", "Select!",class="buttons")),
                                                              conditionalPanel(condition="!is.null(input.htl_allelle_btn)",column(4,br(),br(),br(),br(),br(),br(),
                                                                                                                                  div(
                                                                                                                                    style = "display: flex; justify-content: center;",
                                                                                                                                    dataTableOutput("htl_textarea"))))
                                                              
                            ),br(),br(),br(),br(),HTML("<p> Reference </p>"),br(),p("Reynisson B, Barra C, Kaabinejadian S, Hildebrand WH, Peters B, Nielsen M. Improved Prediction of MHC II Antigen Presentation through Integration and Motif Deconvolution of Mass Spectrometry MHC Eluted Ligand Data. J Proteome Res. 2020 Jun 5;19(6):2304-2315. doi: 10.1021/acs.jproteome.9b00874. Epub 2020 Apr 30. PMID: 32308001.")

                            ),
                            tabPanel(h4("VaxiJen"),column(12,
                                                          selectInput("h_vax_target",h3("VaxiJen Target"),choices=c('bacteria','virus','tumour','parasite','fungal'),selected='bacteria'),br(),
                                                          p("Use VaxiJen Target to specify the nature of the vaccine you are interested in, whether it's for bacteria, tumors, viruses, or other targets, for tailored antigenicity prediction"),br(),
                                                          
                                                          numericInput("h_vax_thres_num",h3("VaxiJen Threshold"),value=0.4,min=0,max=1),
                                                          sliderInput("h_vax_thres",NULL,min=0,max=1,value=0.4),
                                                          p("Set the VaxiJen Threshold to adjust the sensitivity of antigen prediction, controlling the minimum antigenicity score required for protein sequences to be considered as potential vaccine candidates."),
                                                          
                                                          br(),br(),br(),br(),HTML("<p> Reference </p>"),br(),p("Doytchinova IA, Flower DR. VaxiJen: a server for prediction of protective antigens, tumour antigens and subunit vaccines. BMC Bioinformatics. 2007 Jan 5;8:4. doi: 10.1186/1471-2105-8-4. PMID: 17207271; PMCID: PMC1780059.")
                            )),
                            tabPanel(h4("ToxinPred"),column(12,
                                                            radioButtons("h_model_tox",h3("Model"),choices=c('1','2'),selected='2'),br(),
                                                            p("Model 1 stands for Machine Learning Prediction and 2 stands for Hybrid approach (ML + BLAST + MERCI)"),
                                                            br(),
                                                            numericInput("h_tox_thres_num",h3("ToxinPred Threshold"),min=0,max=1,value=0.6),
                                                            sliderInput("h_tox_thres",NULL,min=0,max=1,value=0.6),
                                                            p("Adjust the ToxinPred Threshold to fine-tune the sensitivity of toxin prediction. Higher threshold values will yield more stringent predictions, while lower values may include a broader range of potential toxins"),
                                                            
                                                            br(),br(),br(),br(),HTML("<p> Reference </p>"),br(),p("Sharma N, Naorem LD, Jain S, Raghava GPS. ToxinPred2: an improved method for predicting toxicity of proteins. Brief Bioinform. 2022 Sep 20;23(5):bbac174. doi: 10.1093/bib/bbac174. PMID: 35595541.")
                            )),
                            tabPanel(h4("AlgPred"),column(12,
                                                          radioButtons("h_model_alg",h3("Model"),choices=c('1','2'),selected='2'),br(),
                                                          p("Model 1 stands for Machine Learning Prediction and 2 stands for Hybrid approach (ML + BLAST + MERCI)"),
                                                          br(),
                                                          numericInput("h_alg_thres_num",h3("Threshold AlgPred"),value=0.3,min=0,max=1),
                                                          sliderInput('h_alg_thres',NULL,min=0,max=1,value=0.3),
                                                          p("Adjust the AlgPred Threshold to control the sensitivity of allergenicity prediction. Higher threshold values yield more stringent predictions, while lower values may include a wider range of potential allergenic proteins."),
                                                          br(),br(),br(),br(),HTML("<p> Reference </p>"),br(),p("Sharma N, Patiyal S, Dhall A, Pande A, Arora C, Raghava GPS. AlgPred 2.0: an improved method for predicting allergenic proteins and mapping of IgE epitopes. Brief Bioinform. 2021 Jul 20;22(4):bbaa294. doi: 10.1093/bib/bbaa294. PMID: 33201237.")
                            ))
                            #Parameters here
                            #tabPanel(h4("IFN epitopes"))
                     )))
    ),
    
    ##HTL execute ------------------- "IFN Epitopes" = 4
    tabItem(tabName='H_res',
            column(6,align='center',
                   checkboxGroupInput("h_softwares",h3("Softwares to include "),choices=list("VaxiJen" = 1,"ToxinPred" = 2,"AlgPred" = 3),selected=c(1,2,3)),
                   br(),br(),br(),br(),
                   actionButton("h_exe","Execute!",class="buttons"),br(),br(),textOutput("net_II_text"),textOutput("h_vax_text"),textOutput("h_tox_text"),textOutput("h_alg_text"),textOutput("ifn_text")),
            column(6,align='center',box(width=NULL,status=NULL,solidHeader = TRUE, title = "Sequence Table",dataTableOutput('sel_H_Table')),hr(),br(),
                   HTML("<p><u>Parameters </u></p>"),
                   div(
                     style = "text-align: left;",
                     HTML("<p><u> netMHCIIpan</u>: <p>"),
                     textOutput("h_net_params"),br(),
                     HTML("<p><u> VaxiJen</u>: <p> "),
                     textOutput("h_vaxi_params"),br(),
                     HTML("<p><u>ToxinPred</u>: <p>"),
                     textOutput("h_tox_params"),br(),
                     HTML("<p><u>AlgPred</u>: <p>"),
                     textOutput("h_alg_params")
                   )),br(),br(),br(),
            column(id='h_res_col1',width=12,align='center',HTML('<h4><u>Final Results</u></h4>'),br(),
                   dataTableOutput('res_df_h'),br(),br(),br(),
                   p("Move to Final Results Filter tab for results filtering"))
            ),
    tabItem(tabName = "H_ranked",
            fluidRow(
              column(width = 12, 
                     tabBox(width = NULL,
                            tabPanel(h4("netMHCIIpan"),
                                     column(12,align='center',HTML('<h4><u>HTL epitope sequences</u></h4>'),br(),
                                            dataTableOutput('net_II_df'),br()
                                     )),
                            tabPanel(h4("VaxiJen"),
                                     column(12,align='center',HTML('<h4><u>VaxiJen Results</u></h4>'),br(),
                                            dataTableOutput('h_vax_df'),br()
                                     )
                            ),
                            tabPanel(h4("ToxinPred"),
                                     column(12,align='center',HTML('<h4><u>ToxinPred Results</u></h4>'),br(),
                                            dataTableOutput('h_toxi_df'),br()
                                     )
                            ),
                            tabPanel(h4("AlgPred"),
                                     column(12,align='center',HTML('<h4><u>AlgPred Results</u></h4>'),br(),
                                            dataTableOutput('h_alg_df'),br()
                                     )
                            ),
                            #tabPanel(h4("IFN epitope"),
                            #        column(12,align='center',HTML('<h4><u>IFN epitope prediction Results</u></h4>'),br(),
                            #              dataTableOutput('ifn_df'),br()
                            #      )
                            #),
                            tabPanel(h4('Filter'),
                                     column(12,align='center',HTML('<h4><u>Final Results</u></h4>'),br(),
                                            dataTableOutput('filter_df_h'),br(), uiOutput("sel_h_filter_ui"),actionButton("h_fin_fil","Filter!",class="buttons"),
                                            uiOutput("sel_h_ranking_ui"),actionButton("h_fin_rank","Rank!",class="buttons"),br(),br(),p("Download the results in CSV format for the Vaccine Designer pipeline"),dataTableOutput('filtered_df_h')
                                     )
                            )
                     )
              )
            )
    ),
    
    
    
    #Revacc
    #Upload Epitopes
    tabItem(tabName = "Re_upload",
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Upload B cell epitopes",
                fileInput("B_rank_epitopes", "B cell epitopes (.csv format)",multiple=TRUE,accept='.csv'),
                br(),
                actionButton("b_example_csv","Load example")
            ),
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Upload CTL epitopes",
                fileInput("C_rank_epitopes", "CTL epitopes (.csv format)",multiple=TRUE,accept='.csv'),
                br(),
                actionButton("c_example_csv","Load example")
            ),
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Upload HTL epitopes",
                fileInput("H_rank_epitopes", "HTL epitopes (.csv format)",multiple=TRUE,accept='.csv'),
                br(),
                actionButton("h_example_csv","Load example")
            ),
            box(width=NULL,status=NULL,solidHeader=TRUE,title="B cell epitope candidates",
                dataTableOutput('B_epitopes_table'),dataTableOutput('B_epitopes_table_sel')),
            fluidRow(column(width=4,uiOutput("B_epitope_rank_ui")),column(width=4,actionButton("B_epitope_rank","Rank!",class="buttons")),column(width=4,actionButton("B_epitope_select","Select",class="buttons"))),
            box(width=NULL,status=NULL,solidHeader=TRUE,title="CTL epitope candidates",
                dataTableOutput('C_epitopes_table'),dataTableOutput('C_epitopes_table_sel')),
            fluidRow(column(width=4,uiOutput("C_epitope_rank_ui")),column(width=4,actionButton("C_epitope_rank","Rank!",class="buttons")),column(width=4,actionButton("C_epitope_select","Select",class="buttons"))),
            box(width=NULL,status=NULL,solidHeader=TRUE,title="HTL epitope candidates",
                dataTableOutput('H_epitopes_table'),dataTableOutput('H_epitopes_table_sel')),
            fluidRow(column(width=4,uiOutput("H_epitope_rank_ui")),column(width=4,actionButton("H_epitope_rank","Rank!",class="buttons")),column(width=4,actionButton("H_epitope_select","Select",class="buttons")))
            
    ),
    #Evaluation tab
    tabItem(tabName="Re_eval",
            fluidRow(
              column(width = 12, 
                     tabBox(width = NULL,
                            tabPanel(h5("Vaccine Components"),column(4,
                                                                     numericInput("b_cell_num",h3("Number of B-cell epitopes"),value=2),
                                                                     br(),br(),
                                                                     numericInput("ctl_num",h3("Number of CTL epitopes"),value=2),
                                                                     br(),br(),
                                                                     numericInput("htl_num",h3("Number of HTL epitopes"),value=2),hr()),
                                     column(4,textInput("bcell_linker",h3("B-cell epitope linker"),value="KK"),
                                            textInput("ctl_linker",h3("CTL epitope linker"),value="AAY"),
                                            textInput("htl_linker",h3("HTL epitope linker"),value="GPGPG")),
                                     column(4,textInput("N_term",h3("N terminus Sequence")),
                                            selectInput("order_epitopes",h3("Order for epitope combination"),choices=list("B - C - H" = 1,"B - H - C" = 2,"H - B - C" = 3, "H - C - B" = 4,"C - H - B" = 5,"C - B - H" = 6),selected=1)
                                            ,column(4,align='center',actionButton("generation_epitope","Generate!",class="buttons"),br(),br()),dataTableOutput("Multiepitope_sequences"))
                            ),
                            tabPanel(h5('Filters'),column(4,
                                                          checkboxGroupInput("filters_vac",h3("Filtering for certain properties"),choices=list("VaxiJen" = 1,"ToxinPred" = 2,"AlgPred" = 3,"ProtParam" = 4),selected=c(1,2)),
                                                          br(),br(),
                                                          numericInput("num_vac",h3("Number of final multiepitope sequences"),value=50)),
                                     column(4,
                                            br(),selectInput("gene_seqs",h3("Sequences to filter"),choices=list("Use the generated "= '1',"Upload new sequences"='2'),selected='1'),
                                            hidden(fileInput("seqs_file", "Upload File",multiple=FALSE,accept='.txt')),
                                            box(width=NULL,status=NULL,solidHeader=TRUE,title='Selected Sequences to filter',dataTableOutput("seqs_tab"))),
                                     column(4,align='center',
                                            br(),actionButton("filter_seqs","Filter!",class="buttons"),dataTableOutput("vacc_seqs_df"))
                                     
                            ),
                            tabPanel(h5("VaxiJen"),column(12,
                                                          selectInput("vacc_vax_target",h3("VaxiJen Target"),choices=c('bacteria','virus','tumour','parasite','fungal'),selected='bacteria'),br(),
                                                          p("Use VaxiJen Target to specify the nature of the vaccine you are interested in, whether it's for bacteria, tumors, viruses, or other targets, for tailored antigenicity prediction"),br(),
                                                          
                                                          numericInput("vacc_vax_thres_num",h3("VaxiJen Threshold"),value=0.4,min=0,max=1),
                                                          sliderInput("vacc_vax_thres",NULL,min=0,max=1,value=0.4),
                                                          p("Set the VaxiJen Threshold to adjust the sensitivity of antigen prediction, controlling the minimum antigenicity score required for protein sequences to be considered as potential vaccine candidates."),
                                                          
                                                          br(),br(),br(),br(),HTML("<p> Reference </p>"),br(),p("Sharma N, Patiyal S, Dhall A, Pande A, Arora C, Raghava GPS. AlgPred 2.0: an improved method for predicting allergenic proteins and mapping of IgE epitopes. Brief Bioinform. 2021 Jul 20;22(4):bbaa294. doi: 10.1093/bib/bbaa294. PMID: 33201237.")
                            )),
                            tabPanel(h5("ToxinPred"),column(12,
                                                            radioButtons("vacc_model_tox",h3("Model"),choices=c('1','2'),selected='2'),br(),
                                                            p("Model 1 stands for Machine Learning Prediction and 2 stands for Hybrid approach (ML + BLAST + MERCI)"),
                                                            br(),
                                                            numericInput("vacc_tox_thres_num",h3("Threshold ToxinPred"),value=0.6,min=0,max=1),
                                                            sliderInput('vacc_tox_thres',NULL,min=0,max=1,value=0.6),
                                                            p("Set the VaxiJen Threshold to adjust the sensitivity of antigen prediction, controlling the minimum antigenicity score required for protein sequences to be considered as potential vaccine candidates."),
                                                            
                                                            br(),br(),br(),br(),HTML("<p> Reference </p>"),br(),p("Sharma N, Naorem LD, Jain S, Raghava GPS. ToxinPred2: an improved method for predicting toxicity of proteins. Brief Bioinform. 2022 Sep 20;23(5):bbac174. doi: 10.1093/bib/bbac174. PMID: 35595541.")
    
                            )),
                            tabPanel(h5("AlgPred"),column(12,
                                                          radioButtons("vacc_model_alg",h3("Model"),choices=list('1'= 1,'2' = 2),selected=2),br(),
                                                          p("Model 1 stands for Machine Learning Prediction and 2 stands for Hybrid approach (ML + BLAST + MERCI)"),
                                                          br(),
                                                          numericInput("vacc_alg_thres_num",h3("Threshold AlgPred"),value=0.3,min=0,max=1),
                                                          sliderInput('vacc_alg_thres',NULL,min=0,max=1,value=0.3),
                                                          p("Adjust the AlgPred Threshold to control the sensitivity of allergenicity prediction. Higher threshold values yield more stringent predictions, while lower values may include a wider range of potential allergenic proteins."),
                                                          
                                                          br(),br(),br(),br(),HTML("<p> Reference </p>"),br(),p("Sharma N, Patiyal S, Dhall A, Pande A, Arora C, Raghava GPS. AlgPred 2.0: an improved method for predicting allergenic proteins and mapping of IgE epitopes. Brief Bioinform. 2021 Jul 20;22(4):bbaa294. doi: 10.1093/bib/bbaa294. PMID: 33201237.")
                                                          
                            )),
                            tabPanel(h5('Protparam'),column(12,
                                                            checkboxGroupInput("prot_proper",h3("Properties"),choices=list("Molecular weight" = 1,"Instability" = 2,"Aliphatic" = 3,"GRAVY" = 4),selected=2),
                                                            p("Choose the properties from Protparam software to be included in the Final Results"),
                                                            br(),br(),br(),br(),HTML("<p> Reference </p>"),br(),p("Wilkins MR, Gasteiger E, Bairoch A, Sanchez JC, Williams KL, Appel RD, Hochstrasser DF. Protein identification and analysis tools in the ExPASy server. Methods Mol Biol. 1999;112:531-52. doi: 10.1385/1-59259-584-7:531. PMID: 10027275.")
                                                            
                                                            
                            ))
                            # ,tabPanel(h5('DeepVacPred'),column(12,
                            #                                   checkboxGroupInput("deepvacpred_proper",h3("Properties"),choices=list("Analyze whole sequence" = 1,"Analyze per 30" = 2),selected=1),
                            #                                   numericInput("deepvac_thres_num",h3("DeepVacPred threshold"),value=0.95,min=0,max=1),
                            #                                   sliderInput("deepvac_thres",NULL,value=0.95,min=0,max=1))
                            # )
                     )))
    ),
    tabItem(tabName = "Re_results",
            box(width = NULL, status = "primary", solidHeader = TRUE, title="Final Sequences Results",
                dataTableOutput("final_results_vaccine"),
                br()
            )),
    tabItem(
      tabName = "readme",
      tags$head(
        tags$style(
          HTML("
        /* Define global CSS styles here */
        p {
          font-size: 20px;
          line-height: 1.5; /* Set the desired font size for all <p> elements */
          text-align: justify;
        }
        h2{
         font-size: 20px;
        }
        h3{
          font-size: 20px;
          
        }
        h4{
          font-size: 20px;
          
        }
        h5{
          font-size: 15px;
        }
        h6{
          font-size: 10px;
        }
        h7{
          font-size: 10px;
        }
      ")
        )
      ),
      
      h2("Vaccine Designer guidelines"),
      br(),
      h3("1. B epitope prediction"),
      br(),
      h4("1.1 Fasta protein sequence file upload"),
      br(),
      p("The first step for the B cell epitope prediction is the fasta file upload which contains the protein of interest. The upload procedure begins upon pressing the Browser button indicated by Figure 1. Proteins contained in the uploaded fasta file, are demonstrated to the following table with their respective ID and Sequence. Either the file contains one or more proteins, the user must select the row with the protein of interest as shown in Figure 1."),
      br(),
      imageOutput("img_001",height="50%"),
      h5("Figure 1: Fasta file upload and protein selection for B cell epitope prediction"),
      br(),
      h4("1.2 (Optional) B parameters"),
      br(),
      p("The next step of the process is to set the appropriate parameters in the Parameters tab according to the biological application. The parameters are set to default based on the instructions of VaxiJen, ToxinPred2 and AlgPred2 guidelines. For Bepipred the default parameters for epitope prediction procedure are set to Analysis No 1, which dictates that the epitope formed don't necessary have the same length, and the user can set the minimum length only. The Analysis No. 2 dictates that the epitopes will be formed based on a standard epitope length, moduled by the user. The interface for the parameters panel is shown in Figure 2."),
      br(),
      imageOutput("img_002",height="50%"),
      h5("Figure 2: Parameters Interface for the B cell epitope prediction"),
      br(),
      h4("1.3 B execute."),
      br(),
      p("After setting the parameters or using the default, the next step is the Execute tab panel. The tab panel is depicted in Figure 3. The first choice of the user is to determine which software to include for the analysis. BepiPred software is executed by default and the user can select the VaxiJen, ToxinPred and/or AlgPred software to proceed. The second step of the process is the pressing of the Execute button, depicted in Figure 3. Upon the successful execution of the software, the according message is shown below the Execute button."),
      imageOutput("img_003",height="50%"),
      h5("Figure 3: Execute Tab Panel"),
      br(),
      h4("1.4 B filters/download"),
      br(),
      p("Upon the execution of the pipeline, the Filter Results tab panel contains the results for the B cell epitope of the protein of interest. In the Filter tab, the final table of all scores is demonstrated, as shown in Figure 4. In this graphic interface the user can select the filters that can apply to the predicted epitopes based on the executed software. The selection of the filters is indicated in Figure 4, in Select Filters selection. After selecting the optimal filters for the analysis, the user should press the Filter button, as indicated in Figure 4."),
      imageOutput("img_004",height="50%"),
      h5("Figure 4: Filter application to B cell epitopes."),
      br(),
      p("Upon pressing the Filter button, a new table appears at the bottom of the page, which contains only the epitopes that meet the filtering criteria. As shown in Figure 5, the user should press the CSV button to download the table with the B cell epitopes results. Optionally, the user has the ability to rank the epitopes based on their predicted scores with the Rank button."),
      br(),
      imageOutput("img_005",height="50%"),
      h5("Figure 5: Download the B cell epitope results table."),
      br(),
      
      
      h3("2. CTL epitope prediction"),
      br(),
      h4("2.1 CTL upload"),
      br(),
      p("Like in the B epitope prediction pipeline, the first step of the process is to upload the protein of interest in the CTL epitopes Upload panel."),
      p("Upon uploading, the sequences, just as before, are demonstrated in the arising table, and the user should select the sequence of interest, as depicted in Figure 6."),
      br(),
      imageOutput("img_006",height="50%"),
      h5("Figure 6: Fasta upload and protein sequence selection for CTL epitopes."),
      br(),
      
      h4("2.2 CTL parameters"),
      br(),
      p("The next step of the process is the selection of MHC I alleles in the Parameters panel of CTL epitopes tab. The parameters of VaxiJen, ToxinPred and AlgPred are optional for the user. In the netMHCpan tab in the Parameters section, the user must manually select the alleles of interest. The selection process is indicated in Figure 7."),
      br(),
      imageOutput("img_007",height="50%"), 
      h5("Figure 7: Allele selection and parameter setting for the CTL epitope prediction."),
      br(),
      h4("2.3 CTL execute"),
      br(),
      p("After selecting the alleles and setting the parameters according to the application, the user should proceed to the Execute tab panel and select the software to include and then press the Execute! button."),
      br(),
      imageOutput("img_008",height="50%"),
      h5("Figure 8: Execute panel for the CTL epitope prediction."),
      br(),
      h4("2.4 CTL filters/download"),
      br(),
      p("After the execution of the software, the results of CTL epitopes are stored in the Filter Results tab panel. To view the table with all the results, select the Filter panel, on the page as shown in Figure 9."),
      br(),
      imageOutput("img_009",height="50%"),
      h5("Figure 9: Filter panel for the CTL epitopes."),
      br(),
      p("In the bottom of this page the user can select as before the filters that can apply to the predicted epitopes and press the Filter! button. After filtering button, as new table arises in the bottom of the page as before, and the user can rank the epitopes based on their predicted scores and/or download the table in CSV format, indicated in Figure 10."),
      br(),
      imageOutput("img_010",height="50%"),
      h5("Figure 10: Filter application and Download the CTL epitope results table."),
      br(),
      
      
      h3("3. HTL epitope prediction"),
      br(),
      h4("3.1 HTL upload"),
      br(),
      p("Just as the previous pipelines, the first step is the upload of a fasta protein file. The contained protein(s) appear in a new table. The user must select the protein of interest as indicated with the bottom red arrow in Figure 11."),
      br(),
      imageOutput("img_011",height="50%"),
      h5("Figure 11: Fasta file upload and protein selection for CTL epitope prediction."),
      br(),
      h4("3.2 HTL parameters"),
      br(),
      p("The next step of the process is the Parameters panel in the HTL epitopes section. The parameter setting for VaxiJen, ToxinPred, and AlgPred are optional as before. The necessary parameters need to be set, are the selection of the MHC II alleles. The selection procedure is identical to the procedure in the CTL section. The first step is the manual selection of the alleles based on the application, and second is the Select! button, as shown in Figure 12."),
      br(),
      imageOutput("img_012",height="50%"),
      h5("Figure 12: Allele selection and parameter setting for the HTL epitope prediction."),
      br(),
      h4("3.3 HTL execute"),
      br(),
      p("After setting the necessary parameters for the software, the next step is the execution. In the Execute tab panel of HTL epitopes section, the user, like the previous pipelines, selects the software to include and then needs to press the Execute button as depicted in Figure 13. Upon the successful execution of the software, the message shown in Figure 13 appears."),
      br(),
      imageOutput("img_013",height="50%"),
      h5("Figure 13: Execute panel for the HTL epitope prediction."),
      br(),
      h4("3.4 HTL filter/download"),
      br(),
      p("The final steps are the filter application and download of the results, like previous pipelines. The user selects the necessary filters and then proceeds with the Filter button at the bottom of the page, as shown in Figure 14. Finally, it is necessary to download the filtered results in CSV format indicated with the red arrow in Figure 15."),
      imageOutput("img_014",height="50%"),
      h5("Figure 14: Filter application for HTL epitope results table."),
      br(),
      img("img_015",height="50%"),
      h5("Figure 15: Download of the HTL epitope filtered results table."),
      br(),
      h3("4. Vaccine Designer - Multiepitope sequence generation"),
      br(),
      h4("4.1 Vaccine Designer csv upload and epitope selection"),
      br(),
      p("After acquiring all the necessary CSV result files of epitopes for the protein or proteins of interest, the user needs to proceed to the Vaccine Designer section, to Upload tab panel. This page enables the upload of multiple csv files for B cell, CTL, and HTL epitopes in the respective boxes. The upload procedure begins by pressing the Browse button as indicated in Figure 16."),
      br(),
      imageOutput("img_016",height="50%"),
      h5("Figure 16: Upload of CSV files of the results for the B cell, CTL, and HTL epitopes."),
      br(),
      p("Upon completion of the upload process, the user can select the epitopes that will be part of the final multiepitope sequences. The selection procedure is demonstrated in the Figure 17. Firstly the user needs to rank the predicted epitopes based on their scores by pressing the Rank button. After ranking, the user can manually select the epitope of their choice in their order of preference. After selecting the epitopes, the user needs to press the Select button, and a new table of selected sequences will appear. If the user does not select any epitope, the Select button will appear all the epitopes."),
      br(),
      imageOutput("img_017",height="50%"),
      h5("Figure 17: Epitope selection for the multiepitope sequences"),
      br(),
      h4("4.2 Vaccine Designer multiepitope sequences generation"),
      br(),
      p("After selecting the B cell, CTL, and/or HTL epitopes, the next steps are completed in the Evaluation tab panel in the Vaccine Designer section. In this page, the user can define how many B cell (B), CTL (C), or HTL (H) epitopes they want in every multiepitope sequence. The user also can define which linkers and/or which N terminus sequence will be used. Finally, the user selects the order of the epitope in the final sequences and presses the Generate! Button as indicated with a red arrow in Figure 18. After completion of all combinations of epitopes, the multiepitope sequence library will appear at the bottom of the screen."),
      br(),
      imageOutput("img_018",height="50%"),
      h5("Figure 18: Parameter setting for the generation of multiepitope sequences."),
      br(),
      h4("4.3 Vaccine Designer filter and selection of sequences"),
      br(),
      p("The next step of the process takes place in the Filters panel. The user defines the filters that they want to apply to every multiepitope sequence and then selects how many sequences they want as a result. After deciding the filters and setting the parameters of the filters in the other panels of the Evaluation tab, the user needs to press the Filter button. The whole procedure is depicted in Figure 19."),
      br(),
      imageOutput("img_019",height="50%"),
      h5("Figure 19: Filter application and the number of sequences to select from the multiepitope sequence library."),
      br(),
      h4("4.4 Vaccine Designer results"),
      br(),
      p("Upon successful completion of the procedure, the sequences that passed the filtering steps are stored in the Final Results tab panel."),
      br(),
      imageOutput("img_020",height="50%"),
      h5("Figure 20: Download of the final multiepitope sequences table.")
    ),
    tabItem(
      tabName="about",
      h2("Overview"),
      br(),
      p("Vaccine Designer is an innovative and comprehensive bioinformatics software solution designed to empower researchers and immunologists in the field of vaccine development and immunotherapy. This cutting-edge software offers a suite of advanced tools meticulously crafted to predict B-cell epitopes, CTL (Cytotoxic T Lymphocyte) epitopes, and HTL (Helper T Lymphocyte) epitopes from protein sequences, enabling the precise design and construction of multiepitope vaccines.Vaccine Designer not only simplifies epitope prediction but also provides a user-friendly interface for optimizing vaccine candidates, allowing researchers to harness the power of epitope-based immunization strategies. With its robust features, customizable vaccine design options, and detailed reporting capabilities, Vaccine Designer stands at the forefront of bioinformatics software, facilitating breakthroughs in vaccine development, cancer immunotherapy, and personalized medicine.")
      )
      
    
    
  )
)  
ui <- dashboardPage(
  header,
  sidebar,
  body,
  skin="blue",
  footer = dashboardFooter(
    left = "By D. Trygoniaris",
    right = "Medicine AUTH, 2023"
  )
)