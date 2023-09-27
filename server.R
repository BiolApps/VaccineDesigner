library(shiny)
library(shinydashboard)
library(DT)
library(Biostrings)
library(shinyFiles)
library(reticulate)
library(dplyr)
library(shinyjs)
#Set path
main_path <- getwd()

##UI

#Parameters
allelle_ctl_df <- read.csv('allelenames',header=FALSE)
allelle_ctl_df<- data.frame(Allelles=allelle_ctl_df$V1)


allelle_htl_df <- read.table("allelle_htl.txt")
colnames(allelle_htl_df) <- 'Allelles'

path_bepipred=file.path(main_path,'Local_Tools/BepiPred3_src')

path_res_bepi=file.path(main_path,'temp_data/Bepipred')

path_toxinpred=file.path(main_path,'Local_Tools/toxinpred2/toxinpred2')

path_res_toxinpred=file.path(main_path,'temp_data/Toxinpred')

path_algpred=file.path(main_path,'Local_Tools/algpred2_new/algpred2_new')

path_res_algpred=file.path(main_path,'temp_data/Algpred')

path_immunogen<-file.path(main_path,'Local_Tools/immunogenicity')

path_res_immunogen<-file.path(main_path)

allele_list_ctl <- 'allelenames'
allele_list_ctl <- read.csv(allele_list_ctl,header=FALSE)

allele_list_htl_path <-'allelelist.txt'
allele_list_htl <- read.csv(allele_list_htl_path,header=FALSE,sep="\t")


#Functions
convert_to_linux_path <- function(file_path) {
  # Check if the path contains a Windows-style drive letter (e.g., "C:")
  if (grepl("^[A-Za-z]:", file_path)) {
    # Extract the drive letter
    drive_letter <- substr(file_path, 1, 1)
    # Convert to the corresponding Linux mount point
    linux_path <- gsub(paste0("^", drive_letter, ":"), paste0("/mnt/", tolower(drive_letter)), file_path)
    return(linux_path)
  } else {
    # If the path doesn't contain a drive letter, return it as is
    return(file_path)
  }
}

Toxinpred_fnc <- function(df,path_toxinpred,path_res_toxinpred,model,threshold){
  #get current working directory
  main_path <- getwd()
  #Sequences data from dataframe
  seq=df$Sequences
  #New fasta file with the sequences
  file_path <- "temp.fasta"
  lines <- c()
  for (i in 1:length(seq)){
    name <- paste('>seq',i,sep="")
    seq_i<-seq[i]
    lines <-c(lines,name,seq_i)
  }
  writeLines(lines,file_path,sep='\n')
  #Path for the file
  file_path <- file.path(main_path,file_path)
  #Directory of toxinpred
  setwd(path_toxinpred)
  final_path <-file.path(path_res_toxinpred,'results_toxinpred2.csv')
  #Command for toxinpred
  command<-paste('python toxinpred2.py -i',file_path,'-d 2','-m',as.numeric(model),'-t',as.numeric(threshold),'-o',final_path)
  #Toxinpred run
  system(command)
  #Move to stored data folder
  setwd(path_res_toxinpred)
  res <- read.csv('results_toxinpred2.csv')
  setwd(main_path)
  res$Subject <- df$Sequences
  colnames(res)[colnames(res)=='Subject'] <- 'Sequences'
  return(res)
}

Algpred_fnc <- function(df,path_algpred,path_res_algpred,model,threshold){
  main_path <- getwd()
  seq=df$Sequences
  file_path <- "temp.fasta"
  lines <- c()
  for (i in 1:length(seq)){
    name <- paste('>seq',i,sep="")
    seq_i<-seq[i]
    lines <-c(lines,name,seq_i)
  }
  writeLines(lines,file_path,sep="\n")
  #Path for the file
  file_path <- file.path(main_path,file_path)
  #Directory of Algpred
  setwd(path_algpred)
  final_path <-file.path(path_res_algpred,'results_algpred2.csv')
  #Command for toxinpred
  command<-paste('python algpred2.py -i',file_path,'-d 2','-m',as.numeric(model),'-t',as.numeric(threshold),'-o',final_path)
  #Toxinpred run
  system(command)
  #Move to stored data folder
  setwd(path_res_algpred)
  res <- read.csv('results_algpred2.csv')
  setwd(main_path)
  res$Subject <- df$Sequences
  colnames(res)[colnames(res)=='Subject'] <- 'Sequences'
  return(res)
}

#Combined allel
combination_allel_fnc <- function(dat,dat2){
  seq1 <- dat$Sequences
  seq2<-dat2$Sequences
  for (i in 1:length(seq2)){
    is_in_list <- seq2[i] %in% seq1
    # Find the index of each element in list2 in list1
    if (is_in_list==TRUE){
      index <- match(seq2[i], seq1)
      #Alleles
      temp_alleles_1=dat$Alleles[index]
      #temp_alleles_1<-strsplit(temp_alleles_1,",")[[1]]
      temp_alleles_2=dat2$Alleles[i]
      fin_allele<-append(temp_alleles_1[[1]],temp_alleles_2[[1]])
      dat$Alleles[index][[1]]<-fin_allele
      
      #Weak_Binders
      temp_weak_1=dat$Weak_Binders_Counts[index]
      temp_weak_2=dat2$Weak_Binders_Counts[i]
      fin_weak <- as.numeric(temp_weak_1)+as.numeric(temp_weak_2)
      dat$Weak_Binders_Counts[index]<-fin_weak
      
      #Strong Binders
      temp_strong_1=dat$Strong_Binders_Counts[index]
      temp_strong_2=dat2$Strong_Binders_Counts[i]
      fin_strong<-as.numeric(temp_strong_1)+as.numeric(temp_strong_2)
      dat$Strong_Binders_Counts[index]<-fin_strong
      
      #Affinites
      temp_affin_1=subset(dat,select=c("Affinities(nM)"))[[1]][index]
      temp_affin_2=subset(dat2,select=c('Affinities(nM)'))[[1]][i]
      fin_affin<-append(temp_affin_1[[1]],temp_affin_2[[1]])
      dat$`Affinities(nM)`[index][[1]]<-fin_affin
      
      #score_el
      temp_el_score_1=dat$netMHCpan_Rank_Scores_EL[index]
      temp_el_score_2=dat2$netMHCpan_Rank_Scores_EL[i]
      fin_el<-append(temp_el_score_1[[1]],temp_el_score_2[[1]])
      dat$netMHCpan_Rank_Scores_EL[index][[1]]<-fin_el
      
      #score_ba
      temp_ba_score_1=dat$netMHCpan_Rank_Scores_BA[index]
      temp_ba_score_2=dat2$netMHCpan_Rank_Scores_BA[i]
      fin_ba <- append(temp_ba_score_1[[1]],temp_ba_score_2[[1]])
      dat$netMHCpan_Rank_Scores_BA[index][[1]]<-fin_ba
      
      #bind levels
      temp_bind_levels_1=subset(dat,select=c('Bind-levels'))[[1]][index]
      temp_bind_levels_2=subset(dat2,select=c('Bind-levels'))[[1]][i]
      fin_bind_levels<-append(temp_bind_levels_1[[1]],temp_bind_levels_2[[1]])
      dat$`Bind-levels`[index][[1]]<-fin_bind_levels
      
    }
    else{
      dat<-rbind(dat,dat2[i,])
    }
  }
  return(dat)
}
#Combined allel II
combination_allel_II_fnc <- function(dat,dat2){
  seq1 <- dat$Sequences
  seq2<-dat2$Sequences
  for (i in 1:length(seq2)){
    is_in_list <- seq2[i] %in% seq1
    # Find the index of each element in list2 in list1
    if (is_in_list==TRUE){
      index <- match(seq2[i], seq1)
      #Alleles
      temp_alleles_1=dat$Alleles[index]
      #temp_alleles_1<-strsplit(temp_alleles_1,",")[[1]]
      temp_alleles_2=dat2$Alleles[i]
      fin_allele<-append(temp_alleles_1[[1]],temp_alleles_2[[1]])
      dat$Alleles[index][[1]]<-fin_allele
      
      #Weak_Binders
      temp_weak_1=dat$Weak_Binders_Counts[index]
      temp_weak_2=dat2$Weak_Binders_Counts[i]
      fin_weak <- as.numeric(temp_weak_1)+as.numeric(temp_weak_2)
      dat$Weak_Binders_Counts[index]<-fin_weak
      
      #Strong Binders
      temp_strong_1=dat$Strong_Binders_Counts[index]
      temp_strong_2=dat2$Strong_Binders_Counts[i]
      fin_strong<-as.numeric(temp_strong_1)+as.numeric(temp_strong_2)
      dat$Strong_Binders_Counts[index]<-fin_strong
      
      #Affinites
      temp_affin_1=subset(dat,select=c("Affinities(nM)"))[[1]][index]
      temp_affin_2=subset(dat2,select=c('Affinities(nM)'))[[1]][i]
      fin_affin<-append(temp_affin_1[[1]],temp_affin_2[[1]])
      dat$`Affinities(nM)`[index][[1]]<-fin_affin
      
      #score_el
      temp_el_score_1=dat$netMHCIIpan_Rank_Scores_EL[index]
      temp_el_score_2=dat2$netMHCIIpan_Rank_Scores_EL[i]
      fin_el<-append(temp_el_score_1[[1]],temp_el_score_2[[1]])
      dat$netMHCIIpan_Rank_Scores_EL[index][[1]]<-fin_el
      
      #score_ba
      temp_ba_score_1=dat$netMHCIIpan_Rank_Scores_BA[index]
      temp_ba_score_2=dat2$netMHCIIpan_Rank_Scores_BA[i]
      fin_ba <- append(temp_ba_score_1[[1]],temp_ba_score_2[[1]])
      dat$netMHCIIpan_Rank_Scores_BA[index][[1]]<-fin_ba
      
      #bind levels
      temp_bind_levels_1=subset(dat,select=c('Bind-levels'))[[1]][index]
      temp_bind_levels_2=subset(dat2,select=c('Bind-levels'))[[1]][i]
      fin_bind_levels<-append(temp_bind_levels_1[[1]],temp_bind_levels_2[[1]])
      dat$`Bind-levels`[index][[1]]<-fin_bind_levels
      
    }
    else{
      dat<-rbind(dat,dat2[i,])
    }
  }
  return(dat)
}  


#netMHCpan function
netmhc_fnc <- function(seq,num_peptides,allel,allele_list_ctl){
  #Main path
  #net_path <- "/mnt/c/Users/30697/Desktop/netMHCpan-4.1/Linux_x86_64/bin/netMHCpan"
  net_path <- "netMHCpan"
  main_path <- getwd()
  #Fasta file construction
  file_path <- 'temp1.fasta'
  file_path <- file.path(main_path,file_path)
  lines<-c('>seq1',seq)
  writeLines(lines,file_path,sep="\n")
  file_path=convert_to_linux_path(file_path)
  
  if (length(allel)<=10){
    #Allele format for command from list allele
    allel_ls <- allel
    allel <- paste(allel,collapse=",")
    #Command and run the software
    command=paste(net_path,'-allname',file.path(main_path,'netMHCpan-4.1/Linux_x86_64/data/allelenames'),'-thrfmt',file.path(main_path,'netMHCpan-4.1/Linux_x86_64/data/threshold/%s.thr.%s'),'-version',file.path(main_path,'netMHCpan-4.1/Linux_x86_64/data/version'),'-rdir',file.path(main_path,'netMHCpan-4.1/Linux_x86_64'),'-syn',file.path(main_path,'netMHCpan-4.1/Linux_x86_64/data/synlist.bin'),'-hlapseudo',file.path(main_path,'netMHCpan-4.1/Linux_x86_64/data/MHC_pseudo.dat'),'-s -l',as.character(num_peptides),'-BA -a',allel,file_path,'> res.out')
    
    system(command)
    final_path<-file.path(main_path,'res.out')
    #Mapping allele
    inds <- which(allele_list_ctl$V1 %in% allel_ls)
    map_allele <- allele_list_ctl$V2[inds]
    
    #Script import and run
    net_script <- import("netMHCpan_script")
    res <- net_script$netMHCpan_fnc('res.out',map_allele)
    return(res)
  }
  else{
    # Specify the size of each sublist
    sublist_size <- 10
    # Calculate the number of sublists needed
    num_sublists <- ceiling(length(allel) / sublist_size)
    # Initialize an empty list to store the sublists
    sublists_allel<- list()
    # Split the original list into sublists
    for (i in 1:num_sublists) {
      start_idx <- (i - 1) * sublist_size + 1
      end_idx <- min(i * sublist_size, length(allel))
      sublist <- allel[start_idx:end_idx]
      sublists_allel[[i]] <- sublist
    }
    for (i in 1:num_sublists){
      allel<-sublists_allel[[i]]
      allel_ls <- allel
      allel <- paste(allel,collapse=",")
      #Command and run the software
      command=paste(net_path,'-allname',file.path(main_path,'netMHCpan-4.1/Linux_x86_64/data/allelenames'),'-thrfmt',file.path(main_path,'netMHCpan-4.1/Linux_x86_64/data/threshold/%s.thr.%s'),'-version',file.path(main_path,'netMHCpan-4.1/Linux_x86_64/data/version'),'-rdir',file.path(main_path,'netMHCpan-4.1/Linux_x86_64'),'-syn',file.path(main_path,'netMHCpan-4.1/Linux_x86_64/data/synlist.bin'),'-hlapseudo',file.path(main_path,'netMHCpan-4.1/Linux_x86_64/data/MHC_pseudo.dat'),'-s -l',as.character(num_peptides),'-BA -a',allel,file_path,'> res.out')
      system(command)
      final_path<-file.path(main_path,'res.out')
      #Mapping allele
      inds <- which(allele_list_ctl$V1 %in% allel_ls)
      map_allele <- allele_list_ctl$V2[inds]
      #Script import and run
      net_script <- import("netMHCpan_script")
      res <- net_script$netMHCpan_fnc('res.out',map_allele)
      if (i==1){
        fin_dat <- res
      }
      else{
        fin_dat <- combination_allel_fnc(res,fin_dat)
      }
      
    }
    return(fin_dat)
  }
}

####Immunogenicity function
immunogen_fnc <- function(df,path_immunogen,path_res_immunogen){
  main_path <- getwd()
  seq=df$Sequences
  #Check Sequences if they have X in the end and replace it with N
  matching_indices <- grep("X", seq)
  seq[matching_indices] <- gsub("X", "N", seq[matching_indices])
  #Write fasta file
  file_path <- "temp.fa"
  writeLines(seq,file_path,sep="\n")
  #Path for the file
  file_path <- file.path(main_path,file_path)
  
  setwd(path_immunogen)
  path_res=file.path(path_res_immunogen,'res_immuno.out')
  
  command<-paste('python predict_immunogenicity.py',file_path)
  a<-system(command,intern=TRUE)
  #Output
  setwd(path_res_immunogen)
  lin=a[5:length(a)]
  res_path <-'res_immuno.csv'
  writeLines(lin,res_path,sep="\n")
  res_path<-file.path(path_res_immunogen,'res_immuno.csv')
  data<-read.csv('res_immuno.csv',header=FALSE)
  
  setwd(main_path)
  colnames(data)<-c('Sequences','Peptide_Length','Immunogenicity_score')
  data$Immunogenicity <-ifelse(data$Immunogenicity_score>0,TRUE,FALSE)
  df <- select(data,c('Sequences','Immunogenicity_score','Immunogenicity'))
  return(df)
}

#netMHCIIpan function
netmhcII_fnc <- function(seq,num_peptides,allel,allele_list_htl){
  #Main path
  net_path <- "netMHCIIpan"
  main_path <- getwd()
  #Fasta file construction
  file_path <- 'temp1.fasta'
  file_path <- file.path(main_path,file_path)
  lines<-c('>seq1',seq)
  writeLines(lines,file_path,sep="\n")
  file_path=convert_to_linux_path(file_path)
  
  #Allele format for command from list allele
  if (length(allel)<=10){
    allel_ls <- allel
    allel <- paste(allel,collapse=",")
    #Command and run the software
    command=paste(net_path,'-s -length',as.character(num_peptides),'-BA -a',allel,'-f',file_path,'> res_HTL.out')
    system(command)
    final_path<-file.path(main_path,'res_HTL.out')
    
    #Mapping allele
    inds <- which(allele_list_htl$V1 %in% allel_ls)
    map_allele <- allele_list_htl$V2[inds]
    
    #Script import and run
    net_script <- import("netMHCIIpan_script")
    res <- net_script$netMHCIIpan('res_HTL.out',map_allele)
    
    return(res)
  }else{
    sublist_size <- 10
    # Calculate the number of sublists needed
    num_sublists <- ceiling(length(allel) / sublist_size)
    # Initialize an empty list to store the sublists
    sublists_allel<- list()
    # Split the original list into sublists
    for (i in 1:num_sublists) {
      start_idx <- (i - 1) * sublist_size + 1
      end_idx <- min(i * sublist_size, length(allel))
      sublist <- allel[start_idx:end_idx]
      sublists_allel[[i]] <- sublist
    }
    for (i in 1:num_sublists){
      allel<-sublists_allel[[i]]
      allel_ls <- allel
      allel <- paste(allel,collapse=",")
      #Command and run the software
      command=paste(net_path,'-s -length',as.character(num_peptides),'-BA -a',allel,'-f',file_path,'> res_HTL.out')
      system(command)
      final_path<-file.path(main_path,'res_HTL.out')
      
      #Mapping allele
      inds <- which(allele_list_htl$V1 %in% allel_ls)
      map_allele <- allele_list_htl$V2[inds]
      
      #Script import and run
      net_script <- import("netMHCIIpan_script")
      res <- net_script$netMHCIIpan('res_HTL.out',map_allele)
      
      if (i==1){
        fin_dat <- res
      }
      else{
        fin_dat <- combination_allel_II_fnc(res,fin_dat)
      }
      
    }
    return(fin_dat)
  }
}

server <- function(input,output,session) {
  Sys.setlocale("LC_NUMERIC", "C")
  ###Shinyjs
  shinyjs::hideElement("B_ranked")
  ###Sliders
  observeEvent(input$bepi_thres,{
    updateNumericInput(session, "bepi_thres_num", value = input$bepi_thres)
  })
  observeEvent(input$bepi_thres_num,{
    updateSliderInput(session, "bepi_thres", value = input$bepi_thres_num)
  })
  observeEvent(input$second_thres,{
    updateNumericInput(session,"second_thres_num",value = input$second_thres)
  })
  observeEvent(input$second_thres_num,{
    updateSliderInput(session,"second_thres",value = input$second_thres_num)
  })
  observeEvent(input$vax_thres,{
    updateNumericInput(session,"vax_thres_num",value = input$vax_thres)
  })
  observeEvent(input$vax_thres_num,{
    updateNumericInput(session,"vax_thres",value = input$vax_thres_num)
  })
  observeEvent(input$tox_thres,{
    updateNumericInput(session,"tox_thres_num",value = input$tox_thres)
  })
  observeEvent(input$tox_thres_num,{
    updateNumericInput(session,"tox_thres",value = input$tox_thres_num)
  })
  observeEvent(input$alg_thres,{
    updateNumericInput(session,"alg_thres_num",value = input$alg_thres)
  })
  observeEvent(input$alg_thres_num,{
    updateNumericInput(session,"alg_thres",value=input$alg_thres_num)
  })
  observeEvent(input$c_vax_thres,{
    updateNumericInput(session,"c_vax_thres_num",value=input$c_vax_thres)
  })
  observeEvent(input$c_vax_thres_num,{
    updateNumericInput(session,"c_vax_thres",value=input$c_vax_thres_num)
  })
  observeEvent(input$c_tox_thres,{
    updateNumericInput(session,"c_tox_thres_num",value=input$c_tox_thres)
  })
  observeEvent(input$c_tox_thres_num,{
    updateNumericInput(session,"c_tox_thres",value=input$c_tox_thres_num)
  })
  observeEvent(input$c_alg_thres,{
    updateNumericInput(session,"c_alg_thres_num",value=input$c_alg_thres)
  })
  observeEvent(input$c_alg_thres_num,{
    updateNumericInput(session,"c_alg_thres",value=input$c_alg_thres_num)
  })
  observeEvent(input$h_vax_thres,{
    updateNumericInput(session,"h_vax_thres_num",value=input$h_vax_thres)
  })
  observeEvent(input$h_vax_thres_num,{
    updateNumericInput(session,"h_vax_thres",value=input$h_vax_thres_num)
  })
  observeEvent(input$h_tox_thres,{
    updateNumericInput(session,"h_tox_thres_num",value=input$h_tox_thres)
  })
  observeEvent(input$h_tox_thres_num,{
    updateNumericInput(session,"h_tox_thres",value=input$h_tox_thres_num)
  })
  observeEvent(input$h_alg_thres,{
    updateNumericInput(session,"h_alg_thres_num",value=input$h_alg_thres)
  })
  observeEvent(input$h_alg_thres_num,{
    updateNumericInput(session,"h_alg_thres",value=input$h_alg_thres_num)
  })
  observeEvent(input$vacc_vax_thres,{
    updateNumericInput(session,"vacc_vax_thres_num",value=input$vacc_vax_thres)
  })
  observeEvent(input$vacc_vax_thres_num,{
    updateNumericInput(session,"vacc_vax_thres",value=input$vacc_vax_thres_num)
  })
  observeEvent(input$vacc_tox_thres,{
    updateNumericInput(session,"vacc_tox_thres_num",value=input$vacc_tox_thres)
  })
  observeEvent(input$vacc_tox_thres_num,{
    updateNumericInput(session,"vacc_tox_thres",value=input$vacc_tox_thres_num)
  })
  observeEvent(input$vacc_alg_thres,{
    updateNumericInput(session,"vacc_alg_thres_num",value=input$vacc_alg_thres)
  })
  observeEvent(input$vacc_alg_thres_num,{
    updateNumericInput(session,"vacc_alg_thres",value=input$vacc_alg_thres_num)
  })
  # observeEvent(input$deepvac_thres,{
  #   updateNumericInput(session,"deepvac_thres_num",value=input$deepvac_thres)
  # })
  # observeEvent(input$deepvac_thres_num,{
  #   updateNumericInput(session,"deepvac_thres",value=input$deepvac_thres_num)
  # })
  
  ###Upload protein for B cell epitopes
  b_sequences <- reactive({
    if (is.null(input$B_protein) && is.null(input$b_example_fasta))
      return(NULL)
    if (is.null(input$B_protein)){
      req(input$b_example_fasta)
      fasta <- readAAStringSet('ATL_LNK_conserved.fasta')
      fasta_df <- data.frame(ID = names(fasta), Sequence = as.character(fasta), stringsAsFactors = FALSE)
      return(fasta_df)
    }
    fasta <- readAAStringSet(input$B_protein$datapath)
    fasta_df <- data.frame(ID = names(fasta), Sequence = as.character(fasta), stringsAsFactors = FALSE)
    
    return(fasta_df)
  })
  output$B_Table <- renderDataTable({
    b_sequences() %>% select(-ID)
  },options=list(scrollX=TRUE),selection = 'single')
  
  ###########################Check button sequence --- python script#########################
  python_protein_check <- import("Protein_check")
  #Bcell
  observeEvent(input$B_Table_rows_selected,{
    selectedRow <- input$B_Table_rows_selected
    if (length(selectedRow)>0){
      selectedData <- b_sequences()[selectedRow, ]
      res <- python_protein_check$Protein_check(selectedData)
      output$b_res_check <- renderText({return(res)})}
    else{
      output$b_res_check <-renderText({return("")})
    }
    
  })
  #CTL
  observeEvent(input$C_Table_rows_selected,{
    selectedRow <- input$C_Table_rows_selected
    if (length(selectedRow)>0){
      selectedData <- c_sequences()[selectedRow, ]
      res <- python_protein_check$Protein_check(selectedData)
      output$c_res_check <- renderText({return(res)})}
    else{
      output$c_res_check <-renderText({return("")})
    }
    
  })
  #HTL
  observeEvent(input$H_Table_rows_selected,{

    selectedRow <- input$H_Table_rows_selected
    if (length(selectedRow)>0){
      selectedData <- h_sequences()[selectedRow, ]
      res <- python_protein_check$Protein_check(selectedData)
      output$h_res_check <- renderText({return(res)})}
    else{
      output$h_res_check <-renderText({return("")})
    }
    
    
  })

  ###NEW PART
  sel_BData <- eventReactive(input$B_Table_rows_selected, {
    selectedRow <- input$B_Table_rows_selected
    if (length(selectedRow) > 0) {
      # Get the selected row from table 1
      selectedData <- b_sequences()[selectedRow, ]
      
      # Generate sample data for table 2 based on selected row
      data.frame(
        Sequence = selectedData$Sequence
      )
    } else {
      # If no row is selected, return an empty data frame
      data.frame(Sequence = character(), stringsAsFactors = FALSE)
    }
  })
  
  # Render table --- selection 
  output$sel_B_Table <- DT::renderDataTable({
    sel_BData()
  },options=list(scrollX=TRUE,pageLength=10))
  
  ##Parameters for B cell epitopes #############################################
  
  output$bepi_params <- renderText({
  if(input$bepipred_opt==1){
    if (input$base_neg > 0){
      text <- paste("Bepipred Threshold:",input$bepi_thres_num,"Base to neglect:",input$base_neg,"Secondary Threshold: ",input$second_thres_num,"Minimum Amino Region:",input$amino_region,sep="  ")
    }else{text <- paste("Bepipred Threshold:",input$bepi_thres_num,"Base to neglect:",input$base_neg,"Minimum Amino Region:",input$amino_region,sep="  ")}
  }else{
    text <-paste("Bepipred Threshold:",input$bepi_thres_num,"Epitope length:",input$amino_region,sep="  ")
    }
    text
  })
  output$vaxi_params <- renderText({
    if (1 %in% input$b_softwares){
      text <- paste("Target:",input$vax_target, "Threshold:",input$vax_thres_num,sep=" ")
    }else{ text <- NULL}
    text
  })
  output$tox_params <- renderText({
    if (2 %in% input$b_softwares){
      text <- paste("Model selected:",input$model_tox,"Threshold:",input$tox_thres_num,sep=" ")
    }else{text<-NULL}
    text
  })
  output$alg_params <- renderText({
    if (3 %in% input$b_softwares){
      text <- paste("Model selected:",input$model_alg,"Threshold:",input$alg_thres_num,sep=" ")
    }else{text<-NULL}
    text
  })
  ################################################
  
  ###Upload protein for CTL epitopes
  c_sequences <- reactive({
    if (is.null(input$C_protein) && is.null(input$c_example_fasta))
      return(NULL)
    if (is.null(input$C_protein)){
      req(input$c_example_fasta)
      fasta <- readAAStringSet('ATL_LNK_conserved.fasta')
      fasta_df <- data.frame(ID = names(fasta), Sequence = as.character(fasta), stringsAsFactors = FALSE)
      return(fasta_df)
    }
    fasta <- readAAStringSet(input$C_protein$datapath)
    fasta_df <- data.frame(ID = names(fasta), Sequence = as.character(fasta), stringsAsFactors = FALSE)
    
    return(fasta_df)
  })
  
  output$C_Table <- DT::renderDataTable({
    c_sequences() %>% select(-ID)
  },options=list(scrollX=TRUE),selection='single')
  
  ################## Selection CTL #####################################
  sel_CData <- eventReactive(input$C_Table_rows_selected, {
    selectedRow <- input$C_Table_rows_selected
    if (length(selectedRow) > 0) {
      # Get the selected row from table 1
      selectedData <- c_sequences()[selectedRow, ]
      
      # Generate sample data for table 2 based on selected row
      data.frame(
        Sequence = selectedData$Sequence
      )
    } else {
      # If no row is selected, return an empty data frame
      data.frame(Sequence = character(), stringsAsFactors = FALSE)
    }
  })
  
  # Render table --- selection 
  output$sel_C_Table <- DT::renderDataTable({
    sel_CData()
  },options=list(scrollX=TRUE,pageLength=10))
  
  #########################Parameters for CTL epitopes #############################################
  
  output$c_net_params<- renderText({
    
    text<-paste("Length:",input$c_length,sep=" ")
    
    text
  })
  output$c_vaxi_params <- renderText({
    if (1 %in% input$c_softwares){
      text <- paste("Target:",input$c_vax_target, "Threshold:",input$c_vax_thres_num,sep=" ")
    }else{ text <- NULL}
    text
  })
  output$c_tox_params <- renderText({
    if (2 %in% input$c_softwares){
      text <- paste("Model selected:",input$c_model_tox,"Threshold:",input$c_tox_thres_num,sep=" ")
    }else{text<-NULL}
    text
  })
  output$c_alg_params <- renderText({
    if (3 %in% input$c_softwares){
      text <- paste("Model selected:",input$c_model_alg,"Threshold:",input$c_alg_thres_num,sep=" ")
    }else{text<-NULL}
    text
  })
  
  ##############################################################################################################################
  
  ###Upload protein for HTL epitopes
  h_sequences <- reactive({
    if (is.null(input$H_protein) && is.null(input$h_example_fasta))
      return(NULL)
    if (is.null(input$H_protein)){
      req(input$h_example_fasta)
      fasta <- readAAStringSet('ATL_LNK_conserved.fasta')
      fasta_df <- data.frame(ID = names(fasta), Sequence = as.character(fasta), stringsAsFactors = FALSE)
      
      return(fasta_df)
    }
    fasta <- readAAStringSet(input$H_protein$datapath)
    fasta_df <- data.frame(ID = names(fasta), Sequence = as.character(fasta), stringsAsFactors = FALSE)
    
    return(fasta_df)
  })
  
  output$H_Table <- DT::renderDataTable({
    h_sequences() %>% select(-ID)
  },options=list(scrollX=TRUE),selection='single')
  
  ### Selection
  sel_HData <- eventReactive(input$H_Table_rows_selected, {
    selectedRow <- input$H_Table_rows_selected
    if (length(selectedRow) > 0) {
      # Get the selected row from table 1
      selectedData <- h_sequences()[selectedRow, ]
      
      # Generate sample data for table 2 based on selected row
      data.frame(
        Sequence = selectedData$Sequence
      )
    } else {
      # If no row is selected, return an empty data frame
      data.frame(Sequence = character(), stringsAsFactors = FALSE)
    }
  })
  
  # Render table --- selection 
  output$sel_H_Table <- DT::renderDataTable({
    sel_HData()
  },options=list(scrollX=TRUE,pageLength=10))
  
  ###################################Parameters for HTL epitopes #############################################
  
  output$h_net_params<- renderText({
    
    text<-paste("Length:",input$h_length,sep=" ")
    text
  })
  output$h_vaxi_params <- renderText({
    if (1 %in% input$h_softwares){
      text <- paste("Target:",input$h_vax_target, "Threshold:",input$h_vax_thres_num,sep=" ")
    }else{text<-NULL}
    text
  })
  output$h_tox_params <- renderText({
    if (2 %in% input$h_softwares){
      text <- paste("Model selected:",input$h_model_tox,"Threshold:",input$h_tox_thres_num,sep=" ")
    }else{text<-NULL}
    text
  })
  output$h_alg_params <- renderText({
    if (3 %in% input$h_softwares){
      text <- paste("Model selected:",input$h_model_alg,"Threshold:",input$h_alg_thres_num,sep=" ")
    }else{text<-NULL}
    text
  })
  
  ######################################################################################################
  
  ####Select Allele CTL
  observeEvent(input$allel_opt,{
    data <- allelle_ctl_df
    opt <- input$allel_opt
    if(opt=="all"){
      allelle_f <- data
    }
    else{
    allel_sel <- data[grep(opt,data$Allelles),]
    allelle_f <-data.frame(Allelles=allel_sel)
    }
    output$allelle_ctl_tab <- DT::renderDataTable({allelle_f},options=list(pageLength=20,scrollX=TRUE,selection="multiple"))
  
  })
  # Track selected rows
  ctl_selectedRows <- reactiveVal()
  #Update selected rows on table selection
  observeEvent(input$ctl_allelle_btn, {
    ctl_selectedRows(input$allelle_ctl_tab_rows_selected)
  })
  output$ctl_textarea <- DT::renderDataTable({
    req(input$ctl_allelle_btn)
    rows <- ctl_selectedRows()
    data <- allelle_ctl_df
    opt <- input$allel_opt
    if(opt=="all"){
      allelle_f <- data
    }
    else{
      allel_sel <- data[grep(opt,data$Allelles),]
      allelle_f <-data.frame(Allelles=allel_sel)
    }
    
    data <- allelle_f[rows,]
    df <- data.frame(data)
    colnames(df) <- c('Selected Allelles')
    df
  })
  
  ### Select Allelle HTL
  observeEvent(input$allel_opt_h,{
    data <- allelle_htl_df
    opt <- input$allel_opt_h
    if(opt=="all"){
      allelle_f <- data
    }
    else{
      allel_sel <- data[grep(opt,data$Allelles),]
      allelle_f <-data.frame(Allelles=allel_sel)
    }
    output$allelle_htl_tab <- DT::renderDataTable({allelle_f},options=list(pageLength=20,scrollX=TRUE,selection="multiple"))
    
  })

  # Track selected rows
  htl_selectedRows <- reactiveVal()
  #Update selected rows on table selection
  observeEvent(input$htl_allelle_btn, {
    htl_selectedRows(input$allelle_htl_tab_rows_selected)
  })
  output$htl_textarea <- DT::renderDataTable({
    req(input$htl_allelle_btn)
    rows <- htl_selectedRows()
    data <- allelle_htl_df
    opt <- input$allel_opt_h
    if(opt=="all"){
      allelle_f <- data
    }
    else{
      allel_sel <- data[grep(opt,data$Allelles),]
      allelle_f <-data.frame(Allelles=allel_sel)
    }
    data <- allelle_f[rows,]
    df <- data.frame(data)
    colnames(df) <- c('Selected Allelles')
    df
  })
  
  
  
  ############################################B cell epitope candidates##########################################
  output$b_files <- DT::renderDT({
    DT::datatable(input$B_rank_epitopes, selection = c("single"))
  })
  data_b<-reactiveVal(NULL)
  b_all_files <- reactive({
    if (is.null(input$B_rank_epitopes) && is.null(input$b_example_csv))
      return(NULL)
    if (is.null(input$B_rank_epitopes)){
      req(input$b_example_csv)
      data <- read.csv("B_cell.csv")
      data_b(data)
      return(data)
    }
    files <- input$B_rank_epitopes
    data <- lapply(files$datapath, read.csv)
    concatenatedData <- bind_rows(data)
    data_b(concatenatedData)
    return(concatenatedData)
  })
  
  
  output$B_epitopes_table <- DT::renderDataTable({
    #req(input$B_rank_epitopes)
    b_all_files()},options=list(pageLength=10,scrollX=TRUE))
  
  output$B_epitope_rank_ui <- renderUI({
    #req(input$B_rank_epitopes)
    data <- b_all_files()
    cols <- colnames(data)
    ranking <- c()
    if ("Bepipred.mean.score" %in% cols){
      ranking<-c('Bepipred mean score')
    }
    if ('Prediction_Score' %in% cols){
      ranking <- c(ranking,'Vaxijen score')
    }
    if ('ML.Score_tox' %in% cols){
      ranking <- c(ranking,'Toxinpred score')
    }
    if ('ML.Score_alg' %in% cols){
      ranking<-c(ranking,'Algpred score')
    }
    selectInput("B_epitope_rank_selected","Selected Ranking method for B cell epitopes:",
                choices=ranking)
  })
  ###Rank final epitope table and render It
  
  final_epitope_b<-reactiveVal(NULL)
  observeEvent(input$B_epitope_rank,{
    #req(input$B_rank_epitopes)
    data <- b_all_files()
    if (input$B_epitope_rank_selected == "Bepipred mean score"){
      ranked_indices <- order(data$Bepipred.mean.score,decreasing = TRUE)
      ranked_data <- data[ranked_indices, ]
    }
    if (input$B_epitope_rank_selected == 'Vaxijen score'){
      ranked_indices <- order(data$Prediction_Score,decreasing = TRUE)
      ranked_data <- data[ranked_indices, ]
    }
    if (input$B_epitope_rank_selected == 'Toxinpred score'){
      ranked_indices <- order(data$ML.Score_tox)
      ranked_data <- data[ranked_indices, ]
    }
    if (input$B_epitope_rank_selected == 'Algpred score'){
      ranked_indices <- order(data$ML.Score_alg)
      ranked_data <- data[ranked_indices, ]
    }
    
    data_b(ranked_data)
    output$B_epitopes_table <- DT::renderDataTable({
      ranked_data},options=list(pageLength=10,scrollX=TRUE))
  })
  
  ##Select B epitopes
  observeEvent(input$B_epitope_select,{
    #req(input$B_rank_epitopes)
    dat <-data_b()
    rows <- input$B_epitopes_table_rows_selected
    print(length(rows))
    if (length(rows)>0){
      dat_rows<-dat[rows,]}
    else{
      dat_rows<-dat
    }
    final_epitope_b(dat_rows)
    output$B_epitopes_table_sel <- DT::renderDataTable({
      data.frame("Sequences selected"=dat_rows$Sequences)},options=list(pageLength=10,scrollX=TRUE))
  })
  
  
  ####################################### CTL epitope candidates#############################################################
  output$c_files <- DT::renderDT({
    DT::datatable(input$C_rank_epitopes, selection = c("single"))
  })
  data_c<-reactiveVal(NULL)
  c_all_files <- reactive({
    if (is.null(input$C_rank_epitopes) && is.null(input$c_example_csv))
      return(NULL)
    if (is.null(input$C_rank_epitopes)){
      req(input$c_example_csv)
      data <- read.csv("CTL.csv")
      data_c(data)
      return(data)
    }
    
    
    files <- input$C_rank_epitopes
    data <- lapply(files$datapath, read.csv)
    concatenatedData <- bind_rows(data)
    data_c(concatenatedData)
    return(concatenatedData)
  })
  
  output$C_epitopes_table <- DT::renderDataTable({
    #req(input$C_rank_epitopes)
    c_all_files()},options=list(pageLength=10,scrollX=TRUE))
  
  output$C_epitope_rank_ui <- renderUI({
    #req(input$C_rank_epitopes)
    data <- c_all_files()
    cols <- colnames(data)
    print(cols)
    ranking <- c()
    if ("Strong_Binders_Counts" %in% cols){
      ranking<-c('Strong Binders Counts')
    }
    if ("Weak_Binders_Counts" %in% cols){
      ranking<-c(ranking,'Weak Binders Counts')
    }
    if ("Prediction_Score" %in% cols){
      ranking <- c(ranking,'Vaxijen score')
    }
    if ("ML.Score_tox" %in% cols){
      ranking <- c(ranking,'Toxinpred score')
    }
    if ("ML.Score_alg" %in% cols){
      ranking<-c(ranking,'Algpred score')
    }
    if ("Immunogenicity_score" %in% cols){
      ranking<-c(ranking,'Immunogenicity score')
    }
    selectInput("C_epitope_rank_selected","Selected Ranking method for CTL epitopes:",
                choices=ranking)
  })
  ###Rank final epitope table and render It
  final_epitope_c<-reactiveVal(NULL)
  observeEvent(input$C_epitope_rank,{
    #req(input$C_rank_epitopes)
    data <- c_all_files()
    if (input$C_epitope_rank_selected == "Strong Binders Counts"){
      ranked_indices <- order(data$Strong_Binders_Counts,decreasing = TRUE)
      ranked_data <- data[ranked_indices, ]
    }
    if (input$C_epitope_rank_selected == "Weak Binders Counts"){
      ranked_indices <- order(data$Weak_Binders_Counts,decreasing = TRUE)
      ranked_data <- data[ranked_indices, ]
    }
    if (input$C_epitope_rank_selected == 'Vaxijen score'){
      ranked_indices <- order(data$Prediction_Score,decreasing = TRUE)
      ranked_data <- data[ranked_indices, ]
    }
    if (input$C_epitope_rank_selected == 'Toxinpred score'){
      ranked_indices <- order(data$ML.Score_tox)
      ranked_data <- data[ranked_indices, ]
    }
    if (input$C_epitope_rank_selected == 'Algpred score'){
      ranked_indices <- order(data$ML.Score_alg)
      ranked_data <- data[ranked_indices, ]
    }
    if (input$C_epitope_rank_selected == "Immunogenicity score"){
      ranked_indices <- order(data$Immunogenicity_score,decreasing = TRUE)
      ranked_data <- data[ranked_indices, ]
    }
    data_c(ranked_data)
    
    output$C_epitopes_table <- DT::renderDataTable({
      ranked_data},options=list(pageLength=10,scrollX=TRUE))
  })
  
  ##Select CTL epitopes
  observeEvent(input$C_epitope_select,{
    #req(input$C_rank_epitopes)
    dat <-data_c()
    rows <- input$C_epitopes_table_rows_selected
    if (length(rows)>0){
      dat_rows<-dat[rows,]
    }
    else{
      dat_rows <- dat
    }
    final_epitope_c(dat_rows)
    output$C_epitopes_table_sel <- DT::renderDataTable({
      data.frame("Sequences selected"=dat_rows$Sequences)},options=list(pageLength=10,scrollX=TRUE))
  })
  
  ############################################# HTL epitope candidates####################################################
  output$h_files <- DT::renderDT({
    DT::datatable(input$H_rank_epitopes, selection = c("single"))
  })
  data_h<-reactiveVal(NULL)
  h_all_files <- reactive({
    if (is.null(input$H_rank_epitopes) && is.null(input$h_example_csv))
      return(NULL)
    if (is.null(input$H_rank_epitopes)){
      req(input$h_example_csv)
      data <- read.csv("HTL.csv")
      data_h(data)
      return(data)
    }
    
    
    files <- input$H_rank_epitopes
    data <- lapply(files$datapath, read.csv)
    concatenatedData <- bind_rows(data)
    data_h(concatenatedData)
    return(concatenatedData)
  })
  
  output$H_epitopes_table <- DT::renderDataTable({
    #req(input$H_rank_epitopes)
    h_all_files()},options=list(pageLength=10,scrollX=TRUE))
  
  
  
  output$H_epitope_rank_ui <- renderUI({
    #req(input$H_rank_epitopes)
    data <- h_all_files()
    cols <- colnames(data)
    ranking <- c()
    if ("Strong_Binders_Counts" %in% cols){
      ranking<-c('Strong Binders Counts')
    }
    if ("Weak_Binders_Counts" %in% cols){
      ranking<-c(ranking,'Weak Binders Counts')
    }
    if ("Prediction_Score" %in% cols){
      ranking <- c(ranking,'Vaxijen score')
    }
    if ("ML.Score_tox" %in% cols){
      ranking <- c(ranking,'Toxinpred score')
    }
    if ("ML.Score_alg" %in% cols){
      ranking<-c(ranking,'Algpred score')
    }
    #if ("IFN prediction score" %in% cols){
    # ranking<-c(ranking,'IFN prediction score')
    #}
    selectInput("H_epitope_rank_selected","Selected Ranking method for HTL epitopes:",
                choices=ranking)
  })
  ###Rank final epitope table and render It
  final_epitope_h<-reactiveVal(NULL)
  observeEvent(input$H_epitope_rank,{
    #req(input$H_rank_epitopes)
    data <- h_all_files()
    if (input$H_epitope_rank_selected == "Strong Binders Counts"){
      ranked_indices <- order(data$Strong_Binders_Counts,decreasing = TRUE)
      ranked_data <- data[ranked_indices, ]
    }
    if (input$H_epitope_rank_selected == "Weak Binders Counts"){
      ranked_indices <- order(data$Weak_Binders_Counts,decreasing = TRUE)
      ranked_data <- data[ranked_indices, ]
    }
    if (input$H_epitope_rank_selected == 'Vaxijen score'){
      ranked_indices <- order(data$Prediction_Score,decreasing = TRUE)
      ranked_data <- data[ranked_indices, ]
    }
    if (input$H_epitope_rank_selected == 'Toxinpred score'){
      ranked_indices <- order(data$ML.Score_tox)
      ranked_data <- data[ranked_indices, ]
    }
    if (input$H_epitope_rank_selected == 'Algpred score'){
      ranked_indices <- order(data$ML.Score_alg)
      ranked_data <- data[ranked_indices, ]
    }
    if (input$H_epitope_rank_selected == "IFN prediction score"){
      ranked_indices <- order(data$IFN_score,decreasing = TRUE)
      ranked_data <- data[ranked_indices, ]
    }
    data_h(ranked_data)
    
    output$H_epitopes_table <- DT::renderDataTable({
      ranked_data},options=list(pageLength=10,scrollX=TRUE))
  })
  
  ##Select HTL  epitopes
  observeEvent(input$H_epitope_select,{
    #req(input$H_rank_epitopes)
    dat <-data_h()
    rows <- input$H_epitopes_table_rows_selected
    if(length(rows)>0){
      dat_rows<-dat[rows,]}
    else{
      dat_rows <-dat
    }
    final_epitope_h(dat_rows)
    output$H_epitopes_table_sel <- DT::renderDataTable({
      data.frame("Sequences selected"=dat_rows$Sequences)},options=list(pageLength=10,scrollX=TRUE))
  })
  ##################################################################Execute Button B Cell epitopes#############################################
  shinyjs::hide("B_ranked")
  
  final_data_b <- reactiveVal(NULL)
  filtered_data_b <-reactiveVal(NULL)
  
  shinyjs::hide("b_res_col1")
  observeEvent(input$b_exe,
               {
                 shinyjs::hide("b_res_col1")
                 withProgress(message="B cell epitope computing",detail="Bepipred",value=0,{
                 setwd(main_path)
                 #Python script import
                 bepi_fnc <- import("Bepipred_script")
                 #Parameters
                 ##Bepipred_opt
                 bepipred_opt <- as.numeric(input$bepipred_opt)
                 ###Threshold
                 bepi_thres <- input$bepi_thres_num
                 ###Base neg
                 base_neg <- input$base_neg
                 ###Amino region
                 amino_region <- input$amino_region
                 ###Secondary threshold
                 bepi_sec_thres <- input$second_thres_num
                 ###Dataframe
                 df <- sel_BData()
                 ##Function
                 res_bepi <- bepi_fnc$Bepipred3(df,path_bepipred,path_res_bepi,bepi_thres,base_neg,amino_region,bepi_sec_thres,bepipred_opt)
                 
                 #res_bepi$Bepipred.mean.scores <- as.numeric(res_bepi$Bepipred.mean.scores)
                 final_df <- res_bepi
                 if (nrow(final_df) == 0) {
                   output$bepi_text <- renderText({
                     "No epitopes detected"
                   })
                   
                   filtered_data_b(NULL)
                   final_data_b(NULL)
                   
                   output$bepi_df <- DT::renderDataTable({
                     NULL})
                   output$vaxi_df <- DT::renderDataTable({
                     NULL
                   })
                   output$toxi_df <- DT::renderDataTable({
                     NULL
                   })
                   output$alg_df <- DT::renderDataTable({
                     NULL
                   })
                   output$filter_df_b <- DT::renderDataTable({
                     NULL
                   })
                   output$vax_text <-renderText({
                     ''
                   })
                   output$tox_text <- renderText({
                     ''
                   })
                   output$alg_text <- renderText({
                     ''
                   })
                   output$filtered_df_b <- DT::renderDataTable({
                     NULL
                   })
                   return()  # Exit the observeEvent block
                 }
                 #Outputs tables
                 output$bepi_df <- DT::renderDataTable({
                   res_bepi
                 },extensions = 'Buttons', # Load the 'Buttons' extension
                 options = list(
                   dom = 'Bfrtip', # Specify where to place the buttons
                   scrollX = TRUE, # Enable horizontal scrolling
                   pageLength = 10, # Set the number of rows per page
                   buttons = list(
                     'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                   )))
                 
                 setwd(main_path)
                 
                 output$bepi_csv <- DT::renderDataTable({
                   data <- read.csv("./temp_data/Bepipred/raw_output.csv")
                   
                 },extensions = 'Buttons', # Load the 'Buttons' extension
                 options = list(
                   dom = 'Bfrtip', # Specify where to place the buttons
                   scrollX = TRUE, # Enable horizontal scrolling
                   pageLength = 30, # Set the number of rows per page
                   buttons = list(
                     'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                   )))
                 
                 
                 
                 output$bepi_text <- renderText({
                   text <- 'Bepipred Finished'
                   text
                   
                 })
                 
                 
                 
                 ############Vaxijen
                 if(1 %in% input$b_softwares){
                   incProgress(0.25,detail="Vaxijen")
                   #Python script
                   vaxi_fnc <-import("Vaxijen_script")
                   #Parameters
                   vax_thres <- input$vax_thres_num
                   vax_target <- input$vax_target
                   df_v <- res_bepi
                   #Functions
                   res_vax_b <- vaxi_fnc$Vaxijen(df_v,vax_target,vax_thres)
                   for (i in colnames(res_vax_b)){
                     if(i!='Sequences'){final_df[i]=res_vax_b[i]}
                   }
                   #Outputs
                   output$vaxi_df <- DT::renderDataTable({
                     res_vax_b
                   },extensions = 'Buttons', # Load the 'Buttons' extension
                   options = list(
                     dom = 'Bfrtip', # Specify where to place the buttons
                     scrollX = TRUE, # Enable horizontal scrolling
                     pageLength = 10, # Set the number of rows per page
                     buttons = list(
                       'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                     )))
                   output$vax_text <- renderText({
                     text <- 'Vaxijen Finished'
                     text
                     
                   })
                 }else{output$vaxi_df <- NULL
                 output$vax_text <- NULL}
                 
                 ###########Toxinpred
                 if (2 %in% input$b_softwares){
                   incProgress(0.4,detail="Vaxijen")
                   #Parameters
                   model_tox <- input$model_tox
                   thresh_tox <- input$tox_thres_num
                   df_tox <- res_bepi
                   #Functions
                   res_tox_b <-Toxinpred_fnc(df_tox,path_toxinpred,path_res_toxinpred,model_tox,thresh_tox)
                   for (i in colnames(res_tox_b)){
                     if(i!='Sequences'){j<-paste(i,'_tox',sep="")
                     final_df[j]<-res_tox_b[i]}
                   }
                   #Outputs
                   output$toxi_df <- DT::renderDataTable({
                     res_tox_b
                   },extensions = 'Buttons', # Load the 'Buttons' extension
                   options = list(
                     dom = 'Bfrtip', # Specify where to place the buttons
                     scrollX = TRUE, # Enable horizontal scrolling
                     pageLength = 10, # Set the number of rows per page
                     buttons = list(
                       'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                     )))
                   output$tox_text <- renderText({
                     text <- 'Toxinpred finished'
                     text
                   })
                 }
                 else{output$toxi_df <- NULL
                 output$tox_text <- NULL}
                 
                 ###########Algpred
                 if (3 %in% input$b_softwares){
                   incProgress(0.75,detail="Algpred")
                   #Parameters
                   model_alg <- input$model_alg
                   thresh_alg <- input$alg_thres_num
                   df_alg <- res_bepi
                   #Functions
                   res_alg_b <- Algpred_fnc(df_alg,path_algpred,path_res_algpred,model_alg,thresh_alg)
                   for (i in colnames(res_alg_b)){
                     if (i!='Sequences'){j<-paste(i,'_alg',sep="")
                     final_df[j]<-res_alg_b[i]}
                   }
                   #Outputs
                   output$alg_df <- DT::renderDataTable({
                     res_alg_b
                   },extensions = 'Buttons', # Load the 'Buttons' extension
                   options = list(
                     dom = 'Bfrtip', # Specify where to place the buttons
                     scrollX = TRUE, # Enable horizontal scrolling
                     pageLength = 10, # Set the number of rows per page
                     buttons = list(
                       'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                     )))
                   output$alg_text <- renderText({
                     text <- 'Algpred finished'
                     text
                   })
                 }else{output$alg_df <- NULL
                 output$alg_text <- NULL}
                 ##Filter/Rank table B-cell
                 output$filter_df_b <- DT::renderDataTable({
                   final_df
                 },extensions = 'Buttons', # Load the 'Buttons' extension
                 options = list(
                   dom = 'Bfrtip', # Specify where to place the buttons
                   scrollX = TRUE, # Enable horizontal scrolling
                   pageLength = 10, # Set the number of rows per page
                   buttons = list(
                     'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                   )))
                 
                 
                 
                 shinyjs::show("b_res_col1")
                 output$res_df_b <- DT::renderDataTable({
                   final_df
                 },extensions = 'Buttons', # Load the 'Buttons' extension
                 options = list(
                   dom = 'Bfrtip', # Specify where to place the buttons
                   scrollX = TRUE, # Enable horizontal scrolling
                   pageLength = 10, # Set the number of rows per page
                   buttons = list(
                     'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                   )))
                 
                 output$sel_b_filter_ui <- renderUI({
                   filters <- c()
                   if (1 %in% input$b_softwares){
                     filters <- c(filters,'Vaxijen')
                   }
                   if(2 %in% input$b_softwares){
                     filters <- c(filters,'Toxinpred')
                   }
                   if (3 %in% input$b_softwares){
                     filters <- c(filters,'Algpred')
                   }
                   selectInput("b_filters_selected", "Select Filters :",
                               choices = filters,multiple=TRUE)
                 })
                 output$sel_b_ranking_ui <- renderUI({
                   ranking <- c()
                   if (1 %in% input$b_softwares){
                     ranking <- c(ranking,'Vaxijen score')
                   }
                   if (2 %in% input$b_softwares){
                     ranking <- c(ranking,'Toxinpred score')
                   }
                   if (3 %in% input$b_softwares){
                     ranking <- c(ranking,'Algpred score')
                   }
                   selectInput("b_rank_selected","Selected Ranking method:",
                               choices=ranking)
                 })
                 
                 final_data_b(final_df)
                 incProgress(1) 
                 })
               })
  #Filter button
  filtered_data_b <-reactiveVal(NULL)
  observeEvent(input$b_fin_fil,{
    data <- final_data_b()
    filters <- input$b_filters_selected
    if ('Vaxijen' %in% filters){
      data <- subset(data,Antigen == TRUE)
    }
    if ('Toxinpred' %in% filters){
      data <- subset(data,Prediction_tox != 'Toxin')
    }
    if ('Algpred' %in% filters){
      data <- subset(data,Prediction_alg == 'Non-Allergen')
    }
    filtered_data_b(data)
    output$filtered_df_b <- DT::renderDataTable({
      filtered_data_b()
    },extensions = 'Buttons', # Load the 'Buttons' extension
    options = list(
      dom = 'Bfrtip',
      scrollX = TRUE, # Enable horizontal scrolling
      pageLength = 10,# Set the number of rows per page
      buttons = list('copy', 'csv', 'excel', 'pdf', 'print') # Add the buttons you want
    ))
    
  })
  #Rank button b cell epitopes
  observeEvent(input$b_fin_rank,{
    data1 <- filtered_data_b()
    data2 <- final_data_b()
    if (input$b_rank_selected == 'Bepipred mean score'){
      
      ranked_indices <- order(data1['Bepipred mean score'])
      ranked_data1 <- data1[ranked_indices, ]
      
      ranked_indices <- order(data2['Bepipred mean score'])
      ranked_data2 <- data2[ranked_indices, ]
      
    }
    if (input$b_rank_selected == 'Vaxijen score'){
      ranked_indices <- order(data1$Prediction_Score,decreasing = TRUE)
      ranked_data1 <- data1[ranked_indices, ]
      
      ranked_indices <- order(data2$Prediction_Score)
      ranked_data2 <- data2[ranked_indices, ]
      
    }
    if (input$b_rank_selected == 'Toxinpred score'){
      ranked_indices <- order(data1$ML.Score_tox)
      ranked_data1 <- data1[ranked_indices, ]
      
      ranked_indices <- order(data2$ML.Score_tox)
      ranked_data2 <- data2[ranked_indices, ]
    }
    if (input$b_rank_selected == 'Algpred score'){
      ranked_indices <- order(data1$ML.Score_alg)
      ranked_data1 <- data1[ranked_indices, ]
      
      ranked_indices <- order(data2$ML.Score_alg)
      ranked_data2 <- data2[ranked_indices, ]
    }
    filtered_data_b(ranked_data1)
    final_data_b(ranked_data2)
    
    output$filtered_df_b <- DT::renderDataTable({
      filtered_data_b()
    },extensions = 'Buttons', # Load the 'Buttons' extension
    options = list(
      dom = 'Bfrtip',
      scrollX = TRUE, # Enable horizontal scrolling
      pageLength = 10,# Set the number of rows per page...buttons also
      buttons = list(
        'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
      )))
    
  })
  
  
  
  
  
  ############################################################Execute Button CTL epitopes#################################################
  
  final_data_c <- reactiveVal(NULL)
  shinyjs::hide("c_res_col1")
  observeEvent(input$c_exe,
               {
                 ##############netMHCpan############################
                 shinyjs::hide("c_res_col1")
                 withProgress(message="CTL epitope computing",detail="netMHCpan",value=0,{
                 #Parameters
                 #Number of peptides
                 num_peptides_ctl <- input$c_length
                 #Alleles
                 ctl_rows <- ctl_selectedRows()
                 alleles_ctl <- allelle_ctl_df[ctl_rows,]
                 
                 #Dataframe
                 df <- sel_CData()
                 seq <- df$Sequence
                 ##Function
                 res_df <- netmhc_fnc(seq,num_peptides_ctl,alleles_ctl,allele_list_ctl)
                 
                 final_df_c <- res_df
                 print(nrow(final_df_c))
                 if (nrow(final_df_c) == 0) {
                   output$net_I_text <- renderText({
                     "No epitopes detected"
                   })
                   return()  # Exit the observeEvent block
                 }
                 
                 #Outputs tables
                 output$net_I_df <- DT::renderDataTable({
                   res_df
                 },extensions = 'Buttons', # Load the 'Buttons' extension
                 options = list(
                   dom = 'Bfrtip', # Specify where to place the buttons
                   scrollX = TRUE, # Enable horizontal scrolling
                   pageLength = 10, # Set the number of rows per page
                   buttons = list(
                     'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                   )))
                 
                 setwd(main_path)
                 
                 output$net_I_text <- renderText({
                   text <- 'netMHCpan Finished'
                   text
                 })
                 
                 ############Vaxijen
                 if(1 %in% input$c_softwares){
                   incProgress(0.25,detail="Vaxijen")
                   #Python script
                   vaxi_fnc <-import("Vaxijen_script")
                   #Parameters
                   vax_thres <- input$c_vax_thres_num
                   vax_target <- input$c_vax_target
                   df_v <- res_df
                   #Functions
                   res_vax_c <- vaxi_fnc$Vaxijen(df_v,vax_target,vax_thres)
                   for (i in colnames(res_vax_c)){
                     if(i!='Sequences'){final_df_c[i]=res_vax_c[i]}
                   }
                   #Outputs
                   output$c_vax_df <- DT::renderDataTable({
                     res_vax_c
                   },extensions = 'Buttons', # Load the 'Buttons' extension
                   options = list(
                     dom = 'Bfrtip', # Specify where to place the buttons
                     scrollX = TRUE, # Enable horizontal scrolling
                     pageLength = 10, # Set the number of rows per page
                     buttons = list(
                       'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                     )))
                   output$c_vax_text <- renderText({
                     text <- 'Vaxijen Finished'
                     text
                     
                   })
                 }else{output$c_vax_df <- NULL
                 output$c_vax_text <- NULL}
                 
                 ###########Toxinpred
                 if (2 %in% input$c_softwares){
                   incProgress(0.5,detail="Toxinpred")
                   #Parameters
                   model_tox <- input$c_model_tox
                   thresh_tox <- input$c_tox_thres_num
                   df_tox <- res_df
                   #Functions
                   res_tox_c <-Toxinpred_fnc(df_tox,path_toxinpred,path_res_toxinpred,model_tox,thresh_tox)
                   for (i in colnames(res_tox_c)){
                     if(i!='Sequences'){j<-paste(i,'_tox',sep="")
                     final_df_c[j]<-res_tox_c[i]}
                   }
                   #Outputs
                   output$c_toxi_df <- DT::renderDataTable({
                     res_tox_c
                   },extensions = 'Buttons', # Load the 'Buttons' extension
                   options = list(
                     dom = 'Bfrtip', # Specify where to place the buttons
                     scrollX = TRUE, # Enable horizontal scrolling
                     pageLength = 10, # Set the number of rows per page
                     buttons = list(
                       'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                     )))
                   output$c_tox_text <- renderText({
                     text <- 'Toxinpred finished'
                     text
                   })
                 }
                 else{output$c_toxi_df <- NULL
                 output$c_tox_text <- NULL}
                
                 #########################Algpred
                 if (3 %in% input$c_softwares){
                   incProgress(0.7,detail="Algpred")
                   #Parameters
                   model_alg <- input$c_model_alg
                   thresh_alg <- input$c_alg_thres_num
                   df_alg <- res_df
                   #Functions
                   res_alg_c <- Algpred_fnc(df_alg,path_algpred,path_res_algpred,model_alg,thresh_alg)
                   for (i in colnames(res_alg_c)){
                     if (i!='Sequences'){j<-paste(i,'_alg',sep="")
                     final_df_c[j]<-res_alg_c[i]}
                   }
                   #Outputs
                   output$c_alg_df <- DT::renderDataTable({
                     res_alg_c
                   },extensions = 'Buttons', # Load the 'Buttons' extension
                   options = list(
                     dom = 'Bfrtip', # Specify where to place the buttons
                     scrollX = TRUE, # Enable horizontal scrolling
                     pageLength = 10, # Set the number of rows per page
                     buttons = list(
                       'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                     )))
                   output$c_alg_text <- renderText({
                     text <- 'Algpred finished'
                     text
                   })
                 }else{output$c_alg_df <- NULL
                 output$c_alg_text <- NULL}
                 
                 #################### MHCI Immunogenicity
                 # if (4 %in% input$c_softwares){
                 #   incProgress(0.9,detail="MHCI Immunogenicity")
                 #   df_imm <- res_df
                 #   
                 #   #Function
                 #   res_imm <- immunogen_fnc(df_imm,path_immunogen,path_res_immunogen)
                 #   for (i in colnames(res_imm)){
                 #     if (i!='Sequences'){final_df_c[i]<-res_imm[i]}
                 #   }
                 #   #Render Table
                 #   #Outputs
                 #   output$immuno_df <- DT::renderDataTable({
                 #     res_imm
                 #   },extensions = 'Buttons', # Load the 'Buttons' extension
                 #   options = list(
                 #     dom = 'Bfrtip', # Specify where to place the buttons
                 #     scrollX = TRUE, # Enable horizontal scrolling
                 #     pageLength = 10, # Set the number of rows per page
                 #     buttons = list(
                 #       'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                 #     )))
                 #   
                 #   output$immuno_text <- renderText({
                 #     text <- 'MHCI Immunogenicity finished'
                 #     text
                 #   })
                 # }else{output$immuno_df <- NULL
                 # output$immuno_text <- NULL} 
                 incProgress(1)
                 ######Filter/Rank table CTL
                 output$filter_df_c <- DT::renderDataTable({
                   final_df_c
                 },extensions = 'Buttons', # Load the 'Buttons' extension
                 options = list(
                   dom = 'Bfrtip', # Specify where to place the buttons
                   scrollX = TRUE, # Enable horizontal scrolling
                   pageLength = 10, # Set the number of rows per page
                   buttons = list(
                     'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                   )))
                 
                 
                 output$res_df_c <- DT::renderDataTable({
                   final_df_c
                 },extensions = 'Buttons', # Load the 'Buttons' extension
                 options = list(
                   dom = 'Bfrtip', # Specify where to place the buttons
                   scrollX = TRUE, # Enable horizontal scrolling
                   pageLength = 10, # Set the number of rows per page
                   buttons = list(
                     'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                   )))
                 shinyjs::show("c_res_col1")
                 
                 output$sel_c_filter_ui <- renderUI({
                   filters <- c()
                   if (1 %in% input$c_softwares){
                     filters <- c(filters,'Vaxijen')
                   }
                   if(2 %in% input$c_softwares){
                     filters <- c(filters,'Toxinpred')
                   }
                   if (3 %in% input$c_softwares){
                     filters <- c(filters,'Algpred')
                   }
                   if (4 %in% input$c_softwares){
                     filters <- c(filters,'MHCI Immunogenicity')
                   }
                   selectInput("c_filters_selected", "Select Filters :",
                               choices = filters,multiple=TRUE)
                 })
                 output$sel_c_ranking_ui <- renderUI({
                   ranking <- c('Weak Binders Counts','Strong Binders Counts')
                   if (1 %in% input$c_softwares){
                     ranking <- c(ranking,'Vaxijen score')
                   }
                   if (2 %in% input$c_softwares){
                     ranking <- c(ranking,'Toxinpred score')
                   }
                   if (3 %in% input$c_softwares){
                     ranking <- c(ranking,'Algpred score')
                   }
                   if (4 %in% input$c_softwares){
                     ranking <- c(ranking,'Immunogenicity score')
                   }
                   selectInput("c_rank_selected","Selected Ranking method:",
                               choices=ranking)
                 })
                 final_data_c(final_df_c)
                 
                 })
               })
  
  #############Filter button CTL epitopes
  filtered_data_c <-reactiveVal(NULL)
  observeEvent(input$c_fin_fil,{
    data <- final_data_c()
    filters <- input$c_filters_selected
    if ('Vaxijen' %in% filters){
      data <- subset(data,Antigen == TRUE)
    }
    if ('Toxinpred' %in% filters){
      data <- subset(data,Prediction_tox != 'Toxin')
    }
    if ('Algpred' %in% filters){
      data <- subset(data,Prediction_alg == 'Non-Allergen')
    }
    if ('MHCI Immunogenicity' %in% filters){
      data <- subset(data,Immunogenicity == TRUE)
    }
    filtered_data_c(data)
    output$filtered_df_c <- DT::renderDataTable({
      filtered_data_c()
    },extensions = 'Buttons', # Load the 'Buttons' extension
    options = list(
      dom = 'Bfrtip',
      scrollX = TRUE, # Enable horizontal scrolling
      pageLength = 10,# Set the number of rows per page
      buttons = list('copy', 'csv', 'excel', 'pdf', 'print') # Add the buttons you want
    ))
    
  })
  #Rank button ctl epitopes
  observeEvent(input$c_fin_rank,{
    data1 <- filtered_data_c()
    print(input$c_rank_selected)
    if (input$c_rank_selected=="Strong Binders Counts"){
      ranked_indices <- order(data1$Strong_Binders_Counts,decreasing=TRUE)
      ranked_data1 <- data1[ranked_indices, ]
      filtered_data_c(ranked_data1)
    }
    if (input$c_rank_selected=="Weak Binders Counts"){
      ranked_indices <- order(data1$Weak_Binders_Counts,decreasing=TRUE)
      ranked_data1 <- data1[ranked_indices, ]
      filtered_data_c(ranked_data1)
    }
    if (input$c_rank_selected=="Vaxijen score"){
      ranked_indices <- order(data1$Prediction_Score,decreasing = TRUE)
      ranked_data1 <- data1[ranked_indices, ]
      filtered_data_c(ranked_data1)
    }
    if (input$c_rank_selected=="Toxinpred score"){
      ranked_indices <- order(data1$ML.Score_tox)
      ranked_data1 <- data1[ranked_indices, ]
      filtered_data_c(ranked_data1)
    }
    if (input$c_rank_selected=="Algpred score"){
      ranked_indices <- order(data1$ML.Score_alg)
      ranked_data1 <- data1[ranked_indices, ]
      filtered_data_c(ranked_data1)
    }
    if (input$c_rank_selected=="Immunogenicity score"){
      ranked_indices <- order(data1$Immunogenicity_score)
      ranked_data1<-data1[ranked_indices, ]
      filtered_data_c(ranked_data1)
    }
    output$filtered_df_c <- DT::renderDataTable({
      filtered_data_c()
    },extensions = 'Buttons', # Load the 'Buttons' extension
    options = list(
      dom = 'Bfrtip',
      scrollX = TRUE, # Enable horizontal scrolling
      pageLength = 10,# Set the number of rows per page...buttons also
      buttons = list(
        'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
      )))
    
    
  })
  
  
  
  
  
  
  
  ############################################################Execute Button HTL epitopes#################################################
  
  final_data_h <- reactiveVal(NULL)
  shinyjs::hide("h_res_col1")
  observeEvent(input$h_exe,
               {
                 ##############netMHCpan############################
                 shinyjs::hide("h_res_col1")
                 withProgress(message="HTL epitope computing",detail="netMHCIIpan",value=0,{
                 #Parameters
                 #Number of peptides
                 num_peptides <- input$h_length
                 #Alleles
                 htl_rows <- htl_selectedRows()
                 alleles_htl <- allelle_htl_df[htl_rows,]
                 
                 #Dataframe
                 df <- sel_HData()
                 seq <- df$Sequence
                 ##Function
                 res_df <- netmhcII_fnc(seq,num_peptides,alleles_htl,allele_list_htl)
                 final_df_h <- res_df
                 if (nrow(final_df_h) == 0) {
                   output$net_II_text <- renderText({
                     "No epitopes detected"
                   })
                   return()  # Exit the observeEvent block
                 }
                 #Outputs tables
                 output$net_II_df <- DT::renderDataTable({
                   res_df
                 },extensions = 'Buttons', # Load the 'Buttons' extension
                 options = list(
                   dom = 'Bfrtip', # Specify where to place the buttons
                   scrollX = TRUE, # Enable horizontal scrolling
                   pageLength = 10, # Set the number of rows per page
                   buttons = list(
                     'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                   )))
                 
                 setwd(main_path)
                 output$net_II_text <- renderText({
                   text <- 'netMHCIIpan Finished'
                   text
                 })
                
                 #################Vaxijen
                 if(1 %in% input$h_softwares){
                   incProgress(0.25,detail="Vaxijen")
                   #Python script
                   vaxi_fnc <-import("Vaxijen_script")
                   #Parameters
                   vax_thres <- input$h_vax_thres_num
                   vax_target <- input$h_vax_target
                   df_v <- res_df
                   #Functions
                   res_vax_h <- vaxi_fnc$Vaxijen(df_v,vax_target,vax_thres)
                   for (i in colnames(res_vax_h)){
                     if(i!='Sequences'){final_df_h[i]=res_vax_h[i]}
                   }
                   #Outputs
                   output$h_vax_df <- DT::renderDataTable({
                     res_vax_h
                   },extensions = 'Buttons', # Load the 'Buttons' extension
                   options = list(
                     dom = 'Bfrtip', # Specify where to place the buttons
                     scrollX = TRUE, # Enable horizontal scrolling
                     pageLength = 10, # Set the number of rows per page
                     buttons = list(
                       'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                     )))
                   output$h_vax_text <- renderText({
                     text <- 'Vaxijen Finished'
                     text
                     
                   })
                 }else{output$h_vax_df <- NULL
                 output$h_vax_text <- NULL}
                 
                 
                 ###########################Toxinpred
                 if (2 %in% input$h_softwares){
                   incProgress(0.5,detail="Toxinpred")
                   #Parameters
                   model_tox <- input$h_model_tox
                   thresh_tox <- input$h_tox_thres_num
                   df_tox <- res_df
                   #Functions
                   res_tox_h <-Toxinpred_fnc(df_tox,path_toxinpred,path_res_toxinpred,model_tox,thresh_tox)
                   for (i in colnames(res_tox_h)){
                     if(i!='Sequences'){j<-paste(i,'_tox',sep="")
                     final_df_h[j]<-res_tox_h[i]}
                   }
                   #Outputs
                   output$h_toxi_df <- DT::renderDataTable({
                     res_tox_h
                   },extensions = 'Buttons', # Load the 'Buttons' extension
                   options = list(
                     dom = 'Bfrtip', # Specify where to place the buttons
                     scrollX = TRUE, # Enable horizontal scrolling
                     pageLength = 10, # Set the number of rows per page
                     buttons = list(
                       'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                     )))
                   output$h_tox_text <- renderText({
                     text <- 'Toxinpred finished'
                     text
                   })
                 }
                 else{output$h_toxi_df <- NULL
                 output$h_tox_text <- NULL}
                 
                 #########################Algpred
                 if (3 %in% input$h_softwares){
                   incProgress(0.75,detail="Algpred")
                   #Parameters
                   model_alg <- input$h_model_alg
                   thresh_alg <- input$h_alg_thres_num
                   df_alg <- res_df
                   #Functions
                   res_alg_h <- Algpred_fnc(df_alg,path_algpred,path_res_algpred,model_alg,thresh_alg)
                   for (i in colnames(res_alg_h)){
                     if (i!='Sequences'){j<-paste(i,'_alg',sep="")
                     final_df_h[j]<-res_alg_h[i]}
                   }
                   #Outputs
                   output$h_alg_df <- DT::renderDataTable({
                     res_alg_h
                   },extensions = 'Buttons', # Load the 'Buttons' extension
                   options = list(
                     dom = 'Bfrtip', # Specify where to place the buttons
                     scrollX = TRUE, # Enable horizontal scrolling
                     pageLength = 10, # Set the number of rows per page
                     buttons = list(
                       'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                     )))
                   output$h_alg_text <- renderText({
                     text <- 'Algpred finished'
                     text
                   })
                 }else{output$h_alg_df <- NULL
                 output$h_alg_text <- NULL}
                 incProgress(1)
                 ########################IFN epitope predictions
                 #if(4 %in% input$h_softwares){
                 #
                 #   df_ifn <- res_df
                 #  ifn_fnc <- import("IFNepitope_script")
                 # res_ifn <- ifn_fnc$IFNepitope_fnc(df_ifn$Sequences)
                 #
                 #for (i in colnames(res_ifn)){
                 # if (i!='Sequences'){
                 #  final_df_h[i]<-res_ifn[i]
                 #}
                 #}
                 #Outputs
                 #     output$ifn_df <- DT::renderDataTable({
                 #       res_ifn
                 #     },extensions = 'Buttons', # Load the 'Buttons' extension
                 #     options = list(
                 #       dom = 'Bfrtip', # Specify where to place the buttons
                 #       scrollX = TRUE, # Enable horizontal scrolling
                 #       pageLength = 10, # Set the number of rows per page
                 #       buttons = list(
                 #         'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                 #       )))
                 #     output$ifn_text <- renderText({
                 #       text <- 'IFN epitope prediction finished'
                 #       text
                 #     })
                 # }else{output$ifn_df <- NULL
                 # output$ifn_text <- NULL}
                 
                 ######Filter/Rank table CTL
                 output$filter_df_h <- DT::renderDataTable({
                   final_df_h
                 },extensions = 'Buttons', # Load the 'Buttons' extension
                 options = list(
                   dom = 'Bfrtip', # Specify where to place the buttons
                   scrollX = TRUE, # Enable horizontal scrolling
                   pageLength = 10, # Set the number of rows per page
                   buttons = list(
                     'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                   )))
                 
                 
                 
                 output$res_df_h <- DT::renderDataTable({
                   final_df_h
                 },extensions = 'Buttons', # Load the 'Buttons' extension
                 options = list(
                   dom = 'Bfrtip', # Specify where to place the buttons
                   scrollX = TRUE, # Enable horizontal scrolling
                   pageLength = 10, # Set the number of rows per page
                   buttons = list(
                     'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                   )))
                 shinyjs::show("h_res_col1")
                 output$sel_h_filter_ui <- renderUI({
                   filters <- c()
                   if (1 %in% input$h_softwares){
                     filters <- c(filters,'Vaxijen')
                   }
                   if(2 %in% input$h_softwares){
                     filters <- c(filters,'Toxinpred')
                   }
                   if (3 %in% input$h_softwares){
                     filters <- c(filters,'Algpred')
                   }
                   #if (4 %in% input$h_softwares){
                   # filters <- c(filters,'IFN epitope')
                   #}
                   selectInput("h_filters_selected", "Select Filters :",
                               choices = filters,multiple=TRUE)
                 })
                 output$sel_h_ranking_ui <- renderUI({
                   ranking <- c('Weak Binders Counts','Strong Binders Counts')
                   if (1 %in% input$h_softwares){
                     ranking <- c(ranking,'Vaxijen score')
                   }
                   if (2 %in% input$h_softwares){
                     ranking <- c(ranking,'Toxinpred score')
                   }
                   if (3 %in% input$h_softwares){
                     ranking <- c(ranking,'Algpred score')
                   }
                   #if (4 %in% input$h_softwares){
                   # ranking <- c(ranking,'IFN prediction score')
                   #}
                   selectInput("h_rank_selected","Selected Ranking method:",
                               choices=ranking)
                 })
                 final_data_h(final_df_h)
                 
                 })
               })
  
  #############Filter button HTL epitopes
  filtered_data_h <-reactiveVal(NULL)
  observeEvent(input$h_fin_fil,{
    data <- final_data_h()
    filters <- input$h_filters_selected
    if ('Vaxijen' %in% filters){
      data <- subset(data,Antigen == TRUE)
    }
    if ('Toxinpred' %in% filters){
      data <- subset(data,Prediction_tox != 'Toxin')
    }
    if ('Algpred' %in% filters){
      data <- subset(data,Prediction_alg == 'Non-Allergen')
    }
    if ('IFN epitope' %in% filters){
      data <- subset(data,`IFN predictions` == 'POSITIVE')
    }
    filtered_data_h(data)
    output$filtered_df_h<- DT::renderDataTable({
      filtered_data_h()
    },extensions = 'Buttons', # Load the 'Buttons' extension
    options = list(
      dom = 'Bfrtip',
      scrollX = TRUE, # Enable horizontal scrolling
      pageLength = 10,# Set the number of rows per page
      buttons = list('copy', 'csv', 'excel', 'pdf', 'print') # Add the buttons you want
    ))
    
  })
  #Rank button ctl epitopes
  observeEvent(input$h_fin_rank,{
    data1 <- filtered_data_h()
    print(colnames(data1))
    if (input$h_rank_selected=="Strong Binders Counts"){
      ranked_indices <- order(data1$Strong_Binders_Counts,decreasing=TRUE)
      ranked_data1 <- data1[ranked_indices, ]
      filtered_data_h(ranked_data1)
    }
    if (input$h_rank_selected=="Weak Binders Counts"){
      ranked_indices <- order(data1$Weak_Binders_Counts,decreasing=TRUE)
      ranked_data1 <- data1[ranked_indices, ]
      filtered_data_h(ranked_data1)
    }
    if (input$h_rank_selected=="Vaxijen score"){
      ranked_indices <- order(data1$Prediction_Score,decreasing = TRUE)
      ranked_data1 <- data1[ranked_indices, ]
      filtered_data_h(ranked_data1)
    }
    if (input$h_rank_selected=="Toxinpred score"){
      ranked_indices <- order(data1$ML.Score_tox)
      ranked_data1 <- data1[ranked_indices, ]
      filtered_data_h(ranked_data1)
    }
    if (input$h_rank_selected=="Algpred score"){
      ranked_indices <- order(data1$ML.Score_alg)
      ranked_data1 <- data1[ranked_indices, ]
      filtered_data_h(ranked_data1)
    }
    if (input$h_rank_selected=="IFN prediction score"){
      ranked_indices <- order(data1["IFN scores"])
      ranked_data1<-data1[ranked_indices, ]
      filtered_data_h(ranked_data1)
    }
    output$filtered_df_h <- DT::renderDataTable({
      filtered_data_h()
    },extensions = 'Buttons', # Load the 'Buttons' extension
    options = list(
      dom = 'Bfrtip',
      scrollX = TRUE, # Enable horizontal scrolling
      pageLength = 10,# Set the number of rows per page...buttons also
      buttons = list(
        'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
      )))
    
    
  })
  
  

  
  
  ################Second threshold Bepipred
  observe({
    if(input$base_neg>0){
      shinyjs::show("second_thres")
      shinyjs::show("secondary_thres_message")
    }else{
      shinyjs::hide("second_thres")
      shinyjs::hide("secondary_thres_message")
    }
  })
  
  #########################################Generate multiepitope button#############################################
  python_generate <- import("Generate_multiepitope_fnc")
  res1 <-reactive({
    if(input$b_cell_num>0){
      #req(input$B_rank_epitopes)
      #req(input$B_epitope_rank)
      seq_b <- final_epitope_b()
      df_b<-data.frame(seq_b)
    }else{df_b<-data.frame()}
    if(input$ctl_num>0){
      #req(input$C_rank_epitopes)
      #req(input$C_epitope_rank)
      seq_c <- final_epitope_c()
      df_c <-data.frame(seq_c)
    }else{df_c<-data.frame()}
    if(input$htl_num>0){
      #req(input$H_rank_epitopes)
      #req(input$H_epitope_rank)
      seq_h <- final_epitope_h()
      df_h <- data.frame(seq_h)
    }else{df_h<-data.frame()}
    dfs <- list(df_b,df_c,df_h)
    nums <- list(input$b_cell_num,input$ctl_num,input$htl_num)
    linkers <- list(input$bcell_linker,input$ctl_linker,input$htl_linker)
    N_term <- input$N_term
    order <- input$order_epitopes
    res1 <- python_generate$multiepitope_fnc(dfs,nums,linkers,order,N_term)
    df<-data.frame(res1)
    colnames(df)<-c('Sequences')
    df
    return(df)
  })
  observeEvent(input$generation_epitope,{
    output$Multiepitope_sequences <- DT::renderDataTable({
      res1()
    },options=list(scrollX=TRUE,pageLength=20))
  })
  ###############################Sequences to filter ---- generated ones-----------------
  ############################### Sequences to filter --- upload new ones
  observe({
    if (input$gene_seqs == '2') {
      shinyjs::show("seqs_file")  # Show the fileInput
    } else {
      shinyjs::hide("seqs_file")  # Hide the fileInput
    }
  })
  data_seqs <- reactive({
    req(input$seqs_file)
    df <- read.table(input$seqs_file$datapath)
    return(df)
  })
  
  output$seqs_tab <- DT::renderDataTable({
    if(input$gene_seqs=='2'){
      data_seqs()}else{
        res1()
      }
  },options=list(pageLength=10,scrollX=TRUE))
  
  sequences_final<-reactive({
    if(input$gene_seqs=='2'){
      return(data_seqs())}else{
        return(res1())
      }
  })
  ###############################Filtering multiepitope sequences##################
  final_vacc<-reactiveVal(NULL)
  temp_val<-reactiveVal(NULL)
  observeEvent(input$filter_seqs,{
    #Python script-Vaxijen
    withProgress(message="Multiepitope Sequences filtering",value=0,{
    vaxi_fnc <-import("Vaxijen_script")
    eval_fnc <- import("Prot_DeepVac_script")
    data <- sequences_final()
    #Batch size
    batch=20
    #Number of sequences
    num_sequences<-input$num_vac
    #Final dataframe
    final_dataframe <- c()
    #Temp dataframe which is filtered
    counter = 1
    while (counter < length(data$Sequences)){
      
      counter2 <- counter+batch
      if (counter2>length(data$Sequences)){
        counter2 <- length(data$Sequences)
      }
      temp_seqs<-data.frame(Sequences=data$Sequences[counter:counter2])
      ###Vaxijen
      #Parameters
      vax_thres <- input$vacc_vax_thres_num
      vax_target <- input$vacc_vax_target
      df_v <- temp_seqs
      #Functions
      res_vax_vacc <- vaxi_fnc$Vaxijen(df_v,vax_target,vax_thres)
      for (i in colnames(res_vax_vacc)){
        if(i!='Sequences'){temp_seqs[i]=res_vax_vacc[i]}
      }
      
      ###Toxinpred
      #Parameters
      model_tox <- input$vacc_model_tox
      thresh_tox <- input$vacc_tox_thres_num
      print(thresh_tox)
      df_tox <- temp_seqs
      #Functions
      res_tox_vacc <-Toxinpred_fnc(df_tox,path_toxinpred,path_res_toxinpred,model_tox,thresh_tox)
      for (i in colnames(res_tox_vacc)){
        if(i!='Sequences'){j<-paste(i,'_tox',sep="")
        temp_seqs[j]<-res_tox_vacc[i]}
      }
      
      ###Algpred
      model_alg <- input$vacc_model_alg
      thresh_alg <- input$vacc_alg_thres_num
      df_alg <- temp_seqs
      #Functions
      res_alg_vacc <- Algpred_fnc(df_alg,path_algpred,path_res_algpred,model_alg,thresh_alg)
      for (i in colnames(res_alg_vacc)){
        if (i!='Sequences'){j<-paste(i,'_alg',sep="")
        temp_seqs[j]<-res_alg_vacc[i]}
      }
      
      ###Protparam
      df_param <- temp_seqs
      res_prot <- eval_fnc$Protparam(df_param$Sequences)
      
      if (1 %in% input$prot_proper){
        temp_seqs["Molecular weight"] <- res_prot["Molecular weight"]
      }
      if (2 %in% input$prot_proper){
        temp_seqs["Instability index"] <- res_prot["Instability index"]
        temp_seqs["Instability"] <- res_prot["Instability"]
      }
      if (3 %in% input$prot_proper){
        temp_seqs["Aliphatic index"] <- res_prot["Aliphatic index"]
      }
      if (4 %in% input$prot_proper){
        temp_seqs["GRAVY score"] <- res_prot["GRAVY score"]
      }
      
      # ###DeepVacPred
      # df_vacpred <- temp_seqs
      # thres_vac <- input$deepvac_thres_num
      # res_vacpred <- eval_fnc$Predictions_per_30(df_vacpred,thres_vac)
      # 
      # if (1 %in% input$deepvacpred_proper){
      #   temp_seqs$Prediction_Score_DeepVac <-res_vacpred$Prediction_Score_DeepVac
      #   temp_seqs$Prediction_deepvacpred <-res_vacpred$Prediction_deepvacpred
      # }
      # if (2 %in% input$deepvacpred_proper){
      #   temp_seqs$Prediction_per_30 <- res_vacpred$Prediction_per_30
      # }
      
      #Filters
      if (1 %in% input$filters_vac){
        temp_seqs <- subset(temp_seqs,Antigen == TRUE)
      }
      if (2 %in% input$filters_vac){
        temp_seqs <- subset(temp_seqs,Prediction_tox != 'Toxin')
      }
      if (3 %in% input$filters_vac){
        temp_seqs <- subset(temp_seqs,Prediction_alg == 'Non-Allergen')
      }
      # if (4 %in% input$filters_vac){
      #   #DeepVacPred
      #   temp_seqs<- subset(temp_seqs,Prediction_deepvacpred ==TRUE)
      # }
      if (4 %in% input$filters_vac){
        #Protparam
        temp_seqs <- subset(temp_seqs,Instability=='Stable')
      }
      counter<- counter2
      final_dataframe <- rbind(final_dataframe,temp_seqs)
      
      if(length(final_dataframe$Sequences) <= num_sequences){
        incProgress(length(final_dataframe$Sequences)/num_sequences,detail=paste(length(final_dataframe$Sequences),'/',num_sequences))
      }
      if (length(final_dataframe$Sequences) > num_sequences){
        incProgress(1)
        #final_vacc
        final_vacc(final_dataframe)
        break
      }
      if (counter >= length(data$Sequences)){
        final_vacc(final_dataframe)
      }
      
    }
    
    output$vacc_seqs_df <- DT::renderDataTable({
      final_vaccine <- final_vacc()
      seqs<-data.frame(Sequences=final_vaccine$Sequences[1:input$num_vac])
      seqs
    },extensions = 'Buttons', # Load the 'Buttons' extension
    options = list(
      dom = 'Bfrtip', # Specify where to place the buttons
      scrollX = TRUE, # Enable horizontal scrolling
      pageLength = 10, # Set the number of rows per page
      buttons = list(
        'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
      )))
    
    output$final_results_vaccine <- DT::renderDataTable({
      final_vaccine <- final_vacc()[1:input$num_vac,]
      final_vaccine
    },extensions = 'Buttons', # Load the 'Buttons' extension
    options = list(
      dom = 'Bfrtip', # Specify where to place the buttons
      scrollX = TRUE, # Enable horizontal scrolling
      pageLength = 10, # Set the number of rows per page
      buttons = list(
        'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
      )))
   
  })
  })
  ################IMAGES  MANUAL ###################################################
  output$img_001 <- renderImage({
    return(list(src = "Manual_software.files/image001.png",
         contentType="image/png",alt="Figure 1",width="100%"))
  }, deleteFile = FALSE)
  output$img_002 <- renderImage({
    return(list(src = "Manual_software.files/image003.png",
                contentType="image/png",alt="Figure 2",width="100%"))
  }, deleteFile = FALSE)
  output$img_003 <- renderImage({
    return(list(src = "Manual_software.files/image005.png",
                contentType="image/png",alt="Figure 3",width="100%"))
  }, deleteFile = FALSE)
  
  output$img_004 <- renderImage({
    return(list(src = "Manual_software.files/image007.png",
                contentType="image/png",alt="Figure 4",width="100%"))
  }, deleteFile = FALSE)
  
  output$img_005 <- renderImage({
    return(list(src = "Manual_software.files/image009.png",
                contentType="image/png",alt="Figure 5",width="100%"))
  }, deleteFile = FALSE)
  
  output$img_006 <- renderImage({
    return(list(src = "Manual_software.files/image011.png",
                contentType="image/png",alt="Figure 5",width="100%"))
  }, deleteFile = FALSE)
  
  output$img_007 <- renderImage({
    return(list(src = "Manual_software.files/image013.png",
                contentType="image/png",alt="Figure 5",width="100%"))
  }, deleteFile = FALSE)
  
  output$img_008 <- renderImage({
    return(list(src = "Manual_software.files/image015.png",
                contentType="image/png",alt="Figure 5",width="100%"))
  }, deleteFile = FALSE)
  
  output$img_009 <- renderImage({
    return(list(src = "Manual_software.files/image017.png",
                contentType="image/png",alt="Figure 5",width="100%"))
  }, deleteFile = FALSE)
  
  output$img_010 <- renderImage({
    return(list(src = "Manual_software.files/image019.png",
                contentType="image/png",alt="Figure 5",width="100%"))
  }, deleteFile = FALSE)
  
  output$img_011 <- renderImage({
    return(list(src = "Manual_software.files/image021.png",
                contentType="image/png",alt="Figure 5",width="100%"))
  }, deleteFile = FALSE)
  
  output$img_012 <- renderImage({
    return(list(src = "Manual_software.files/image023.png",
                contentType="image/png",alt="Figure 5",width="100%"))
  }, deleteFile = FALSE)
  
  output$img_013 <- renderImage({
    return(list(src = "Manual_software.files/image025.png",
                contentType="image/png",alt="Figure 5",width="100%"))
  }, deleteFile = FALSE)
  
  output$img_014 <- renderImage({
    return(list(src = "Manual_software.files/image027.png",
                contentType="image/png",alt="Figure 5",width="100%"))
  }, deleteFile = FALSE)
  
  output$img_015 <- renderImage({
    return(list(src = "Manual_software.files/image029.png",
                contentType="image/png",alt="Figure 5",width="100%"))
  }, deleteFile = FALSE)
  
  output$img_016 <- renderImage({
    return(list(src = "Manual_software.files/image031.png",
                contentType="image/png",alt="Figure 5",width="100%"))
  }, deleteFile = FALSE)
  
  output$img_017 <- renderImage({
    return(list(src = "Manual_software.files/image033.png",
                contentType="image/png",alt="Figure 5",width="100%"))
  }, deleteFile = FALSE)
  
  output$img_018 <- renderImage({
    return(list(src = "Manual_software.files/image035.png",
                contentType="image/png",alt="Figure 5",width="100%"))
  }, deleteFile = FALSE)
  
  output$img_019 <- renderImage({
    return(list(src = "Manual_software.files/image037.png",
                contentType="image/png",alt="Figure 5",width="100%"))
  }, deleteFile = FALSE)
  
  output$img_020 <- renderImage({
    return(list(src = "Manual_software.files/image039.png",
                contentType="image/png",alt="Figure 5",width="100%"))
  }, deleteFile = FALSE)
  
  
}