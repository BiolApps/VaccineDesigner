library(shiny)
library(shinydashboard)
library(DT)
library(Biostrings)
library(shinyFiles)
library(reticulate)
library(dplyr)
library(shinyjs)
library(shinyWidgets)
#Set path
main_path <- "/home/guest/Revacc/webserver"
setwd(main_path)



##UI
add_colon <- function(string) {
  # Split the string by "*"
  parts <- strsplit(string, "\\*")[[1]]
  
  # Check if ":" exists in parts[2]
  if (length(parts) == 2 && !grepl(":", parts[2])) {
    # Add ":" before the last two digits
    parts[2] <- paste0(substr(parts[2], 1, nchar(parts[2]) - 2), ":", substr(parts[2], nchar(parts[2]) - 1, nchar(parts[2])))
  }
  
  # Combine the parts back together
  paste(parts, collapse = "*")
}

#Parameters
allelle_ctl_df <- read.csv('allelenames',header=FALSE)
allelle_ctl_df<- data.frame(Allelles=allelle_ctl_df$V1)

allele_ctl_df_con <-read.csv('final_alleles.csv',header=TRUE)
allele_ctl_df_con <-data.frame(Alleles=unique(allele_ctl_df_con$MHC))

allele_htl_df_con <-read.csv('final_alleles_htl.csv',header=TRUE)
allele_htl_df_con<-data.frame(Alleles=allele_htl_df_con$Allele)

allelle_htl_df <- read.table("allelle_htl.txt")
colnames(allelle_htl_df) <- 'Allelles'


population_alleles_I <-read.table('Class_1_coded.csv',sep=",",header=TRUE)
population_alleles_I<-data.frame(Alleles=population_alleles_I$Codes)
population_alleles_I<-data.frame(Alleles=sapply(population_alleles_I$Alleles,add_colon))
aa<-data.frame(Alleles=population_alleles_I$Alleles)
population_alleles_I<-aa


population_alleles_II<-read.table('Class_2_coded.csv',sep=",",header=FALSE)
population_alleles_II<-data.frame(Alleles=population_alleles_II$V2)

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


########################################################Functions #####################################################################
blast_function <- function(sequence,database){
  main_path<-getwd()
  file_path <-'test.fasta'
  file_path<-file.path(main_path,file_path)
  
  lines<-c('>protein1',sequence)
  writeLines(lines,file_path,sep="\n")
  file_path=convert_to_linux_path(file_path)
  if (database==0){
    db <- "/home/guest/Desktop/ncbi-blast-2.15.0+/db/swissprot/swissprot"
    script_opt <- 0
  }else{
    db <- "/home/guest/Desktop/ncbi-blast-2.15.0+/db/env_nr/env_nr"
    script_opt <-1
  }
  #if (taxid==0){
  command<-paste('/home/guest/Desktop/ncbi-blast-2.15.0+/bin/blastp -db',db,'-query',file_path,'-out results_blast.out')
  #else{
    #command<-paste('/home/guest/Desktop/ncbi-blast-2.15.0+/bin/blastp -db',db,'-taxids',taxid,'-query',file_path,'-out results_blast.out')}
  
  system(command)
  blast_script <-import("blast_parsing_script")
  if(script_opt==0){
  res_df <- blast_script$blast_parsing_fnc("results_blast.out")}
  else{
    res_df<-blast_script$blast_parsing_fnc_meta("results_blast.out")
  }
  return(res_df)
  }

find_epitope_on_vaccine <- function(vaccine,col_name,data){
  seq <- vaccine$Sequence
  epitopes <-data[[col_name]][[1]]
  start<-c()
  end<-c()
  
  for(i in seq_along(epitopes)){
    ss <- gregexpr(epitopes[i], seq)[[1]]
    start<-c(start,ss[1])
    end<-c(end,ss[1]+nchar(epitopes[i])-1)
  }
  dd<-data.frame(Start=start,End=end)
  return(dd)
}

Epitope_comparison <- function(data,new_epitope_list,type_of_epitope){
  sequence <- data$Sequences
  b_cell <- data$B_cell[[1]]
  ctl <- data$CTL[[1]]
  htl <- data$HTL[[1]]
  
  new_epitope_all <- new_epitope_list
  epitopes_wanted <- c()
  mismatch <- c()
  epitopes_found <- c()
  if(type_of_epitope=='B-cell'){
    sequences_of_interest <-b_cell
  }
  else if (type_of_epitope=='CTL'){
    sequences_of_interest<-ctl
  }
  else if(type_of_epitope=='HTL'){
    sequences_of_interest <-htl
  }
  
  #### Mismatch 0
  seqs <- sequences_of_interest
  
  for(i in sequences_of_interest){
    for(j in new_epitope_all){
      res <- vcountPattern(i,j,max.mismatch=0)
      
      if(res>0){
        epitopes_wanted<-c(epitopes_wanted,j)
        mismatch<-c(mismatch,0)
        epitopes_found <-c(epitopes_found,i)
        seqs<-seqs[!seqs%in%i]
      }
    }
  }
  #### Mismatch 1 if sequences_of_interest >1
  if(length(seqs)>0){
    ###Mismatch 1
    for(i in sequences_of_interest){
      for(j in new_epitope_all){
        res <- vcountPattern(i,j,max.mismatch=1)
        if(res>0){
          epitopes_wanted<-c(epitopes_wanted,j)
          mismatch<-c(mismatch,1)
          epitopes_found <-c(epitopes_found,i)
          sequences_of_interest<-sequences_of_interest[!sequences_of_interest%in%i]
        }
      }
    }
  }else{ res <-data.frame(Initial_Epitopes=epitopes_found,Predicted_Epitopes=epitopes_wanted,Mismatches=mismatch)
    return(res)}
  #### Mismatch 2 if sequences_of_interest >1
  if(length(seqs)>0){
    ###Mismatch 1
    for(i in sequences_of_interest){
      for(j in new_epitope_all){
        res <- vcountPattern(i,j,max.mismatch=2)
        if(res>0){
          epitopes_wanted<-c(epitopes_wanted,j)
          mismatch<-c(mismatch,2)
          epitopes_found <-c(epitopes_found,i)
          sequences_of_interest<-sequences_of_interest[!sequences_of_interest%in%i]
        }
      }
    }
  }else{ res <-data.frame(Initial_Epitopes=epitopes_found,Predicted_Epitopes=epitopes_wanted,Mismatches=mismatch)
  return(res)}
  
  ff <-data.frame(Initial_Epitopes=epitopes_found,Predicted_Epitopes=epitopes_wanted,Mismatches=mismatch)
  return(ff)
  
}


cons_I_fnc <- function(seq,num_peptides,allel){
  #Main path
  main_path <- getwd()
  #Fasta file construction
  file_path <- 'temp1.fasta'
  file_path <- file.path(main_path,file_path)
  lines<-c('>seq1',seq)
  writeLines(lines,file_path,sep="\n")
  file_path=convert_to_linux_path(file_path)
  
  
  #Allele format for command from list allele
  allel_list<-allel
  num_peptides_all <-c()
  for(i in allel_list){
    num_peptides_all<-c(num_peptides_all,num_peptides)
  }
  allel <- paste(allel,collapse=",")
  
  num_peptides_all<-paste(num_peptides_all,collapse=",")
  #Command and run the software
  command<-paste("python ./Local_Tools/mhc_i/src/predict_binding.py","consensus",allel,num_peptides_all,file_path,"> res_con.out")
  command
  #print(command)
  system(command,timeout=20)
  #args <- c("./Local_Tools/mhc_i/src/predict_binding.py",'consensus',allel,num_peptides_all,file_path)#,'>','res_con.out')
  #system2(command,args=args,stdout='res_con.out')
  final_path<-file.path(main_path,'res_con.out')
  #Mapping allele
  #inds <- which(allele_list_ctl$V1 %in% allel_ls)
  #map_allele <- allele_list_ctl$V2[inds]
  return(final_path)
  
}

cons_II_fnc <- function(seq,num_peptides,allel){
  #Main path
  main_path <- getwd()
  #Fasta file construction
  file_path <- 'temp2.fasta'
  file_path <- file.path(main_path,file_path)
  lines<-c('>seq1',seq)
  writeLines(lines,file_path,sep="\n")
  file_path=convert_to_linux_path(file_path)
  
  
  #Allele format for command from list allele
  allel_list<-allel
  #num_peptides_all <-c()
  #for(i in allel_list){
   # num_peptides_all<-c(num_peptides_all,num_peptides)
  #}
  allel <- paste(allel,collapse=",")
  
  #num_peptides_all<-paste(num_peptides_all,collapse=",")
  #Command and run the software
  command<-paste("python ./Local_Tools/mhc_ii/mhc_II_binding.py","consensus3",allel,file_path,num_peptides,"> res_con_II.out")
  
  #print(command)
  system(command,timeout=20)
  
  final_path<-file.path(main_path,'res_con_II.out')
 
  return(final_path)
  
}

consensus_prediction_fnc_v2 <- function(path, cutoff) {
  # Read data
  data <- read.table(path, sep="\t", header=TRUE)
  
  # Filter data based on cutoff
  data_interest <- subset(data, consensus_percentile_rank < cutoff)
  
  # Return empty if no data_interest
  if (nrow(data_interest) == 0) {
    return(data_interest)
  }
  grouped_data <- aggregate(data_interest, 
                            by = list(data_interest$start), 
                            FUN = function(x) x)
  
  data_new<-grouped_data[,!names(grouped_data)%in% c('Group.1','seq_num')]
  data_new$start <- sapply(data_new$start,function(x) x[[1]])
  data_new$end<- sapply(data_new$end,function(x) x[[1]])
  data_new$Start<-data_new$start
  data_new$End<-data_new$end
  data_new<-data_new[,!names(data_new)%in%c('start','end')]
  data_new$length<- sapply(data_new$length,function(x) x[[1]])
  data_new$peptide<-sapply(data_new$peptide,function(x) x[[1]])
  data_new$Sequences<-data_new$peptide
  data_new<-data_new[,!names(data_new)%in%c('peptide')]
  data_new$number_of_alleles<-sapply(data_new$allele,function(x) length(x))
  data_new$Alleles<-data_new$allele
  data_new<-data_new[,!names(data_new)%in%c('allele')]
  data_new$minimum_percentile_score <-sapply(data_new$consensus_percentile_rank,min)

  
  df <- data_new[order(-data_new$number_of_alleles,data_new$minimum_percentile_score),]
  df <- df[, c("Sequences", "number_of_alleles", "minimum_percentile_score","Alleles","Start","End", setdiff(names(df), c("Sequences", "number_of_alleles", "minimum_percentile_score","Alleles","Start","End")))]
  
  
  return(df)
}

consensus3_prediction_fnc_v2 <- function(path, cutoff) {
  # Read data
  data <- read.table(path, sep="\t", header=TRUE)
  
  # Filter data based on cutoff
  data_interest <- subset(data, consensus_percentile_rank < cutoff)
  
  # Return empty if no data_interest
  if (nrow(data_interest) == 0) {
    return(data_interest)
  }
  grouped_data <- aggregate(data_interest, 
                            by = list(data_interest$start), 
                            FUN = function(x) x)
  
  data_new<-grouped_data[,!names(grouped_data)%in% c('Group.1','seq_num')]
  data_new$start <- sapply(data_new$start,function(x) x[[1]])
  data_new$end<- sapply(data_new$end,function(x) x[[1]])
  data_new$Start<-data_new$start
  data_new$End<-data_new$end
  data_new<data_new[,!names(data_new)%in%c('start','end')]
  data_new$length<- sapply(data_new$length,function(x) x[[1]])
  data_new$peptide<-sapply(data_new$peptide,function(x) x[[1]])
  data_new$Sequences<-data_new$peptide
  data_new<-data_new[,!names(data_new)%in%c('peptide')]
  data_new$number_of_alleles<-sapply(data_new$allele,function(x) length(x))
  data_new$Alleles<-data_new$allele
  data_new<-data_new[,!names(data_new)%in%c('allele')]
  data_new$minimum_percentile_score <-sapply(data_new$consensus_percentile_rank,min)
  
  
  df <- data_new[order(-data_new$number_of_alleles,data_new$minimum_percentile_score),]
  df <- df[, c("Sequences", "number_of_alleles", "minimum_percentile_score","Alleles","Start","End", setdiff(names(df), c("Sequences", "number_of_alleles", "minimum_percentile_score","Alleles","Start","End")))]
  
  
  return(df)
}


netchop3_fnc <- function(main_path,vaccine_seq,method,threshold){
  setwd(main_path)
  line1<-'>seq1'
  line2<-vaccine_seq
  file_path <- 'temp_netchop.fasta'
  file_path <- file.path(main_path,file_path)
  
  lines<-c('>seq1',vaccine_seq$Sequence)
  writeLines(lines,file_path,sep="\n")
  
  
  
  file_path=convert_to_linux_path(file_path)
  command='./Local_Tools/netchop-3.1/netchop'
  args <- c(file_path,'-v',method,'-t',threshold,'>','res_netchop.out')
  
  res_path <- 'res_netchop.out'
  system2(command,args=args)
  
  net_script <- import("netchop_script")
  res <- net_script$netchop_fnc('res_netchop.out',vaccine_seq)
  df<-data.frame(res)
  
  df_f <- df[df$C=="S",]
  return(df_f)
}

ranking_results <- function(data){
  cols <- c("Prediction_Score","Instability index")
  
  if("Hybrid.Score_tox" %in% colnames(data)) {
    # Drop the specified columns
    cols<-c(cols,"Hybrid.Score_tox")
    columns_to_drop <- c("ML.Score_tox", "MERCI.Score_tox", "BLAST.Score_tox")
    data <- data[, !(names(data) %in% columns_to_drop)]
  }else{
    cols <-c(cols,"ML.Score_tox")
  }
  
  if("Hybrid.Score_alg" %in% colnames(data)) {
    # Drop the specified columns
    cols<-c(cols,"Hybrid.Score_alg")
    columns_to_drop <- c("ML.Score_alg", "MERCI.Score_alg", "BLAST.Score_alg")
    data <- data[, !(names(data) %in% columns_to_drop)]
    
  }else{
    cols<-c(cols,"ML.Score_alg")
  }
  
  if("X" %in% colnames(data)) {
    # Drop the specified columns
    data <- data[, !(names(data) %in% "X")]
  }
  
  numeric_cols <-sapply(data,is.numeric)
  numeric_cols_indexes <- which(numeric_cols)
  numeric_cols_names <-colnames(data)[(numeric_cols_indexes)]
  
 
  
  
  
  ranked_data <- sapply(data[,cols],order)
  
  if("Prediction_Score" %in% numeric_cols_names){
    ranked_values <-order(data[["Prediction_Score"]],decreasing=TRUE)
    ranked<-data.frame(ranked_data)
    colnames(ranked)<-colnames(ranked_data)
    ranked$Prediction_Score<-ranked_values
  }
  
  if('Prediction_Score' %in% colnames(ranked)){
    ranked$Antigenicity_ranking <- ranked$Prediction_Score
  }
  if('Instability index' %in% colnames(ranked)){
    ranked$Stability_ranking <-ranked[,"Instability index"]
  }
  if('ML.Score_tox' %in% colnames(ranked)){
    ranked$Toxicity_ranking <-ranked[,"ML.Score_tox"]
  }
  if('ML.Score_alg'%in% colnames(ranked)){
    ranked$Allergenicity_ranking <- ranked[,"ML.Score_alg"]
  }
  
  if('Hybrid.Score_tox' %in% colnames(ranked)){
    ranked$Toxicity_ranking <- ranked[,'Hybrid.Score_tox']
  }
  if('Hybrid.Score_alg' %in% colnames(ranked)){
    ranked$Allergenicity_ranking <- ranked[,"Hybrid.Score_alg"]
  }
  
  ranked<- ranked[, !(names(ranked) %in% cols)]
  
  
  sums <- rowSums(ranked)
  
  return(list(order(sums),ranked))
}
to_snake_case <- function(string) {
  words <- strsplit(string, " ")[[1]]  # Split string into words
  if (length(words) == 2) {
    return(tolower(paste(words, collapse = "_")))  # Join words with underscore
  } else {
    return(tolower(string))  # Return original string if only one word
  }
}
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
  print(file_path)
  #Directory of toxinpred
  setwd(path_toxinpred)
  final_path <-file.path(path_res_toxinpred,'results_toxinpred2.csv')
  #Command for toxinpred
  command<-paste('python toxinpred2.py -i',file_path,'-d 2','-m',as.numeric(model),'-t',as.numeric(threshold),'-o',final_path)
  print(command)
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
  net_path <- "/usr/local/bin/netMHCpan"
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
    res<-data.frame(res)
    res$Start<-as.numeric(as.character(res$Position))
    res$End<-res$Start+as.numeric(num_peptides)-1
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
      print(i)
      print('Error in allel?')
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
      print('Error before function?')
      net_script <- import("netMHCpan_script")
      res <- net_script$netMHCpan_fnc('res.out',map_allele)
      if (i==1){
        fin_dat <- res
      }
      else{
        print('Error in combination_allel_fnc?')
        fin_dat <- combination_allel_fnc(res,fin_dat)
      }
    }
      fin_dat<-data.frame(fin_dat)
      fin_dat$Start<-as.numeric(as.character(fin_dat$Position))
      fin_dat$End<-fin_dat$Start+as.numeric(num_peptides)-1
    
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
    res<-data.frame(res)
    res$Start<-as.numeric(as.character(res$Position))
    res$End <-res$Start+as.numeric(num_peptides-1)
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
      fin_dat<-data.frame(fin_dat)
      fin_dat$Start<-as.numeric(as.character(fin_dat$Position))
      fin_dat$End<-fin_dat$Start+as.numeric(num_peptides)-1
    }
    return(fin_dat)
  }
}

############################################################ Server ################################################################
server <- function(input,output,session) {
  Sys.setlocale("LC_NUMERIC", "C")
  ###Shinyjs
  shinyjs::hideElement("B_ranked")
  ########################################## B TEXT AREA ################################################
  b_sequences<-reactiveVal(NULL)
  
  
  observeEvent(input$b_text_update,{
    req(input$B_protein_text)
    
    temp_file <-tempfile()
    writeLines(input$B_protein_text,temp_file)
    fasta <-readAAStringSet(temp_file)
    fasta_df <- data.frame(ID = names(fasta), Sequence = as.character(fasta), stringsAsFactors = FALSE)
    
    output$B_Table <- renderDataTable({
      fasta_df
    },options=list(scrollX=TRUE),selection = 'single')
    
    b_sequences(fasta_df)
  })
  observeEvent(input$B_protein,{
    req(input$B_protein)
    fasta <- readAAStringSet(input$B_protein$datapath)
    fasta_df <- data.frame(ID = names(fasta), Sequence = as.character(fasta), stringsAsFactors = FALSE)
    output$B_Table <- renderDataTable({
      fasta_df
    },options=list(scrollX=TRUE),selection = 'single')
    b_sequences(fasta_df)
  })
  observeEvent(input$b_example_fasta,{
    req(input$b_example_fasta)
    fasta <- readAAStringSet('ATL_LNK_conserved.fasta')
    fasta_df <- data.frame(ID = names(fasta), Sequence = as.character(fasta), stringsAsFactors = FALSE)
    output$B_Table <- renderDataTable({
      fasta_df
    },options=list(scrollX=TRUE),selection = 'single')
    b_sequences(fasta_df)
  })
  
  ####################################################################################
  ####################################################### C TEXT AREA #################################################
  c_sequences<-reactiveVal(NULL)
  
  
  observeEvent(input$c_text_update,{
    req(input$C_protein_text)
    
    temp_file <-tempfile()
    writeLines(input$C_protein_text,temp_file)
    fasta <-readAAStringSet(temp_file)
    fasta_df <- data.frame(ID = names(fasta), Sequence = as.character(fasta), stringsAsFactors = FALSE)
    
    output$C_Table <- renderDataTable({
      fasta_df
    },options=list(scrollX=TRUE),selection = 'single')
    
    c_sequences(fasta_df)
  })
  observeEvent(input$C_protein,{
    req(input$C_protein)
    fasta <- readAAStringSet(input$C_protein$datapath)
    fasta_df <- data.frame(ID = names(fasta), Sequence = as.character(fasta), stringsAsFactors = FALSE)
    output$C_Table <- renderDataTable({
      fasta_df
    },options=list(scrollX=TRUE),selection = 'single')
    c_sequences(fasta_df)
  })
  observeEvent(input$c_example_fasta,{
    req(input$c_example_fasta)
    fasta <- readAAStringSet('ATL_LNK_conserved.fasta')
    fasta_df <- data.frame(ID = names(fasta), Sequence = as.character(fasta), stringsAsFactors = FALSE)
    output$C_Table <- renderDataTable({
      fasta_df
    },options=list(scrollX=TRUE),selection = 'single')
    c_sequences(fasta_df)
  })
  
  observeEvent(input$ctl_allelle_clear,{
    table_proxy <- dataTableProxy("allelle_ctl_tab")
    table_proxy %>% selectRows(selected = NULL)
    
    output$ctl_textarea <- DT::renderDataTable({
      NULL
    })
    
  })
  
  
  
  ################################################# H TEXT AREA #######################################################
  h_sequences<-reactiveVal(NULL)
  
  
  observeEvent(input$h_text_update,{
    req(input$H_protein_text)
    
    temp_file <-tempfile()
    writeLines(input$H_protein_text,temp_file)
    fasta <-readAAStringSet(temp_file)
    fasta_df <- data.frame(ID = names(fasta), Sequence = as.character(fasta), stringsAsFactors = FALSE)
    
    output$H_Table <- renderDataTable({
      fasta_df
    },options=list(scrollX=TRUE),selection = 'single')
    
    h_sequences(fasta_df)
  })
  observeEvent(input$H_protein,{
    req(input$H_protein)
    fasta <- readAAStringSet(input$H_protein$datapath)
    fasta_df <- data.frame(ID = names(fasta), Sequence = as.character(fasta), stringsAsFactors = FALSE)
    output$H_Table <- renderDataTable({
      fasta_df
    },options=list(scrollX=TRUE),selection = 'single')
    h_sequences(fasta_df)
  })
  observeEvent(input$h_example_fasta,{
    req(input$h_example_fasta)
    fasta <- readAAStringSet('ATL_LNK_conserved.fasta')
    fasta_df <- data.frame(ID = names(fasta), Sequence = as.character(fasta), stringsAsFactors = FALSE)
    output$H_Table <- renderDataTable({
      datatable(fasta_df)
    },options=list(scrollX=TRUE),selection = 'single')
    h_sequences(fasta_df)
  })
  
  observeEvent(input$htl_allelle_clear,{
    table_proxy <- dataTableProxy("allelle_htl_tab")
    table_proxy %>% selectRows(selected = NULL)
    
    output$htl_textarea <- DT::renderDataTable({
      NULL
    })
    
  })
  
  ##########################################################################################################  
  ###################################################################################
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
  

  ###########################Check button sequence --- python script#########################
  python_protein_check <- import("Protein_check")
  #Bcell
  observeEvent(input$B_Table_rows_selected,{
    selectedRow <- input$B_Table_rows_selected
    if (length(selectedRow)>0){
      selectedData <- b_sequences()[selectedRow, ]
      res <- python_protein_check$Protein_check(selectedData,opt="B_cell")
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
    text <-paste("Bepipred Threshold:",input$bepi_thres_num,"Epitope length:",input$amino_region_ii,sep="  ")
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
    if(input$c_software_opt==0){
      text<-paste("Method : NetMHCpan","Length:",input$c_length,sep=" ")
    }
    else{
      text<-paste("Method : Consensus (TepiTool)","Length:",input$c_length_con,"Cutoff value:",input$c_cutoff_val,sep=" ")
    }
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
    if(input$h_software_opt==0){
      text<-paste("Method : NetMHCIIpan","Length:",input$h_length,sep=" ")
    }
    else{
      text<-paste("Method : Consensus (TepiTool)","Length:",input$h_length_con,"Cutoff value:",input$h_cutoff_val,sep=" ")
    }
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
  # observeEvent(input$allel_opt,{
  #   data <- allelle_ctl_df
  #   opt <- input$allel_opt
  #   if(opt=="all"){
  #     allelle_f <- data
  #   }
  #   else{
  #   allel_sel <- data[grep(opt,data$Allelles),]
  #   allelle_f <-data.frame(Allelles=allel_sel)
  #   }
  #   output$allelle_ctl_tab <- DT::renderDataTable({allelle_f},options=list(pageLength=20,scrollX=TRUE,selection="multiple"))
  # 
  # })
  output$allelle_ctl_tab <- DT::renderDataTable({
    data <-allelle_ctl_df
    allelle_f<-data
    allelle_f
  },options=list(pageLength=20,scrollX=TRUE,selection="multiple"))
  
  
  # Track selected rows
  ctl_selectedRows <- reactiveVal()
  #Update selected rows on table selection
  observeEvent(input$ctl_allelle_btn, {
    ctl_selectedRows(input$allelle_ctl_tab_rows_selected)
    rows <- ctl_selectedRows()
    data <- allelle_ctl_df
    allelle_f <- data
    
    data <- allelle_f[rows,]
    df <- data.frame(data)
    colnames(df) <- c('Selected Allelles')
    output$ctl_textarea <- DT::renderDataTable({
      df
    })
    
  })
  
  ### Select Allelle HTL
  # observeEvent(input$allel_opt_h,{
  #   data <- allelle_htl_df
  #   opt <- input$allel_opt_h
  #   if(opt=="all"){
  #     allelle_f <- data
  #   }
  #   else{
  #     allel_sel <- data[grep(opt,data$Allelles),]
  #     allelle_f <-data.frame(Allelles=allel_sel)
  #   }
  #   
  #})
  output$allelle_htl_tab <- DT::renderDataTable({
    data <- allelle_htl_df
    allelle_f <- data
    allelle_f
  },options=list(pageLength=20,scrollX=TRUE,selection="multiple"))
  
  htl_selectedRows <- reactiveVal()
  #Update selected rows on table selection
  observeEvent(input$htl_allelle_btn, {
    htl_selectedRows(input$allelle_htl_tab_rows_selected)
    rows <- htl_selectedRows()
    data <- allelle_htl_df
    allelle_f <- data
    
    data <- allelle_f[rows,]
    df <- data.frame(data)
    colnames(df) <- c('Selected Allelles')
    output$htl_textarea <- DT::renderDataTable({
      df
    })
    
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
    b_all_files()},options=list(pageLength=100,scrollX=TRUE))
  
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
       
      ranked_data},options=list(pageLength=100,scrollX=TRUE))
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
      data.frame("Sequences selected"=dat_rows$Sequences)},options=list(pageLength=100,scrollX=TRUE))
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
    c_all_files()},options=list(pageLength=100,scrollX=TRUE))
  
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
      ranked_data},options=list(pageLength=100,scrollX=TRUE))
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
      data.frame("Sequences selected"=dat_rows$Sequences)},options=list(pageLength=100,scrollX=TRUE))
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
    h_all_files()},options=list(pageLength=100,scrollX=TRUE))
  
  
  
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
      ranked_data},options=list(pageLength=100,scrollX=TRUE))
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
      data.frame("Sequences selected"=dat_rows$Sequences)},options=list(pageLength=100,scrollX=TRUE))
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
                 if(input$bepipred_opt==1){
                 amino_region <- input$amino_region}
                 else{
                   amino_region<-input$amino_region_ii
                 }
                 ###Secondary threshold
                 bepi_sec_thres <- input$second_thres_num
                 ###Dataframe
                 df <- sel_BData()
                 ##Functio
                 res_bepi <- bepi_fnc$Bepipred3(df,path_bepipred,path_res_bepi,bepi_thres,base_neg,amino_region,bepi_sec_thres,bepipred_opt)
                 ranked_indices <- order(res_bepi[["Bepipred mean score"]],decreasing = TRUE)
                 
                 res_bepi <- res_bepi[ranked_indices, ]
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
                   res_bepi <- lapply(res_bepi, function(col) {
                     if (is.numeric(col)&&all(col%%1 != 0)){
                        formatC(col, format = "f", digits = 4)}
                      else {
                       col
                      }
                 })
                   res_bepi <- data.frame(res_bepi)
                   res_bepi
                 },extensions = 'Buttons', # Load the 'Buttons' extension
                 options = list(
                   dom = 'Bfrtip', # Specify where to place the buttons
                   scrollX = TRUE, # Enable horizontal scrolling
                   buttons = list(
                     'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                   )))
                 
                 setwd(main_path)
                 
                 output$bepi_csv <- DT::renderDataTable({
                   data <- read.csv("./temp_data/Bepipred/raw_output.csv")
                   data <-lapply(data, function(col) {
                     if (is.numeric(col)&&all(col%%1 != 0)){
                       formatC(col, format = "f", digits = 4)}
                     else {
                       col
                     }
                   }
                   )
                   data <- data.frame(data)
                   data
                 },extensions = 'Buttons', # Load the 'Buttons' extension
                 options = list(
                   dom = 'Bfrtip', # Specify where to place the buttons
                   scrollX = TRUE, # Enable horizontal scrolling
                   pageLength=100,
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
                     res_vax_b <-lapply(res_vax_b, function(col) {
                       if (is.numeric(col)&&all(col%%1 != 0)){
                         formatC(col, format = "f", digits = 4)}
                       else {
                         col
                       }
                     }
                     )
                     res_vax_b <- data.frame(res_vax_b)
                     res_vax_b
                   },extensions = 'Buttons', # Load the 'Buttons' extension
                   options = list(
                     dom = 'Bfrtip', # Specify where to place the buttons
                     scrollX = TRUE, # Enable horizontal scrolling
                     pageLength=100,
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
                     pageLength=100,
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
                     pageLength=100,
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
                   final_df<-lapply(final_df, function(col) {
                     if (is.numeric(col)&&all(col%%1 != 0)){
                       formatC(col, format = "f", digits = 4)}
                     else {
                       col
                     }
                   }
                   )
                   final_df <- data.frame(final_df)
                   final_df
                 },extensions = 'Buttons', # Load the 'Buttons' extension
                 options = list(
                   dom = 'Bfrtip', # Specify where to place the buttons
                   scrollX = TRUE, # Enable horizontal scrolling
                   pageLength=100,
                   buttons = list(
                     'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                   )))
                 
                 
                 
                 shinyjs::show("b_res_col1")
                 output$res_df_b <- DT::renderDataTable({
                   final_df<-lapply(final_df, function(col) {
                     if (is.numeric(col)&&all(col%%1 != 0)){
                       formatC(col, format = "f", digits = 4)}
                     else {
                       col
                     }
                   }
                   )
                   final_df <- data.frame(final_df)
                   final_df
                 },extensions = 'Buttons', # Load the 'Buttons' extension
                 options = list(
                   dom = 'Bfrtip', # Specify where to place the buttons
                   scrollX = TRUE, # Enable horizontal scrolling
                   pageLength=100,
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
      pageLength=100,
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
      pageLength=100,
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
                 if(input$c_software_opt==0){
                   detail <-'netMHCpan'
                 }
                 else{
                   detail<-'Consensus'
                 }
                 withProgress(message="CTL epitope computing",detail=detail,value=0,{
                 #Parameters
                 #Number of peptides
                if (input$c_software_opt==0){
                 num_peptides_ctl <- input$c_length
                 #Alleles
                 ctl_rows <- ctl_selectedRows()
                 alleles_ctl <- allelle_ctl_df[ctl_rows,]
                 
                 #Dataframe
                 df <- sel_CData()
                 seq <- df$Sequence
                 ##Function
                 res_df <- netmhc_fnc(seq,num_peptides_ctl,alleles_ctl,allele_list_ctl)
                 
                 final_df_c <- res_df}
                else{
                  num_peptides<-input$c_length_con
                  cutoff<-input$c_cutoff_val
                  
                  rows<-ctl_selectedRows_con()
                  allele<-allele_ctl_df_con[rows,]
                  
                  df<-sel_CData()
                  seq<-df$Sequence
                  res_path<-cons_I_fnc(seq,num_peptides,allel=allele)
                  

                  res_df <- consensus_prediction_fnc_v2(res_path,cutoff)
                  final_df_c<-res_df
                }
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
                   pageLength=100,
                   buttons = list(
                     'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                   )))
                 
                 setwd(main_path)
                 
                 output$net_I_text <- renderText({
                   text <- 'Epitopes predicted'
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
                     pageLength=100,
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
                     pageLength=100,
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
                     pageLength=100,
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
                   pageLength=100,
                   buttons = list(
                     'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                   )))
                 
                 
                 output$res_df_c <- DT::renderDataTable({
                   final_df_c
                 },extensions = 'Buttons', # Load the 'Buttons' extension
                 options = list(
                   dom = 'Bfrtip', # Specify where to place the buttons
                   scrollX = TRUE, # Enable horizontal scrolling
                   pageLength=100,
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
                   if(input$c_software_opt==0){
                   ranking <- c('Weak Binders Counts','Strong Binders Counts')}
                   else{
                     ranking<-c()
                   }
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
      pageLength=100,
      scrollX = TRUE, # Enable horizontal scrolling
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
      scrollX = TRUE,
      pageLength=100,# Enable horizontal scrolling
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
                   if(input$h_software_opt==0){
                     detail <-'netMHCIIpan'
                   }
                   else{
                     detail<-'Consensus'
                   }
                   withProgress(message="HTL epitope computing",detail=detail,value=0,{
                 #Parameters
                 #Number of peptides
                  if(input$h_software_opt==0){
                 num_peptides <- input$h_length
                 #Alleles
                 htl_rows <- htl_selectedRows()
                 alleles_htl <- allelle_htl_df[htl_rows,]
                 
                 #Dataframe
                 df <- sel_HData()
                 seq <- df$Sequence
                 ##Function
                 res_df <- netmhcII_fnc(seq,num_peptides,alleles_htl,allele_list_htl)
                 final_df_h <- res_df}
                  else{
                    num_peptides<-input$h_length_con
                    cutoff<-input$h_cutoff_val
                    
                    rows<-htl_selectedRows_con()
                    allele<-allele_htl_df_con[rows,]
                    
                    df<-sel_HData()
                    seq<-df$Sequence
                    res_path<-cons_II_fnc(seq,num_peptides,allel=allele)
                    
                    res_df <- consensus3_prediction_fnc_v2(res_path,cutoff)
                    
                    final_df_h<-res_df
                  }
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
                   pageLength=100,
                   buttons = list(
                     'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                   )))
                 
                 setwd(main_path)
                 output$net_II_text <- renderText({
                   text <- 'HTL Epitopes predicted'
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
                     pageLength=100,
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
                     pageLength=100,
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
                     pageLength=100,
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
                   pageLength=100,
                   buttons = list(
                     'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
                   )))
                 
                 
                 
                 output$res_df_h <- DT::renderDataTable({
                   final_df_h
                 },extensions = 'Buttons', # Load the 'Buttons' extension
                 options = list(
                   dom = 'Bfrtip', # Specify where to place the buttons
                   scrollX = TRUE, # Enable horizontal scrolling
                   pageLength=100,
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
                   if(input$h_software_opt==0){
                   ranking <- c('Weak Binders Counts','Strong Binders Counts')}
                  else{
                    ranking<-c()
                  }
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
      pageLength=100,
      scrollX = TRUE, # Enable horizontal scrolling
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
      pageLength=100,
      scrollX = TRUE, # Enable horizontal scrolling
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
    print(df)
    df
    return(df)
  })
  observeEvent(input$generation_epitope,{
    output$Multiepitope_sequences <- DT::renderDataTable({
      res1()
    },options=list(scrollX=TRUE,pageLength=100))
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
  filtered_vaccine<-reactiveVal(NULL)
  temp_val<-reactiveVal(NULL)
  observeEvent(input$filter_seqs,{
    #Python script-Vaxijen
    withProgress(message="Multiepitope Sequences filtering",value=0,{
    vaxi_fnc <-import("Vaxijen_script")
    eval_fnc <- import("Prot_DeepVac_script")
    data <- sequences_final()
    #Batch size
    if(length(data$Sequences)<20){
      batch<-length(data$Sequences)-1
    }
    else{
    batch=20}
    #Number of sequences
    num_sequences<-input$num_vac
    #Final dataframe
    final_dataframe <- c()
    #Temp dataframe which is filtered
    counter = 1
    while (counter < length(data$Sequences)){
      print(batch)
      counter2 <- counter+batch
      if (counter2>length(data$Sequences)){
        counter2 <- length(data$Sequences)
      }
      
      temp_seqs<-data[counter:counter2,]
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
      print('--------------------Toxinpred----------------------')
      #Parameters
      model_tox <- input$vacc_model_tox
      thresh_tox <- input$vacc_tox_thres_num
      df_tox <- temp_seqs
      print(df_tox)
      #Functions
      res_tox_vacc <-Toxinpred_fnc(df_tox,path_toxinpred,path_res_toxinpred,model_tox,thresh_tox)
      for (i in colnames(res_tox_vacc)){
        if(i!='Sequences'){j<-paste(i,'_tox',sep="")
        temp_seqs[j]<-res_tox_vacc[i]}
      }
      print('---------------------Algpred---------------------')
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
      print('------------------------------Protparam-----------------')
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
      
     
      print('---------------------------Filters--------------------------')
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
      print('-----------------------------After filters------------------------')
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
      pageLength = 100, # Set the number of rows per page
      buttons = list(
        'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
      )),selection='single')
    
    output$final_results_vaccine <- DT::renderDataTable({
      final_vaccine <- final_vacc()[1:input$num_vac,]
      ranked_all <- ranking_results(final_vaccine)
      ranked_indices<-ranked_all[[1]]
      ranking_df <-ranked_all[[2]]
      
      final_vaccine$Ranking <- ranked_indices
      for(i in colnames(ranking_df)){
        final_vaccine[,i]<-ranking_df[,i]
      }
      f<-final_vaccine[order(final_vaccine$Ranking),]
      f <- f[, c("Sequences","Ranking","Stability_ranking","Antigenicity_ranking","Toxicity_ranking","Allergenicity_ranking", setdiff(names(f), 
    c("Sequences","Ranking","Stability_ranking","Antigenicity_ranking","Toxicity_ranking","Allergenicity_ranking","Toxicity_ranking")))]
      
      filtered_vaccine(f)
      f
    },extensions = 'Buttons', # Load the 'Buttons' extension
    options = list(
      dom = 'Bfrtip', # Specify where to place the buttons
      scrollX = TRUE, # Enable horizontal scrolling
      pageLength = 100, # Set the number of rows per page
      buttons = list(
        'copy', 'csv', 'excel', 'pdf', 'print' # Add the buttons you want
      )),selection='single')
   
  })
  })
  
  ################################################################### AFTER REVISION ################################
  sel_vaccine_sequence <- eventReactive(input$final_results_vaccine_rows_selected, {
    selectedRow <- input$final_results_vaccine_rows_selected
    if (length(selectedRow) > 0) {
      # Get the selected row from table 1
      selectedData <- filtered_vaccine()[selectedRow, ]
      
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
  output$sel_vaccine_seq <- DT::renderDataTable({
    sel_vaccine_sequence()
  },options=list(scrollX=TRUE,pageLength=10,selection='single'))
  
  ############################################################################## Population coverage ---- Actual
  output$population_allele_I_table <- DT::renderDataTable({
    data <-population_alleles_I
    allelle_f<-data
    allelle_f
  },options=list(pageLength=10,scrollX=TRUE,selection="multiple"))
  
  output$population_allele_II_table <- DT::renderDataTable({
    data <-population_alleles_II
    allelle_f<-data
    allelle_f
  },options=list(pageLength=10,scrollX=TRUE,selection="multiple"))
  observe({
    if(input$population_sel){
  output$tables_popul_ui<-renderUI({
    column(width=12,
          column(width=6,h3("Select Alleles for Class I"),dataTableOutput("population_allele_I_table"),align='center'),
          column(width=6,h3("Select Alleles for Class II"),dataTableOutput("population_allele_II_table"),align='center')
    )
  })
  }
    else{
      output$tables_popul_ui<-renderUI({NULL})
    }
  
  })
  observeEvent(input$calc_popul,{
    withProgress(message="Population coverage",value=0.1,{
    if(input$population_sel){
      
      setwd(main_path)
      #Step 1 --- Find all the epitopes sequences and alleles class I class II
      ctl_epitopes_table <- c_all_files()
      ctl_epitope_sequences<-ctl_epitopes_table$Sequences
      #ctl_epitope_allelle<-ctl_epitopes_table$Alleles
      
      htl_epitopes_table <- h_all_files()
      htl_epitope_sequences<-htl_epitopes_table$Sequences
      #htl_epitope_allelle<-htl_epitopes_table$Alleles
      
      #Step 2 ---- Find all the CTL , HTL epitopes presented on the vaccine sequence
      vaccine_seq <- sel_vaccine_sequence()
      
      #CTL sequences and alleles
      ctl_seqs <-c()
      ctl_indexes<-c()
      for (i in seq_along(ctl_epitope_sequences)){
        if (grepl(ctl_epitope_sequences[i],vaccine_seq)){
          ctl_indexes<-c(ctl_indexes,i)
          ctl_seqs<-c(ctl_seqs,ctl_epitope_sequences[i])
        }
        
      }
      #ctl_alleles <-ctl_epitope_allelle[ctl_indexes]
      #HTL sequences and alleles
      htl_seqs <-c()
      htl_indexes<-c()
      for (i in seq_along(htl_epitope_sequences)){
        if(grepl(htl_epitope_sequences[i],vaccine_seq)){
          htl_indexes<-c(htl_indexes,i)
          htl_seqs<-c(htl_seqs,htl_epitope_sequences[i])
        }
      }
      
      all_htl_seqs <- population_alleles_II
      all_ctl_seqs <-population_alleles_I
      htl_rows <-input$population_allele_II_table_rows_selected
      ctl_rows <-input$population_allele_I_table_rows_selected
      #### Here from Class_2_coded.csv$V2
      final_htl_alleles<-all_htl_seqs[htl_rows,]
      #### Here from Class_1_coded.csv$Codes
      ctl_alleles <- all_ctl_seqs[ctl_rows,]
      combine_i<-c()
      combine_ii<-c()
      
      
      for(i in ctl_seqs){
          combine_i<-c(combine_i,paste(i,paste(ctl_alleles,collapse=', '),sep="\t"))
      }
      
      for(i in htl_seqs){
        combine_ii<-c(combine_ii,paste(i,paste(final_htl_alleles,collapse=', '),sep="\t"))
      }
      
      writeLines(combine_i,"input_i.txt")
      writeLines(combine_ii,"input_ii.txt")
      
      #### 2 command and a script for mining the information
      population_selection <- input$population
      command <- "python"
      args_i<-c("./Local_Tools/population_coverage/calculate_population_coverage.py", "-p", sprintf("\"%s\"", population_selection), "-c", "I", "-f", "input_i.txt", "--plot", ".")
      args_ii<-c("./Local_Tools/population_coverage/calculate_population_coverage.py", "-p", sprintf("\"%s\"", population_selection), "-c", "II", "-f", "input_ii.txt", "--plot", ".")
      
      system2(command,args=args_i,stdout="res_i.txt")
      system2(command,args=args_ii,stdout="res_ii.txt")
      incProgress(0.5)

      #Step 3 --- Code for rendering the tables and graphs
      check_var <- read.table('res_i.txt',header = TRUE, skip = 1, stringsAsFactors = FALSE, sep = "\t", fill = TRUE)[1,1]
      pattern<-'data for following combinations are not available'
      if(grepl(pattern,check_var)){
        var<-0
      }
      else{
        var<-1
      }
      if(var==1){
      output$table_1_i <- renderDataTable({
        data_i <-read.table('res_i.txt', header = TRUE, skip = 1, stringsAsFactors = FALSE, sep = "\t", fill = TRUE)
        split_index <-which(data_i[,1]=="population/area")
        table_1 <- data_i[1:(split_index-1),]
        data.frame(table_1)
      },options=list(scrollX=TRUE,pageLength=10,selection='single'))
      output$table_2_i <- renderDataTable({
        data_i <-read.table('res_i.txt', header = TRUE, skip = 1, stringsAsFactors = FALSE, sep = "\t", fill = TRUE)
        split_index <-which(data_i[,1]=="population/area")
        table_2 <- data_i[(split_index+1):(nrow(data_i)-1),]
        data.frame(table_2)
      })
      output$pop_text_i<-renderText({
        NULL
      })
      output$popcov_i <- renderImage({
        plot_i <-paste0("./popcov_",to_snake_case(population_selection),"_i.png")
        list(src = plot_i)  # Adjust the height as needed
      }, deleteFile = FALSE)
      
      }
      else{
        
        output$table_1_i<-renderDataTable({
          NULL
        })
        output$table_2_i<-renderDataTable({
          NULL
        })
        output$pop_text_i<-renderText({
          'Combinations not found'
        })
        output$popcov_i<-NULL
      }
      output$pop_cov_ui <-renderUI({
        box(width=NULL,status="primary",solidHeader=TRUE,title='Class I Population Coverage',textOutput("pop_text_i"),column(width=6,dataTableOutput("table_1_i")),column(width=6,dataTableOutput('table_2_i')),br(),
        column(width=12,imageOutput('popcov_i'),align='center')   
        )
      })
      
      
      check_var_ii<-read.table('res_ii.txt', header = TRUE, skip = 1, stringsAsFactors = FALSE, sep = "\t", fill = TRUE)[1,1]
      if(grepl(pattern,check_var_ii)){
        var_ii<-0
      }
      else{
        var_ii<-1
      }
      if(var_ii==1){
      output$table_1_ii <- renderDataTable({
        data_i <-read.table('res_ii.txt', header = TRUE, skip = 1, stringsAsFactors = FALSE, sep = "\t", fill = TRUE)
        split_index <-which(data_i[,1]=="population/area")
        table_1 <- data_i[1:(split_index-1),]
        data.frame(table_1)
      },options=list(scrollX=TRUE,pageLength=10,selection='single'))
      
      output$table_2_ii <- renderDataTable({
        data_i <-read.table('res_ii.txt', header = TRUE, skip = 1, stringsAsFactors = FALSE, sep = "\t", fill = TRUE)
        split_index <-which(data_i[,1]=="population/area")
        table_2 <- data_i[(split_index+1):(nrow(data_i)-1),]
        data.frame(table_2)
      })
      output$popcov_ii <- renderImage({
        plot_ii <-paste0("./popcov_",to_snake_case(population_selection),"_ii.png")
        list(src = plot_ii)  # Adjust the height as needed
      }, deleteFile = FALSE)
      output$pop_text_ii<-renderText({
        NULL
      })
      }
      
      else{
        output$table_1_ii<-renderDataTable({
          NULL
        })
        output$table_2_ii<-renderDataTable({
          NULL
        })
        output$popcov_ii<-NULL
        output$pop_text_ii<-renderText({
          "Combinations not found"
        })
      }
      output$pop_cov_ui_ii <-renderUI({
        box(width=NULL,status="primary",solidHeader=TRUE,title='Class II Population Coverage',textOutput("pop_text_ii"),column(width=6,dataTableOutput("table_1_ii")),column(width=6,dataTableOutput('table_2_ii')),br(),
            column(width=12,imageOutput('popcov_ii'),align='center') 
        )
      })
      
      
    }
    })
  })
  
  ############################################################################## Show-hide elements
  observe({
    
    if(input$population_sel){
      shinyjs::show("population")
      shinyjs::show("calc_popul")
      shinyjs::show("pop_cov_ui")
      shinyjs::show('pop_cov_ui_ii')
    }
    else{
      shinyjs::hide("population")
      shinyjs::hide("calc_popul")
      shinyjs::hide("pop_cov_ui")
      shinyjs::hide('pop_cov_ui_ii')
    }
  })
  
  ############################################################################## 2. Protein adjuvants #######################################
  observe({
    # Update the text input based on the selected value
    if(input$adjuvant==0){
      text <- 'AKFVAAWTLKAAA'
    }
    else if(input$adjuvant==1){
      text <-'QQKFQFQFEQQ'
    }
    else if(input$adjuvant==2){
      text <-'QYIKANSKFIGTEL'
    }
    else if(input$adjuvant==3){
      text <-'QRLSTGSRINSAKDDAAGLQIA'
    }
    else if(input$adjuvant==4){
      text <-'CSARDILVKGVDESGASRFTVLFPSGTPLPEEVD'
    }
    else if(input$adjuvant==5){
      text <-'GFKGLVDDADIGNEY'
    }
    else{
      text <- ''
    }
    updateTextInput(session, "N_term", value = text)
  })
  
  
  ############################################################################## 3. Re-analyze - Epitopes ################################################
  ##hide button if checkbox uncheck
  observe({
    if(is.null(input$c_length)){
      updateNumericInput(session, "c_length", value = 9)
    }
    
    if(input$re_analyze_selection){
      shinyjs::show("re_analyze_btn")
      shinyjs::show("re_analyze_params")
      shinyjs::show('bepi_re_analyze')
      shinyjs::show('ctl_re_analyze')
      shinyjs::show('htl_re_analyze')
      shinyjs::show('re_analyze_params_ui')
      shinyjs::show('vaccine_diagram')
      shinyjs::show('b_cell_compared_table')
      shinyjs::show('ctl_compared_table')
      shinyjs::show('htl_compared_table')
      shinyjs::show('bepi_ui_1')
      shinyjs::show('ctl_ui_1')
      shinyjs::show('htl_ui_1')
      shinyjs::show('bepi_ui_2')
      shinyjs::show('ctl_ui_2')
      shinyjs::show('htl_ui_2')
    }
    else{
      shinyjs::hide("re_analyze_btn")
      shinyjs::hide("re_analyze_params")
      shinyjs::hide('bepi_re_analyze')
      shinyjs::hide('ctl_re_analyze')
      shinyjs::hide('htl_re_analyze')
      shinyjs::hide('re_analyze_params_ui')
      shinyjs::hide('vaccine_diagram')
      shinyjs::hide('b_cell_compared_table')
      shinyjs::hide('ctl_compared_table')
      shinyjs::hide('htl_compared_table')
      shinyjs::hide('bepi_ui_1')
      shinyjs::hide('ctl_ui_1')
      shinyjs::hide('htl_ui_1')
      shinyjs::hide('bepi_ui_2')
      shinyjs::hide('ctl_ui_2')
      shinyjs::hide('htl_ui_2')
    }
  })
  ##Bepi/c-net/h-net
  bepi_text<-reactive({if(input$bepipred_opt==1){
    if (input$base_neg > 0){
      text <- paste("Bepipred Threshold:",input$bepi_thres_num,"Base to neglect:",input$base_neg,"Secondary Threshold: ",input$second_thres_num,"Minimum Amino Region:",input$amino_region,sep="  ")
    }else{text <- paste("Bepipred Threshold:",input$bepi_thres_num,"Base to neglect:",input$base_neg,"Minimum Amino Region:",input$amino_region,sep="  ")}
  }else{
    text <-paste("Bepipred Threshold:",input$bepi_thres_num,"Epitope length:",input$amino_region,sep="  ")
  }
    
    text})
  c_net_text<-reactive({
    
    if(input$c_software_opt==0){
      if(is.null(input$c_length)){
        text<-paste("Method : NetMHCpan","Length:",9,sep=" ")
      }
      else{
      text<-paste("Method : NetMHCpan","Length:",input$c_length,sep=" ")}
    }
    else{
      
      text<-paste("Method : Consensus (TepiTool)","Length:",input$c_length_con,"Cutoff value:",input$c_cutoff_val,sep=" ")
    }
    text
  })
  h_net_text <-reactive({
    if(input$h_software_opt==0){
      if(is.null(input$h_length)){
        text<-paste("Method : NetMHCIIpan","Length:",15,sep=" ")
      }
      else{text<-paste("Method : NetMHCIIpan","Length:",input$h_length,sep=" ")}
    }
    else{
      text<-paste("Method : Consensus (TepiTool)","Length:",input$h_length_con,"Cutoff value:",input$h_cutoff_val,sep=" ")
    }
    text
  })

  ## ui element re_analyze_params
  output$re_analyze_params_ui <- renderUI({
    bepi_text_temp <- bepi_text()
    column(width=12,column(width=4,h3('B cell epitopes parameters'),
                           br(),h4(bepi_text_temp),br(),align='center'
    ),
    column(width=4,h3('CTL epitopes parameters'),
           br(),h4(c_net_text()),br(),h3('Select MHC I alleles'),br(),uiOutput('allele_table_ui_I'),align='center'
    ),
    column(width=4,h3('HTL epitopes parameters'),
           br(),h4(h_net_text()),br(),h3('Select MHC II alleles'),br(),uiOutput('allele_table_ui_II'),align='center')
    )
    
  })
  ## re-analyze btn all the steps here
  
  #### epitope selection for C,H
  #### CTL allele netMHCpan
  output$allele_ctl_tab_re <- DT::renderDataTable({
    data <-allelle_ctl_df
    allelle_f<-data
    allelle_f
  },options=list(pageLength=10,scrollX=TRUE,selection="multiple"))
  #### HTL allele netMHCIIpan
  output$allele_htl_tab_re <- DT::renderDataTable({
    data <- allelle_htl_df
    allelle_f <- data
    allelle_f
  },options=list(pageLength=10,scrollX=TRUE,selection="multiple"))
  ###### CTL allele Consensus
  output$allele_ctl_tab_con_re <- DT::renderDataTable({
    data <-allele_ctl_df_con
    allelle_f<-data
    allelle_f
  },options=list(pageLength=10,scrollX=TRUE,selection="multiple"))
  ###### HTL allele Consensus
  output$allele_htl_tab_con_re <- DT::renderDataTable({
    data <-allele_htl_df_con
    allelle_f<-data
    allelle_f
  },options=list(pageLength=10,scrollX=TRUE,selection="multiple"))
  ###### UI for CTL allele table
  output$allele_table_ui_I <-renderUI({
    if(input$c_software_opt==0){
      dataTableOutput("allele_ctl_tab_re")
    }else{
      dataTableOutput("allele_ctl_tab_con_re")
    }
  })
  output$allele_table_ui_II<-renderUI({
    if(input$h_software_opt==0){
      dataTableOutput('allele_htl_tab_re')
    }else{
      dataTableOutput('allele_htl_tab_con_re')
    }
  })
  
  observeEvent(input$re_analyze_btn,{
    withProgress(message='Re-analyze',detail="Prepare",value=0.1,{
######################
    #1. Find all the epitopes on vaccine sequence analyzed and corresponding alleles for the evaluation
    ##B epitope sequences - all##
    b_epitopes_table <- b_all_files()
    b_epitopes_sequences <-b_epitopes_table$Sequences
    ##CTL epitope sequences -all and alleles##
    ctl_epitopes_table <- c_all_files()
    ctl_epitope_sequences<-ctl_epitopes_table$Sequences
    ctl_epitope_allelle<-ctl_epitopes_table$Alleles
    ##HTL epitope sequences - all and alleles
    htl_epitopes_table <- h_all_files()
    htl_epitope_sequences<-htl_epitopes_table$Sequences
    htl_epitope_allelle<-htl_epitopes_table$Alleles
    #### Vaccine sequence to be analyzed
    vaccine_seq <- sel_vaccine_sequence()
    
    #B epitope sequences presented on vaccine
    b_seqs <- c()
    for (i in seq_along(b_epitopes_sequences)){
      if(grepl(b_epitopes_sequences[i],vaccine_seq)){
        b_seqs<-c(b_seqs,b_epitopes_sequences[i])
      }
    }
    #CTL sequences and alleles
    ctl_seqs <-c()
    ctl_indexes<-c()
    for (i in seq_along(ctl_epitope_sequences)){
      if (grepl(ctl_epitope_sequences[i],vaccine_seq)){
        ctl_indexes<-c(ctl_indexes,i)
        ctl_seqs<-c(ctl_seqs,ctl_epitope_sequences[i])
      }
      
    }
    
    if(input$c_software_opt==0){
    ctl_alleles <-allelle_ctl_df[input$allele_ctl_tab_re_rows_selected,]
    
    }else{
      ctl_alleles<-allele_ctl_df_con[input$allele_ctl_tab_con_re_rows_selected,]
      
    }
    
    #HTL sequences and alleles
    htl_seqs <-c()
    htl_indexes<-c()
    for (i in seq_along(htl_epitope_sequences)){
      if(grepl(htl_epitope_sequences[i],vaccine_seq)){
        htl_indexes<-c(htl_indexes,i)
        htl_seqs<-c(htl_seqs,htl_epitope_sequences[i])
      }
    }
    
    if(input$h_software_opt==0){
      htl_alleles <-allelle_htl_df[input$allele_htl_tab_re_rows_selected,]
      
    }else{
      htl_alleles<-allele_htl_df_con[input$allele_htl_tab_con_re_rows_selected,]
      
    }
    
    
    
    incProgress(0.4,detail='B cell epitopes')
######## Run bepipred#################
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
    if(input$bepipred_opt==1){
      amino_region <- input$amino_region}
    else{
      amino_region<-input$amino_region_ii
    }
    ###Secondary threshold
    bepi_sec_thres <- input$second_thres_num
    ###Dataframe
    df <- data.frame(Sequence=vaccine_seq)
    ##Function
    res_bepi <- bepi_fnc$Bepipred3(df,path_bepipred,path_res_bepi,bepi_thres,base_neg,amino_region,bepi_sec_thres,bepipred_opt)
    output$bepi_re_analyze<-renderDataTable({
      res_bepi
    },options=list(scrollX=TRUE,pageLength=5))
    
    output$bepi_ui_1 <-renderUI({
      column(width=4,h3('B cell epitopes predicted'),br(),dataTableOutput('bepi_re_analyze'),align='center')
    })
    
    
    incProgress(0.7,detail='CTL epitopes')
############# Run netMHCpan - CTL#########################
    
    setwd(main_path)
    if(input$c_software_opt==0){
    if(is.null(input$c_length)){
      num_peptides_ctl<-9
    }else{
    num_peptides_ctl <- input$c_length}
    #Alleles
    alleles_ctl <- ctl_alleles
    print(alleles_ctl)
    seq <- df$Sequence
    ##Function
    res_c <- netmhc_fnc(seq,num_peptides_ctl,alleles_ctl,allele_list_ctl)
  }else{
      num_peptides<-input$c_length_con
      allele<-ctl_alleles
      cutoff<-input$c_cutoff_val
      seq<-df$Sequence
      res_path<-cons_I_fnc(seq,num_peptides,allel=allele)
      res_c <- consensus_prediction_fnc_v2(res_path,cutoff)
  } 
    output$ctl_re_analyze<-renderDataTable({
      res_c
    },options=list(scrollX=TRUE,pageLength=5))
   output$ctl_ui_1 <-renderUI({
     column(width=4,h3('CTL epitopes predicted'),br(),dataTableOutput('ctl_re_analyze'),align='center')
   })
    
    incProgress(0.9,detail='HTL epitopes')
########### RUn netMHCIIpan - HTL####################################
    setwd(main_path)
    if(input$h_software_opt==0){
      if(is.null(input$h_length)){
        num_peptides<-15
      }
      else{num_peptides <- input$h_length}
    alleles_htl <-htl_alleles
    seq <- df$Sequence
    res_h<- netmhcII_fnc(seq,num_peptides,alleles_htl,allele_list_htl)
    }else{
      num_peptides<-input$h_length_con
      cutoff<-input$h_cutoff_val
      allele<-htl_alleles
      
      seq<-df$Sequence
      res_path<-cons_II_fnc(seq,num_peptides,allel=allele)
      res_h <- consensus3_prediction_fnc_v2(res_path,cutoff)
     }
    output$htl_re_analyze<-renderDataTable({
      res_h
    },options=list(scrollX=TRUE,pageLength=5))
    
    output$htl_ui_1 <- renderUI({
      column(width=4,h3('HTL epitopes predicted'),br(),dataTableOutput('htl_re_analyze'),align='center')
    })
    
    
    ### HERE
    
    output$vaccine_diagram <- renderPlot({
    vaccine_seq <-df$Sequence
    df <- data.frame(
      position = 1:nchar(vaccine_seq),
      amino_acid = strsplit(vaccine_seq, "")[[1]]
    )
    all_start <- c(res_bepi$Start,res_c$Start,res_h$Start)
    all_end<-c(res_bepi$End,res_c$End,res_h$End)
    all_types <- c(rep("B-cell Epitopes", length(res_bepi$Start)),
                    rep("CTL Epitopes", length(res_c$Start)),
                    rep("HTL Epitopes", length(res_h$Start)))
    all_yy <-c(rep(0.02, length(res_bepi$Start)),
               rep(0.03, length(res_c$Start)),
               rep(0.04, length(res_h$Start)))
    all_data <- data.frame(Start = all_start,End=all_end,type = all_types,y_val=all_yy)
    
    # Plot vaccine sequence with horizontal red line and other colored lines
    p<-ggplot(df, aes(x = position)) +
      geom_text(aes(y = 0, label = amino_acid), vjust = -1) +
      geom_segment(data = all_data, aes(x = Start, xend = End, y = y_val, yend = y_val, color = type), linewidth = 1) +
      labs(x = "Position", y = "", title = "Vaccine Epitope Visualization") +
      scale_color_manual(values = c("red", "blue", "green"), 
                         breaks = c("B-cell Epitopes", "CTL Epitopes", "HTL Epitopes"),
                         labels = c("B-cell Epitopes", "CTL Epitopes", "HTL Epitopes")) +
      theme_minimal() +
      theme(axis.text.y = element_blank(),
            axis.title = element_text(size = 14, face = "bold"),
            plot.title = element_text(size = 16, face = "bold"))
    
    p + coord_cartesian(ylim = c(-0.1, 0.1))
    
    })
    
    
    
    
    
    ###################### TABLES COMPARISON HERE
    selectedRow <- input$final_results_vaccine_rows_selected
    
    # Get the selected row from table 1
    selectedData <- filtered_vaccine()[selectedRow, ]
    print(selectedData)
    b_epitope_list <-res_bepi$Sequences
    print(b_epitope_list)
    c_epitope_list <-res_c$Sequences
    h_epitope_list <-res_h$Sequences
    compared_df_b  <- Epitope_comparison(selectedData,b_epitope_list,'B-cell')
    compared_df_c  <- Epitope_comparison(selectedData,c_epitope_list,'CTL')
    compared_df_h  <- Epitope_comparison(selectedData,h_epitope_list,'HTL')
    
    print(compared_df_c)
    print(compared_df_h)
    #Epitope_comparison <- function(data,new_epitope_list,type_of_epitope)
    
    output$b_cell_compared_table <-DT::renderDataTable({
      compared_df_b
    },options=list(scrollX=TRUE))
    
    output$ctl_compared_table <-DT::renderDataTable({
      compared_df_c
      
    },options=list(scrollX=TRUE))
    output$htl_compared_table <-DT::renderDataTable({
      compared_df_h
      
    },options=list(scrollX=TRUE))
    
    
    output$bepi_ui_2 <-renderUI({ column(width=4,h3('B epitopes'),dataTableOutput('b_cell_compared_table'),align='center')})
    output$ctl_ui_2<-renderUI({column(width=4,h3('CTL epitopes'),dataTableOutput('ctl_compared_table'),align='center')
      })
    output$htl_ui_2<-renderUI({column(width=4,h3('HTL epitopes'),dataTableOutput('htl_compared_table'),align='center')})
    
    
    
    
##########################################
    })
  })
  
  ############################################################################ 4.BLAST #############################################
  observe({
    if(input$blast_selection){
      shinyjs::show('blast_btn')
      shinyjs::show('blast_database')
      shinyjs::show('taxid_opt')
    }
    else{
      shinyjs::hide('blast_btn')
      shinyjs::hide('blast_database')
      shinyjs::hide('taxid_blast_opt')
      shinyjs::hide('taxid_opt')
    }
  })
  

  
  
  
  observeEvent(input$blast_btn,{
    withProgress(message="BLAST",value=0.3,{
    vaccine_seq <- sel_vaccine_sequence()
    seq <- vaccine_seq$Sequence
    if(!is.null(input$blast_database)){
    database <-input$blast_database}
    else{
      database <- 0
    }
    #if(input$taxid_opt==0){
      #taxid <- 1
    #}else{
    #  taxid <-input$taxid_selection
  #  }
    res_table <-blast_function(seq,database)
    
    incProgress(0.7)
    output$blast_table<-DT::renderDataTable({
      res_table
    },options=list(scrollX=TRUE,pageLength=5))
    
    })
  })
  

  
  
  ############################################################################ 5. Netchop ##############################################
  ####################hide elements
  observe({
    if(input$netchop_selection){
      shinyjs::show('netchop_btn')
      shinyjs::show('netchop_thres')
      shinyjs::show('netchop_method')
      shinyjs::show('netchop_res')
      shinyjs::show('vaccine_diagram_chop')}
    
      
    else{
      shinyjs::hide('netchop_btn')
      shinyjs::hide('netchop_thres')
      shinyjs::hide('netchop_method')
      shinyjs::hide('netchop_res')
      shinyjs::hide('vaccine_diagram_chop')
    }
    
  })
  #################### Action button 
  observeEvent(input$netchop_btn,{
    withProgress(message="Netchop",value=0.5,{
    method <- input$netchop_method
    threshold <-input$netchop_thres
    vaccine_seq <-sel_vaccine_sequence()
    
    res <- netchop3_fnc(main_path,vaccine_seq,method,threshold)
    output$netchop_res <- renderDataTable({
      data.frame(res)
    })
    selectedRow <- input$final_results_vaccine_rows_selected
    
     # Get the selected row from table 1
    selectedData <- filtered_vaccine()[selectedRow, ]

    ###vaccine-seq a dataframe with B-cell start,end
    res_bepi <-find_epitope_on_vaccine(vaccine_seq,'B_cell',selectedData)
    #print(res_bepi)
    res_c <- find_epitope_on_vaccine(vaccine_seq,'CTL',selectedData)
    #print(res_c)
    res_h <- find_epitope_on_vaccine(vaccine_seq,'HTL',selectedData)
    #print(res_h)
    
    
    output$vaccine_diagram_chop <- renderPlot({
      vaccine_seq <-vaccine_seq$Sequence
      df <- data.frame(
        position = 1:nchar(vaccine_seq),
        amino_acid = strsplit(vaccine_seq, "")[[1]]
      )
      all_start <- c(res_bepi$Start,res_c$Start,res_h$Start)
      all_end<-c(res_bepi$End,res_c$End,res_h$End)
      all_types <- c(rep("B-cell Epitopes", length(res_bepi$Start)),
                     rep("CTL Epitopes", length(res_c$Start)),
                     rep("HTL Epitopes", length(res_h$Start)))
      all_yy <-c(rep(0.02, length(res_bepi$Start)),
                 rep(0.02, length(res_c$Start)),
                 rep(0.02, length(res_h$Start)))
      all_data <- data.frame(Start = all_start,End=all_end,type = all_types,y_val=all_yy)
      
      res_chop <-data.frame(res)
      res_chop$position <-as.numeric(res_chop$pos)
      res_chop$Cleavages <-c(rep('Proteasome cleavages',length(res_chop$position)))
      
      p<-ggplot(df, aes(x = position)) +
        geom_text(aes(y = 0, label = amino_acid), vjust = -1) +
        
        
        
        geom_segment(data=res_chop,aes(x = position, xend = position,y = 0.01,yend = 0.03,color=Cleavages),linewidth=1)+
        #scale_color_manual(values=c("red"),breaks=c("Proteasome cleavages"),labels=c('Proteasome cleavages'))+
        geom_segment(data = all_data, aes(x = Start, xend = End, y = y_val, yend = y_val, color = type), linewidth = 1) +
        labs(x = "Position", y = "", title = "Proteasome Cleavages on multiepitope Vaccine Visualization") +
        scale_color_manual(values = c("yellow", "blue", "green","red"), 
                           breaks = c("B-cell Epitopes", "CTL Epitopes", "HTL Epitopes","Proteasome cleavages"),
                           labels = c("B-cell Epitopes", "CTL Epitopes", "HTL Epitopes","Proteasome cleavages")) +
        theme_minimal() +
        theme(axis.text.y = element_blank(),
              axis.title = element_text(size = 14, face = "bold"),
              plot.title = element_text(size = 16, face = "bold"))
      
      p + coord_cartesian(ylim = c(-0.1, 0.1))
      
    })
    
  
    })
    })
  
  ############################################################################# 6. ui---- ctl,htl params #########################
  observe({
    if(input$c_software_opt==0){
      output$ctl_params_ui <-renderUI({
        box(width=NULL,status='primary',title='netmhcpan',column(12,column(6,align='center',
                         numericInput("c_length",h3("Length of CTL epitopes"),value= 9,min=8,max=12),
                         actionButton("ctl_allelle_btn", "Select!",class="buttons"),br(),
                         br(),p("Choose the MHC class I alleles of interest for peptide binding prediction by clicking one by one"),
                         h3("Select Allelles"),dataTableOutput("allelle_ctl_tab")),
                        conditionalPanel(condition="!is.null(input.ctl_allelle_btn)",column(6,br(),br(),br(),br(),br(),br(),
                      actionButton("ctl_allelle_clear","Clear",class="buttons"),br(),br(),
                      div(style = "display: flex; justify-content: center;",dataTableOutput("ctl_textarea"))))
                      ),
                    column(width=12,
        br(),br(),br(),br(),HTML("<h4> Reference </h4>"),br(),h5("Birkir Reynisson, Bruno Alvarez, Sinu Paul, Bjoern Peters, Morten Nielsen, NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data, Nucleic Acids Research, Volume 48, Issue W1, 02 July 2020, Pages W449–W454, https://doi.org/10.1093/nar/gkaa379")
        ))
      })
    }
    else{
      output$ctl_params_ui <-renderUI({
        box(width=NULL,status='primary',title='Consensus TepiTool',column(width=12,column(6,align='center',
            numericInput("c_length_con",h3("Length of CTL epitopes"),value=9,min=8,max=12),br(),
            
            numericInput("c_cutoff_val",h3("Select cutoff value for percentile ranking of Consensus prediction"),value=3.0,min=0.1,max=30.0,step=0.1),
            p('Cutoff values range 0.1-30, smaller the threshold, stronger the binders'),br(),br(),
            actionButton("ctl_allele_btn_con","Select!",class="buttons"),br(),br(),
            p("Choose the MHC class I alleles of interest for peptide binding prediction by clicking one by one"),br(),br(),
            h3("Select Alleles"),dataTableOutput("allele_ctl_tab_con")),
            conditionalPanel(condition='is.null(input.ctl_allele_btn_con)',column(6,br(),br(),br(),br(),br(),br(),
            actionButton("ctl_allele_clear_con","Clear",class="buttons"),br(),br(),
            div(style = "display: flex; justify-content: center;",dataTableOutput("ctl_textarea_con"))))
          ),
          column(width=12,
        br(),br(),br(),br(),HTML("<h4> Reference </h4>"),br(),h5("Birkir Reynisson, Bruno Alvarez, Sinu Paul, Bjoern Peters, Morten Nielsen, NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data, Nucleic Acids Research, Volume 48, Issue W1, 02 July 2020, Pages W449–W454, https://doi.org/10.1093/nar/gkaa379")
        )
        )
      })
    }
    
    ###############HERE NEW PLOT
      
    
  })
  
  ######### HTL ui
  observe({
    if(input$h_software_opt==0){
      output$htl_params_ui <-renderUI({
        box(width=NULL,status='primary',title='netMHCIIpan',column(12,column(6,align='center',
                                                                           numericInput("h_length",h3("Length of HTL epitopes"),value= 15,min=12,max=30),
                                                                           actionButton("htl_allelle_btn", "Select!",class="buttons"),br(),
                                                                           br(),p("Choose the MHC class II alleles of interest for peptide binding prediction by clicking one by one"),
                                                                           h3("Select Allelles"),dataTableOutput("allelle_htl_tab")),
                                                                 conditionalPanel(condition="!is.null(input.htl_allelle_btn)",column(6,br(),br(),br(),br(),br(),br(),
                                                                                                                                     actionButton("htl_allelle_clear","Clear",class="buttons"),br(),br(),
                                                                                                                                     div(style = "display: flex; justify-content: center;",dataTableOutput("htl_textarea"))))
        ),
        column(width=12,
               br(),br(),br(),br(),HTML("<h4> Reference </h4>"),br(),h5("Birkir Reynisson, Bruno Alvarez, Sinu Paul, Bjoern Peters, Morten Nielsen, NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data, Nucleic Acids Research, Volume 48, Issue W1, 02 July 2020, Pages W449–W454, https://doi.org/10.1093/nar/gkaa379")
        ))
      })
    }
    else{
      output$htl_params_ui <-renderUI({
        box(width=NULL,status='primary',title='Consensus TepiTool',column(width=12,column(6,align='center',
                                                                                          numericInput("h_length_con",h3("Length of HTL epitopes"),value=15,min=12,max=30),
                                                                                          br(),
                                                                                          numericInput("h_cutoff_val",h3("Select cutoff value for percentile ranking of Consensus prediction"),value=3.0,min=0.1,max=30.0,step=0.1),
                                                                                          p('Cutoff values range 0.1-30, smaller the threshold, stronger the binders'),br(),br(),
                                                                                          
                                                                                          actionButton("htl_allele_btn_con","Select!",class="buttons"),br(),br(),
                                                                                          p("Choose the MHC class II alleles of interest for peptide binding prediction by clicking one by one"),br(),
                                                                                          h3("Select Alleles"),dataTableOutput("allele_htl_tab_con")),
                                                                          conditionalPanel(condition='is.null(input.htl_allele_btn_con)',column(6,br(),br(),br(),br(),br(),br(),
                                                                                                                                                actionButton("htl_allele_clear_con","Clear",class="buttons"),br(),br(),
                                                                                                                                                div(style = "display: flex; justify-content: center;",dataTableOutput("htl_textarea_con"))))
        ),
        column(width=12,
               br(),br(),br(),br(),HTML("<h4> Reference </h4>"),br(),h5("Birkir Reynisson, Bruno Alvarez, Sinu Paul, Bjoern Peters, Morten Nielsen, NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data, Nucleic Acids Research, Volume 48, Issue W1, 02 July 2020, Pages W449–W454, https://doi.org/10.1093/nar/gkaa379")
        )
        )
      })
    }
  })
  
  
  
  ######allele_ctl_tab_con
  output$allele_ctl_tab_con <- DT::renderDataTable({
    data <-allele_ctl_df_con
    allelle_f<-data
    allelle_f
  },options=list(pageLength=20,scrollX=TRUE,selection="multiple"))
  ######allele_htl_tab_con
  output$allele_htl_tab_con <- DT::renderDataTable({
    data <-allele_htl_df_con
    allelle_f<-data
    allelle_f
  },options=list(pageLength=20,scrollX=TRUE,selection="multiple"))
  ######buttons
  ctl_selectedRows_con<-reactiveVal(NULL)
  observeEvent(input$ctl_allele_btn_con, {
    ctl_selectedRows_con(input$allele_ctl_tab_con_rows_selected)
    rows <- ctl_selectedRows_con()
    data <- allele_ctl_df_con
    allelle_f <- data
    
    data <- allelle_f[rows,]
    df <- data.frame(data)
    colnames(df) <- c('Selected Allelles')
    output$ctl_textarea_con <- DT::renderDataTable({
      df
    })
    
  })
  observeEvent(input$ctl_allele_clear_con,{
    table_proxy <- dataTableProxy("allele_ctl_tab_con")
    table_proxy %>% selectRows(selected = NULL)
    
    output$ctl_textarea_con <- DT::renderDataTable({
      NULL
    })
    
  })
  
  ######htl-buttons
  htl_selectedRows_con<-reactiveVal(NULL)
  observeEvent(input$htl_allele_btn_con, {
    htl_selectedRows_con(input$allele_htl_tab_con_rows_selected)
    rows <- htl_selectedRows_con()
    data <- allele_htl_df_con
    allelle_f <- data
    
    data <- allelle_f[rows,]
    df <- data.frame(data)
    colnames(df) <- c('Selected Allelles')
    output$htl_textarea_con <- DT::renderDataTable({
      df
    })
    
  })
  observeEvent(input$htl_allele_clear_con,{
    table_proxy <- dataTableProxy("allele_htl_tab_con")
    table_proxy %>% selectRows(selected = NULL)
    
    output$htl_textarea_con <- DT::renderDataTable({
      NULL
    })
    
  })
  
  
  ################################### Future steps
  observe({
    if(input$Next_steps_opt){
  
      shinyjs::show('next_steps_text_cimmsim')
      shinyjs::show('next_steps_text_3dstab')}
    
    
    else{
      
      shinyjs::hide('next_steps_text_cimmsim')
      shinyjs::hide('next_steps_text_3dstab')
    
    }
    
  })
  output$next_steps_text_cimmsim<-renderText({
    text1 <- "1. Evaluation of multiepitope sequence for Immune Simulation with C-ImmSIm (https://kraken.iac.rm.cnr.it/C-IMMSIM/)"
    text2<-'2. Evaluation of sequence stability and 3D interactions through molecular docking and dynamics '
    text1
  })
  output$next_steps_text_3dstab<-renderText({
    text1 <- "1. Evaluation of multiepitope sequence for Immune Simulation with C-ImmSIm (https://kraken.iac.rm.cnr.it/C-IMMSIM/)"
    text2<-'2. Evaluation of sequence stability and 3D interactions through molecular docking and dynamics '
    text2
  })
  
  
  
  ################IMAGES  MANUAL ###################################################
  output$fig1 <- renderImage({
    return(list(src = "Manual_software.files/Figure_1.png",
                contentType="image/png",alt="Figure 1"))
  }, deleteFile = FALSE)
  
  output$img_001 <- renderImage({
    return(list(src = "Manual_software.files/img_01.png",
         contentType="image/png",alt="Figure 1"))
  }, deleteFile = FALSE)
  output$img_002 <- renderImage({
    return(list(src = "Manual_software.files/img_02.png",
                contentType="image/png",alt="Figure 2"))
  }, deleteFile = FALSE)
  output$img_003 <- renderImage({
    return(list(src = "Manual_software.files/img_03.png",
                contentType="image/png",alt="Figure 3"))
  }, deleteFile = FALSE)
  
  output$img_004 <- renderImage({
    return(list(src = "Manual_software.files/img_04.png",
                contentType="image/png",alt="Figure 4"))
  }, deleteFile = FALSE)
  
  output$img_005 <- renderImage({
    return(list(src = "Manual_software.files/img_05.png",
                contentType="image/png",alt="Figure 5"))
  }, deleteFile = FALSE)
  
  output$img_006 <- renderImage({
    return(list(src = "Manual_software.files/img_06.png",
                contentType="image/png",alt="Figure 5"))
  }, deleteFile = FALSE)
  
  output$img_007 <- renderImage({
    return(list(src = "Manual_software.files/img_07.png",
                contentType="image/png",alt="Figure 5"))
  }, deleteFile = FALSE)
  
  output$img_008 <- renderImage({
    return(list(src = "Manual_software.files/img_08.png",
                contentType="image/png",alt="Figure 5"))
  }, deleteFile = FALSE)
  
  output$img_009 <- renderImage({
    return(list(src = "Manual_software.files/img_09.png",
                contentType="image/png",alt="Figure 5"))
  }, deleteFile = FALSE)
  
  output$img_010 <- renderImage({
    return(list(src = "Manual_software.files/img_10.png",
                contentType="image/png",alt="Figure 5"))
  }, deleteFile = FALSE)
  
  output$img_011 <- renderImage({
    return(list(src = "Manual_software.files/img_11.png",
                contentType="image/png",alt="Figure 5"))
  }, deleteFile = FALSE)
  
  output$img_012 <- renderImage({
    return(list(src = "Manual_software.files/img_12.png",
                contentType="image/png",alt="Figure 5"))
  }, deleteFile = FALSE)
  
  output$img_013 <- renderImage({
    return(list(src = "Manual_software.files/img_13.png",
                contentType="image/png",alt="Figure 5"))
  }, deleteFile = FALSE)
  
  output$img_014 <- renderImage({
    return(list(src = "Manual_software.files/img_14.png",
                contentType="image/png",alt="Figure 5"))
  }, deleteFile = FALSE)
  
  output$img_015 <- renderImage({
    return(list(src = "Manual_software.files/img_15.png",
                contentType="image/png",alt="Figure 5"))
  }, deleteFile = FALSE)
  
  output$img_016 <- renderImage({
    return(list(src = "Manual_software.files/img_16.png",
                contentType="image/png",alt="Figure 5"))
  }, deleteFile = FALSE)
  
  output$img_017 <- renderImage({
    return(list(src = "Manual_software.files/img_17.png",
                contentType="image/png",alt="Figure 5"))
  }, deleteFile = FALSE)
  
  output$img_018 <- renderImage({
    return(list(src = "Manual_software.files/img_18.png",
                contentType="image/png",alt="Figure 5"))
  }, deleteFile = FALSE)
  
  output$img_019 <- renderImage({
    return(list(src = "Manual_software.files/img_19.png",
                contentType="image/png",alt="Figure 5"))
  }, deleteFile = FALSE)
  
  output$img_epitopes<-renderImage({
    return(list(src="./figure_epitopes.png",contentType="image/png",alt="Figure epitopes Los Alamos"))
  },deleteFile=FALSE)
  
  output$table_epitopes<-renderDataTable({
    data <-read.csv('epitopes.csv')
    data.frame(data)
  })
  
  output$table_epitopes_properties<-renderDataTable({
    data<-read.csv('epitope_results.csv')
    data.frame(data)
  })
  #output$img_020 <- renderImage({
   # return(list(src = "Manual_software.files/image039.png",
    #            contentType="image/png",alt="Figure 5",width="100%"))
  #}, deleteFile = FALSE)
  
  
}