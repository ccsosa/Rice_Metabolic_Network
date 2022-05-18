require(dplyr)
require(data.table)
require(readxl)
require(stringr)
require(stringi)
require(xlsx)
require(BBmisc)
require(riceidconverter)
require(plyr)
require(data.table)

################################################################################
dir_o <- "D:/TESIS_PHD/CHAPTER2"
dir_2 <-
  "D:/PROGRAMAS/Dropbox/shared/Metabolic_network_manual_curation/Aluminio/GO_ANNOTATION"
dir_3 <-
  "D:/PROGRAMAS/Dropbox/shared/Metabolic_network_manual_curation/Aluminio"
dir_in <- "D:/PROGRAMAS/Dropbox/shared/Metabolic_network_manual_curation/Aluminio/INPUTS"
################################################################################

iOS2164 <-
  as.data.frame(readxl::read_xls(paste0(dir_o, "/", "iOS2164_GSM.xls"), "S2A"))
iOS_name <- colnames(iOS2164)


#loading OryzaCyc file
ORYZACYC <-
  as.data.frame(fread(paste0(dir_in,"/","All_reactions_of_O._sativa_Japonica_Group_2022_05_10_CC.tsv")))
#  as.data.frame(readxl::read_xlsx(paste0(dir_o, "/", "ORYZACYC_2022_03_04.xlsx")))
################################################################################
#genes to query
genes_to_explore <-
  as.data.frame(readxl::read_xls(paste0(dir_2, "/", "per_gen_summary_GO.xls"), sheet = "SOURCE"))
genes_to_explore$ORYZACYC <- NA
genes_to_explore_NPM <- genes_to_explore[which(genes_to_explore$AL_IOS==FALSE),]

################################################################################
#In OryzaCyc
#oryacyc_genes <- unlist(unique(str_split(ORYZACYC$`Genes of a reaction`," // ")))
ORYZACYC_REP <- lapply(1:nrow(ORYZACYC),function(i){
x_i <- ORYZACYC[i,]
x_genes_i <- trimws(unlist(str_split(x_i$`Genes of a reaction`," // ")))
x_i_sub <- data.frame(matrix(ncol=ncol(x_i)+1,nrow = length(x_genes_i)))
x_i_sub[,1] <- i
x_i_sub[,2] <- x_i$Reaction
x_i_sub[,3] <- x_i$Reaction_id
x_i_sub[,4] <- x_i$`EC-Number`
x_i_sub[,5] <- x_genes_i
x_i_sub[,6] <- x_i$`Gibbs-0`
x_i_sub[,7] <- x_i$`In-Pathway`
x_i_sub[,8] <- x_i$`In-Pathway_id`
x_i_sub[,9] <- x_i$`Common-Name`
colnames(x_i_sub) <- c("i",colnames(x_i))
return(x_i_sub)
})

ORYZACYC_REP <- do.call(rbind,ORYZACYC_REP)
#oryzacyc rxns with information for no mapped genes
ORYZACYC_REP2 <- 
  ORYZACYC_REP[ ORYZACYC_REP$`Genes of a reaction` %in% genes_to_explore_NPM$Ensembl_gene_id,]
#ORYZACYC_REP[genes_to_explore_NPM$Ensembl_gene_id %in% ORYZACYC_REP$`Genes of a reaction`,]
################################################################################

un_or_rxn <- unique(ORYZACYC_REP2$Reaction)
un_or_rep2_genes <- unique(na.omit(ORYZACYC_REP2$`Genes of a reaction`))
un_or_rep2_genes <- un_or_rep2_genes[which(un_or_rep2_genes!="")]
################################################################################
genes_to_explore$ORYZACYC[genes_to_explore$Ensembl_gene_id %in% un_or_rep2_genes] <- TRUE
genes_to_explore$ORYZACYC[!genes_to_explore$Ensembl_gene_id %in% un_or_rep2_genes] <- FALSE
length(genes_to_explore$ORYZACYC[which(genes_to_explore$ORYZACYC==TRUE)])
################################################################################

write.csv(genes_to_explore,
          paste0(dir_3, "/MAPPINGS/", "genes_to_explore.csv"))


################################################################################
################################################################################
################################################################################


x_sub_mapping <- lapply(1:length(un_or_rxn),function(i){
  
  #message(i)
  sub_NMAP <- ORYZACYC_REP2[which(ORYZACYC_REP2$Reaction==un_or_rxn[[i]]),]
  x_sub <- data.frame(matrix(ncol = length(iOS_name),nrow=1))
  colnames(x_sub) <- iOS_name
  x_sub$Abbreviation <- NA
  x_sub$Name <- sub_NMAP$`Common-Name`[[1]]
  x_sub$`Reaction Equation` <- sub_NMAP$Reaction[[1]]
  x_sub$`GPR (MSU)` <- NA
  x_sub$`GPR (RAP)` <- paste0("(",paste(sub_NMAP$`Genes of a reaction`,collapse = " or "),")")
  x_sub$status <- NA
  x_sub$DeltaG <- sub_NMAP$`Gibbs-0`[[1]]
  x_sub$RXN_ID <- sub_NMAP$Reaction_id[[1]]
  x_sub$Subsystem <- sub_NMAP$`In-Pathway`[[1]]
  x_sub$`EC Number` <- sub_NMAP$`EC-Number`[[1]]
  return(x_sub)
})



x_sub_mapping <- do.call(rbind,x_sub_mapping)

write.xlsx(x_sub_mapping,
           paste0(dir_3, "/MAPPINGS/", "ORYZACYC_NO_MAPPED.xlsx"),
           showNA = F)

