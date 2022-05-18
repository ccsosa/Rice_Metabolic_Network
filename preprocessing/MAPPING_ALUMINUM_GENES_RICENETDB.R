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
################################################################################

#loading ricenetdb files
ricenetdb_data <-
  as.data.frame(fread(paste0(dir_in,"/","ricenetdb_2022_05_10.tsv")))

ricenetdb_data_rxn <-
  as.data.frame(fread(paste0(dir_in,"/","ricenetdb_2022_05_10_rxn.tsv")))
#  as.data.frame(readxl::read_xlsx(paste0(dir_o, "/", "ORYZACYC_2022_03_04.xlsx")))
################################################################################
#genes to query
genes_to_explore <-
  as.data.frame(readxl::read_xls(paste0(dir_2, "/", "per_gen_summary_GO.xls"), sheet = "SOURCE"))
genes_to_explore_NMAPPED <- genes_to_explore[which(genes_to_explore$AL_IOS==FALSE & genes_to_explore$ORYZACYC==FALSE),]

genes_to_explore$RICENETDB <- NA

################################################################################
#In OryzaCyc
#oryacyc_genes <- unlist(unique(str_split(ORYZACYC$`Genes of a reaction`," // ")))
ricenet_rxn_REP <- lapply(1:nrow(ricenetdb_data_rxn),function(i){

  x_i <- ricenetdb_data_rxn[i,]
  x_genes_i <- trimws(unlist(str_split(x_i$protein,";")))
  x_i_sub <- data.frame(matrix(ncol=6,nrow = length(x_genes_i)))
  x_i_sub[,1] <- i
  x_i_sub[,2] <- x_i$id
  x_i_sub[,3] <- x_i$formula
  x_i_sub[,4] <- x_i$enzyme
  x_i_sub[,5] <- x_genes_i
  x_i_sub[,6] <- x_i$crosslink
  colnames(x_i_sub) <- c("i","rxn_id","formula","enzyme","gene","crosslink")
  return(x_i_sub)
})

ricenet_rxn_REP <- do.call(rbind,ricenet_rxn_REP)
ricenet_rxn_REP$enzyme <- gsub("[||]","/",ricenet_rxn_REP$enzyme)
#oryzacyc rxns with information for no mapped genes
ricenet_rxn_REP2 <- 
  ricenet_rxn_REP[ ricenet_rxn_REP$gene %in% genes_to_explore_NMAPPED$MSU,]
#ORYZACYC_REP[genes_to_explore_NPM$Ensembl_gene_id %in% ORYZACYC_REP$`Genes of a reaction`,]
################################################################################
#check genes id match in rxns
identical(sort(unique(ricenetdb_data$id_fix)),#,unique(ricenet_rxn_REP2$gene))
sort(unique(ricenet_rxn_REP2$gene)))
################################################################################
un_or_rxn <- unique(ricenet_rxn_REP2$rxn_id)

un_or_rep2_genes <- unique(na.omit(ricenet_rxn_REP2$gene))
un_or_rep2_genes <- un_or_rep2_genes[which(un_or_rep2_genes!="")]
################################################################################
genes_to_explore$RICENETDB[genes_to_explore$MSU %in% un_or_rep2_genes] <- TRUE
genes_to_explore$RICENETDB[!genes_to_explore$MSU %in% un_or_rep2_genes] <- FALSE
length(genes_to_explore$RICENETDB[which(genes_to_explore$RICENETDB==TRUE)])
################################################################################

write.csv(genes_to_explore,
          paste0(dir_3, "/MAPPINGS/", "genes_to_explore_rice_net.csv"))


################################################################################
################################################################################
################################################################################


x_sub_mapping <- lapply(1:length(un_or_rxn),function(i){
  #message(i)
  sub_NMAP <- ricenet_rxn_REP2[which(ricenet_rxn_REP2$rxn_id==un_or_rxn[[i]]),]
  x_sub <- data.frame(matrix(ncol = length(iOS_name),nrow=1))
  colnames(x_sub) <- iOS_name
  x_sub$Abbreviation <- NA
  x_sub$Name <- NA
  x_sub$`Reaction Equation` <- sub_NMAP$formula[[1]]
  x_sub$`GPR (MSU)` <- paste0("(",paste(sub_NMAP$gene,collapse = " or "),")")
  x_sub$`GPR (RAP)` <- NA
  x_sub$status <- NA
  x_sub$DeltaG <- NA
  x_sub$RXN_ID <- sub_NMAP$rxn_id[[1]]
  x_sub$Subsystem <- sub_NMAP$crosslink[[1]]
  x_sub$`EC Number` <- sub_NMAP$enzyme[[1]]
  return(x_sub)
})



x_sub_mapping <- do.call(rbind,x_sub_mapping)

write.xlsx(x_sub_mapping,
           paste0(dir_3, "/MAPPINGS/", "RICENETDB_NO_MAPPED.xlsx"),
           showNA = F)

