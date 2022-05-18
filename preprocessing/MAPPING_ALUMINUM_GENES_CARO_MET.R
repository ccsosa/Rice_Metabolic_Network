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
dir_file <- "D:/PROGRAMAS/Dropbox/shared/Metabolic_network_manual_curation/Caro_metabolism"
################################################################################


iOS2164 <-
  as.data.frame(readxl::read_xls(paste0(dir_o, "/", "iOS2164_GSM.xls"), "S2A"))
iOS_name <- colnames(iOS2164)
################################################################################

genes_to_explore <-
  as.data.frame(readxl::read_xls(paste0(dir_2, "/", "per_gen_summary_GO.xls"), sheet = "SOURCE"))
genes_to_explore$CARO_FILE <- NA

CARO_FILE <- genes_to_explore <-
  as.data.frame(readxl::read_xls(paste0(dir_2, "/", "per_gen_summary_GO.xls"), sheet = "Hoja1"))

NM_C <- genes_to_explore$CARO_RAP[which(genes_to_explore$`CHECK_MANUALLY?`=="CHECK BY HAND")]
  
  
NM_C <-  
  trimws(unique(na.omit(NM_C)))

CARO_FILE <- CARO_FILE[which(CARO_FILE$rap...1 %in% NM_C),]

unique_rxns <- unique(trimws(unlist(strsplit(CARO_FILE$`Reaction Formula`,"//"))))

################################################################################

ricenet_rxn_REP <- lapply(1:nrow(CARO_FILE),function(i){
#i <- 1
  x_i <- CARO_FILE[i,]
  x_genes_i <- trimws(unlist(str_split(x_i$`Reaction Formula`,"//")))
  #x_genes_i_id <- trimws(unlist(str_split(x_i$Reaction,"//")))
  gibbs_i_id <- trimws(unlist(str_split(x_i$`Gibbs-0`,"//")))
  x_i_sub <- data.frame(matrix(ncol=7,nrow = length(x_genes_i)))
  x_i_sub[,1] <- i
  x_i_sub[,2] <- x_genes_i
  x_i_sub[,3] <- x_i$Reaction#_genes_i_id
  x_i_sub[,4] <- x_i$`EC-Number`
  x_i_sub[,5] <- gibbs_i_id#x_i$`Gibbs-0`
  x_i_sub[,6] <- x_i$`Metabolic Subsystem`
  x_i_sub[,7] <- x_i$rap...1
  colnames(x_i_sub) <- c("i","rxn","rxn_id","EC","gibbs","subsystem","rap_or")
  return(x_i_sub)
})

ricenet_rxn_REP <- do.call(rbind,ricenet_rxn_REP)
################################################################################
NM_C <- sub("G","g",NM_C)
NM_C <- sub("OS","Os",NM_C)
genes_to_explore$CARO_FILE[genes_to_explore$Ensembl_gene_id %in% NM_C] <- TRUE
genes_to_explore$CARO_FILE[!genes_to_explore$Ensembl_gene_id %in% NM_C] <- FALSE
length(genes_to_explore$CARO_FILE[which(genes_to_explore$CARO_FILE==TRUE)])

write.csv(genes_to_explore,
          paste0(dir_3, "/MAPPINGS/", "genes_to_explore_Caro.csv"))
# ricenet_rxn_REP2 <- 
#   ricenet_rxn_REP[ ricenet_rxn_REP$gene %in% NM_C,]




x_sub_mapping <- lapply(1:length(unique_rxns),function(i){
  #message(i)
  sub_NMAP <- ricenet_rxn_REP[which(ricenet_rxn_REP$rxn==unique_rxns[[i]]),]
  x_sub <- data.frame(matrix(ncol = length(iOS_name),nrow=1))
  colnames(x_sub) <- iOS_name
  x_sub$Abbreviation <- NA
  x_sub$Name <- NA
  x_sub$`Reaction Equation` <- sub_NMAP$rxn[[1]]
  x_sub$`GPR (MSU)` <- paste0("(",paste(sub_NMAP$rap_or,collapse = " or "),")")
  x_sub$`GPR (RAP)` <- NA
  x_sub$status <- NA
  x_sub$DeltaG <- sub_NMAP$gibbs
  x_sub$RXN_ID <- sub_NMAP$rxn_id[[1]]
  x_sub$Subsystem <- sub_NMAP$subsystem[[1]]
  x_sub$`EC Number` <- sub_NMAP$EC[[1]]
  return(x_sub)
})



x_sub_mapping <- do.call(rbind,x_sub_mapping)

write.xlsx(x_sub_mapping,
           paste0(dir_3, "/MAPPINGS/", "CARO_NO_MAPPED.xlsx"),
           showNA = F)


