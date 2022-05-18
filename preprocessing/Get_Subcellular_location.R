#install.packages("UniprotR")
require(UniprotR);require(biomaRt);require(riceidconverter)
################################################################################
# ensembl_pl = biomaRt::useMart(biomart="plants_mart",host="plants.ensembl.org")
# # avail_datasets <- listDatasets(ensembl2)
# ensembl_os = biomaRt::useDataset("osativa_eg_gene", mart = ensembl_pl)

################################################################################
dir_o <- "D:/TESIS_PHD/CHAPTER2"
dir_2 <-
  "D:/PROGRAMAS/Dropbox/shared/Metabolic_network_manual_curation/Aluminio/GO_ANNOTATION"
dir_3 <-
  "D:/PROGRAMAS/Dropbox/shared/Metabolic_network_manual_curation/Aluminio"
dir_in <- "D:/PROGRAMAS/Dropbox/shared/Metabolic_network_manual_curation/Aluminio/INPUTS"


################################################################################
genes_to_explore <-
  as.data.frame(readxl::read_xls(paste0(dir_2, "/", "per_gen_summary_GO.xls"), sheet = "SOURCE"))

uniprot_file <-
  as.data.frame(readxl::read_xlsx(paste0(dir_2, "/INPUTS/", "uniprot-yourlist_M202205134ABAA9BC7178C81CEBC9459510EDDEA34D3D641.xlsx"), sheet = "Sheet0"))

xx <- GetSubcellular_location(uniprot_file$Entry)
#xx2 <- GetSubcellular_location(uniprot_file$Entry[1001:19100])




xx_uni <- lapply(1:nrow(genes_to_explore),function(i){
  #i <- 1
  x_uni_i <- uniprot_file[which(uniprot_file$`yourlist:M202205134ABAA9BC7178C81CEBC9459510EDDEA34D3D641`==genes_to_explore$Ensembl_gene_id[[i]]),]
  xx_uni_i_uni <- xx[row.names(xx) %in% x_uni_i$Entry,]
  
  xx_uni_i_uni$Subcellular.location..CC.[which(is.na(xx_uni_i_uni$Subcellular.location..CC.))] <- ""
  xx_uni_i_uni$Intramembrane[which(is.na(xx_uni_i_uni$Intramembrane))] <- ""
  xx_uni_i_uni$Topological.domain[which(is.na(xx_uni_i_uni$Topological.domain))] <- ""
  xx_uni_i_uni$Transmembrane[which(is.na(xx_uni_i_uni$Transmembrane))] <- ""


  xx_i <- data.frame(Ensembl_gene_id=genes_to_explore$Ensembl_gene_id[[1]],
                   UNIPROT= paste0(x_uni_i$Entry,collapse = "//"),
                   UNIPROT_STATUS= paste0(x_uni_i$Status,collapse = "//"),
                   Subcellular.location=paste0(xx_uni_i_uni$Subcellular.location..CC.,collapse = "//",recycle0 = T),
                   Intramembrane=paste0(xx_uni_i_uni$Intramembrane,collapse = "//",recycle0 = T),
                   Topological.domain=paste0(xx_uni_i_uni$Topological.domain,collapse = "//",recycle0 = T),
                   Transmembrane= paste0(xx_uni_i_uni$Transmembrane,collapse = "//",recycle0 = T)
                   )
  
  
  
  xx_i$Subcellular.location[which(xx_i$Subcellular.location=="//")] <- NA
  xx_i$Subcellular.location[which(xx_i$Subcellular.location=="////")] <- NA
  xx_i$Subcellular.location[which(xx_i$Subcellular.location=="//////")] <- NA
  xx_i$Subcellular.location[which(xx_i$Subcellular.location=="////////")] <- NA
  xx_i$Subcellular.location[which(xx_i$Subcellular.location=="//////////")] <- NA
  xx_i$Subcellular.location[which(xx_i$Subcellular.location=="///////////")] <- NA
  xx_i$Subcellular.location[which(xx_i$Subcellular.location=="////////////")] <- NA

  xx_i$Intramembrane[which(xx_i$Intramembrane=="//")] <- NA
  xx_i$Intramembrane[which(xx_i$Intramembrane=="////")] <- NA
  xx_i$Intramembrane[which(xx_i$Intramembrane=="//////")] <- NA
  xx_i$Intramembrane[which(xx_i$Intramembrane=="////////")] <- NA
  xx_i$Intramembrane[which(xx_i$Intramembrane=="//////////")] <- NA
  xx_i$Intramembrane[which(xx_i$Intramembrane=="///////////")] <- NA
  xx_i$Intramembrane[which(xx_i$Intramembrane=="////////////")] <- NA
  
  
  xx_i$Topological.domain[which(xx_i$Topological.domain=="//")] <- NA
  xx_i$Topological.domain[which(xx_i$Topological.domain=="////")] <- NA
  xx_i$Topological.domain[which(xx_i$Topological.domain=="//////")] <- NA
  xx_i$Topological.domain[which(xx_i$Topological.domain=="////////")] <- NA
  xx_i$Topological.domain[which(xx_i$Topological.domain=="//////////")] <- NA
  xx_i$Topological.domain[which(xx_i$Topological.domain=="///////////")] <- NA
  xx_i$Topological.domain[which(xx_i$Topological.domain=="////////////")] <- NA
  
  
  xx_i$Transmembrane[which(xx_i$Transmembrane=="//")] <- NA
  xx_i$Transmembrane[which(xx_i$Transmembrane=="////")] <- NA
  xx_i$Transmembrane[which(xx_i$Transmembrane=="//////")] <- NA
  xx_i$Transmembrane[which(xx_i$Transmembrane=="////////")] <- NA
  xx_i$Transmembrane[which(xx_i$Transmembrane=="//////////")] <- NA
  xx_i$Transmembrane[which(xx_i$Transmembrane=="///////////")] <- NA
  xx_i$Transmembrane[which(xx_i$Transmembrane=="////////////")] <- NA
  
  return(xx_i)
  
})

xx_uni <- do.call(rbind, xx_uni)


write.csv(xx_uni,
          paste0(dir_3, "/MAPPINGS/", "genes_to_explore_sub.csv"))
#View(biomaRt::listAttributes(ensembl_os)) #get attributes

################################################################################

# go_ids= biomaRt::getBM(attributes=c("external_gene_name","ensembl_gene_id","uniprot_gn_id"),
#                        filters="ensembl_gene_id",#list(c(en='external_gene_name',en_id="ensembl_gene_id")),
#                        values=genes_to_explore$Ensembl_gene_id, 
#                        mart=ensembl_os)
################################################################################



################################################################################