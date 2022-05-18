options(java.parameters = "-Xmx8192m")
require(data.table);require(dplyr);require(plyr);require(stringr);require(xlsx);require(rBiopaxParser)

dir <- "D:/MET_BIOCYC_DBs/Ricecyc/v_3_3/3.3"
out_dir <- "D:/TESIS_PHD/CHAPTER2/ricecyc"
# enzymes <- fread(paste0(dir,"/data/","enzymes.col"),sep = "\t",skip = 22)
# 
# genes<- fread(paste0(dir,"/data/","genes.col"),sep = "\t",skip = 22)
# 
# pathways <- fread(paste0(dir,"/data/","pathways.col"),sep = "\t",skip = 22,fill = T)
# 
# enzymes <- fread(paste0(dir,"/data/","enzymes.col"),sep = "\t",skip = 22)

#############################################################################################
enzrxn <- fread(paste0(dir,"/data/","reactions.dat"),sep = "\\",skip = 58,header=F,fill=T) #reactions.dat #enzrxns.dat
enzrxn <- strsplit(enzrxn$V1," - ")
#############################################################################################
i_i <- enzrxn %in% "//"
i_i_n <- 1:length(i_i)
i_i_n <- i_i_n[i_i]
i_i_n <-  c(1,i_i_n)


df_total <- lapply(1:(length(i_i_n)-1),function(i){
  
  x <- enzrxn[i_i_n[i]:i_i_n[i+1]]
  xx <- x[!x %in% "//$"]
  #x <- x[-length(x)]
  xx <- t(as.data.frame((do.call(cbind, x))))
  xx[,1] <- trimws(xx[,1])
  xx[,2] <- trimws(xx[,2])
  title <- unique(xx[,1])
  title <- title[!str_detect(title, "/")]
  title <- title[!title %in% "COMMENT"]
  
  df <- data.frame(matrix(ncol=length(title)+1,nrow=1))
  colnames(df) <- c("i",title)
  df[1,1] <- i
  for(j in 2:ncol(df)){
    df[1,j] <- paste(xx[,2][which(xx[,1]==colnames(df)[j])],collapse = "//")
  }
  rm(x,xx,title)
  return(df)
  # setTxtProgressBar(pb, i)
  
  
})


df_total2 <- rbind.fill(df_total)
#  write.xlsx(df_total2,paste0(out_dir,"/","reactions_dat.xlsx"),row.names = F,showNA = F)

write.table(df_total2,paste0(out_dir,"/","reactions_dat.tsv"),row.names = F,na = "",sep = "\t")


#############################################################################################
#############################################################################################
enzrxn <- fread(paste0(dir,"/data/","enzrxns.dat"),sep = "\\",skip = 46,header=F) #reactions.dat #enzrxns.dat
enzrxn <- strsplit(enzrxn$V1," - ")
#############################################################################################
i_i <- enzrxn %in% "//"
i_i_n <- 1:length(i_i)
i_i_n <- i_i_n[i_i]
i_i_n <-  c(1,i_i_n)


df_total <- lapply(1:(length(i_i_n)-1),function(i){
  
  x <- enzrxn[i_i_n[i]:i_i_n[i+1]]
  xx <- x[!x %in% "//$"]
  #x <- x[-length(x)]
  xx <- t(as.data.frame((do.call(cbind, x))))
  xx[,1] <- trimws(xx[,1])
  xx[,2] <- trimws(xx[,2])
  title <- unique(xx[,1])
  title <- title[!str_detect(title, "/")]
  title <- title[!title %in% "COMMENT"]
  
  df <- data.frame(matrix(ncol=length(title)+1,nrow=1))
  colnames(df) <- c("i",title)
  df[1,1] <- i
  for(j in 2:ncol(df)){
    df[1,j] <- paste(xx[,2][which(xx[,1]==colnames(df)[j])],collapse = "//")
  }
  rm(x,xx,title)
  return(df)
  # setTxtProgressBar(pb, i)
  
  
})

df_total2 <- rbind.fill(df_total)

#write.xlsx(df_total2,paste0(out_dir,"/","enzrxns_dat.xlsx"),row.names = F,showNA = F)
write.table(df_total2,paste0(out_dir,"/","enzrxns_dat.tsv"),row.names = F,na = "",sep = "\t")
#############################################################################################
#############################################################################################
#############################################################################################
enzrxn <- fread(paste0(dir,"/data/","proteins.dat"),sep = "\\",skip = 82,header=F,fill=T) #reactions.dat #enzrxns.dat
enzrxn <- strsplit(enzrxn$V1," - ")
#############################################################################################
i_i <- enzrxn %in% "//"
i_i_n <- 1:length(i_i)
i_i_n <- i_i_n[i_i]
i_i_n <-  c(1,i_i_n)


df_total <- lapply(1:(length(i_i_n)-1),function(i){
  
  x <- enzrxn[i_i_n[i]:i_i_n[i+1]]
  xx <- x[!x %in% "//$"]
  #x <- x[-length(x)]
  xx <- t(as.data.frame((do.call(cbind, x))))
  xx[,1] <- trimws(xx[,1])
  xx[,2] <- trimws(xx[,2])
  title <- unique(xx[,1])
  title <- title[!str_detect(title, "/")]
  #title <- title[!title %in% "COMMENT"]
  
  df <- data.frame(matrix(ncol=length(title)+1,nrow=1))
  colnames(df) <- c("i",title)
  df[1,1] <- i
  for(j in 2:ncol(df)){
    df[1,j] <- paste(xx[,2][which(xx[,1]==colnames(df)[j])],collapse = "//")
  }
  rm(x,xx,title)
  return(df)
  # setTxtProgressBar(pb, i)
  
  
})

df_total2 <- rbind.fill(df_total)
#write.xlsx(df_total2,paste0(out_dir,"/","proteins_dat.xlsx"),row.names = F,showNA = F)
write.table(df_total2,paste0(out_dir,"/","proteins_dat.tsv"),row.names = F,na = "",sep = "\t")

#############################################################################################
#############################################################################################

protein <- fread(paste0(out_dir,"/","proteins_dat.tsv"))
protein <- protein[which(protein$CATALYZES!=""),]
#protein[which(!is.na(protein$CATALYZES) |protein$CATALYZES!=""),]

enzrxns_dat <- fread(paste0(out_dir,"/","enzrxns_dat.tsv"))
#reactions_dat <- fread(paste0(out_dir,"/","reactions_dat.tsv"))
reactions_dat <- fread(paste0(dir,"/data/","enzymes.col"),skip = 16)
reactions_dat$enz2 <- gsub(".1-MONOMER","",reactions_dat$`SUBUNIT-COMPOSITION`)
reactions_dat$enz2 <- gsub("[1*]","",reactions_dat$enz2)

reactions_dat_dat <- fread(paste0(out_dir,"/","reactions_dat.tsv"))

#write.table(reactions_dat,paste0(out_dir,"/","enzymes_col.tsv"),row.names = F,na = "",sep = "\t")
  #############################################################################################
xx <- left_join(protein,reactions_dat,by=c("CATALYZES"="UNIQUE-ID"))

reactions_dat

require(data.table)
dir2 <- "D:/TESIS_PHD/CHAPTER2/RICENETDB"
data <- fread(paste0(dir2,"/","batchresult_genes.txt"))
data_filter <- data[which(data$reaction!=""),]
data_rxn <- data.frame(rxn=trimws(unique(unlist(strsplit(data$reaction,",")))))

write.table(data_rxn,paste0(dir2,"/","rxn_ids_ricenetdb.tsv"),row.names = F,na = "",sep = "\t")



#########load rxns to match genes

data_rxns_ricecyc <- fread(paste0(dir2,"/","batchresult_rxns.txt"))
data_rxns_ricecyc$rap <- NA
data_rxns_ricecyc$msu <- NA 
data_rxns_ricecyc$source <- "RICENETDB"

genes_to <- fread(paste0(dir2,"/","genes_to_map.csv"),sep = ",")

for(i in 1:nrow(genes_to)){
  message(i)
  msu_i <- genes_to$MSU_list[genes_to$MSU_list %in% unlist(strsplit(data_rxns_ricecyc$protein[i],";"))]
  rap_i <- genes_to$ci_gen[genes_to$MSU_list %in% unlist(strsplit(data_rxns_ricecyc$protein[i],";"))]
  if(length(msu_i)>0){
    m_i <- paste0(msu_i,collapse = " or ")
    m_i <- paste0("(",m_i,")")
    data_rxns_ricecyc$msu[[i]] <- m_i
    r_i <- paste0(rap_i,collapse = " or ")
    r_i <- paste0("(",r_i,")")
  data_rxns_ricecyc$rap[[i]] <- r_i
  }
  rm(msu_i,rap_i,m_i,r_i)
};rm(i)

write.xlsx(data_rxns_ricecyc,paste0(dir2,"/","ricenetdb_rxns.xlsx"),row.names = F,showNA = F)

split_genes_i <- sub("[(]", "", data_rxns_ricecyc$msu)
split_genes_i <- sub("[)]", "", split_genes_i)
split_genes_i <- strsplit(split_genes_i, "or")
split_genes_i <- trimws(unlist(split_genes_i))
split_genes_i <- na.omit(split_genes_i)
split_genes_i <- unique(split_genes_i)
