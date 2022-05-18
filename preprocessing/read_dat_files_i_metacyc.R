options(java.parameters = "-Xmx8192m")
require(data.table);require(dplyr);require(plyr);require(stringr);require(xlsx);require(rBiopaxParser)
  
  dir <- "D:/MET_BIOCYC_DBs/METACYC/25.1"
  out_dir <- "D:/TESIS_PHD/CHAPTER2/metacyc"
  # enzymes <- fread(paste0(dir,"/data/","enzymes.col"),sep = "\t",skip = 22)
  # 
  # genes<- fread(paste0(dir,"/data/","genes.col"),sep = "\t",skip = 22)
  # 
  # pathways <- fread(paste0(dir,"/data/","pathways.col"),sep = "\t",skip = 22,fill = T)
  # 
  # enzymes <- fread(paste0(dir,"/data/","enzymes.col"),sep = "\t",skip = 22)
  
  #############################################################################################
  enzrxn <- fread(paste0(dir,"/data/","reactions.dat"),sep = "\\",skip = 68,header=F,fill=T) #reactions.dat #enzrxns.dat
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
  enzrxn <- fread(paste0(dir,"/data/","enzrxns.dat"),sep = "\\",skip = 62,header=F) #reactions.dat #enzrxns.dat
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
  enzrxn <- fread(paste0(dir,"/data/","proteins.dat"),sep = "\\",skip = 110,header=F,fill=T) #reactions.dat #enzrxns.dat
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
  #write.xlsx(df_total2,paste0(out_dir,"/","proteins_dat.xlsx"),row.names = F,showNA = F)
  write.table(df_total2,paste0(out_dir,"/","proteins_dat.tsv"),row.names = F,na = "",sep = "\t")
  
  #############################################################################################
  #############################################################################################
