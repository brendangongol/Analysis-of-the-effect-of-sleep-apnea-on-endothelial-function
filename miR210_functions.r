##################################
#### Write out table function ####
##################################
writab <- function(DT, fname){
   write.table(DT, fname, row.names=FALSE, quote=FALSE, sep="\t")
}

#### frequency counting function ####
#####################################
freqsdt <- function(DTstr, xstr, percent=TRUE)
{    
    if (percent)
        return(eval(parse(text=sprintf('%s[,.(frequency=.N),.(%s)]', DTstr, xstr)))[
            order(-frequency)][,percent:=100*frequency/sum(frequency)]) 
    else
        return(eval(parse(text=sprintf('%s[,.(frequency=.N),.(%s)]', DTstr, xstr)))[
            order(-frequency)]) 
}  

#### miR calculation functions ####
###################################
meanna <- function(var) mean(var,na.rm=TRUE)

norm <- function(d210 = NA, d92 = NA, d21 = NA, source){
    if(length(d210) > 1){
    meang1 <- ifelse(source == "UCSD", meangdt[[1]][1]$meanave, ifelse(source == "taiwan", meangdt[[1]][2]$meanave, meangdt[[1]][3]$meanave))  
    fc1 <- d210/meang1
    }else if(length(d92) > 1){
    meang1 <- ifelse(source == "UCSD", meangdt[[2]][1]$meanave, ifelse(source == "taiwan", meangdt[[2]][2]$meanave, meangdt[[2]][3]$meanave))  
    fc1 <- d92/meang1
    } else{
    meang1 <- ifelse(source == "UCSD", meangdt[[3]][1]$meanave, ifelse(source == "taiwan", meangdt[[3]][2]$meanave, meangdt[[3]][3]$meanave))  
    fc1 <- d21/meang1
    }
  return(list(fc1))}

#### function that runs a logistic regression model on each independent variable in a data table ####
#### The key is a numeric string indicating the columns to use ######################################
#### e.g. key <- 2:length(colnames(dt_sub))-1 #######################################################
univariantglmR <- function(DT, key, significant){
  compiled_dt <- data.table()
  pb <- txtProgressBar(min = 0, max = length(key), style = 3)
  for(m in 2:length(key)){
  #print(m)
    #### subset data table to contain only required columns
    newdt <- DT[,c(as.numeric(key[m])), with = FALSE]
    #### add outcome variable to last column 
    newdt$outcome <- DT$outcome
    
    #### generate formula 
    nam <- colnames(newdt[, -(ncol(newdt)), with = FALSE])
    if(length(nam) > 1){
      plus <- "+"
      p1 <- paste(nam[-length(nam)], plus, sep = " ")
      p2 <- paste(p1, nam[length(nam)])
      fmla <- as.formula(paste0(colnames(newdt)[ncol(newdt)], " ~ ", p2))
    }else{
      fmla <- as.formula(paste0(colnames(newdt)[ncol(newdt)], " ~ ", nam))
      }
    
    #### fit glm model
    glm.fit <- glm(fmla, data=newdt, family = binomial)
    
    #### create data table 
    dtf <- data.table(coef(summary(glm.fit)))#[,"Pr(>|z|)"]
    dtf$names <- rownames(summary(glm.fit)$coefficients)
    chr <- as.character(glm.fit$formula)
    dtf$formula <- paste(chr[2], chr[1], chr[3:length(chr)])
    if(significant == "T"){
        #### subset data table to retain only significant indepenent variables 
        dtf <- dtf[dtf$`Pr(>|z|)` < 0.05 & (!dtf$names == "(Intercept)"),]
    }
    if(nrow(dtf) > 0){
      compiled_dt <- rbind(compiled_dt, dtf)
      }
    setTxtProgressBar(pb, m)
  }
  close(pb)
  compiled_dt <- compiled_dt[!duplicated(compiled_dt[,c(1:6)]),]
  compiled_dt <- compiled_dt[(!compiled_dt$names == "(Intercept)"),]
  
  compiled_dt$Significance <- "NA"
  compiled_dt[`Pr(>|z|)` < 0.05,]$Significance <- "Significant"
  compiled_dt[`Pr(>|z|)` > 0.05,]$Significance <- "Not-Significant"
  return(compiled_dt)
}

#### function that runs a linear regression model on each independent variable in a data table ####
#### The key is a numeric string indicating the columns to use ######################################
#### e.g. key <- 2:length(colnames(dt_sub))-1 #######################################################
univariantlmR <- function(DT, key, significant){
  compiled_dt <- data.table()
  pb <- txtProgressBar(min = 0, max = length(key), style = 3)
  for(m in 2:length(key)){
    #print(m)
    #### subset data table to contain only required columns
    newdt <- DT[,c(as.numeric(key[m])), with = FALSE]
    #### add outcome variable to last column 
    newdt$outcome <- DT$outcome
    
    #### generate formula 
    nam <- colnames(newdt[, -(ncol(newdt)), with = FALSE])
    if(length(nam) > 1){
      plus <- "+"
      p1 <- paste(nam[-length(nam)], plus, sep = " ")
      p2 <- paste(p1, nam[length(nam)])
      fmla <- as.formula(paste0(colnames(newdt)[ncol(newdt)], " ~ ", p2))
    }else{
      fmla <- as.formula(paste0(colnames(newdt)[ncol(newdt)], " ~ ", nam))
    }
    
    #### fit glm model
    glm.fit <- lm(fmla, data=newdt)
    
    #### create data table 
    dtf <- data.table(coef(summary(glm.fit)))#[,"Pr(>|z|)"]
    dtf$names <- rownames(summary(glm.fit)$coefficients)
    chr <- as.character(fmla)
    dtf$formula <- paste(chr[2], chr[1], chr[3:length(chr)])
    if(significant == "T"){
      #### subset data table to retain only significant indepenent variables 
      dtf <- dtf[dtf$`Pr(>|t|)` < 0.05 & (!dtf$names == "(Intercept)"),]
    }
    if(nrow(dtf) > 0){
      compiled_dt <- rbind(compiled_dt, dtf)
    }
    setTxtProgressBar(pb, m)
  }
  close(pb)
  compiled_dt <- compiled_dt[!duplicated(compiled_dt[,c(1:6)]),]
  compiled_dt <- compiled_dt[(!compiled_dt$names == "(Intercept)"),]
  
  compiled_dt$Significance <- "NA"
  compiled_dt[compiled_dt$`Pr(>|t|)` < 0.05,]$Significance <- "Significant"
  compiled_dt[compiled_dt$`Pr(>|t|)` > 0.05,]$Significance <- "Not-Significant"
  return(compiled_dt)
}

#### Function that generates all combinations or models for a partucilar set of columns #########################
#### The data table must contain only columns one wants to run and the dependent variable needs to ##############
#### be in a separate column called "outcome" ###################################################################
#### If significant = T, only components of a model that are significant will be returned #######################
glmcompileR <- function(DT, key, significant){
  compiled_dt <- data.table()
  pb <- txtProgressBar(min = 0, max = length(key), style = 3)
  for(m in 2:length(key)){
    #### subset data table to contain only required columns
    spl <- strsplit(key[m], split = "_")[[1]]
    spl <- spl[!spl == "0"]
    newdt <- data.table()
    for(i in 1:length(spl)){
      temp <- DT[,c(as.numeric(spl[i])), with = FALSE]
      newdt <- cbind(newdt, temp)
      }
    #### add outcome variable to last column 
    newdt$outcome <- DT$outcome
    
    #### generate formula 
    nam <- colnames(newdt[, -(ncol(newdt)), with = FALSE])
    if(length(nam) > 1){
      plus <- "+"
      p1 <- paste(nam[-length(nam)], plus, sep = " ")
      p2 <- paste(p1, nam[length(nam)])
      fmla <- as.formula(paste0(colnames(newdt)[ncol(newdt)], " ~ ", p2))
    }else{
      fmla <- as.formula(paste0(colnames(newdt)[ncol(newdt)], " ~ ", nam))
      }
    
    #### fit glm model
    glm.fit <- glm(fmla, data=newdt, family = binomial)
    
    #### create data table 
    dtf <- data.table(coef(summary(glm.fit)))#[,"Pr(>|z|)"]
    dtf$names <- rownames(summary(glm.fit)$coefficients)
    chr <- as.character(glm.fit$formula)
    dtf$formula <- paste(chr[2], chr[1], chr[3:length(chr)])
    if(significant == "T"){
        #### subset data table to retain only significant indepenent variables 
        dtf <- dtf[dtf$`Pr(>|z|)` < 0.05 & (!dtf$names == "(Intercept)"),]
    }
    if(nrow(dtf) > 0){
      compiled_dt <- rbind(compiled_dt, dtf)
      }
    setTxtProgressBar(pb, m)
  }
  close(pb)
  compiled_dt <- compiled_dt[!duplicated(compiled_dt[,c(1:6)]),]
  compiled_dt <- compiled_dt[(!compiled_dt$names == "(Intercept)"),]
 
  compiled_dt$Significance <- "NA"
  compiled_dt[`Pr(>|z|)` < 0.05,]$Significance <- "Significant"
  compiled_dt[`Pr(>|z|)` > 0.05,]$Significance <- "Not-Significant"
  return(compiled_dt)
}

################################
#### OSA severity annotater ####
################################
# Normal: AHI<5
# Mild sleep apnea: AHI < 15
# Moderate sleep apnea: AHI < 30
# Severe sleep apnea: AHI > 30
OSAannotate <- function(DT){
  ahi <- DT$ahi
  severity <- as.character(vector())
  condition <- c("normal","mild","moderate","severe","unknown")
  for(i in 1:length(ahi)){
    #### annotate severity for pediatric patients
    if(DT$source[i] == namestr[['location']][1]){
      if(ahi[i] < 2) severity[i] <- condition[1]
      else if(ahi[i] >= 2 & ahi[i] < 5) severity[i] <- condition[2]
      else if(ahi[i] >= 5 & ahi[i] < 10) severity[i] <- condition[3]
      else if(ahi[i] >= 10) severity[i] <- condition[4]
    }
    #### annotate severity for adult patients
    else if(DT$source[i] == namestr[['location']][2]){
      if(ahi[i] < 5) severity[i] <- condition[1]
      else if(ahi[i] >= 5 & ahi[i] < 15) severity[i] <- condition[2]
      else if(ahi[i] >= 15 & ahi[i] < 30) severity[i] <- condition[3]
      else if(ahi[i] >= 30) severity[i] <- condition[4]
    }
    else if(DT$source[i] == namestr[['location']][3]){
      if(ahi[i] < 5) severity[i] <- condition[1]
      else if(ahi[i] >= 5 & ahi[i] < 15) severity[i] <- condition[2]
      else if(ahi[i] >= 15 & ahi[i] < 30) severity[i] <- condition[3]
      else if(ahi[i] >= 30) severity[i] <- condition[4]
    }}
  return(severity)}

#############################
#### Taiwan med combiner ####
#############################
TaiwanMeds <- function(DT ){
  meds <- NULL
  combined_meds <- NULL
  for(i in 1:nrow(DT)){
    meds <- NULL
    if(is.na(DT[i,])[1] == TRUE){
      combined_meds[i] <- NA
    }else if(is.na(DT[i,])[1] == FALSE){
      if(DT$asa[i] == 1){one <- "aspirin"; meds <- c(meds, one)}
      if(DT$plavix[i] == 1){two <-  "plavix" ; meds <- c(meds, two)}
      if(DT$acei[i] == 1){three <-  "ACE inhibitor"; meds <- c(meds, three)}
      if(DT$arb[i] == 1){four <-  "Angiotensin II receptor blocker"; meds <- c(meds, four)}
      if(DT$bb[i] == 1){five <-  "beta blocker"; meds <- c(meds, five)}
      if(DT$ccb[i] == 1){six <-  "calcium channel blocker"; meds <- c(meds, six)}
      if(DT$ntg[i] == 1){seven <-  "Nitroglycerin" ; meds <- c(meds, seven)}
      if(DT$diuretics[i] == 1){eight <- "diuretics"; meds <- c(meds, eight)}
      if(DT$statin[i] == 1){nine <- "statin"; meds <- c(meds, nine)}
      # meds
      if(length(meds > 0)){meds <- paste(meds, collapse = ", "); combined_meds[i] <- meds
      }else{
        combined_meds[i] <- "None"
      }}}
  return(combined_meds)}

####################################
#### flatten correlation matrix ####
####################################
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )}

####################################
#### scatterplot graph function ####
####################################
scatterplotgraph <- function(dataset, xaxis, yaxis, title, xlabel, ylabel){
  p <- ggplot(dataset, aes(x=xaxis, y=yaxis))+
  geom_point(shape = 21, size = 3, colour = "black", fill = "#08519C")+
  scale_color_manual(values=c("#E31A1C")) +
  geom_smooth(method = "lm")+
  ggtitle(title)+ xlab(xlabel) + ylab(ylabel) +
  theme_pubr()
return(p)
}

################################
#### dotplot graph function ####
################################
dotplotgraph <- function(dataset, xaxis, yaxis, title, dodge, bin, xlabel, ylabel){
  p <- ggplot(dataset, aes(x= xaxis, y= yaxis, fill = source)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(dodge), binwidth = bin) + 
  ylab(ylabel) + xlab(xlabel) + ggtitle(title) + theme_pubr()
return(p)
}

#### This function runs a mann whitney u test to compare continuous variables between two data tables ####
#### Compair = a character name annotating what variable the data tables are partitioned on. #############
#### DT1 & DT2 = The two partitioned data tables to compare ##############################################
#### keys = a character string indicating the overlap between the two columns to compare #################
continuous_compare <- function(DT1, DT2, compair, key){
test_dt <- NULL
comparisons <- compair
for(i in 1:length(key)){
  test <- suppressWarnings(wilcox.test(DT1[[key[i]]], DT2[[key[i]]], conf.int = TRUE, paired = FALSE, formula = "lhs"))
  if(test$p.value < 0.05){sig <- "TRUE"
  }else{sig <- "FALSE"
  }
test_dt <- rbind(test_dt, data.table(comparison = comparisons, parameter = key[i], Mann_Whitney_U_p_value = test$p.value, significant = sig))}
return(test_dt)
}

#### This function runs a chai squared test to compare categorical variables between two data tables ####
#### Compair = a character name annotating what variable the data tables are partitioned on. ############
#### compair_column = a character indicating the name of the column to partition data on ################
#### DT1 = The data table to partition the data on ######################################################
#### keys = a character string indicating the overlap between the two columns to compare ################
categorical_compare <- function(DT, compair, compair_column, keys){

test_dt <- NULL
comparisons <- compair
for(i in 1:length(keys)){
    mytab <- with(DT,table(DT[[compair_column]],DT[[keys[i]]]))
    test <- chisq.test(mytab)
  if(test$p.value < 0.05){sig <- "TRUE"
  }else{sig <- "FALSE"
  }
test_dt <- rbind(test_dt, data.table(comparison = comparisons, parameter = keys[i], Chai_squared_p_value = test$p.value, significant = sig))}
return(test_dt)
}

#### function to download gene ontology information from BioMart ####
#####################################################################
BioMartR <- function(value, mart, species, attribute, filter_by){
Values <- value
#listMarts()
myMart <- useMart(mart)
# listDatasets(myMart)
myMart <- useMart(mart, dataset=species)
# listAttributes(myMart)[1:200,] # Choose data types you want to download
# myMart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast") # use if biomart sever is down.
go <- getBM(attributes=attribute, mart=myMart, values = Values, filters = filter_by)
return (go)
}


#### function to plot ChIP peaks ####
#####################################
ChiPseqPeakPlotter_single <- function(bw_file, chr, start, end, average_dist, fill_dist, type="all", peaks, x_lab, color){
  
  if(type == "all"){
    ##############################################
    ######### Code to get score from bw files ####
    ##############################################
    which <- GRanges(c(chr), IRanges(c(start), c(end)))
    bw_tr <- import(bw_file, which = which)
    
    #### format treatment data ####
    ###############################
    scores <- data.table()
    for(i in 1:length(start(ranges(bw_tr)))){
      num <- start(ranges(bw_tr))[i]:end(ranges(bw_tr))[i]
      sco <- score(bw_tr)[i]
      dt <- data.table(num, sco)
      scores <- rbind(scores, dt)
    }
    scores
    len <- data.table(start:end); setnames(len, colnames(len), "num")
    zz <- merge(len, scores, all = TRUE)
    zz[is.na(zz)] <- 0
    zz$key <- 1:nrow(zz)
    
    #### create averages across specified windows ####
    ##################################################
    mea <- suppressWarnings(colMeans(matrix(zz$sco, average_dist)))######################
    dt <- data.table(mea)
    #### fill in lost rows ####
    data_tx <- data.table()
    for(i in 1:nrow(dt)){
      d <- data.table(rep(dt$mea[i], fill_dist)); setnames(d, colnames(d), "mea")
      data_tx <- rbind(data_tx, d)
    }
    data_tx$num <- 1:nrow(data_tx)
    data_tx$fill <- "tx"
    dt <- data_tx
    
    #### plot data ####
    ###################
    p <- ggplot(aes(x=num, y=mea), data = dt)+
      geom_line(aes(color = fill))+   #color = c("#1B9E77"))+
      #scale_color_brewer(palette = "Dark2")+ # "Set1"
      
      #### add fill below line and color ####
    geom_area(aes(fill=fill))+
      #scale_fill_brewer(palette="Dark2")+ # "Dark2"  "Set1"
      scale_fill_manual(values=c( color))+ # Manually fill colors #"#1B9E77", "#7570B3"
      scale_color_manual(values=c(color)) + 
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # axis.title.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.x=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.y=element_blank())+
      xlab(x_lab)
    return(p)
  }
  
  else if(type == "called"){
    
    #########################################
    ### plot only peaks that were called ####
    #########################################
    
    #### get peaks for treatment and control files ####
    ###################################################
    df_tx <- read.delim(peaks, comment = "#")
    df_tx$seqnames <- gsub("chr", "", df_tx$seqnames, ignore.case = TRUE)
    peaks_tx <- as(df_tx, "GRanges")
    

    #### subset out peaks for region of interest ####
    #################################################
    which <- GRanges(c(chr), IRanges(c(start), c(end)))
    
    which_tx <- subsetByOverlaps(peaks_tx, which)
    
    #### get score for selected regions ####
    ########################################
    bw_tr <- suppressWarnings(import(treatment_bw, which = which_tx))########################
    
    #### format treatment data ####
    ###############################
    scores <- data.table()
    for(i in 1:length(start(ranges(bw_tr)))){
      num <- start(ranges(bw_tr))[i]:end(ranges(bw_tr))[i]
      sco <- score(bw_tr)[i]
      dt <- data.table(num, sco)
      scores <- rbind(scores, dt)
    }
    scores
    len <- data.table(start:end); setnames(len, colnames(len), "num")
    zz <- merge(len, scores, all = TRUE)
    zz[is.na(zz)] <- 0
    zz$key <- 1:nrow(zz)
    
    #### create averages across specified windows ####
    ##################################################
    mea <- suppressWarnings(colMeans(matrix(zz$sco, average_dist)))#####################
    dt <- data.table(mea)
    #### fill in lost rows ####
    data_tx <- data.table()
    for(i in 1:nrow(dt)){
      d <- data.table(rep(dt$mea[i], fill_dist)); setnames(d, colnames(d), "mea")
      data_tx <- rbind(data_tx, d)
    }
    data_tx$num <- 1:nrow(data_tx)
	data_tx$fill <- "tx"
    dt <- data_tx
    
    #### plot data ####
    ###################
    p <- ggplot(aes(x=num, y=mea), data = dt)+
      geom_line(color = "#1B9E77")+
      geom_area(aes(fill=fill))+
      #scale_fill_brewer(palette="Dark2")+	  
      scale_fill_manual(values=c( color))+ # Manually fill colors #"#1B9E77", "#7570B3"
      scale_color_manual(values=c(color)) + 
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # axis.title.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.x=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.y=element_blank())+
      xlab(x_lab)
    return(p)
  }
  
}







