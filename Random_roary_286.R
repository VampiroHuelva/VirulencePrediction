#### Definitive Wrangling for 286 strains roary
library(xlsx)
library(tidyr)
library(dplyr)
library(tibble)
library(knitr)
library(rJava)
library(rpart)
#library(janitor)
library(party)
library(partykit)
library(ggplot2)
library(RWeka)
library(randomForest)
library(iRF)
library(AUC)
library(ggpubr)

### github instalation of plot package
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")


rm(list = ls())


setwd("~/Documentos/TFM.LoadingData/Other_Input")
binary_BC <- read.table("Pres_ausen_286.Rtab",header=TRUE)

### Criva de genes
# Primero eliminamos los genes que esten en menos de 
cepa_sum <- apply(binary_BC[,-1],1,sum)
binary_BC <- binary_BC[which(cepa_sum >= 15),]


### Transpose it
genes_name <- binary_BC[,1]
satrins_name <- colnames(binary_BC)
binary_BC <- binary_BC[,-1]
binary_BC_t <- as.data.frame(t(binary_BC))
colnames(binary_BC_t) <- genes_name


### Eliminate the genes that are equal so it would not give information to the algoritmh
bad_genes <- c()
for (n in 1:length(binary_BC_t[1,])) {
  if(length(unique(binary_BC_t[,n]))==1){
    bad_genes <- append(bad_genes,n)
  }
}
binary_BC_t <- binary_BC_t[,-bad_genes]


### Add the class column
# From "Cross-Validation Scheme.Rmd" take the class column in "BC_Balanced"
# check if the order of the strains is the same so it posible to add it

binary_BC_t[1:10,1:3]
# Get the strains in the same order as the class column
binary_BC_t <- binary_BC_t[rownames(BC_Balanced),]
# insert the class column
binary_BC_DEF <- cbind(BC_Balanced[,1],binary_BC_t)
colnames(binary_BC_DEF)[1] <- "Class"

genes_name <- colnames(binary_BC_DEF)[c(-1)]





###### Repeat the Random Forest Analysis

### store it in other variable so you can call BC_train later for the class_column
data_RF_1 <- binary_BC_DEF

### Error rate variable
error_rate_1 <-  data.frame(n_var= numeric(0), OOB_error= numeric(0))


### RF ANALYSIS Functions
## Create the variables to store the number of variables and the number of iteration
#################### RF Nºof VARIABLES  Function################
# Create a function that
#n_var; nº of variebles from rh
#min_var; Minimun number of variables to stop the iteration
#reps; Number of times the same experiment is executed on the second loop. 
var_selection <- function(n_var, min_var, reps, data){
   
  names(data) <- make.names(names(data))
  n_iter <- 0
  while (n_var > min_var) {
    n_iter<- n_iter+1
    print(n_iter)
  
    if(n_iter>1){
      n_var <- round(n_var/2,digits = 0)
    ### Get the most important features
      import_table <- as.data.frame(importance(RF_1))
      import_table <- import_table[order(import_table$MeanDecreaseAccuracy,decreasing = TRUE),]
    ### Make the step of reducing the dataframe
      reduced_table <- data[,rownames(import_table)[1:n_var]]# change 300 by n_var
      reduced_table <- as.data.frame(cbind(binary_BC_DEF$Class,reduced_table))
      colnames(reduced_table)[1] <- "Class"
      data <- reduced_table
      names(data) <- make.names(names(data))
    
    ### Store the variable names
      nam <- paste("Genes_Names", n_var, sep = "_")
      assign(nam, colnames(reduced_table[,-1]))
    }
  
  
    for(j in 1:8){
      if (j > 1){RF_1 <- randomForest(Class~.,data= data ,importance=TRUE, proximity= TRUE, mtry.select.prob= sel.prob)}
      else {RF_1 <- randomForest(Class~.,data= data ,importance=TRUE, proximity= TRUE)}
    ### Calculate store and print the error form the confusion matrix
      err <- sum(RF_1$confusion[,3])/2
      error_rate_1 <- rbind(error_rate_1,c(n_var,err))
      print(error_rate_1)
      sel.prob <- RF_1$importance
    ### Extract the var_importance and the intersect from the previous one
    
    
    }
  
  }

  colnames(error_rate_1) <- c("Number_Variables","OOB_error")
}

##### PLOTS.1#####

### Making plots 
mtry_plot <- ggplot(error_rate_1, aes(x=Number_Variables, y=OOB_error)) + geom_point() +stat_summary(aes(y = OOB_error), fun.y=mean, colour="red", geom="line") 
mtry_plot_reduced <- ggplot(error_rate_1, aes(x=Number_Variables, y=OOB_error)) + geom_point() +stat_summary(aes(y = OOB_error), fun.y=mean, colour="red", geom="line")+coord_cartesian(xlim = c(0, 2100))
mtry_plot_super_reduced <- ggplot(error_rate_1, aes(x=Number_Variables, y=OOB_error)) + geom_point() +stat_summary(aes(y = OOB_error), fun.y=mean, colour="red", geom="line")+coord_cartesian(xlim = c(0, 200))

mtry_plot_arranged <- ggarrange(mtry_plot, mtry_plot_super_reduced, 
                                    labels = c("A", "B"),
                                    ncol = 2, nrow = 1)

### Store the image


### Store the plots


######### ~~ RF NºVariables 2 ###########
### let do it from 20 to 632
## There is a problem with the paralogs cause the code assign a different puntuation, 
names(binary_BC_DEF) <- make.names(names(binary_BC_DEF))

data_RF_2 <- binary_BC_DEF[, Genes_Names_632]



data_RF_2 <- as.data.frame(cbind(binary_BC_DEF$Class,data_RF_2))
colnames(data_RF_2)[1] <- "Class"

error_rate_2 <-  data.frame(n_var= numeric(0), OOB_error= numeric(0))

n_var <- 632
n_iter <- 0
while (n_var > 4) {
  n_iter<- n_iter+1
  n_var <- n_var-10
  
  if(n_iter>1){
    ### Get the most important features
    import_table <- as.data.frame(importance(RF_2))
    import_table <- import_table[order(import_table$MeanDecreaseAccuracy,decreasing = TRUE),]
    ### Make the step of reducing the dataframe
    reduced_table <- data_RF_2[,rownames(import_table)[1:n_var]]
    reduced_table <- as.data.frame(cbind(binary_BC_DEF$Class,reduced_table))
    colnames(reduced_table)[1] <- "Class"
    data_RF_2 <- reduced_table
    names(data_RF_2) <- make.names(names(data_RF_2))
    
    ### Store the variable names
    nam <- paste("Last_Genes_Names", n_var, sep = "_")
    assign(nam, colnames(reduced_table[,-1]))
  }
  
  
  for(j in 1:8){
    
    RF_2 <- randomForest(Class~.,data= data_RF_2 ,importance=TRUE, proximity= TRUE)
    ### Calculate store and print the error form the confusion matrix
    err <- sum(RF_2$confusion[,3])/2
    error_rate_2 <- rbind(error_rate_2,c(n_var,err))
    print(error_rate_2)
    
    ### Extract the var_importance and the intersect from the previous one
    
    
  }
  
}

colnames(error_rate_2) <- c("Number_Variables","OOB_error")
mtry_plot_2 <- ggplot(error_rate_2, aes(x=Number_Variables, y=OOB_error)) + geom_point() +stat_summary(aes(y = OOB_error), fun.y=mean, colour="red", geom="line")
mtry_plot_zoom_2 <- ggplot(error_rate_2, aes(x=Number_Variables, y=OOB_error)) + geom_point() +stat_summary(aes(y = OOB_error), fun.y=mean, colour="red", geom="line")+coord_cartesian(xlim = c(0, 150)) 

mtry_plot_arranged_2 <- ggarrange(mtry_plot_2, mtry_plot_zoom_2, 
                                labels = c("A", "B"),
                                ncol = 2, nrow = 1)




######### ~~ RF NºVariables 3 ###########
### let do it from 20 to 158
## There is a problem with the paralogs cause the code assign a different puntuation, 
names(binary_BC_DEF) <- make.names(names(binary_BC_DEF))

data_RF_3 <- binary_BC_DEF[, Genes_Names_158]



data_RF_3 <- as.data.frame(cbind(binary_BC_DEF$Class,data_RF_3))
colnames(data_RF_3)[1] <- "Class"

error_rate_3 <-  data.frame(n_var= numeric(0), OOB_error= numeric(0))

n_var <- 158
n_iter <- 0
while (n_var > 10) {
  n_iter<- n_iter+1
  n_var <- n_var-1
  print(n_iter)
  if(n_iter>1){
    ### Get the most important features
    import_table <- as.data.frame(importance(RF_3))
    import_table <- import_table[order(import_table$MeanDecreaseAccuracy,decreasing = TRUE),]
    ### Make the step of reducing the dataframe
    reduced_table <- data_RF_3[,rownames(import_table)[1:n_var]]
    reduced_table <- as.data.frame(cbind(binary_BC_DEF$Class,reduced_table))
    colnames(reduced_table)[1] <- "Class"
    data_RF_3 <- reduced_table
    names(data_RF_3) <- make.names(names(data_RF_3))
    
    ### Store the variable names
    nam <- paste("Last_Genes_Names", n_var, sep = "_")
    assign(nam, colnames(reduced_table[,-1]))
  }
  
  
  for(j in 1:8){
    
    RF_3 <- randomForest(Class~.,data= data_RF_3 ,importance=TRUE, proximity= TRUE)
    ### Calculate store and print the error form the confusion matrix
    err <- sum(RF_3$confusion[,3])/2
    error_rate_3 <- rbind(error_rate_3,c(n_var,err))
    
    
    ### Extract the var_importance and the intersect from the previous one
    
    
  }
  
}

colnames(error_rate_3) <- c("Number_Variables","OOB_error")
mtry_plot_3 <- ggplot(error_rate_3, aes(x=Number_Variables, y=OOB_error)) + geom_point() +stat_summary(aes(y = OOB_error), fun.y=mean, colour="red", geom="line")
mtry_plot_zoom_3 <- ggplot(errsor_rate_3, aes(x=Number_Variables, y=OOB_error)) + geom_point() +stat_summary(aes(y = OOB_error), fun.y=mean, colour="red", geom="line")+coord_cartesian(xlim = c(10,47)) 


mtry_plot_arranged_3 <- ggarrange(mtry_plot_3, mtry_plot_zoom_3, 
                                labels = c("A", "B"),
                                ncol = 2, nrow = 1)



### Store the list of genes that score the minimun error in building the RF object
setwd("~/Documentos/TFM.LoadingData/Results")
saveRDS(Last_Genes_Names_47,file = "Genes_47.RData")



######### Final Random Forest 47#########
names(binary_BC_DEF) <- make.names(names(binary_BC_DEF))
Genes_47 <- readRDS("Genes_47.RData")
data_RF_47 <- binary_BC_DEF[, Genes_47]

data_RF_47 <- as.data.frame(cbind(binary_BC_DEF$Class,data_RF_47))
colnames(data_RF_47)[1] <- "Class"
RF_47 <- randomForest(Class~.,data= data_RF_47 ,importance=TRUE, proximity= TRUE)

####### There is another minimun but all the methodology has to be done again when the terminology is ready.

### Trying something on the iteration
names(binary_BC_DEF) <- make.names(names(binary_BC_DEF))
Genes_47 <- readRDS("Genes_47.RData")
data_RF_47 <- binary_BC_DEF[, Genes_47]

data_RF_47 <- as.data.frame(cbind(binary_BC_DEF$Class,data_RF_47))
colnames(data_RF_47)[1] <- "Class"
RF_47 <- randomForest(Class~.,data= data_RF_47 ,importance=TRUE, proximity= TRUE)
sel.prob <- RF_47$importance
RF_47_sel <- randomForest(Class~.,data= data_RF_47 ,importance=TRUE, proximity= TRUE,mtry.select.prob=sel.prob)








