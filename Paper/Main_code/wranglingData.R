### Dataset DEF. Wrangling and Class Balance
### Export files created; "Bac_Col_Full.csv"
## Libraries


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

rm(list = ls())
setwd("~/Documentos/TFM.LoadingData/Other_Input")
#Load the data; Presence/Ausence & Metadata
binary <- read.table("gene_presence_absence.Rtab",header=TRUE)
metadata <- read.xlsx(file = "Last_Clases.xlsx", header = TRUE,sheetIndex = 1)




#### Criva ####
cepa_sum <- apply(BC_Balanced[,-c(1)],1,function(x){
  as.integer(x)
  sum(x)
})
binary <- binary[which(cepa_sum >= 122),]

apply(BC_Balanced, 2, function(x){if})


Gene_npresen <- vector(mode="integer", length=0)
for (i in 2:19494){
  x <- as.numeric(as.character(BC_Balanced[,i]))
  y <- sum(x)
  print(y)
  Gene_npresen <- c(Gene_npresen,y)
}


cepa_sum <- apply(binary[,-1],1,sum)
binary <- binary[which(cepa_sum >= 122),]

########### Metadata wrangling
metadata_class <- metadata[,c(1,18:21)]
rownames(metadata_class)<- c()

metadata_class <- metadata_class[which(metadata_class$HOST_DISEASE_RECOD=="Bacteriemia"| metadata_class$HOST_DISEASE_RECOD == "Colonización"
                                       | metadata_class$HOST_DISEASE_RECOD=="IPPB" | metadata_class$HOST_DISEASE_RECOD == "ITU" 
                                       | metadata_class$HOST_DISEASE_RECOD == "Respiratorio"),]

### Ambiguous Column
metadata_class$Ambiguous <- FALSE

## Respiratory

metadata_class[grepl("utum",metadata_class$ISOLATION_SOURCE),"Ambiguous"] <-TRUE 
metadata_class[grepl("eumon",metadata_class$HOST_DISEASE),"Ambiguous"] <-TRUE 

# Ambigua ya que mezcla IPPT con Colonizacion
metadata_class[which(metadata_class$Nº==2294),"Ambiguous"]<- TRUE

## Colonización/Perirectal

## IPPB; Infeccion de piel y partes blandas

# Crees usted que las muestras de esta clase donde hay mucho unknowm en ISOLATION_SOURCE tienen alta probabilidad de venir de la piel o de partes blandes y que
#significa partes blandas.
## ITU class

# Cree usted que  puede haber cepas probenientes de muestras de orina que vengan de la sangre de una bateremia? Supogno que si es asi es mas probable que sean las que
# el Host_disease esta sin completar. O cree usted que si esta sin completad es porque seguramente no tenga bacteremia ya que analisis de sangre es una prueba muy comun 
# y se la harian a todos los pacientes?
metadata_class[1222:1247,"HOST_DISEASE"] <- as.factor("Unknown")
metadata_class[which(metadata_class$HOST_DISEASE=="Unknown" & metadata_class$HOST_DISEASE_RECOD=="ITU"),"Ambiguous"] <- TRUE
metadata_class[which(metadata_class$HOST_DISEASE=="unknown" & metadata_class$HOST_DISEASE_RECOD=="ITU"),"Ambiguous"] <- TRUE

## Bacteremia
metadata_class[which(metadata_class$HOST_DISEASE_RECOD == "Bacteriemia" & metadata_class$ISOLATION_SOURCE_RECOD== "Desconocido"),"Ambiguous"] <- TRUE
metadata_class[which(metadata_class$HOST_DISEASE_RECOD=="Bacteriemia" & metadata_class$ISOLATION_SOURCE_RECOD == "Respiratorio"),"Ambiguous"]<- TRUE
### It is necesary to change the column names to a nummber so it can be used from the metadata


### Export metadata class
write.csv(metadata_class, file = "metadata_class_def.csv",row.names = TRUE,col.names = TRUE)

####################Option 1. Pick all the class labels#############################
class_labels <- c("Bacteremia","Colonización", "ITU","IPPB")


###################Option2.  Pick only Bactermia and Colonization####################
class_labels <- c("Bacteremia", "Colonization")


bac_N <- metadata_class[which(metadata_class$HOST_DISEASE_RECOD == "Bacteriemia" & metadata_class$Ambiguous== FALSE),1]  # 259 strains
col_N <- metadata_class[which(metadata_class$HOST_DISEASE_RECOD == "Colonización" & metadata_class$Ambiguous== FALSE),1]  # 147 strains

###### Binary Wrangling 
## Change the colnames of the data to select the correct one

gene_names <- binary$Gene
gath <- binary[,-1]%>%
  t()%>%
  as.data.frame()
colnames(gath) <- gene_names
Strains_ID <- as.integer(gsub("ab","",rownames(gath)))
gath<- as.data.frame(cbind(Strains_ID,gath))

### Select only the rows for bac_N and col_N class
Bac_Col <- gath[gath$Strains_ID %in% c(bac_N,col_N),]

#Cheking if the nummbers are correct
bac_N2 <- intersect(Strains_ID,bac_N)
col_N2 <- intersect(Strains_ID,col_N)
class_column <- replicate(397,"")
Bac_Col <- as.data.frame(cbind(class_column,Bac_Col))
Bac_Col$class_column <- as.character(Bac_Col$class_column)
Bac_Col[Bac_Col$Strains_ID %in% bac_N2,"class_column"] <- class_labels[1]
Bac_Col[Bac_Col$Strains_ID %in% col_N2,"class_column"] <- class_labels[2]
Bac_Col$class_column <- as.factor(Bac_Col$class_column)

### Export the full dataset
write.csv(Bac_Col, file = "Bac_Col_Full.csv",row.names = TRUE,col.names = TRUE)

## Read the full dataset 
sfiug <- read.csv("Bac_Col_Full.csv",header = TRUE)

###create the training and the test set

bac_Random<- sample(1:254,100,replace = FALSE)
bac_Ntrain <- bac_N2[bac_Random]
bac_Nest <- bac_N2[-bac_Random]

col_Random <- sample(1:143,100,replace = FALSE)
col_Ntrain <- col_N2[col_Random]
col_Nest <- col_N2[-col_Random]

Bac_Col_Train <- Bac_Col[Bac_Col$Strains_ID %in% c(bac_Ntrain,col_Ntrain),]
write.csv(Bac_Col_Train, file = "Bac_Col_Train.csv",row.names = TRUE,col.names = TRUE)

Bac_Col_Rest <- Bac_Col[Bac_Col$Strains_ID %in% c(bac_Nest, col_Nest),]
write.csv(Bac_Col_Rest, file = "Bac_Col_Rest.csv",row.names = TRUE,col.names = TRUE)

### Use 100 of each one and 47 as training and test set
############ Make the  test set ###########

BC_rest <- read.csv("Bac_Col_Rest.csv",header = TRUE)
col_test <- BC_rest[BC_rest$class_column== class_labels[1],]
col_test <- col_test[sample(1:154, 43,replace = FALSE),]

BC_test <- rbind(col_test,BC_rest[which(BC_rest$class_column== class_labels[2]),])
rownames(BC_test) <- BC_test$X 
BC_test <- BC_test[,-1]



prove <- metadata_class[metadata_class$Nº %in% bac_N,]



intersect(metadata_class$Nº,bac_N)



#### Input for Phandango###
# Right now, this is not going to be used
#BC_Balanced is obtained runing the code of the first part of "Coss-Validation Scheme.RMD"
BC_pres_Phan <- BC_Balanced[,-1]

write.csv(BC_pres_Phan, file = "Presencia_Phandango.csv",row.names = TRUE,col.names = TRUE)

