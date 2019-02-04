##### RF Functions ######

### RF ANALYSIS Functions
## Create the variables to store the number of variables and the number of iteration
#################### RF Nºof VARIABLES  Function: var_selection()################
### Create a function that calculate a error rate dataframe with steps of variables reduction to elucidate the correct number of variables 
### that should be given as input. 
#n_var; nº of variebles for the input(if mising,all the variables are selected)
#min_var; Minimun number of variables to stop the iteration
#reps; Number of times the same experiment is executed on the second loop. 
# data: dataframe with a the variables as columns and the samples as rows plus a class column for each sample
#rest_var: number of variables reduced in each iteration(if msising select half of the best scored variables)
#List_N: Give the string "List_N" with N being the number of genes to have on the output.
####################

var_selection <- function(n_var, min_var, reps, data, rest_var, List_N){
  error_rate <-  data.frame(n_var= numeric(0), OOB_error= numeric(0))
  names(data) <- make.names(names(data))
  #Store the original dataframe so nothing is lost with the iterations
  orig_data <- data
  #if n_var is not given, use the maximun number of variables.
  if(missing(n_var)){n_var <- length(colnames(data))}
  n_iter <- 0
  
  while (n_var > min_var) {
    n_iter<- n_iter+1
    print(n_iter)
    
    if(n_iter>1){
      if(missing(rest_var)){n_var <- round(n_var/2,digits = 0)} else{n_var <- n_var + rest_var}
      ### Get the most important features
      import_table <- as.data.frame(importance(RF_1))
      import_table <- import_table[order(import_table$MeanDecreaseAccuracy,decreasing = TRUE),]
      ### Make the step of reducing the dataframe
      reduced_table <- data[,rownames(import_table)[1:n_var]]# change 300 by n_var
      reduced_table <- as.data.frame(cbind(orig_data$Class,reduced_table))
      colnames(reduced_table)[1] <- "Class"
      data <- reduced_table
      names(data) <- make.names(names(data))
      
      ### Store the variable names
      nam <- paste("List", n_var, sep = "_")
      assign(nam, colnames(reduced_table[,-1]))
    }
    
    
    for(j in 1:reps){
      RF_1 <- randomForest(Class~.,data= data,importance=TRUE)
      ### Calculate store and print the error form the confusion matrix
      err <- sum(RF_1$confusion[,3])/2
      error_rate <- rbind(error_rate,c(n_var,err))
      ### Extract the var_importance and the intersect from the previous one
      
      
    }
    
  }
  
  colnames(error_rate) <- c("Number_Variables","OOB_error")
  if(missing(List_N)){return(error_rate)} else {return(list(error_rate,Genes_Names_50))}
}



