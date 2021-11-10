core_abundance_RA33 <- function(seq_list, unique_core_list,ag_data){
 
    core_index <- sapply(unique_core_list,grep,seq_list) 

    df <- plyr::ldply(core_index, rbind)

    abundances = data.frame()
    length_of_cores = length(unique_core_list)
    
    for (i in 1:length_of_cores){
        abundances[i,1] = unique_core_list[i]
        abundances[i,2] = 0
        abundances[i,3] = 0
        abundances[i,4] = 0
    }
    
    colnames(abundances) <- c("Core","Native","PAD2","PAD4")
    
    for (i in 1:length_of_cores){ # For each core...
        total_abundance = 0
        num_peptides = length(core_index[[i]]) # 2 if there are 2 peptides containing given core
        if (num_peptides !=0){
            for (j in 1:num_peptides){
                peptide = core_index[[i]][[j]] # This should give us 143 for the first pass
                abundance = sum(c(ag_data$`Native_Replicate_1`[[peptide]],ag_data$`Native_Replicate_2`[[peptide]], ag_data$`Native_Replicate_3`[[peptide]], ag_data$`Native_Replicate_4`[[peptide]],ag_data$`Native_Replicate_5`[[peptide]],ag_data$`Native_Replicate_6`[[peptide]]),na.rm = TRUE)/num_peptides
                
                total_abundance = total_abundance + abundance
            }
            abundances[[i,2]] = total_abundance
        }
    }
    
    for (i in 1:length_of_cores){ # For each core...
        total_abundance = 0
        num_peptides = length(core_index[[i]]) # 2 if there are 2 peptides containing given core
        if (num_peptides !=0){
            for (j in 1:num_peptides){
                peptide = core_index[[i]][[j]] # This should give us 143 for the first pass
                abundance = sum(c(ag_data$`PAD2-Citrullinated_Replicate_1`[[peptide]],ag_data$`PAD2-Citrullinated_Replicate_2`[[peptide]], ag_data$`PAD2-Citrullinated_Replicate_3`[[peptide]], ag_data$`PAD2-Citrullinated_Replicate_4`[[peptide]],ag_data$`PAD2-Citrullinated_Replicate_5`[[peptide]],ag_data$`PAD2-Citrullinated_Replicate_6`[[peptide]]),na.rm = TRUE)/num_peptides
            
                total_abundance = total_abundance + abundance
            }
            abundances[[i,3]] = total_abundance
        }
    }
    
    for (i in 1:length_of_cores){ # For each core...
        total_abundance = 0
        num_peptides = length(core_index[[i]]) # 2 if there are 2 peptides containing given core
        if (num_peptides !=0){
            for (j in 1:num_peptides){
                peptide = core_index[[i]][[j]] # This should give us 143 for the first pass
                abundance = sum(c(ag_data$`PAD4-Citrullinated_Replicate_1`[[peptide]],ag_data$`PAD4-Citrullinated_Replicate_2`[[peptide]], ag_data$`PAD4-Citrullinated_Replicate_3`[[peptide]], ag_data$`PAD4-Citrullinated_Replicate_4`[[peptide]],ag_data$`PAD4-Citrullinated_Replicate_5`[[peptide]],ag_data$`PAD4-Citrullinated_Replicate_6`[[peptide]]),na.rm = TRUE)/num_peptides
                
                total_abundance = total_abundance + abundance
            }
            abundances[[i,4]] = total_abundance
        }
    }
    return(abundances)
}