
citrullinate <- function(seq_list, mod_pos){
  #In each sequence of seq_list, replace R with Q
  #at position in mod_pos
  
  num_seq = length(seq_list)
  seq_char <- strsplit(seq_list, split = character(0))
  out_seq <- NA
  maxcit <- dim(mod_pos)[2]
   for (i in 1:num_seq){
  #First check whether mod_pos in seq is truly an arginine
  R_idx <- which(seq_char[[i]] == "R")
  #Figure out how many replacements we're making
  ncit <- maxcit - sum(is.na(mod_pos[i,]))
  if (!all((mod_pos[i,1:ncit])%in% R_idx)){
    print(paste0("Warning: Arg does not exist at a position marked for modification in sequence #",i))
    print(paste0("Aborting modification of sequence #",i))
    }
 else{
  #Then replace with Q
  seq_char[[i]][mod_pos[i,1:ncit]] = "Q"}
  
  #Reduce single character string back into sequence
   out_seq[i] <- paste(seq_char[[i]], collapse = "")
   }
  
  return (out_seq)
}


