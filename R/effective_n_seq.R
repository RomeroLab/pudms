#' Number of effective sequences 
#'
#'@param n_reads number of reads
#'@param n_characterized number of function characterizations
effective_n_seq = function (n_reads, n_characterized){

  
  # Draw, with replacement, n_reads sequences from a set of nc:= n_characterized sequences (seq1,seq2,..., seq nc)
  # What proportion of the n_characterized sequences are in the sampled set?
  # Effective sequences = an expectated number of the present sequences in the sampled set 
  # Note: the number of effective sequences cannot exceed n_characterized
  
  # effective seq = E[ sum_i 1[seq i is present] ] where i = 1,...,nc.
  # P(seq i is present) = 1- P(seq i is never drawn) = 1- (1- 1/n_characterized)^n_reads
  
  n_reads*(1- (1- 1/n_characterized)^n_reads)
}
