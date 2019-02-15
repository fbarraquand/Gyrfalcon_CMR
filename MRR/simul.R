### Function from Kéry & Schaub BPA

# Define function to simulate multistate capture-recapture data
simul.ms <- function(PSI.STATE, PSI.OBS, marked, unobservable = NA){
   # Unobservable: number of state that is unobservable
   n.occasions <- dim(PSI.STATE)[4] + 1
   CH <- CH.TRUE <- matrix(NA, ncol = n.occasions, nrow = sum(marked))
   # Define a vector with the occasion of marking
   mark.occ <- matrix(0, ncol = dim(PSI.STATE)[1], nrow = sum(marked))
   g <- colSums(marked)
   for (s in 1:dim(PSI.STATE)[1]){
      if (g[s]==0) next  # To avoid error message if nothing to replace
      mark.occ[(cumsum(g[1:s])-g[s]+1)[s]:cumsum(g[1:s])[s],s] <-
      rep(1:n.occasions, marked[1:n.occasions,s])
      } #s
   for (i in 1:sum(marked)){
      for (s in 1:dim(PSI.STATE)[1]){
         if (mark.occ[i,s]==0) next
         first <- mark.occ[i,s]
         CH[i,first] <- s
         CH.TRUE[i,first] <- s
         } #s
      for (t in (first+1):n.occasions){
         # Multinomial trials for state transitions
         if (first==n.occasions) next
         state <- which(rmultinom(1, 1, PSI.STATE[CH.TRUE[i,t-1],,i,t-1])==1)
         CH.TRUE[i,t] <- state
         # Multinomial trials for observation process
         event <- which(rmultinom(1, 1, PSI.OBS[CH.TRUE[i,t],,i,t-1])==1)
         CH[i,t] <- event
         } #t
      } #i
   # Replace the NA and the highest state number (dead) in the file by 0
   CH[is.na(CH)] <- 0
   CH[CH==dim(PSI.STATE)[1]] <- 0
   CH[CH==unobservable] <- 0
   id <- numeric(0)
   for (i in 1:dim(CH)[1]){
      z <- min(which(CH[i,]!=0))
      ifelse(z==dim(CH)[2], id <- c(id,i), id <- c(id))
      }
   return(list(CH=CH[-id,], CH.TRUE=CH.TRUE[-id,]))
   # CH: capture histories to be used
   # CH.TRUE: capture histories with perfect observation
   }

