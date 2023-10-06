#' Relabel clusters 
#' @keywords internal
#' @param ix vector of cluster labels
.relabell <- function(ix) 
{
  if(min(ix)==0) ix <- ix+1
  while( max(ix) != length(unique(ix)) )
  {
    missing_label1 <- (1:max(ix))[sapply(1:max(ix), function(x) !(x %in% ix) )][1]
    ix[ix>missing_label1] <- ix[ix>missing_label1]-1
  }
  return(ix)
}
