nmi_wrapper = function(z, y) {
  NMI::NMI(cbind(seq_along(z), z), cbind(seq_along(y), y))$value
}


find_modes <- function(x) {
  # find the modes of a vector (the values with largest frequency)
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}

sort_labels = function(z) {
  z_new = z
  sorted_comms = sort(table(z), dec = T)
  old_labels = as.integer(names(sorted_comms))
  old_labels
  for (i in seq_along(old_labels)) {
    z_new[z == old_labels[i]] = i
  }
  z_new
}

get_seq_nmi = function(zout) {
  sapply(1:(ncol(zout)-1), function(itr) nett::compute_mutual_info(zout[,itr], zout[,itr+1]))
}

############### Computes the MAP estimated of the labels ##################
#' @export 
get_map_labels <- function(z_mat, burnin = NULL){
  niter = ncol(z_mat)
  if (is.null(burnin)) {
      burnin = round(niter/2)
  }
  
  z_mat = z_mat[,(burnin+1):niter]
  n = nrow(z_mat)
  z_map = sapply(1:n, function(i) find_modes(z_mat[i,])[1])
  confidence = sapply(1:n, function(i) sum(z_mat[i,] == z_map[i])/length(z_mat[i,]) )

  list(z_map = z_map, confidence = confidence)
}