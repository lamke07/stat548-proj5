rescale <- function(x, original_from = 0, original_to = 1, target_from = 0, target_to = 100){
  # Function to rescale values from original range to range of interest
  return(target_from + x*(target_to - target_from)/(original_to - original_from))
}
