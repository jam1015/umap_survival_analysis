all_dup <- function (value){
  duplicated(value) | duplicated(value, fromLast = TRUE)
}