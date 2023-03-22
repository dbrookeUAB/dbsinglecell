#' Internal Function Messages
#'
#' @param text
#'
#' @return
#' @export
#'
#' @examples
message_section <- function(text){
  n <- ceiling(options()$width)

  z <- paste0("\033[1m\033[33m[%s]:\033[39m\033[22m",  # bold and  yellow time
              '\033[1m\033[32m %2.', # bold and green start
              n-22,          # max number of characters to print
              's \033[39m\033[22m\n')  # bold and green finish

  cat("\n",rep('-',22), "\n",sep = '')
  cat(sprintf(z,Sys.time(), text))
  cat(rep('-',22), "\n\n",sep = '')
}

#' Internal Function Messages
#'
#' @param text
#'
#' @return
#' @export
#'
#' @examples
message_task <- function(text){
  n <- ceiling(options()$width)

  z <- paste0("\033[1m\033[33m[%s]:\033[39m\033[22m",  # bold and  yellow time
              '\033[37m %s', # white start
              # n-22,          # max number of characters to print
              '\033[39m\n')  # white finish

  cat(sprintf(z,Sys.time(), substr(text, 1,n)))
}

#' Internal Function Messages
#'
#' @param text
#'
#' @return
#' @export
#'
#' @examples
message_append <- function(text){
  n <- ceiling(options()$width)-1

  n.lines <- ceiling(nchar(text)/(n-23))


  z0 <- paste0("\033[90m - %.",  # bold and  yellow time,
               n-23,
               's \033[39m\n')  # white finish

  zn <- paste0("\033[90m %.",  # bold and  yellow time,
               n-23,
               's \033[39m\n')  # white finish

  for(i in 1:n.lines){
    q <- i-1
    start <- q*65+1
    stop <- q*65+65
    if(q==0){
      cat(rep(' ',22), sep = '')
      cat(sprintf(z0,substr(text, start,stop)))
    } else {
      cat(rep(' ',24), sep = '')
      cat(sprintf(zn,substr(text, start,stop)))
    }

  }

}
