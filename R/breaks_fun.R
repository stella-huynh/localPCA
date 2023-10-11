#' Define x-axis break intervals based on chromosome size
#'
#' This function defines break intervals on ggplot x-axis proportional to chromosome size.
#' If :
#'    - chromosome length < 100 Mb, set breaks every 10 Mb
#'    - chromosome length > 100 Mb but < 200 Mb, set breaks every 20 Mb
#'    - chromosome length > 200 Mb but < 350 Mb, set breaks every 50 Mb
#'    - chromosome length > 350 Mb, set breaks every 100 Mb.
#'
#' @param x A numeric value
#'
#' @return Returns a numeric vector with breaks between 0 and x.
#'
#' @export
#'


breaks_fun <- function(x) {

    if (max(x) < 100000000) { # if chr < 100 Mb : breaks every 10 Mb.
    b = seq(0, signif(max(x),2), 10000000)
  } else if (max(x) >= 100000000 & max(x) <= 200000000) { # if 100Mb < chr < 200 Mb : breaks every 20 Mb.
    b = seq(0, signif(max(x),2), 20000000)
  } else if (max(x) >= 200000000 & max(x) <= 350000000) { # if 200Mb < chr < 350 Mb : breaks every 50 Mb.
    b = seq(0, signif(max(x),2), 50000000)
  } else { # if chr > 350 Mb : breaks every 100 Mb.
    b = seq(0, signif(max(x),1), 100000000)
  }

  return(b)
}

