#' Rename Parameters in an Analysis Object
#'
#' Renames parameters in an output object. This might be helpful if you used different
#' parameter names for the same parameter and want to fix this without recomputing.
#'
#' @param x An object of class \code{multiverseMPT}.
#' @param rename. A named character vector of the form
#'   \code{c(oldname1 = "newname", oldname2 = "newname2", ...)}.

#' @rdname rename_parameters
#' @export

rename_parameters <- function(x, rename) {
  UseMethod("rename_parameters")
}

#' @rdname rename_parameters
#' @export

rename_parameters.multiverseMPT <- function(x, rename = c()) {

  cols <- sapply(X = x, FUN = is.list)

  for (i in names(cols)[cols]) {
    x[[i]] <- lapply(X = x[[i]], FUN = function(x) {
      if("parameter" %in% colnames(x)) {
        x$parameter[x$parameter %in% names(rename)] <- rename[x$parameter[x$parameter %in% names(rename)]]
      }
      x
    })
  }
  x
}
