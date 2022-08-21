do.call <- function (what, args, quote = FALSE, envir = parent.frame())
{
    if (!is.list(args))
        stop("second argument must be a list")
    if (quote)
        args <- lapply(args, enquote)
    .Internal(do.call(what, args, envir))
}

