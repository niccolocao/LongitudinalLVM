native_symbol <- function(name) {
  get(paste0("native_", name), envir = asNamespace("LongitudinalLVM"),
      inherits = FALSE)
}

native_has <- function(name) {
  exists(paste0("native_", name), envir = asNamespace("LongitudinalLVM"),
         inherits = FALSE)
}

native_call <- function(name, ...) {
  symbol <- if (is.character(name)) native_symbol(name) else name
  .C(symbol, ...)
}
