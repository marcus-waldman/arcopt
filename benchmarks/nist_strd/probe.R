# Probe NISTnls structure: how do certified values, starting points, and
# the model formula get exposed?

suppressWarnings(library(NISTnls))

datasets <- c("Thurber", "Lanczos1", "Lanczos2", "Lanczos3",
              "MGH09", "MGH10", "MGH17", "Bennett5", "BoxBOD",
              "Eckerle4")

for (nm in datasets) {
  cat("=========================================\n")
  cat(sprintf("Dataset: %s\n", nm))
  cat("=========================================\n")
  if (exists(nm, envir = asNamespace("NISTnls"))) {
    e <- new.env()
    do.call(data, list(nm, package = "NISTnls", envir = e))
    if (exists(nm, envir = e)) {
      d <- get(nm, envir = e)
      cat(sprintf("  rows: %d  cols: %s\n", nrow(d),
                  paste(colnames(d), collapse = ", ")))
      cat("  attributes:", paste(names(attributes(d)), collapse = ", "), "\n")
    } else {
      cat("  NOT FOUND in package\n")
    }
  } else {
    cat("  not exported\n")
  }
}

# Look at one help page for full info on starting values + certified params
cat("\n\n--- attempting to find start values & certified params ---\n")

# NISTnls help pages embed starting values and certified values in
# Examples; the package object itself is just the raw data frame. Use
# tools::Rd_db to extract them.
rd_db <- tools::Rd_db("NISTnls")
for (nm in datasets) {
  rd_key <- paste0(nm, ".Rd")
  if (rd_key %in% names(rd_db)) {
    cat(sprintf("\n--- %s.Rd Examples ---\n", nm))
    rd <- rd_db[[rd_key]]
    # Find Examples section
    for (sec in rd) {
      if (!is.null(attr(sec, "Rd_tag")) &&
          attr(sec, "Rd_tag") == "\\examples") {
        cat(paste(unlist(sec), collapse = ""), "\n")
        break
      }
    }
  }
}
