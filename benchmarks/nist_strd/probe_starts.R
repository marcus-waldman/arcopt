# Dump start values from each NIST Rd file Examples section so we can
# verify against the canonical NIST Start I / Start II.

rd_db <- tools::Rd_db("NISTnls")
problems <- c("Bennett5", "Eckerle4", "Hahn1", "Lanczos1", "Lanczos2",
              "Lanczos3", "MGH09", "MGH10", "MGH17", "Thurber")

for (nm in problems) {
  cat("\n===== ", nm, " =====\n", sep = "")
  rd <- rd_db[[paste0(nm, ".Rd")]]
  for (sec in rd) {
    tag <- attr(sec, "Rd_tag")
    if (!is.null(tag) && tag == "\\examples") {
      txt <- paste(unlist(sec), collapse = "")
      # Print only lines containing "start ="
      lines <- strsplit(txt, "\n", fixed = TRUE)[[1L]]
      for (ln in lines) {
        if (grepl("start\\s*=\\s*c\\(", ln) ||
            grepl("^\\s*start", ln, ignore.case = TRUE)) {
          cat(ln, "\n")
        }
      }
    }
  }
}
