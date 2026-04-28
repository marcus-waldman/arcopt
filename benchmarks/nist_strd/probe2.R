# Look at the full Rd content for one problem (not just Examples)

rd_db <- tools::Rd_db("NISTnls")
rd <- rd_db[["Thurber.Rd"]]

# Print all section tags
for (sec in rd) {
  tag <- attr(sec, "Rd_tag")
  if (!is.null(tag)) {
    cat("==== TAG:", tag, "====\n")
    txt <- paste(unlist(sec), collapse = "")
    cat(substr(txt, 1, 800), "\n\n")
  }
}
