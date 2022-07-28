
read.nrrd.header<-function (file, Verbose = FALSE) 
{
  nrrdspec = list()
  if (!inherits(file, "connection")) {
    con <- file(file, open = "rt")
    attr(nrrdspec, "path") = file
    on.exit(close(con))
  } else {con = file}
  headerLines = readLines(con, 1)
  NRRDMAGIC = "NRRD000"
  if (substring(headerLines, 1, nchar(NRRDMAGIC)) != NRRDMAGIC) 
    stop("This does not appear to be a NRRD file: ", summary(con)$description)
  nrrdkeyvals = vector("character")
  while (length(l <- readLines(con, 1)) > 0 && l != "") {
    headerLines = c(headerLines, l)
    if (substring(l, 1, 1) == "#") 
      next
    if (length(grep(": ", l)) > 0) {
      #https://stackoverflow.com/questions/41717781/warning-input-string-not-available-in-this-locale
      #hingepos = regexpr(": ", l, fixed = TRUE,useBytes = TRUE)
      # as I already set 'Sys.setlocale("LC_ALL", "C")' according to
      # https://stackoverflow.com/questions/4993837/r-invalid-multibyte-string
      # when use nchar(l) report error
      hingepos = regexpr(": ", l, fixed = TRUE)
      fieldname = substring(l, 1, hingepos - 1)
      if (!fieldname %in% c("space dimension", "space units", 
                            "space origin", "space directions", "measurement frame")) 
      {fieldname = gsub(" ", "", fieldname, fixed = TRUE)}
      fieldval = substring(l, hingepos + 2, nchar(l))
      if (fieldname == "content") {
      }
      else if (substring(fieldval, 1, 1) == "(") {
        fieldval = gsub(" ", "", fieldval)
        fieldval = substring(fieldval, 2, nchar(fieldval) - 
                               1)
        vectorstring = unlist(strsplit(fieldval, ")(", 
                                       fixed = TRUE))
        tc = textConnection(vectorstring)
        fieldval = scan(tc, sep = ",", quiet = TRUE)
        if (length(vectorstring) > 1) 
          fieldval = matrix(fieldval, byrow = TRUE, nrow = length(vectorstring))
        close(tc)
      }
      else if (!fieldname %in% c("type", "datafile")) {
        if (length(grep("^[\\-+]{0,1}[0-9.]+", fieldval, 
                        perl = T)) > 0) 
          what = 0
        else what = ""
        tc = textConnection(fieldval)
        fieldval = scan(tc, quiet = TRUE, what = what)
        close(tc)
      }
      else if (fieldname == "datafile") {
        if (substring(fieldval, 1, 4) == "LIST") {
          fieldval = c(fieldval, readLines(con))
        }
      }
      nrrdspec[[fieldname]] = fieldval
    }
    else if (length(grep(":=", l)) > 0) {
      hingepos = regexpr(":=", l, fixed = TRUE)
      nrrdkeyvals[substring(l, 1, hingepos - 1)] = substring(l, 
                                                             hingepos + 2, nchar(l))
    }
    else {
      warning("Skipping malformed line #", length(headerLines), 
              " in NRRD header\n")
    }
  }
  attr(nrrdspec, "headertext") = headerLines
  attr(nrrdspec, "keyvals") = nrrdkeyvals
  nrrdspec
}