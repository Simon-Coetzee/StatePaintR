
#' @importFrom data.table fread
#' @importFrom stringr str_replace
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom GenomicRanges GRangesList GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols mcols<-
GetBioFeatures <- function(bio.features.loc = NULL, search.term = NULL) {
  if (!is.null(bio.features.loc)) {
    bio.features.loc <- str_replace(bio.features.loc, "/$", "")
    peakfiles.pat <- ".bed$|.broadPeak$|.narrowPeak$|.gappedPeak$"
    get.files <- list.files(path = bio.features.loc, pattern = peakfiles.pat, full.names = TRUE)
    if(!is.null(search.term)) {
      get.files <- get.files[grepl(search.term, get.files, ignore.case = TRUE)]
    }
    if (length(get.files) >= 1) {
      hg19 <- Seqinfo(genome="hg19")
      bed.list <- lapply(get.files,
                         function(x, s.info) {
                           xf <- fread(input = x, sep = "\t", header = FALSE,
                                       select = c(1:3), skip = "chr",
                                       col.names = c("chr", "start", "end"),
                                       encoding = "UTF-8",
                                       stringsAsFactors = FALSE,
                                       data.table = FALSE, showProgress = FALSE)
                           name <- sub("(.*)\\..*", "\\1", basename(x))
                           xf <- GRanges(seqnames = xf$chr,
                                         ranges = IRanges(start = xf$start + 1L,
                                                          end = xf$end),
                                         strand = "*",
                                         feature = base::rep.int(name, nrow(xf)),
                                         seqinfo = s.info)
                           return(xf)
                         }, s.info = hg19)
      names(bed.list) <- sub("(.*)\\..*", "\\1", basename(get.files))
    } else {
      bed.list <- NULL
      stop("No files matching search term, try a new search term\n",
           "and make sure that the directory contains files ending in\n",
           ".bed, .broadPeak, .narrowPeak, or .gappedPeak")
    }
    if (is.null(bed.list)) {
      return(GRangesList())
    } else {
      return(GRangesList(bed.list))
    }
  } else {
    return(NULL)
  }
}

footprintlookup <- function(query, definition) {
  q <- t(query)
  qr <- matrix(nrow = dim(q)[2])
  for(i in 1:nrow(definition)) {
    d <- definition[i, ]
    qt <- q[!is.na(d), ]
    d <- d[!is.na(d)]
    mn <- dim(qt)
    ## need to find a way to prioritise rules where MSB == 1
    x <- .colSums((bitwAnd(d, 1L) == bitwAnd(qt, 1L)) |
                    (bitwAnd(d, 2L) == 0L &
                       qt == 0L),
                  m = mn[1],
                  n = mn[2],
                  na.rm = FALSE) == length(d)
    qr[x, 1] <- rownames(definition)[i]
  }
  return(qr)
}

#' getFootprint
#'
#' @param directory
#' @param search
#' @param
#' @param translation_layer
#'
#' @return
#' @importFrom GenomicRanges disjoin findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom stringr str_detect
#' @export
#'
#' @examples
PaintStates <- function(directory = NULL, search = NULL, chrome_states = "default",
                         translation_layer = "default") {
  message(paste("searching for", search, Sys.time(), sep = " "))
  if(translation_layer == "default") {
    deflookup <- list(h3k4me1 = "Regulatory",
                      h3k4me3 = "Promoter",
                      h3k27ac = "Active",
                      h3k27me3 = "Polycomb",
                      h3k36me3 = "Transcription",
                      h3k9me3 = "Heterochromatin",
                      h3k9ac = "Active",
                      dnase = "Core",
                      atac = "Core",
                      tf = "Core",
                      ctcf = "CTCF")
  } else {
    if(class(translation_layer) != "list") {
      stop("translation layer must be a list")
    }
    deflookup <- translation_layer
  }
  if(is.null(directory)) {stop("Please supply a directory with your genomic intervals")}
  if(is.null(search)) {warning(paste0("Using all files in the directory: ", directory, "\n",
                                      " to construct a single footprint. Use the 'search'\n",
                                      " argument to restrict")); Sys.sleep(5)}
  if(chrome_states == "default") {
    d <- Cedars.BFG.states
  } else {
    d <- as.matrix(read.table(chrome_states, row.names = 1))
  }
  x <- GetBioFeatures(directory, search)
  # create query set
  x.f <- disjoin(unlist(x))
  if(length(x) == 1) {
    x.f <- x.f[[1]]
  }
  overlaps <- findOverlaps(x.f, x)
  qover <- queryHits(overlaps)
  sover <- subjectHits(overlaps)
  feature.names <- names(x)
  overlaps.info <- data.frame(qhits = qover, feats = feature.names[sover], stringsAsFactors = FALSE)
  overlap.table <- table(overlaps.info)

  # create rule set
  inputset <- tolower(colnames(overlap.table))
  inputset <- sapply(inputset, function(x, y = deflookup) { unname(y[str_detect(x, names(y))]) } )
  inputset.c <- names(inputset)
  names(inputset.c) <- inputset
  nameorder <- inputset.c[colnames(d)[colnames(d) %in% names(inputset.c)]]
  resmatrix <- matrix(overlap.table, nrow = nrow(overlap.table),
                      dimnames = list(rownames(overlap.table),
                                      tolower(colnames(overlap.table))))[, nameorder]
  if(!inherits(resmatrix, "matrix")) {
    dim(resmatrix) <- c(length(resmatrix), 1)
    dimnames(resmatrix) <- list(rownames(overlap.table),
                                tolower(colnames(overlap.table)))
  }
  colnames(resmatrix) <- names(nameorder)
  resmatrix[resmatrix == 0L] <- 2L
  resmatrix[resmatrix == 1L] <- 3L
  missing.data <- colnames(d)[!(colnames(d) %in% names(inputset.c))]
  d.m <- d[, missing.data] > 1
  d.m[is.na(d.m)] <- FALSE
  d <- d[rowSums(d.m) < 1, colnames(d) %in% names(inputset.c)]
  d <- d[order(rowSums(d, na.rm = TRUE), decreasing = FALSE), ]
  mcols(x.f)$footprint <- footprintlookup(resmatrix, d)[, 1]
  return(x.f)
}
