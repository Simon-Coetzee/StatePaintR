
#' @importFrom data.table fread
#' @importFrom stringr str_replace
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom GenomicRanges GRangesList GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom readr read_tsv
GetBioFeatures <- function(manifest = NULL, forcemerge = FALSE) {
  gzipfiles <- manifest[grepl(".gz$", manifest$FILE), ]
  good.files <- sapply(split(manifest, 1:nrow(manifest)), function(x) { file.exists(x$FILE) })
  if(sum(good.files) < nrow(manifest)) {
    stop(paste0("cannot find file: ", manifest[!good.files, "FILE"], "\n"))
  }
  if(length(gzipfiles) >= 1) {
    manifest[grepl(".gz$", manifest$FILE), "FILE"] <- paste(gzipfiles$OPENWITH,
                                                            "<", gzipfiles$FILE)
  }
  get.files <- split(manifest, 1:nrow(manifest))
  names(get.files) <- paste(manifest$SAMPLE, manifest$MARK, sep = "_")
  #if(anyDuplicated(names(get.files))) {
  #  stop("the manifest describes multiple files with the same sample name and mark \n",
  #       "merge these files before running StatePaintR")
  #}
  if (sum(good.files) >= 1) {
    my.seqinfo <- Seqinfo(genome = manifest[1, "BUILD"])
    bed.list <- lapply(get.files,
                       function(x, s.info) {
                         xf <- fread(input = x$FILE, sep = "\t", header = FALSE,
                                     select = c(1:3), skip = "chr",
                                     col.names = c("chr", "start", "end"),
                                     encoding = "UTF-8",
                                     stringsAsFactors = FALSE,
                                     data.table = FALSE, showProgress = FALSE)
                         name <- paste(x$SAMPLE, x$MARK, sep = "_")
                         xf <- GRanges(seqnames = xf$chr,
                                       ranges = IRanges(start = xf$start + 1L,
                                                        end = xf$end),
                                       strand = "*",
                                       feature = base::rep.int(name, nrow(xf)),
                                       seqinfo = s.info)
                         return(xf)
                       }, s.info = my.seqinfo)
  } else {
    stop("no files")
  }
  if (is.null(bed.list)) {
    return(GRangesList())
  } else {
    return(GRangesList(bed.list))
  }
}

parse.manifest <- function(manifest = NULL) {
  manifest.df <- read.table(manifest, stringsAsFactors = FALSE, header = TRUE)
  valid.columns <- colnames(manifest.df) == c("SAMPLE",
                                              "MARK",
                                              "SRC",
                                              "BUILD",
                                              "FILE")
  if(any(!valid.columns)) {
    stop("manifest file must contain the columns ",
         "'SAMPLE', 'MARK', 'SRC', 'BUILD', and 'FILE'")
  }
  manifest.df$OPENWITH <- "DEFAULT"
  if(any(grepl(".gz$", manifest.df$FILE))) {
    gzipper <- switch(Sys.info()[["sysname"]],
                      Linux  = Sys.which("zcat"),
                      Darwin = Sys.which("gzcat"),
                      ...    = c(zcat = ""))
    if(gzipper == "" | (gzipper == "zcat" | gzipper == "gzcat")) {
        stop("some of your files are gzipped, and the (g)zcat command cannot \n",
             "be found on your system")
    } else {
      manifest.df[grepl(".gz$", manifest.df$FILE), "OPENWITH"] <- gzipper
    }
  }
  if(length(unique(sort(manifest.df$BUILD))) > 1) {
    stop("each celltype must be described by one and only one genome build")
  }
  alls <- manifest.df[manifest.df$SAMPLE == "ALL", ]
  manifest.df <- manifest.df[manifest.df$SAMPLE != "ALL", ]
  by.sample <- split(manifest.df, manifest.df$SAMPLE)
  by.sample <- lapply(by.sample, function(out, add = alls) {
    out <- rbind(out, add)
    return(out)
  })
  names(by.sample) <- tolower(names(by.sample))
  return(by.sample)
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

#from toupper documentation
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

#' getFootprint
#'
#' @param directory
#' @param search
#' @param
#' @param translation_layer
#'
#' @return
#' @importFrom GenomicRanges disjoin findOverlaps reduce
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom stringr str_detect coll
#' @export
#'
#' @examples
PaintStates <- function(manifest = NULL, chrome_states = "default",
                        translation_layer = "default") {
  if(inherits(translation_layer, "list")) {
    deflookup <- translation_layer
    names(deflookup) <- tolower(names(deflookup))
  } else {
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
      stop("translation layer must be a list with translations or the string 'default'")
    }
  }
  if(is.null(manifest)) {stop("provide a manifest describing the location of your files \n",
                              "and the mark that was ChIPed")}
  if(chrome_states == "default") {
    d <- Cedars.BFG.states
  } else {
    d <- as.matrix(read.table(chrome_states, row.names = 1))
  }
  if (!is.null(manifest)) {
    samples <- parse.manifest(manifest)
  }
  output <- list()
  for(cell.sample in seq_along(samples)) {
    cell.sample <- samples[[cell.sample]]
    x <- GetBioFeatures(manifest = cell.sample)
    names(x) <- tolower(names(x))
    inputset <- names(x)
    inputset <- sapply(inputset, function(x, y = deflookup) {
      out <- unlist(y[str_detect(x, paste0(names(y), "$"))], use.names = FALSE)
      return(out)
      })
    names(x) <- inputset
    # multiple marks with the same meaning
    to.merge <- inputset[duplicated(inputset)]
    if(length(to.merge) >= 1) {
      to.merge <- unique(unlist(to.merge))
      for(mark in to.merge) {
        x.merge <- unlist(x[names(x) %in% mark])
        x.names <- unique(x.merge$feature)
        x <- x[!(names(x) %in% mark)]
        names(x.merge) <- NULL
        x.merge <- reduce(x.merge)
        mcols(x.merge)$feature <- paste(x.names, collapse = "|")
        x.merge <- GRangesList(x.merge)
        names(x.merge) <- mark
        x <- append(x, x.merge)
        rm(x.merge)
      }
    }
    inputset.c <- names(inputset)
    names(inputset.c) <- inputset
    nameorder <- inputset.c[colnames(d)[colnames(d) %in% names(inputset.c)]]
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
    col.order <- colnames(d)[colnames(d) %in% colnames(overlap.table)]
    resmatrix <- matrix(overlap.table, nrow = nrow(overlap.table),
                        dimnames = list(rownames(overlap.table),
                                        capwords(colnames(overlap.table))))[, col.order]
    if(!inherits(resmatrix, "matrix")) {
      dim(resmatrix) <- c(length(resmatrix), 1)
      dimnames(resmatrix) <- list(rownames(overlap.table),
                                  tolower(colnames(overlap.table)))
    }
    resmatrix[resmatrix == 0L] <- 2L
    resmatrix[resmatrix == 1L] <- 3L
    missing.data <- colnames(d)[!(colnames(d) %in% colnames(resmatrix))]
    if(length(missing.data) > 0){
      d.m <- d[, missing.data] > 1
      d.m[is.na(d.m)] <- FALSE
      if(!inherits(d.m, "matrix")) d.m <- matrix(d.m, ncol = 1, dimnames = list(c(names(d.m)), c("SAMPLE")))
      d <- d[rowSums(d.m) < 1, colnames(d) %in% colnames(resmatrix)]
    }
    d <- d[order(rowSums(d, na.rm = TRUE), decreasing = FALSE), ]
    mcols(x.f)$name <- cell.sample[1, "SAMPLE"]
    mcols(x.f)$state <- footprintlookup(resmatrix, d)[, 1]
    output <- c(output, x.f)
  }
  if(length(samples) > 1) {
    output <- GRangesList(output)
  }
  names(output) <- names(samples)
  attributes(output)$manifest <- samples
  return(output)
}

#' @importFrom dplyr left_join
write.state <- function(x, y, color, file = stdout()) {
  manifest <- y
  meta <- list(software = "StatePaintR",
               version = packageVersion("StatePaintR"),
               files = paste("#",
                             manifest$SAMPLE,
                             manifest$MARK,
                             manifest$FILE,
                             sep = " "))
  if(all(names(mcols(x)) == c("name", "state"))) {
    my.cols <- data.frame(mcols(x), stringsAsFactors = FALSE)
    my.cols <- left_join(my.cols, color, by = c("state" = "STATE"))
    colnames(my.cols) <- c("sample", "name", "itemRgb")
    mcols(x) <- my.cols
  } else {
    # stop("non default column names in input data, meta data may be name and sample")
  }
  file <- file(file, "w+")
  writeLines(paste("# this file was produced by", meta$software), file)
  writeLines(paste("# version number:", meta$version), file)
  writeLines("# it is the chromatin segmentation of the following files: ", file)
  writeLines(meta$files, file)
  my.track <- new("BasicTrackLine",
                  itemRgb = TRUE,
                  db = x@seqinfo@genome[[1]],
                  name = manifest$SAMPLE[[1]],
                  description = paste0("StatePaintR Segmentation for ",
                                       manifest$SAMPLE[[1]]))
  rtracklayer::export.bed(x, file,
                          trackLine = my.track)
  close(file)
}

#' Write StatePaintR object to bedfiles
#' @importFrom rtracklayer export.bed
#' @return
#' @export
#' @example
ExportStatePaintR <- function(states, color.key = Cedars.BFG.colors, output.dir = tempdir()) {
  m.data <- attributes(states)$manifest
  color.value <- strsplit(color.key$COLOR, ",")
  color.key$COLOR <- sapply(color.value, function(x)
    rgb(x[1], x[2], x[3], maxColorValue=255))
  for(state in seq_along(states)) {
    s.name <- names(m.data)[state]
    s.m.data <- m.data[[state]]
    state <- states[[state]]
    write.state(state, s.m.data, color.key,
                file.path(output.dir,
                          paste0(s.name, ".segmentation.bed")))
  }
}