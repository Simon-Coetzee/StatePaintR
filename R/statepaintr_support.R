#' @importFrom data.table fread
#' @importFrom stringr str_replace
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom GenomicRanges GRangesList GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom utils read.table
#' @importFrom readr read_tsv cols_only col_character col_integer col_number
GetBioFeatures <- function(manifest, my.seqinfo) {
  good.files <- sapply(split(manifest, 1:nrow(manifest)), function(x) { file.exists(x$FILE) })
  if(sum(good.files) < nrow(manifest)) {
    stop(paste0("cannot find file: ", manifest[!good.files, "FILE"], "\n"))
  }
  get.files <- split(manifest, 1:nrow(manifest))
  names(get.files) <- paste(manifest$SAMPLE, manifest$MARK, sep = "_")
  if(anyDuplicated(names(get.files))) {
    stop("the manifest describes multiple files with the same sample name and mark \n",
         "merge these files before running StatePaintR")
  }
  if (sum(good.files) >= 1) {
    bed.list <- lapply(get.files,
                       function(x, s.info) {
                         if(grepl(".narrowPeak", x$FILE, ignore.case = TRUE)) {
                           col.types <- cols_only(chr = col_character(),
                                                  start = col_integer(),
                                                  end = col_integer(),
                                                  name = col_character(),
                                                  score = col_integer(),
                                                  strand = col_character(),
                                                  signalValue = col_number())
                           col.names <- c("chr", "start", "end", "name", "score", "strand", "signalValue")
                           col.numbers <- c(1:7)
                         } else {
                           col.types <- cols_only(chr = col_character(),
                                                  start = col_integer(),
                                                  end = col_integer())
                           col.names <- c("chr", "start", "end")
                           col.numbers <- c(1:3)
                         }
                         if(any(grepl(".gz$", x$FILE))) {
                           xf <- suppressWarnings(read_tsv(file = x$FILE, col_names = col.names,
                                                           col_types = col.types,
                                                           progress = FALSE))
                         } else {
                           xf <- fread(input = x$FILE, sep = "\t", header = FALSE,
                                       select = col.numbers, skip = "chr",
                                       col.names = col.names,
                                       encoding = "UTF-8",
                                       stringsAsFactors = FALSE,
                                       data.table = FALSE, showProgress = FALSE)
                         }
                         name <- paste(x$SAMPLE, x$MARK, sep = "_")
                         if(ncol(xf) < 7) xf$signalValue <- NA
                         xf <- GRanges(seqnames = xf$chr,
                                       ranges = IRanges(start = xf$start + 1L,
                                                        end = xf$end),
                                       strand = "*",
                                       feature = base::rep.int(name, nrow(xf)),
                                       signalValue = xf$signalValue,
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

parse.manifest <- function(manifest) {
  if(!file.exists(manifest)) stop("manifest file: ", manifest, " does not exist")
  manifest.df <- read.table(manifest, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  valid.columns <- colnames(manifest.df) == c("SAMPLE",
                                              "MARK",
                                              "SRC",
                                              "BUILD",
                                              "FILE")
  if(any(!valid.columns)) {
    stop("manifest file must contain the columns ",
         "'SAMPLE', 'MARK', 'SRC', 'BUILD', and 'FILE'")
  }
  manifest.path <- dirname(manifest)
  files.good <- file.exists(manifest.df$FILE)
  if(any(!files.good)) {
    test.files <- manifest.df[!files.good, "FILE"]
    test.files <- file.path(manifest.path, test.files)
    new.files.good <- file.exists(test.files)
    if(any(!new.files.good)) {
      stop("cannot find files listed in manifest")
    } else {
      manifest.df[!files.good, "FILE"] <- test.files
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
  qr <- matrix(nrow = dim(query)[2])
  if(any(query == 0)) {
    for(i in 1:nrow(definition)) {
      d <- definition[i, ]
      qt <- query[!is.na(d), ]
      d <- d[!is.na(d)]
      x <- base::colSums((bitwAnd(d, 1L) == bitwAnd(qt, 1L)) |
                           (bitwAnd(d, 2L) == 0L &
                              qt == 0L),
                         na.rm = FALSE) == length(d)
      qr[x, 1] <- rownames(definition)[i]
    }
  } else {
    for(i in 1:nrow(definition)) {
      d <- definition[i, ]
      qt <- query[!is.na(d)]
      d <- d[!is.na(d)]
      x <- base::.colSums(bitwAnd(d, 1L) == bitwAnd(qt, 1L),
                          m = length(d), n = length(qt)/length(d),
                          na.rm = FALSE) == length(d)
      qr[x, 1] <- rownames(definition)[i]
    }
  }
  return(qr)
}

## from toupper documentation
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

produce.state <- function(state.obj, position = 1) {
  state <- state.obj[[position]]
  order <- sapply(state$features, function(x) {x$order})
  value <- sapply(state$features, function(x) {x$score})
  name <- sapply(state$features, function(x) {x$name})
  output <- value[order]
  output[output == -1] <- NA
  names(output) <- name[order]
  return(output)
}

zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

make.overlap.matrix <- function(query.over, subject.over, samples) {
  overlaps <- data.frame(qhits = query.over,
                         samples = samples[subject.over],
                         stringsAsFactors = FALSE)
  overlap.table <- table(overlaps)
  output.matrix <- matrix(overlap.table,
                          nrow = nrow(overlap.table),
                          dimnames = list(rownames(overlap.table),
                                          colnames(overlap.table)))
  if(!is(output.matrix, "matrix")) {
    dim(output.matrix) <- c(length(output.matrix), 1)
    dimnames(output.matrix) <- list(rownames(overlap.table),
                                    colnames(overlap.table))
  }
  return(output.matrix)
}

#' @importFrom utils packageVersion
write.state <- function(x, y, color, hub.id, file = stdout()) {
  manifest <- y
  meta <- list(software = "StatePaintR",
               version = packageVersion("StatePaintR"),
               files = paste("#",
                             manifest$SAMPLE,
                             manifest$MARK,
                             manifest$FILE,
                             sep = " "))
  if(all(names(mcols(x)) == c("name", "state", "score"))) {
    if(!is.null(color)) {
      my.cols <- data.frame(mcols(x), stringsAsFactors = FALSE)
      my.cols <- left_join(my.cols, color, by = c("state" = "STATE"))
      colnames(my.cols) <- c("sample", "name", "score", "itemRgb")
      mcols(x) <- my.cols
    } else {
      names(mcols(x)) <- c("sample", "name", "score")
    }
  } else {
    stop("non default column names in input data, meta data may be name and sample")
  }
  file <- file(file, "w+")
  writeLines(paste("# this file was produced by", meta$software), file)
  writeLines(paste("# version number:", meta$version), file)
  writeLines(paste("# StateHub Model ID:", hub.id), file)
  writeLines("# it is the chromatin segmentation of the following files: ", file)
  writeLines(meta$files, file)
  if(!is.null(color)) {
    my.track <- new("BasicTrackLine",
                    itemRgb = TRUE,
                    db = x@seqinfo@genome[[1]],
                    name = manifest$SAMPLE[[1]],
                    description = paste0("StatePaintR Segmentation for ",
                                         manifest$SAMPLE[[1]]))
  } else {
    my.track <- new("BasicTrackLine",
                    db = x@seqinfo@genome[[1]],
                    name = manifest$SAMPLE[[1]],
                    description = paste0("StatePaintR Segmentation for ",
                                         manifest$SAMPLE[[1]]))
  }
  rtracklayer::export.bed(x, file,
                          trackLine = my.track)
  close(file)
}

reverse_tl <- function(tl) {
  values <- sapply(tl, length)
  values <- rep(names(values), values)
  keys <- unlist(tl, use.names = FALSE)
  names(values) <- keys
  return(as.list(values))
}


setMethod("show",
          signature = signature(object = "decisionMatrix"),
          function(object) {
            cat(" ID:", object@id, "\n",
                "Name:", object@name, "\n",
                "Author:", object@author, "\n",
                "Revision:", object@revision, "\n",
                "Description:", object@description, "\n")
          })

