
#' @importFrom data.table fread
#' @importFrom stringr str_replace
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom GenomicRanges GRangesList GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom readr read_tsv cols_only col_character col_integer
GetBioFeatures <- function(manifest, my.seqinfo, forcemerge = FALSE) {
  # gzipfiles <- manifest[grepl(".gz$", manifest$FILE), ]
  good.files <- sapply(split(manifest, 1:nrow(manifest)), function(x) { file.exists(x$FILE) })
  if(sum(good.files) < nrow(manifest)) {
    stop(paste0("cannot find file: ", manifest[!good.files, "FILE"], "\n"))
  }
  # if(length(gzipfiles) >= 1) {
  #   manifest[grepl(".gz$", manifest$FILE), "FILE"] <- paste(gzipfiles$OPENWITH,
  #                                                           "<", gzipfiles$FILE)
  # }
  get.files <- split(manifest, 1:nrow(manifest))
  names(get.files) <- paste(manifest$SAMPLE, manifest$MARK, sep = "_")
  #if(anyDuplicated(names(get.files))) {
  #  stop("the manifest describes multiple files with the same sample name and mark \n",
  #       "merge these files before running StatePaintR")
  #}
  if (sum(good.files) >= 1) {
    bed.list <- lapply(get.files,
                       function(x, s.info) {
                         if(any(grepl(".gz$", x$FILE))) {
                           xf <- suppressWarnings(read_tsv(file = x$FILE, col_names = c("chr", "start", "end"),
                                                           col_types = cols_only(chr = col_character(), start = col_integer(), end = col_integer()),
                                                           progress = FALSE))
                         } else {
                           xf <- fread(input = x$FILE, sep = "\t", header = FALSE,
                                       select = c(1:3), skip = "chr",
                                       col.names = c("chr", "start", "end"),
                                       encoding = "UTF-8",
                                       stringsAsFactors = FALSE,
                                       data.table = FALSE, showProgress = FALSE)
                         }
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
  # manifest.df$OPENWITH <- "DEFAULT"
  # if(any(grepl(".gz$", manifest.df$FILE))) {
  #   gzipper <- switch(Sys.info()[["sysname"]],
  #                     Linux  = Sys.which("zcat"),
  #                     Darwin = Sys.which("gzcat"),
  #                     ...    = c(zcat = ""))
  #   if(gzipper == "" | (gzipper == "zcat" | gzipper == "gzcat")) {
  #       stop("some of your files are gzipped, and the (g)zcat command cannot \n",
  #            "be found on your system")
  #   } else {
  #     manifest.df[grepl(".gz$", manifest.df$FILE), "OPENWITH"] <- gzipper
  #   }
  # }
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

#' getFootprint
#'
#' @param manifest
#' @param decisionMatrix
#'
#' @return
#' @importFrom GenomicRanges disjoin findOverlaps reduce
#' @importFrom S4Vectors from to
#' @importFrom stringr str_detect coll str_replace
#' @export
#'
#' @examples
PaintStates <- function(manifest, decisionMatrix, progress = TRUE) {
  start.time <- Sys.time()
  if(missing(manifest)) {stop("provide a manifest describing the location of your files \n",
                              "and the mark that was ChIPed")}
  if(missing(decisionMatrix)) {stop("provide a decisionMatrix object")}
  if(is(decisionMatrix, "decisionMatrix")) {
    deflookup <- reverse_tl(abstractionLayer(decisionMatrix))
    names(deflookup) <- tolower(names(deflookup))
    d <- decisionMatrix(decisionMatrix)
  } else {
    stop("arg: decisionMatrix must be object of class decisionMatrix")
  }
  samples <- parse.manifest(manifest)
  sample.genomes <- as.list(unique(unlist(sapply(samples, "[", "BUILD"))))
  sample.genomes.names <- sample.genomes
  sample.genomes <- lapply(sample.genomes.names, function(x) Seqinfo(genome = x))
  names(sample.genomes) <- sample.genomes.names
  output <- list()
  total <- 20
  pb <- txtProgressBar(min = 0, max = length(samples), style = 3)
  for(cell.sample in seq_along(samples)) {
    setTxtProgressBar(pb, cell.sample)
    cell.sample <- samples[[cell.sample]]
    x <- GetBioFeatures(manifest = cell.sample, my.seqinfo = sample.genomes[[cell.sample[1, "BUILD"]]])
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
    dontbreak <- str_replace(inputset[grep(pattern = "^\\*", inputset)], "^\\*", "")
    inputset <- str_replace(inputset, "^\\*", "")
    if(length(dontbreak) > 0) {
      x.f <- x.f[x.f %outside% x[dontbreak]]
      unfrag <- reduce(sort(do.call(c, list(unlist(x[dontbreak]), ignore.mcols = TRUE))))
      x.f <- c(x.f, unfrag, ignore.mcols = TRUE)
      x.f <- sort(x.f)
    }
    overlaps <- findOverlaps(x.f, x)
    resmatrix <- make.overlap.matrix(from(overlaps), to(overlaps), names(x))
    resmatrix <- resmatrix[, colnames(d)[colnames(d) %in% colnames(resmatrix)]]
    resmatrix[resmatrix == 0L] <- 2L
    resmatrix[resmatrix == 1L] <- 3L
    missing.data <- colnames(d)[!(colnames(d) %in% colnames(resmatrix))]
    if(length(missing.data) > 0){
      d.m <- d[, missing.data] > 1
      d.m[is.na(d.m)] <- FALSE
      if(!inherits(d.m, "matrix")) d.m <- matrix(d.m, ncol = 1, dimnames = list(c(names(d.m)), c("SAMPLE")))
      d <- d[rowSums(d.m) < 1, colnames(d) %in% colnames(resmatrix)]
    }
    dl <- d[order(rowSums(d, na.rm = TRUE), decreasing = FALSE), ]
    mcols(x.f)$name <- cell.sample[1, "SAMPLE"]
    resmatrix <- t(resmatrix);
    mcols(x.f)$state <- footprintlookup(resmatrix, dl)[, 1]
    x.f.l <- split(x.f, x.f$state)
    x.f.l <- lapply(x.f.l, reduce)
    for(state.name in names(x.f.l)) {
      mcols(x.f.l[[state.name]])$name <- cell.sample[1, "SAMPLE"]
      mcols(x.f.l[[state.name]])$state <- state.name
    }
    x.f <- sort(do.call(c, unlist(x.f.l, use.names = FALSE)))
    output <- c(output, x.f)
  }
  if(length(samples) > 1) {
    output <- GRangesList(output)
  }
  close(pb)
  names(output) <- names(samples)
  attributes(output)$manifest <- samples
  done.time <- Sys.time() - start.time
  message("processed ", length(samples), " in ", round(done.time, digits = 2), " ", attr(done.time, "units"))
  return(output)
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
    if(!is.null(color)) {
      my.cols <- data.frame(mcols(x), stringsAsFactors = FALSE)
      my.cols <- left_join(my.cols, color, by = c("state" = "STATE"))
      colnames(my.cols) <- c("sample", "name", "itemRgb")
      mcols(x) <- my.cols
    } else {
      names(mcols(x)) <- c("sample", "name")
    }
  } else {
    stop("non default column names in input data, meta data may be name and sample")
  }
  file <- file(file, "w+")
  writeLines(paste("# this file was produced by", meta$software), file)
  writeLines(paste("# version number:", meta$version), file)
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

#' Write StatePaintR object to bedfiles
#'
#' @param states
#' @param decisionMatrix
#' @param output.dir
#'
#' @importFrom rtracklayer export.bed
#' @importFrom dplyr left_join
#' @return
#' @export
#' @example
ExportStatePaintR <- function(states, decisionMatrix, output.dir) {
  if(missing(output.dir)) { stop("please indicate output directory") }
  if(missing(decisionMatrix)) { stop("please include decisionMatrix") }
  if(!dir.exists(output.dir)) { dir.create(output.dir) }
  m.data <- attributes(states)$manifest
  color.key <- stateColors(decisionMatrix)
  pb <- txtProgressBar(min = 0, max = length(states), style = 3)
  for(state in seq_along(states)) {
    setTxtProgressBar(pb, state)
    s.name <- names(m.data)[state]
    s.m.data <- m.data[[state]]
    state <- states[[state]]
    write.state(state, s.m.data, color.key,
                file.path(output.dir,
                          paste0(s.name, ".", nrow(s.m.data), "mark.segmentation.bed")))
  }
  close(pb)
  message("segmentation files written to: ", output.dir)
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

#' Retrieve Decision Matrix from StateHub
#'
#' @param search character string of either unique id, or search term.
#'
#' @return decisionMatrix object
#' @importFrom httr content GET http_error modify_url
#' @importFrom jsonlite toJSON fromJSON
#' @importFrom stringr str_extract
#' @export
#'
#' @examples
get.decision.matrix <- function(search) {
  if(missing(search)) { stop("missing search argument") }
  stateHub <- modify_url("http://statehub.org/statehub", path = "statehub/getmodel.jsp")
  query <- GET(stateHub, query = list(id = search))
  if(http_error(query)) {
    query <- GET(stateHub, query = list(search = search))
    if(http_error(query)) { stop("site not availible") }
    if(length(query) < 1) { stop("no search results found") }
  }
  query <- content(query)
  for(query.i in seq_along(query)) {
    query.result <- query[[query.i]]
    state.names <- sapply(query.result$states, function(x) {x$name})
    for(state in seq_along(query.result$states)) {
      if(state == 1) {
        decision.matrix <- produce.state(query.result$states, state)
      } else {
        decision.matrix <- rbind(decision.matrix, produce.state(query.result$states, state))
      }
    }
    rownames(decision.matrix) <- state.names
    state.colors <- sapply(query.result$states, "[[", "format")
    state.colors <- str_extract(state.colors, "#([a-f]|[A-F]|[0-9]){3}(([a-f]|[A-F]|[0-9]){3})?\\b")
    names(state.colors) <- state.names
    decision.matrix <- new("decisionMatrix",
                           id = query.result[["id"]],
                           name = query.result[["name"]],
                           author = query.result[["author"]],
                           revision = query.result[["revision"]],
                           description = query.result[["description"]],
                           abstraction.layer = fromJSON(toJSON(query.result[["translation"]], auto_unbox = TRUE)),
                           decision.matrix = decision.matrix,
                           state.colors = state.colors)
    if(query.i > 1) {
      output <- c(output, decision.matrix)
    } else {
      output <- decision.matrix
    }
  }
  return(output)
}

reverse_tl <- function(tl) {
  values <- sapply(tl, length)
  values <- rep(names(values), values)
  keys <- unlist(tl, use.names = F)
  names(values) <- keys
  return(as.list(values))
}

setClass(Class = "decisionMatrix",
         slots = list(id = "character",
                      name = "character",
                      author = "character",
                      revision = "character",
                      description = "character",
                      abstraction.layer = "list",
                      decision.matrix = "matrix",
                      state.colors = "character"))
setMethod("show",
          signature = signature(object = "decisionMatrix"),
          function(object) {
            cat(" ID:", object@id, "\n",
                "Name:", object@name, "\n",
                "Author:", object@author, "\n",
                "Revision:", object@revision, "\n",
                "Description:", object@description, "\n")
          })
setGeneric("decisionMatrix", function(object) standardGeneric("decisionMatrix"))
#' Extract Descision Matrix from descisionMatrix object
#'
#' @param object decisionMatrix.
#'
#' @return matrix descisionMatrix stripped of metadata
#' @export
#'
#' @examples
setMethod("decisionMatrix",
          signature = signature(object = "decisionMatrix"),
          function(object) {
            object@decision.matrix
          })
setGeneric("abstractionLayer", function(object) standardGeneric("abstractionLayer"))
#' Extract abstraction layer from descisionMatrix object
#'
#' @param object decisionMatrix.
#'
#' @return list describing relationship between chromatin marks and functional categories
#' @export
#'
#' @examples
setMethod("abstractionLayer",
          signature = signature(object = "decisionMatrix"),
          function(object) {
            object@abstraction.layer
          })
setGeneric("stateColors", function(object) standardGeneric("stateColors"))
#' Extract color information from descisionMatrix object
#'
#' @param object decisionMatrix.
#'
#' @return named vector of colors for chromatin states
#' @export
#'
#' @examples
setMethod("stateColors",
          signature = signature(object = "decisionMatrix"),
          function(object) {
            data.frame(STATE = names(object@state.colors), COLOR = object@state.colors, stringsAsFactors = FALSE)
          })