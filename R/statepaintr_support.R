#' @importFrom data.table fread
#' @importFrom stringr str_replace
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom GenomicRanges GRangesList GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom utils read.table
#' @importFrom readr read_tsv cols_only col_character col_integer col_number
GetBioFeatures <- function(manifest, dm, my.seqinfo) {
  dm@abstraction.layer <- lapply(abstractionLayer(dm), tolower)
  manifest <- manifest[tolower(manifest$MARK) %in% unlist(abstractionLayer(dm), use.names = FALSE), , drop = FALSE]
  if(nrow(manifest) < 1) return(NULL)
  good.files <- sapply(split(manifest, 1:nrow(manifest)), function(x) file.exists(x$FILE))
  if (sum(good.files) < nrow(manifest)) {
    stop(paste0("cannot find file: ", manifest[!good.files, "FILE"], "\n"))
  }
  get.files <- split(manifest, 1:nrow(manifest))
  names(get.files) <- paste(manifest$SAMPLE, manifest$MARK, sep = "_")
  # if (anyDuplicated(names(get.files))) {
  #   stop("the manifest describes multiple files with the same sample name and mark \n",
  #        "merge these files before running StatePaintR")
  # }
  if (sum(good.files) >= 1) {
    bed.list <- lapply(get.files,
                       function(x, s.info) {
                         if (grepl(".narrowPeak", x$FILE, ignore.case = TRUE)) {
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
                         if (any(grepl(".gz$", x$FILE))) {
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
                         if (ncol(xf) < 7) xf$signalValue <- NA
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

parse.manifest <- function(manifest, check = TRUE) {
  if (is(manifest, "character")) {
    if (check == TRUE) {
      if (!file.exists(manifest)) stop("manifest file: ", manifest, " does not exist")
    }
    col.names <- c("SAMPLE", "MARK", "SRC", "BUILD", "FILE")
    manifest.df <- fread(file = manifest,
                         col.names = col.names,
                         stringsAsFactors = FALSE,
                         data.table = FALSE,
                         showProgress = FALSE)
    manifest.path <- dirname(manifest)
  } else {
    if (!is(manifest, "data.frame")) stop("manifest must be a character vector containing one path name, or a data.frame")
    manifest.df <- manifest
    manifest.path <- ""
  }
  files.good <- file.exists(manifest.df$FILE)
  if (check == TRUE) {
    if (any(!files.good)) {
      test.files <- manifest.df[!files.good, "FILE"]
      test.files <- file.path(manifest.path, test.files)
      new.files.good <- file.exists(test.files)
      if (any(!new.files.good)) {
        stop("cannot find these files listed in manifest:\n",
             paste0(manifest.df[!new.files.good, "FILE"], collapse = "\n"))
      } else {
        manifest.df[!files.good, "FILE"] <- test.files
      }
    }
  }
  if (length(unique(sort(manifest.df$BUILD))) > 1) {
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

flookup <- function(query, definition) {
  qr <- matrix(nrow = dim(query)[2])
  for (i in 1:nrow(definition)) {
    d <- definition[i, , drop = FALSE]
    keepnames <- which(d != 1L)
    d <- d[keepnames]
    qt <- query[keepnames, , drop = FALSE]
    x <- colSums((d == qt | (d == 0L & qt != 3))) == length(d)
    qr[x, 1] <- rownames(definition)[i]
  }
  return(qr)
}

footprintlookup <- function(query, definition) {
  qr <- matrix(nrow = dim(query)[2])
  if (any(query == 0)) {
    for (i in 1:nrow(definition)) {
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
    for (i in 1:nrow(definition)) {
      d <- definition[i, ]
      qt <- query[!is.na(d)]
      d <- d[!is.na(d)]
      x <- base::.colSums(bitwAnd(d, 1L) == bitwAnd(qt, 1L),
                          m = length(d), n = length(qt) / length(d),
                          na.rm = FALSE) == length(d)
      qr[x, 1] <- rownames(definition)[i]
    }
  }
  return(qr)
}

## from toupper documentation
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if (strict) tolower(s) else s},
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
                         stringsAsFactors = TRUE)
  overlap.table <- table(overlaps)
  output.matrix <- matrix(overlap.table,
                          nrow = nrow(overlap.table),
                          dimnames = list(rownames(overlap.table),
                                          colnames(overlap.table)))
  if (!is(output.matrix, "matrix")) {
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
  if (all(names(mcols(x)) == c("name", "state", "score"))) {
    if (!is.null(color)) {
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
  if (!is.null(color)) {
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

#' @importFrom jsonlite toJSON fromJSON unbox
ExportStateHub <- function(states, decisionMatrix, output.dir, as = "json") {
  if (missing(output.dir)) { stop("please indicate output directory") }
  if (missing(decisionMatrix)) { stop("please include decisionMatrix") }
  if (!dir.exists(output.dir)) { dir.create(output.dir) }
  if(inherits(states, "GRangesList")) {
    m.data <- attributes(states)$manifest
  } else {
    m.data <- parse.manifest(states)
  }
  decisionMatrix@abstraction.layer <- lapply(abstractionLayer(decisionMatrix), tolower)
  m.data <- lapply(m.data, function(manifest, dm = decisionMatrix) {
    manifest <- manifest[tolower(manifest$MARK) %in% unlist(abstractionLayer(dm), use.names = FALSE), , drop = FALSE]
    if(nrow(manifest) < 1) {
      return(NULL)
    } else {
      return(manifest)
    }
  })
  m.data <- m.data[!sapply(m.data, is.null)]
  df <- data.frame(name = names(m.data),
                   description = names(m.data),
                   project = sapply(m.data, function(x) {x[1, "SRC"]}),
                   genome = sapply(m.data, function(x) {x[1, "BUILD"]}),
                   marks = NA,
                   bedFileName = paste0(names(m.data), ".", sapply(m.data, length), "mark.segmentation.bed"),
                   bigBedFileName = paste0(names(m.data), ".", sapply(m.data, length), "mark.segmentation.bb"),
                   "statePaintRVersion" = NA,
                   "modelID" = NA,
                   baseURL = NA,
                   Order = NA,
                   check.names = FALSE,
                   stringsAsFactors = FALSE)
  df$marks <- sapply(m.data, function(x) {x[, "MARK"]})
  df$statePaintRVersion <- as.character(packageVersion("StatePaintR"))
  df$modelID <- decisionMatrix@id
  df$baseURL <- "http://s3-us-west-2.amazonaws.com/statehub-trackhub/tracks/"
  df$order <- 0
  row.names(df) <- NULL
  if (as == "json") {
    jlist <- list(modelID = unbox(decisionMatrix@id),
                  tracks = df)
    df.json <- toJSON(jlist, pretty = TRUE)
    writeLines(df.json, con = file.path(output.dir, "manifest.json"))
    return(invisible(jlist))
  } else if (as == "autosql") {
	  for () {
	  }
  }
}

setMethod("show",
          signature = signature(object = "decisionMatrix"),
          function(object) {
            cat(" ID:", object@id, "\n",
                "Name:", object@name, "\n",
                "Author:", object@author, "\n",
                "Revision:", object@revision, "\n",
                "Description:", object@description, "\n")
            if (any(grepl("\\*$", names(abstractionLayer(object))))) {
              cat(" Not Scoring:", grep(pattern = "\\*$",
                                        x = names(abstractionLayer(object)),
                                        value = TRUE), "\n")
            }
            if (any(grepl("^\\[.*\\]$|^\\[.*\\]\\*$", names(abstractionLayer(object))))) {
              cat(" Not Splitting:", grep(pattern = "^\\[.*\\]$|^\\[.*\\]\\*$",
                                          x = names(abstractionLayer(object)),
                                          value = TRUE), "\n")
            }
          })


PaintStatesBenchmark <- function(manifest, decisionMatrix, scoreStates = FALSE, progress = TRUE) {
  start.time <- Sys.time()
  if (missing(manifest)) {stop("provide a manifest describing the location of your files \n",
                               "and the mark that was investigated")}
  if (missing(decisionMatrix)) {stop("provide a decisionMatrix object")}
  if (is(decisionMatrix, "decisionMatrix")) {
    deflookup <- reverse_tl(abstractionLayer(decisionMatrix))
    names(deflookup) <- tolower(names(deflookup))
    d <- decisionMatrix(decisionMatrix)
  } else {
    stop("arg: decisionMatrix must be object of class decisionMatrix")
  }
  samples <- parse.manifest(manifest)
  sample.genomes <- as.list(unique(unlist(sapply(samples, "[", "BUILD"))))
  sample.genomes.names <- sample.genomes
  sample.genomes.names <- lapply(sample.genomes.names, function(x) {x <- str_replace(tolower(x), "grch38", "hg38")})
  do.seqinfo <- try(Seqinfo(genome = sample.genomes.names[[1]]))
  if (inherits(do.seqinfo, "try-error")) {
    sample.genomes <- lapply(sample.genomes.names, function(x) NULL)
    warning("could not download seqinfo for genome")
  } else {
    sample.genomes <- lapply(sample.genomes.names, function(x) Seqinfo(genome = x))
  }
  names(sample.genomes) <- sample.genomes.names
  output <- list()
  skipped <- list()
  lc <- 0
  if (progress) pb <- txtProgressBar(min = 0, max = length(samples) * 3, style = 3)
  for (cell.sample in seq_along(samples)) {
    lc <- lc + 1
    if (progress) setTxtProgressBar(pb, lc)
    skipname <- cell.sample
    cell.sample <- samples[[cell.sample]]
    x <- GetBioFeatures(manifest = cell.sample, my.seqinfo = sample.genomes[[cell.sample[1, "BUILD"]]])
    names(x) <- tolower(names(x))
    inputset <- names(x)
    inputset <- sapply(inputset, function(x, y = deflookup) {
      out <- unlist(y[str_detect(x, paste0(names(y), "$"))], use.names = FALSE)
      return(out)
    })
    notInTranslationLayer <- sapply(inputset, is.null)
    if(all(sapply(inputset, is.null))) {
      warning(cell.sample[1, "SAMPLE"], " has no MARKs that are present in the translation layer \n",
              " MARKs included are ", paste(cell.sample[, "MARK"], collapse = " "))
      skipped <- c(skipped, skipname)
      next()
    }
    x <- x[!notInTranslationLayer]
    inputset <- inputset[!notInTranslationLayer]
    names(x) <- inputset
    to.merge <- inputset[duplicated(inputset)]
    if (length(to.merge) >= 1) {
      to.merge <- unique(unlist(to.merge))
      for (mark in to.merge) {
        x.merge <- unlist(x[names(x) %in% mark])
        x.names <- unique(x.merge$feature)
        x <- x[!(names(x) %in% mark)]
        names(x.merge) <- NULL
        if (scoreStates) {
          x.merge <- binnedAverage(reduce(x.merge), coverage(x.merge, weight = "signalValue"), varname = "signalValue")
          mcols(x.merge)$feature <- paste(x.names, collapse = "|")
          mcols(x.merge) <- mcols(x.merge)[, c(2,1)]
        } else {
          x.merge <- reduce(x.merge)
          mcols(x.merge)$feature <- paste(x.names, collapse = "|")
          mcols(x.merge)$signalValue <- NA
        }
        x.merge <- GRangesList(x.merge)
        names(x.merge) <- mark
        x <- append(x, x.merge)
        rm(x.merge)
        inputset <- inputset[inputset != mark]
        inputset <- c(inputset, merged = mark)
      }
    }

    inputset.c <- names(inputset)
    names(inputset.c) <- inputset
    x.f <- disjoin(unlist(x))

    dontuseScore.i <- grep(pattern = "\\*$", inputset)
    inputset <- str_replace(inputset, "\\*$", "")
    dontbreak.i <- grep(pattern = "^\\[.*\\]$", inputset)
    inputset <- str_replace(inputset, "^\\[", "")
    inputset <- str_replace(inputset, "\\]$", "")

    dontbreak <- inputset[dontbreak.i]
    dontuseScore <- inputset[dontuseScore.i]

    names(x) <- inputset
    if (length(dontbreak) > 0) {
      x.f <- x.f[x.f %outside% x[dontbreak]]
      unfrag <- reduce(sort(do.call(c, list(unlist(x[dontbreak]), ignore.mcols = TRUE))))
      x.f <- c(x.f, unfrag, ignore.mcols = TRUE)
      x.f <- sort(x.f)
    }
    overlaps <- findOverlaps(x.f, x)
    resmatrix <- make.overlap.matrix(from(overlaps), to(overlaps), names(x))
    resmatrix <- resmatrix[, colnames(d)[colnames(d) %in% colnames(resmatrix)], drop = FALSE]
    if (scoreStates) {
      scorematrix <- resmatrix
      signalCol <- sapply(x, function(x) !is.na(mcols(x)[1, "signalValue"]))
      signalCol <- signalCol[!(names(signalCol) %in% dontuseScore)]
      signalCol <- which(signalCol[colnames(scorematrix)])
      for (feature in names(signalCol)) {
        scoreOverlaps <- findOverlaps(x.f, x[[feature]])
        scorematrix[from(scoreOverlaps), feature] <- mcols(x[[feature]])[to(scoreOverlaps), "signalValue"]
        scorematrix[, feature] <- frankv(scorematrix[, feature], ties.method = "average")
      }
      scorematrix <- scorematrix[, names(signalCol), drop = FALSE]
    }

    resmatrix <- resmatrix + 2L

    missing.data <- colnames(d)[!(colnames(d) %in% colnames(resmatrix))]
    if (any(is.na(d))) {
      d[is.na(d)] <- 1L
    }

    lmd <- length(missing.data)
    if (lmd > 0) {
      dl <- d[-which(matrix(bitwAnd(d[, missing.data, drop = FALSE], 2L) == 2L, ncol = lmd), arr.ind = TRUE)[,1], -which(colnames(d) %in% missing.data), drop = FALSE]
    } else {
      dl <- d
    }
    d.order <- dl
    d.order[dl == 1] <- 0
    d.order[dl == 0] <- 1
    dl <- dl[order(rowSums(d.order, na.rm = TRUE), decreasing = FALSE), , drop = FALSE]
    resmatrix <- t(resmatrix)
    mcols(x.f)$name <- cell.sample[1, "SAMPLE"]
    lc <- lc + 1
    if (progress) setTxtProgressBar(pb, lc)
    if (scoreStates) {
      dl.score <- dl[, names(signalCol)]
      score.cells <- which(dl.score == 3L, arr.ind = TRUE)
      if (length(signalCol) == 1) {
        score.cells <- matrix(score.cells, ncol = 1, dimnames = list(names(score.cells), c("row")))
        score.cells <- cbind(score.cells, col = 1)
      }
      segments <- data.frame(state = flookup(resmatrix, dl)[, 1], median = NA, mean = NA, max = NA, stringsAsFactors = FALSE)
      score.features <- unique(segments$state)[(unique(segments$state) %in% rownames(score.cells))]
      for (score.feature in score.features) {
        scorematrix.rows <- segments$state == score.feature
        feature.scores <- scorematrix[segments$state == score.feature,
                                      score.cells[rownames(score.cells) %in% score.feature, "col"],
                                      drop = FALSE]
        feature.scores.df <- data.frame(median = numeric(), mean = numeric(), max = numeric())
        if (is(feature.scores, "matrix")) {
          feature.scores.df <- rbind.data.frame(feature.scores.df, data.frame(median = rowMedians(feature.scores), mean = rowMeans(feature.scores), max = rowMaxs(feature.scores)))
        } else {
          feature.scores.df <- rbind.data.frame(feature.scores.df, data.frame(median = feature.scores / 2, mean = feature.scores / 2, max = feature.scores))
        }
        segments[segments$state == score.feature, c("median", "mean", "max")] <- feature.scores.df
      }
      mcols(x.f)$state <- segments$state
      mcols(x.f)[, c("median", "mean", "max")] <- DataFrame(segments[, c("median", "mean", "max")])
    } else {
      mcols(x.f)$state <- flookup(resmatrix, dl)[, 1]
    }
    x.f.l <- split(x.f, x.f$state)
    if (scoreStates) {
      x.f.ll <- lapply(x.f.l, reduce, with.revmap = TRUE)
      for (state.name in names(x.f.ll)) {
        if (state.name %in% score.features) {
          score.feature <- state.name
          revmap <- mcols(x.f.ll[[score.feature]])$revmap
          revmap <- sapply(revmap, function(x) {x[1]})
          my.scores <- mcols(x.f.l[[score.feature]])[revmap, c("median", "mean", "max")]
        } else {
          my.scores <- DataFrame(median = rep.int(0, length(x.f.ll[[state.name]])), mean = rep.int(0, length(x.f.ll[[state.name]])), max = rep.int(0, length(x.f.ll[[state.name]])))
        }
        mcols(x.f.ll[[state.name]])$revmap <- NULL
        mcols(x.f.ll[[state.name]])$name <- cell.sample[1, "SAMPLE"]
        mcols(x.f.ll[[state.name]])$state <- state.name

        mcols(x.f.ll[[state.name]])[, c("median", "mean", "max")] <- lapply(my.scores, function(x) {
          if(max(x) > 0) {
            ceiling((x/max(x)) * 1000)
          } else {
            0
          }
        })
      }
      x.f.l <- x.f.ll
    } else {
      x.f.l <- lapply(x.f.l, reduce)
      for (state.name in names(x.f.l)) {
        mcols(x.f.l[[state.name]])$name <- cell.sample[1, "SAMPLE"]
        mcols(x.f.l[[state.name]])$state <- state.name
        mcols(x.f.l[[state.name]])$score <- 0
      }
    }
    x.f <- sort(do.call(c, unlist(x.f.l, use.names = FALSE)))
    output <- c(output, x.f)
    lc <- lc + 1
    if (progress) setTxtProgressBar(pb, lc)
  }
  if (length(samples) > 1) {
    output <- GRangesList(output)
  }
  if (progress) close(pb)

  skipped <- unlist(skipped)
  if (is.null(skipped)) {
    names(output) <- names(samples)
  } else {
    names(output) <- names(samples)[-skipped]
  }

  attributes(output)$manifest <- samples
  done.time <- Sys.time() - start.time
  if (progress) message("processed ", length(samples), " in ", round(done.time, digits = 2), " ", attr(done.time, "units"))
  return(output)
}

