#' Run a StatePaintR decisionMatrix against a group of files listed in a manifest.
#'
#' @param manifest A character vector containing the filename of the manifest file. See details for the expected format
#' @param decisionMatrix An object of class \code{\linkS4class{decisionMatrix}}.
#' @param scoreStates logical; if scores have been specified in the decisionMatrix, should the states be scored?
#' @param progress logical; show progress bar and timing details?
#'
#' @return A \code{\link{GRangesList}}, each \code{\link{GRanges}} object describing the states of all
#' segements for a sample. Each range is described by the fields \code{name} indicating the sample,
#' \code{state} indicating the state, and optionally \code{score} indicating the score of the state
#' @details The manifest is a tab delimited file containing five fields; SAMPLE, MARK, SRC, BUILD, and FILE. \cr
#' `SAMPLE` refers to the sample to which the marks are related, like a cell line, or tissue. \cr
#' `MARK` refers to the chromatin mark, or other feature described in the \code{\linkS4class{decisionMatrix}} \code{abstractionLayer} \cr
#' `SRC` refers to the source of the data, retained for documentation, but not used in the functions. \cr
#' `BUILD` refers to the genome build of the data tracks. All tracks under a single sample must be of the same genome build. \cr
#' `FILE` refers to the location of the file describing the mark. Can be .bed, .narrowPeak, .gappedPeak, etc. Only the first three
#' columns are used, `chromosome`, `start`, and `end`, unless \code{scoreStates = TRUE}, in which case a narrowPeak file is required.
#' @examples
#' manifest <- system.file("extdata", "manifest.hmec.txt", package = "StatePaintR")
#' load(system.file("extdata", "poised.promoter.model.rda", package = "StatePaintR"))
#' states <- PaintStates(manifest = manifest,
#'                       decisionMatrix = poised.promoter.model,
#'                       scoreStates = FALSE, progress = FALSE)
#' states
#'
#' @importFrom GenomicRanges disjoin findOverlaps reduce binnedAverage coverage
#' @importFrom S4Vectors from to DataFrame
#' @importFrom stringr str_detect coll str_replace
#' @importFrom data.table frankv
#' @importFrom matrixStats rowMedians rowMaxs
#' @importFrom Matrix rowMeans
#' @importFrom IRanges %outside%
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
PaintStates <- function(manifest, decisionMatrix, scoreStates = FALSE, progress = TRUE) {
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
    x <- GetBioFeatures(manifest = cell.sample, dm = decisionMatrix, my.seqinfo = sample.genomes[[cell.sample[1, "BUILD"]]])
    if(is.null(x)) {
      warning(cell.sample[1, "SAMPLE"], " has no MARKs that are present in the translation layer \n",
              " MARKs included are ", paste(cell.sample[, "MARK"], collapse = " "))
      skipped <- c(skipped, skipname)
      lc <- lc + 2
      next()
    }
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
      lc <- lc + 2
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
      segments <- data.frame(state = flookup(resmatrix, dl)[, 1], score = NA, stringsAsFactors = FALSE)
      #segments <- data.frame(state = flookup(resmatrix, dl)[, 1], median = NA, mean = NA, max = NA, stringsAsFactors = FALSE)
      score.features <- unique(segments$state)[(unique(segments$state) %in% rownames(score.cells))]
      for (score.feature in score.features) {
        scorematrix.rows <- segments$state == score.feature
        feature.scores <- scorematrix[segments$state == score.feature,
                                      score.cells[rownames(score.cells) %in% score.feature, "col"],
                                      drop = FALSE]
        ## debug
        # if (score.feature == score.features[[1]]) {
        #   feature.score.tmp1 <- scorematrix[segments$state == score.feature, ]
        #   feature.score.tmp1 <- cbind(state = segments[segments$state == score.feature, "state"], feature.score.tmp1)
        #   if (is(feature.scores, "matrix")) {
        #     feature.score.tmp1 <- cbind(feature.score.tmp1, totalscore = rowMedians(feature.scores))
        #   } else {
        #     feature.score.tmp1 <- cbind(feature.score.tmp1, totalscore = feature.scores / 2)
        #   }
        # } else {
        #   feature.score.tmp2 <- scorematrix[segments$state == score.feature, ]
        #   feature.score.tmp2 <- cbind(state = segments[segments$state == score.feature, "state"], feature.score.tmp2)
        #   if (is(feature.scores, "matrix")) {
        #     feature.score.tmp2 <- cbind(feature.score.tmp2, totalscore = rowMedians(feature.scores))
        #   } else {
        #     feature.score.tmp2 <- cbind(feature.score.tmp2, totalscore = feature.scores / 2)
        #   }
        #   feature.score.tmp1 <- rbind(feature.score.tmp1, feature.score.tmp2)
        # }
        ## end debug
        #feature.scores.df <- data.frame(median = numeric(), mean = numeric(), max = numeric())
        if (is(feature.scores, "matrix")) {
          feature.scores <- rowMaxs(feature.scores)
          #feature.scores.df <- rbind.data.frame(feature.scores.df, data.frame(median = rowMedians(feature.scores), mean = rowMeans(feature.scores), max = rowMaxs(feature.scores)))
        } else {
          feature.scores <- feature.scores
          #feature.scores <- feature.scores / 2
          #feature.scores.df <- rbind.data.frame(feature.scores.df, data.frame(median = feature.scores / 2, mean = feature.scores / 2, max = feature.scores))
        }
        #if(cell.sample[1,1] == "thyroid gland" & score.feature == "PPR") browser()
        segments[segments$state == score.feature, "score"] <- feature.scores
        #segments[segments$state == score.feature, c("median", "mean", "max")] <- feature.scores.df
      }
      # feature.score.bm2 <<- feature.score.tmp1
      mcols(x.f)$state <- segments$state
      mcols(x.f)$score <- segments$score
      #mcols(x.f)[, c("median", "mean", "max")] <- DataFrame(segments[, c("median", "mean", "max")])
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
          my.scores <- sapply(revmap, function(my.splits, orig.gr.cols) {
            if(length(my.splits) > 1) {
              if(zero_range(orig.gr.cols[my.splits])) {
                my.splits[1]
              } else {
                my.splits[which(orig.gr.cols[my.splits] == max(orig.gr.cols[my.splits]))[[1]]]
              }
            } else {
              my.splits
            }
          }, orig.gr.cols = mcols(x.f.l[[score.feature]])$score)
          my.scores <- mcols(x.f.l[[score.feature]])[my.scores, "score"]
          #my.scores <- mcols(x.f.l[[score.feature]])[revmap, "score"]
          #my.scores <- mcols(x.f.l[[score.feature]])[revmap, c("median", "mean", "max")]
        } else {
          my.scores <- 0
          #my.scores <- DataFrame(median = rep.int(0, length(x.f.ll[[state.name]])), mean = rep.int(0, length(x.f.ll[[state.name]])), max = rep.int(0, length(x.f.ll[[state.name]])))
        }
        mcols(x.f.ll[[state.name]])$revmap <- NULL
        mcols(x.f.ll[[state.name]])$name <- cell.sample[1, "SAMPLE"]
        mcols(x.f.ll[[state.name]])$state <- state.name
        if(max(my.scores) > 0) {
          mcols(x.f.ll[[state.name]])$score <- ceiling((my.scores/max(my.scores)) * 1000)
        } else {
          mcols(x.f.ll[[state.name]])$score <- 0
        }
        #mcols(x.f.ll[[state.name]])[, c("median", "mean", "max")] <- my.scores
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
    #x.f$score <- (x.f$score/max(x.f$score)) * 1000
    # mcols(x.f)[, c("median", "mean", "max")] <- lapply(mcols(x.f)[, c("median", "mean", "max")], function(x) {
    #   x/max(x) * 1000
    # })
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

#' Write StatePaintR object to bedfiles
#'
#' @param states A GRangesList produced by \code{\link{PaintStates}}
#' @param decisionMatrix The \code{\linkS4class{decisionMatrix}} object used to produce the states
#' @param output.dir A character string indicating the directory to save the exported states.
#' The directory will be created if it does not exist.
#' @param progress logical; show progress bar and timing details?
#'
#' @importFrom rtracklayer export.bed
#' @importFrom dplyr left_join
#' @return Invisibly returns the states object.
#' @export
#' @examples
#' \dontrun{
#' ExportStatePaintR(states = states,
#'                   decisionMatrix = poised.promoter.model,
#'                   output.dir = tempdir())
#' }
ExportStatePaintR <- function(states, decisionMatrix, output.dir, progress = TRUE) {
  if (missing(output.dir)) { stop("please indicate output directory") }
  if (missing(decisionMatrix)) { stop("please include decisionMatrix") }
  if (!dir.exists(output.dir)) { dir.create(output.dir) }
  m.data <- attributes(states)$manifest
  color.key <- stateColors(decisionMatrix)
  hub.id <- decisionMatrix@id
  if (progress) pb <- txtProgressBar(min = 0, max = length(states), style = 3)
  for (state in seq_along(states)) {
    if (progress) setTxtProgressBar(pb, state)
    s.name <- str_replace(names(m.data)[state], "/", "-")
    s.m.data <- m.data[[state]]
    state <- states[[state]]
    write.state(state, s.m.data, color.key, hub.id,
                file.path(output.dir,
                          paste0(s.name, ".", nrow(s.m.data), "mark.segmentation.bed")))
  }
  if (progress) close(pb)
  message("segmentation files written to: ", output.dir)
  return(invisible(states))
}

#' Retrieve Decision Matrix from StateHub
#'
#' @param search character string of either unique id, or search term.
#'
#' @return \code{\linkS4class{decisionMatrix}} object
#' @importFrom httr content GET http_error modify_url
#' @importFrom jsonlite toJSON fromJSON
#' @importFrom stringr str_extract
#' @export
#'
#' @examples
#' get.decision.matrix(search = "Cedars-Sinai")
#' poised.promoter.model <- get.decision.matrix("5813b67f46e0fb06b493ceb0")
#' poised.promoter.model
get.decision.matrix <- function(search) {
  if (missing(search)) { stop("missing search argument") }
  stateHub <- modify_url("http://statehub.org/statehub", path = "statehub/getmodel.jsp")
  query <- GET(stateHub, query = list(id = search))
  if (length(content(query)) < 1) {
    query <- GET(stateHub, query = list(search = search))
    if (http_error(query)) { stop("site not availible") }
    if (length(content(query)) < 1) { stop("no search results found") }
  }
  query <- content(query)
  for (query.i in seq_along(query)) {
    query.result <- query[[query.i]]
    state.names <- sapply(query.result$states, function(x) {x$name})
    for (state in seq_along(query.result$states)) {
      if (state == 1) {
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
    if (query.i <= 1) {
      output <- decision.matrix
    } else {
      output <- c(output, decision.matrix)
    }
  }
  return(output)
}

#' decisionMatrix class
#'
#' @slot id character. a unique ID that refers to a model revision
#' @slot name character. The name of the model on StateHub
#' @slot author character. The author of the model on StateHub
#' @slot revision character. The revision of the model. A time stamp
#' @slot description character. A short description of the model.
#' @slot abstraction.layer list. A description of the relationship between the precise data type (e.g. a chromatin mark like H3K27ac)
#' and feature describing the functional category in the decision.matrix (e.g. Regulatory).
#' @slot decision.matrix matrix. The description of what combination of features are required to call a state.
#' @slot state.colors character. The color to assign to each state.
#'
#' @import methods
#' @return an object of class decisionMatrix, normally made in response to a query to StateHub with \code{\link{get.decision.matrix}}
#' @export
#'
setClass(Class = "decisionMatrix",
         slots = list(id = "character",
                      name = "character",
                      author = "character",
                      revision = "character",
                      description = "character",
                      abstraction.layer = "list",
                      decision.matrix = "matrix",
                      state.colors = "character"))




#' decisionMatrix
#'
#' @param object of class decisionMatrix
#'
#' @return a matrix; The description of what combination of features are required to call a state.
#' @export
#'
#' @examples
#' load(system.file("extdata", "poised.promoter.model.rda", package = "StatePaintR"))
#' poised.promoter.model
#' decisionMatrix(poised.promoter.model)
setGeneric("decisionMatrix", function(object) standardGeneric("decisionMatrix"))
#' @describeIn decisionMatrix Extract Descision Matrix from descisionMatrix object
#' @export
#' @examples
#' load(system.file("extdata", "poised.promoter.model.rda", package = "StatePaintR"))
#' poised.promoter.model
#' decisionMatrix(poised.promoter.model)
setMethod("decisionMatrix",
          signature = signature(object = "decisionMatrix"),
          function(object) {
            object@decision.matrix
          })


#' abstractionLayer
#'
#' @param object of class decisionMatrix
#'
#' @return a list; A description of the relationship between the precise data type (e.g. a chromatin mark like H3K27ac)
#' and feature describing the functional category in the decision.matrix (e.g. Regulatory).
#' @details In the abstraction layer one may indicate if a functional category should not be used for scoring states.
#' This is done by placing an asterisk at the end of it's name (e.g. Regulatory*). One may also specify that
#' a certain type of functional category never be split into smaller units by its overlapping with other features.
#' This is done by placing the name of the functional category in square brackets (e.g. [Core])
#' @export
#'
#' @examples
#' load(system.file("extdata", "poised.promoter.model.rda", package = "StatePaintR"))
#' poised.promoter.model
#' abstractionLayer(poised.promoter.model)
setGeneric("abstractionLayer", function(object) standardGeneric("abstractionLayer"))
#' @describeIn decisionMatrix Extract abstraction layer from descisionMatrix object
#' @export
#' @examples
#' abstractionLayer(poised.promoter.model)
setMethod("abstractionLayer",
          signature = signature(object = "decisionMatrix"),
          function(object) {
            object@abstraction.layer
          })

#' doNotScore
#'
#' @param decisionMatrix of class decisionMatrix
#' @param functionalCategory a character vector indicating the functional category to be
#' excluded from scoring.
#' @return a decisionMatrix; the abstractionLayer of the decision matrix will be altered
#' to indicate that the functional category indicated must not be used in scoring calculations.
#' @details Sometimes one may wish to use a functional category for calling a state, but not use
#' the functional category for any scoring. This allows the user that option.
#' @export
#'
#' @examples
#' dm <- get.decision.matrix("5813b67f46e0fb06b493ceb0")
#' dm <- doNotScore(decisionMatrix = dm, functionalCategory = "Regulatory")
#' dm
doNotScore <- function(decisionMatrix, functionalCategory) {
  if (missing(decisionMatrix)) {
    stop("provide a decisionMatrix object")
  }
  if (missing(functionalCategory)) {
    stop("provide a functionalCategory to modify")
  }
  if (inherits(decisionMatrix, "decisionMatrix")) {
    categories <- grep(pattern = functionalCategory, x = names(abstractionLayer(decisionMatrix)))
    if (length(colnames) > 0) {
      for (category.i in categories) {
        category <- names(abstractionLayer(decisionMatrix))[category.i]
        names(decisionMatrix@abstraction.layer)[category.i] <- paste0(category, "*")
      }
    }
  } else {
    stop("decisionMatrix must be an object of class decisionMatrix")
  }
  return(decisionMatrix)
}

#' doNotSplit
#'
#' @param decisionMatrix of class decisionMatrix
#' @param functionalCategory a character vector indicating intervals described by the
#'  functional category will never be split into smaller intervals
#' @return a decisionMatrix; the abstractionLayer of the decision matrix will be altered
#' to indicate that the functional category indicated must not be split.
#' @details Some functional categories described in an abstraction layer may perform better
#' or more like one expects if they are never split into smaller intervals due to being
#' overlapped by other functional categories. One may indicate if this is the case with this
#' function
#' @export
#'
#' @examples
#' dm <- get.decision.matrix("5813b67f46e0fb06b493ceb0")
#' dm <- doNotSplit(decisionMatrix = dm, functionalCategory = "Core")
#' dm
doNotSplit <- function(decisionMatrix, functionalCategory) {
  if (missing(decisionMatrix)) {
    stop("provide a decisionMatrix object")
  }
  if (missing(functionalCategory)) {
    stop("provide a functionalCategory to modify")
  }
  if (inherits(decisionMatrix, "decisionMatrix")) {
    categories <- grep(pattern = functionalCategory, x = names(abstractionLayer(decisionMatrix)))
    if (length(colnames) > 0) {
      for (category.i in categories) {
        category <- names(abstractionLayer(decisionMatrix))[category.i]
        if (grepl("\\*$", category)) {
          category <- str_replace(category, "\\*$", "")
          category <- paste0("[", category, "]*")
        } else {
          category <- paste0("[", category, "]")
        }
        names(decisionMatrix@abstraction.layer)[category.i] <- category
      }
    }
  } else {
    stop("decisionMatrix must be an object of class decisionMatrix")
  }
  return(decisionMatrix)
}

#' stateColors
#'
#' @param object of class decisionMatrix
#'
#' @return character. The color to assign to each state.
#' @export
#'
#' @examples
#' load(system.file("extdata", "poised.promoter.model.rda", package = "StatePaintR"))
#' poised.promoter.model
#' abstractionLayer(poised.promoter.model)
setGeneric("stateColors", function(object) standardGeneric("stateColors"))
#' @describeIn decisionMatrix Extract color information from descisionMatrix object
#' @param object decisionMatrix.
#' @export
#' @examples
#' stateColors(poised.promoter.model)
setMethod("stateColors",
          signature = signature(object = "decisionMatrix"),
          function(object) {
            data.frame(STATE = names(object@state.colors), COLOR = object@state.colors, stringsAsFactors = FALSE)
          })
