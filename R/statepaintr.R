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
#' @importFrom GenomicRanges disjoin findOverlaps reduce
#' @importFrom S4Vectors from to
#' @importFrom stringr str_detect coll str_replace
#' @importFrom data.table frankv
#' @importFrom matrixStats rowMedians
#' @importFrom IRanges %outside%
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
PaintStates <- function(manifest, decisionMatrix, scoreStates = FALSE, progress = TRUE) {
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
  do.seqinfo <- try(Seqinfo(genome = sample.genomes.names[[1]]))
  if(inherits(do.seqinfo, "try-error")) {
    sample.genomes <- lapply(sample.genomes.names, function(x) Seqinfo(genome = NA))
    warning("could not download seqinfo for genome")
  } else {
    sample.genomes <- lapply(sample.genomes.names, function(x) Seqinfo(genome = x))
  }
  names(sample.genomes) <- sample.genomes.names
  output <- list()
  lc <- 0
  if(progress) pb <- txtProgressBar(min = 0, max = length(samples)*3, style = 3)
  for(cell.sample in seq_along(samples)) {
    lc <- lc + 1
    if(progress) setTxtProgressBar(pb, lc)
    cell.sample <- samples[[cell.sample]]
    x <- GetBioFeatures(manifest = cell.sample, my.seqinfo = sample.genomes[[cell.sample[1, "BUILD"]]])
    names(x) <- tolower(names(x))
    inputset <- names(x)
    inputset <- sapply(inputset, function(x, y = deflookup) {
      out <- unlist(y[str_detect(x, paste0(names(y), "$"))], use.names = FALSE)
      return(out)
      })
    names(x) <- inputset
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
    x.f <- disjoin(unlist(x))
    if(length(x) == 1) {
      x.f <- x.f[[1]]
    }
    dontbreak <- str_replace(inputset[grep(pattern = "^\\*", inputset)], "^\\*", "")
    dontbreak <- str_replace(dontbreak, "\\*$", "")

    useScore <- str_replace(inputset[grep(pattern = "\\*$", inputset)], "\\*$", "")
    useScore <- str_replace(useScore, "^\\*", "")

    inputset <- str_replace(inputset, "^\\*", "")
    inputset <- str_replace(inputset, "\\*$", "")

    names(x) <- str_replace(inputset, "^\\*", "")
    names(x) <- str_replace(inputset, "\\*$", "")
    if(length(dontbreak) > 0) {
      x.f <- x.f[x.f %outside% x[dontbreak]]
      unfrag <- reduce(sort(do.call(c, list(unlist(x[dontbreak]), ignore.mcols = TRUE))))
      x.f <- c(x.f, unfrag, ignore.mcols = TRUE)
      x.f <- sort(x.f)
    }
    overlaps <- findOverlaps(x.f, x)
    resmatrix <- make.overlap.matrix(from(overlaps), to(overlaps), names(x))
    resmatrix <- resmatrix[, colnames(d)[colnames(d) %in% colnames(resmatrix)]]
    if(scoreStates) {
      scorematrix <- resmatrix
      signalCol <- sapply(x, function(x) !is.na(mcols(x)[1, "signalValue"]))
      if(length(useScore) < 1) useScore <- names(signalCol)
      signalCol <- signalCol[names(signalCol) %in% useScore]
      signalCol <- which(signalCol[colnames(scorematrix)])
      for(feature in names(signalCol)) {
        scoreOverlaps <- findOverlaps(x.f, x[[feature]])
        scorematrix[from(scoreOverlaps), feature] <- mcols(x[[feature]])[to(scoreOverlaps), "signalValue"]
        scorematrix[, feature] <- frankv(scorematrix[, feature], ties.method = "average")
      }
      scorematrix <- scorematrix[, names(signalCol)]
    }
    resmatrix[resmatrix == 0L] <- 2L
    resmatrix[resmatrix == 1L] <- 3L
    missing.data <- colnames(d)[!(colnames(d) %in% colnames(resmatrix))]
    if(length(missing.data) > 0){
      d.m <- d[, missing.data] > 1
      d.m[is.na(d.m)] <- FALSE
      if(!inherits(d.m, "matrix")) d.m <- matrix(d.m, ncol = 1, dimnames = list(c(names(d.m)), c("SAMPLE")))
      dl <- d[rowSums(d.m) < 1, colnames(d) %in% colnames(resmatrix)]
    }
    dl <- dl[order(rowSums(dl, na.rm = TRUE), decreasing = FALSE), ]
    mcols(x.f)$name <- cell.sample[1, "SAMPLE"]
    resmatrix <- t(resmatrix);
    lc <- lc + 1
    if(progress) setTxtProgressBar(pb, lc)
    if(scoreStates) {
      dl.score <- dl[, names(signalCol)]
      score.cells <- which(dl.score == 3L, arr.ind = TRUE)
      segments <- data.frame(state = footprintlookup(resmatrix, dl)[, 1], score = NA, stringsAsFactors = FALSE)
      score.features <- unique(segments$state)[(unique(segments$state) %in% rownames(score.cells))]
      for(score.feature in score.features) {
        feature.scores <- scorematrix[segments$state == score.feature, score.cells[rownames(score.cells) %in% score.feature, "col"]]
        if(is(feature.scores, "matrix")) {
          feature.scores <- rowMedians(feature.scores)
        } else {
          feature.scores <- feature.scores/2
        }
        segments[segments$state == score.feature, "score"] <- feature.scores
      }
      mcols(x.f)$state <- segments$state
      mcols(x.f)$score <- segments$score
    } else {
      mcols(x.f)$state <- footprintlookup(resmatrix, dl)[, 1]
    }
    x.f.l <- split(x.f, x.f$state)
    if(scoreStates) {
      x.f.ll <- lapply(x.f.l, reduce, with.revmap=TRUE)
      for(state.name in names(x.f.ll)) {
        if(state.name %in% score.features) {
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
        } else {
          my.scores <- 0
        }
        mcols(x.f.ll[[state.name]])$revmap <- NULL
        mcols(x.f.ll[[state.name]])$name <- cell.sample[1, "SAMPLE"]
        mcols(x.f.ll[[state.name]])$state <- state.name
        mcols(x.f.ll[[state.name]])$score <- my.scores
      }
      x.f.l <- x.f.ll
    } else {
      x.f.l <- lapply(x.f.l, reduce)
      for(state.name in names(x.f.l)) {
        mcols(x.f.l[[state.name]])$name <- cell.sample[1, "SAMPLE"]
        mcols(x.f.l[[state.name]])$state <- state.name
        mcols(x.f.l[[state.name]])$score <- 0
      }
    }
    x.f <- sort(do.call(c, unlist(x.f.l, use.names = FALSE)))
    x.f$score <- (x.f$score/max(x.f$score)) * 1000
    output <- c(output, x.f)
    lc <- lc + 1
    if(progress) setTxtProgressBar(pb, lc)
  }
  if(length(samples) > 1) {
    output <- GRangesList(output)
  }
  if(progress) close(pb)
  names(output) <- names(samples)
  attributes(output)$manifest <- samples
  done.time <- Sys.time() - start.time
  if(progress) message("processed ", length(samples), " in ", round(done.time, digits = 2), " ", attr(done.time, "units"))
  return(output)
}

#' Write StatePaintR object to bedfiles
#'
#' @param states A GRangesList produced by \code{\link{PaintStates}}
#' @param decisionMatrix The \code{\linkS4class{decisionMatrix}} object used to produce the states
#' @param output.dir A character string indicating the directory to save the exported states.
#' The directory will be created if it does not exist.
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
ExportStatePaintR <- function(states, decisionMatrix, output.dir) {
  if(missing(output.dir)) { stop("please indicate output directory") }
  if(missing(decisionMatrix)) { stop("please include decisionMatrix") }
  if(!dir.exists(output.dir)) { dir.create(output.dir) }
  m.data <- attributes(states)$manifest
  color.key <- stateColors(decisionMatrix)
  hub.id <- decisionMatrix@id
  pb <- txtProgressBar(min = 0, max = length(states), style = 3)
  for(state in seq_along(states)) {
    setTxtProgressBar(pb, state)
    s.name <- names(m.data)[state]
    s.m.data <- m.data[[state]]
    state <- states[[state]]
    write.state(state, s.m.data, color.key, hub.id,
                file.path(output.dir,
                          paste0(s.name, ".", nrow(s.m.data), "mark.segmentation.bed")))
  }
  close(pb)
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
  if(missing(search)) { stop("missing search argument") }
  stateHub <- modify_url("http://statehub.org/statehub", path = "statehub/getmodel.jsp")
  query <- GET(stateHub, query = list(id = search))
  if(length(content(query)) < 1) {
    query <- GET(stateHub, query = list(search = search))
    if(http_error(query)) { stop("site not availible") }
    if(length(content(query)) < 1) { stop("no search results found") }
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
#' @details In the abstraction layer one may indicate if a functional category should be used for scoring states.
#' This is done by placing an asterisk at the start of it's name (e.g. *Regulatory). One may also specify that
#' a certain type of functional category never be split into smaller units by its overlapping with other features.
#' This is done by placing an asterisk at the end of the name of the functional category (e.g. Core*)
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
