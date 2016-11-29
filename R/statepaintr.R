#' Run a StatePaintR decisionMatrix against a group of files listed in a manifest.
#'
#' @param manifest
#' @param decisionMatrix
#' @param scoreStates
#' @param progress
#'
#' @return
#' @importFrom GenomicRanges disjoin findOverlaps reduce
#' @importFrom S4Vectors from to
#' @importFrom stringr str_detect coll str_replace
#' @importFrom data.table frankv
#' @importFrom matrixStats rowMedians
#' @importFrom IRanges %outside%
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#'
#' @examples
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
      score.cells <- which(dl.score == 3L, arr.ind = T)
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
          #my.scores <- (my.scores/max(my.scores)) * 1000
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
  close(pb)
  names(output) <- names(samples)
  attributes(output)$manifest <- samples
  done.time <- Sys.time() - start.time
  if(progress) message("processed ", length(samples), " in ", round(done.time, digits = 2), " ", attr(done.time, "units"))
  return(output)
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
  hub.id <- x@id
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

#' Title
#'
#' @slot id character.
#' @slot name character.
#' @slot author character.
#' @slot revision character.
#' @slot description character.
#' @slot abstraction.layer list.
#' @slot decision.matrix matrix.
#' @slot state.colors character.
#'
#' @return
#' @export
#'
#' @examples
setClass(Class = "decisionMatrix",
         slots = list(id = "character",
                      name = "character",
                      author = "character",
                      revision = "character",
                      description = "character",
                      abstraction.layer = "list",
                      decision.matrix = "matrix",
                      state.colors = "character"))

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
