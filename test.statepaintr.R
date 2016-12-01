## Vista invivo tested is https://www.encodeproject.org/annotations/ENCSR964TTF/ mm10 vista validations
## VISTA Browser in vivo tested elements summary tracks as of Oct 3, 2016

## ENCSR215ZYV is Reference epigenome for embryonic mouse neural tube (11.5 day).
## ENCSR843IAS is Reference epigenome for mouse midbrain (embryonic 11.5 day).
## ENCSR501OPC is Reference epigenome for mouse hindbrain (embryonic 11.5 day).
## ENCSR283NCE is Reference epigenome for mouse limb (embryonic 11.5 day).

## ENCFF786KUB is the Enhancer-like regions using DNase and H3K27ac for neural tube 11.5 day (embryonic)
## ENCFF520EGD is the Enhancer-like regions using DNase and H3K27ac for limb 11.5 day (embryonic)
## ENCFF733UJT is the Enhancer-like regions using DNase and H3K27ac for midbrain 11.5 day (embryonic)
## ENCFF324INM is the Enhancer-like regions using DNase and H3K27ac for hindbrain 11.5 day (embryonic)

library(prg)
vista.e <- data.table::fread("~/Downloads/VISTA_invivo_tested.bed", data.table=F)
#vista.e <- data.table::fread("~/Downloads/VISTA_invivo_tested_hg19.bed", data.table=F)
vista.validation <- stringr::str_replace(vista.e$V10, "ganglion, cranial", "ganglion-cranial")
vista.validation <- stringr::str_split(vista.validation, ",")

vista.validation <- lapply(vista.validation, function(x) {
  x <- stringr::str_replace(x, " \\(.*\\)", "")
  x <- stringr::str_trim(x, side = "both")
  return(x)
})

vista.validation.no.score <- lapply(vista.validation, function(x) {
  x <- stringr::str_replace(x, "\\[.*\\]", "")
})

vista.e.m <- matrix(0L, nrow = nrow(vista.e), ncol = length(unique(unlist(vista.validation.no.score))))
colnames(vista.e.m) <- unique(unlist(vista.validation.no.score))

for(enhancer in seq_along(vista.validation.no.score)) {
  vista.e.m[enhancer, vista.validation.no.score[[enhancer]]] <- 1L
}

library(GenomicRanges)
vista.e <- GRanges(seqnames = vista.e$V1,
                   ranges = IRanges(start = vista.e$V2,
                                    end = vista.e$V3),
                   name = vista.e$V4,
                   validated = vista.e$V5)

mcols(vista.e) <- cbind(as.data.frame(mcols(vista.e)), vista.e.m)

x <- get.decision.matrix("5813b67f46e0fb06b493ceb0")
x.scores.splits@abstraction.layer
x.splits@abstraction.layer
mouse.embryo.states <- PaintStates(manifest = "/Volumes/single_cell/mouse_mm10_11.5days/IDR.Manifest.txt", decisionMatrix = x, scoreStates = TRUE)

plot.data <- list()
result <- data.frame(tissue = character(), state = character(), precision = numeric(), recall = numeric())
for(tissue in names(mouse.embryo.states)) {
  neural.tube <- vista.e
  tissue <- stringr::str_replace(tissue, " embryo", "")
  mcols(neural.tube) <- data.frame(FOUND=mcols(neural.tube)[, tissue])
  neural.tube.enc <- neural.tube
  neural.tube.scores <- mouse.embryo.states[[paste0(tissue, " embryo")]]

  for(chosen.state in unique(neural.tube.scores$state)) {
    neural.tube.enhancer.scores <- neural.tube.scores[neural.tube.scores$state %in% chosen.state, ]
  #neural.tube.enhancer.scores <- neural.tube.scores[neural.tube.scores$state %in% c("EARC", "ARC", "RPS"), ]
    olaps <- findOverlaps(neural.tube, neural.tube.enhancer.scores)
    mcols(neural.tube)$score <- 0
    mcols(neural.tube)[from(olaps), "score"]<- mcols(neural.tube.enhancer.scores[to(olaps)])$score

    prg_curve <- create_prg_curve(mcols(neural.tube)$FOUND, mcols(neural.tube)$score)
    l.result <- data.frame(tissue = tissue, state = chosen.state, precision = prg_curve$precision_gain, recall = prg_curve$recall_gain)
    result <- rbind(result, l.result)
  }
  #auprg = calc_auprg(prg_curve)
  #convex_hull = prg_convex_hull(prg_curve)
  #plot.tissue <- list(list(tissue = tissue, curve=prg_curve, auprg = auprg, hull = convex_hull))
  #names(plot.tissue) <- tissue
  #plot.data <- c(plot.data, plot.tissue)
}
#median.tests.subscore.nobreak4 <- lapply(plot.data, function(x) {x$auprg})
#unlist(median.tests.subscore.nobreak4) <= unlist(median.tests.subscore.nobreak3)

#median.tests <- lapply(plot.data, function(x) {x$auprg})

#mean.tests <- lapply(plot.data, function(x) {x$auprg})
gg.tissue <- lapply(plot.data, function(x) {
  x <- data.frame(x$tissue, x$curve$precision_gain, x$curve$recall_gain)
  return(x)
})
gg.tissue <- do.call("rbind", gg.tissue)
library(ggplot2)
ggplot(gg.tissue, aes(y = x.curve.precision_gain, x = x.curve.recall_gain, group = x.tissue)) +
  geom_line(aes(color = x.tissue)) +
  geom_point(aes(color=x.tissue)) +
  coord_cartesian(xlim=c(0,1), ylim = c(0,1))

x <- mouse.embryo.states[["limb embryo"]]
x <- mcols(x)[, c("state", "score")]
x <- x[order(x$score), ]
x$order <- nrow(x):1
ggplot(as.data.frame(x), aes(x = order, y = score, color = state)) + geom_point() + facet_wrap( ~ state)

ggplot(result, aes(y = precision, x = recall, group = tissue)) +
  geom_line(aes(color = tissue)) +
  geom_point(aes(color= tissue)) +
  coord_cartesian(xlim=c(0,1), ylim = c(0,1)) + facet_grid(tissue ~ state)


#fig = plot_prg(prg_curve)
#print(prg_curve)
print(auprg)
}
#print(convex_hull)
print(fig)
##              OURS     ENCODE
## limb auprg = 0.82150, 0.83241
## hindbrain  = 0.72679, 0.77427
## midbrain   = 0.70325, 0.79657
## neuraltube = 0.73256, 0.76932

pred <- ROCR::prediction(mcols(neural.tube)$score, mcols(neural.tube)$FOUND)
perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
library(ROCR)
plot(perf, col=rainbow(10))

ef.enhancer <- read.delim("~/Downloads/10.1371%2Fjournal.pcbi.1003677.s018/enhancerfinder_limb_hg19.bed")
ef.enhancer <- GRanges(seqnames = ef.enhancer$X.chromosome,
                       ranges = IRanges(start = ef.enhancer$start,
                                        end = ef.enhancer$end),
                       score = ef.enhancer$MKL_scores)
neural.tube.ef <- vista.e
mcols(neural.tube.ef) <- data.frame(FOUND=neural.tube.ef$limb)
olaps.ef <- findOverlaps(neural.tube.ef, ef.enhancer)
mcols(neural.tube.ef)$score <- 0
mcols(neural.tube.ef)[from(olaps.ef), "score"] <- mcols(ef.enhancer)[to(olaps.ef), "score"]
prg_curve.ef <- create_prg_curve(mcols(neural.tube.ef)$FOUND, mcols(neural.tube.ef)$score)
auprg.ef <- calc_auprg(prg_curve.ef)
convex_hull.ef <- prg_convex_hull(prg_curve.ef)
fig.ef <- plot_prg(prg_curve.ef)
#print(prg_curve.ef)
print(auprg.ef)
print(convex_hull.ef)
print(fig.ef)

pred <- ROCR::prediction(mcols(neural.tube.ef)$score, mcols(neural.tube.ef)$FOUND)
perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
library(ROCR)
plot(perf, col=rainbow(10))

## ENCFF786KUB is the Enhancer-like regions using DNase and H3K27ac for neural tube 11.5 day (embryonic)
## ENCFF520EGD is the Enhancer-like regions using DNase and H3K27ac for limb 11.5 day (embryonic)
## ENCFF733UJT is the Enhancer-like regions using DNase and H3K27ac for midbrain 11.5 day (embryonic)
## ENCFF324INM is the Enhancer-like regions using DNase and H3K27ac for hindbrain 11.5 day (embryonic)

neural.tube.enc <- vista.e
mcols(neural.tube.enc) <- data.frame(FOUND=neural.tube.enc$`midbrain`)

encode.enhancer <- read.delim("~/Downloads/ENCFF520EGD.bed.gz", header = F) ## limb
encode.enhancer <- read.delim("~/Downloads/ENCFF324INM.bed.gz", header = F) ## hindbrain
encode.enhancer <- read.delim("~/Downloads/ENCFF733UJT.bed.gz", header = F) ## midbrain
encode.enhancer <- read.delim("~/Downloads/ENCFF786KUB.bed.gz", header = F) ## neuraltube
encode.enhancer$V5 <- nrow(encode.enhancer):1
encode.enhancer <- GRanges(seqnames = encode.enhancer$V1,
                           ranges = IRanges(start = encode.enhancer$V2,
                                            end = encode.enhancer$V3),
                           score = encode.enhancer$V5)

olaps.enc <- findOverlaps(neural.tube.enc, encode.enhancer)
mcols(neural.tube.enc)$score <- 0
mcols(neural.tube.enc)[from(olaps.enc), "score"] <- mcols(encode.enhancer)[to(olaps.enc), "score"]
prg_curve.enc <- create_prg_curve(mcols(neural.tube.enc)$FOUND, mcols(neural.tube.enc)$score)
auprg.enc <- calc_auprg(prg_curve.enc)
convex_hull.enc <- prg_convex_hull(prg_curve.enc)
fig.enc <- plot_prg(prg_curve.enc)
#print(prg_curve.enc)
print(auprg.enc)
#print(convex_hull.enc)
print(fig.enc)

pred <- ROCR::prediction(mcols(neural.tube.enc)$score, mcols(neural.tube.enc)$FOUND)
perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
library(ROCR)
plot(perf, col=rainbow(10))
