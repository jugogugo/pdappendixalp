require("annotatr")
require("GenomicRanges")


myFisherTest <- function(t) {
  require(data.table)
  require(dplyr)
  exp <- sum(t[2,])/sum(t)
  obs <- t[2,2] / sum(t[,2])
  f <- t %>% fisher.test
  data.table(
    `Expected, %` = exp,
    `Obseverd, %` = obs,
    OR            = f$estimate, 
    Lo            = f$conf.int[1], 
    Hi            = f$conf.int[2], 
    P             = f$p.value,
    `P < 0.05`    = gtools::stars.pval(f$p.value)
  )
}

loadHumanAnnotations <- function() {
	genome <- "hg19"
	annots <- c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic')
	build_annotations(genome = genome, annotations = annots)
}

getAnnotatrEnrichments <- function(fits, annotations, qThresh = 0.05) {
  fitsgr <- with(fits, GRanges(seqnames=Chr, IRanges(SNP, SNP+1)))
  if (!is.null(fits$str))
    strand(fitsgr) <- fits$str
  values(fitsgr) <- fits[, list(ID)]
  annot <- annotate_regions(
      regions = fitsgr,
      annotations = annotations,
      ignore.strand = TRUE,
      quiet = FALSE)
  annot <- data.frame(annot) %>% setDT
  rm(fitsgr)


  # get enrichments
  categories <- annot[, unique(annot.type)]
  pd <- 
    foreach (cat = categories, .combine = rbind) %do% {
      dt <- fits[, list(ID, Sig = adj.P.Val < qThresh, Dir = sign(logFC))]
      dt[, CAT := ID %in% annot[annot.type == cat, ID]]
      A <- dt[, table(CAT, Sig)] %>%
        myFisherTest %>%
        .[, Category := cat %>% as.character] %>%
        .[, Type     := "All"]

      B <- dt[, table(CAT, Sig & (Dir == 1))] %>%
        myFisherTest %>%
        .[, Category := cat %>% as.character] %>%
        .[, Type     := "Hyper-"]

      C <- dt[, table(CAT, Sig & (Dir == -1))] %>%
        myFisherTest %>%
        .[, Category := cat %>% as.character] %>%
        .[, Type     := "Hypo-"]
      rbind(A, B, C)
    }
  pd
}

plotAnnotatrEnrichment <- function(dt) {
  dt <- copy(dt)
  dt[, Element := strsplit(Category, split="_") %>% sapply(., `[[`, 2)]
  dt[, Label := strsplit(Category, split="_") %>% sapply(., `[[`, 3)]
  map <- list(
    inter      = "Open sea", 
    shelves    = "Shelves", 
    shores     = "Shores", 
    islands    = "Islands",
    `5UTRs`    = "3.5' UTRs",
    `3UTRs`    = "6.3' UTRs",
    exons      = "4.Exons", 
    introns    = "5.Introns",
    intergenic = "0.Intergenic", 
    `1to5kb`   = "1.1-5 kb",
    promoters  = "2.Promoters",
    Promoter   = "Active Promoter",
    PoisedEnhancer = "Poised Enhancers",
    ActiveEnhancer = "Active Enhancers"
  )
  old <- dt$Label
  new <- map[as.character(old)]
  dt[, Label := new %>% unlist]


  # CpG islands
  pd <- dt[Element == "cpg"] %>%
        .[, Label := Label %>% factor(., levels = unique(.))]
  p1 <- 
  ggplot(pd, 
    aes(Label, log2(OR), ymin=log2(Lo), ymax=log2(Hi), 
      group=Type, shape=Type)) + 
  geom_point(aes(color=P < 0.05),
    position=position_dodge(width=.5), size=4) + 
  geom_errorbar(width=0.25, position=position_dodge(width=.5), size=0.25) +
  geom_line(aes(linetype=Type),
    position=position_dodge(width=.5), color="grey", size=0.25) +
  geom_hline(yintercept=0, size=0.5) + 
  ggtitle("CpG islands") + 
  xlab("") + ylab("Odds ratio, log2") +
  scale_color_brewer(palette="Dark2") + 
  theme_bw(base_size=14) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1))

  pd <- dt[Element == "genes"] %>%
        .[, Label := Label %>% 
            factor(., labels = sort(unique(.)) %>% 
              gsub("^[0-9]\\.", "", .))]
  p2 <- 
  ggplot(pd, 
    aes(Label, log2(OR), ymin=log2(Lo), ymax=log2(Hi), 
      group=Type, shape=Type)) + 
  geom_point(aes(color=P < 0.05),
    position=position_dodge(width=.5), size=4) + 
  geom_errorbar(width=0.25, position=position_dodge(width=.5), size=0.25) +
  geom_line(aes(linetype=Type),
    position=position_dodge(width=.5), color="grey", size=0.25) +
  geom_hline(yintercept=0, size=0.5) + 
  ggtitle("Genomic features") + 
  xlab("") + ylab("Odds ratio, log2") +
  scale_color_brewer(palette="Dark2") + 
  theme_bw(base_size=14) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1))

  return(list(p1, p2))
}

getPsychEncodeEnrichments <- function(fits, qThresh = 0.05) {
  # Annotate with PsychEncode 
  files <- Sys.glob("../input/Enh_Prom/E07*.bed")
  pd <- foreach(fname = files, .combine=rbind) %do% {
    message(fname)
    f <- fread(fname, skip=1)
    f <- f[, list(annot.seqnames=V1, annot.start = V2, annot.end = V3, 
                  annot.width = V3 - V2, 
                  annot.strand = "*", 
                  annot.type = V4)] %>%
          setkey(annot.seqnames, annot.start, annot.end)


    # For each type of region get overlaps with fits and then produce the
    # enrichment 
    categories <- f[, unique(annot.type)]
    res <- 
    foreach (cat = categories, .combine = rbind) %do% {
      message(cat)
      dt <- fits[, 
              list(Chr, Start=SNP, End=SNP, 
                    Sig = adj.P.Val < qThresh, 
                    Dir = sign(logFC)
                  )
            ] %>%
            setkey(Chr, Start, End)
      dt <- foverlaps(dt, f[annot.type == cat])
      A <-  dt[, table(!is.na(annot.type), Sig)] %>%
            myFisherTest %>%
            .[, Category := cat %>% as.character] %>%
            .[, Type     := "All"]
      B <-  dt[, table(!is.na(annot.type), Sig & (Dir == 1))] %>%
            myFisherTest %>%
            .[, Category := cat %>% as.character] %>%
            .[, Type     := "Hyper-"]
      C <-  dt[, table(!is.na(annot.type), Sig & (Dir == -1))] %>%
            myFisherTest %>%
            .[, Category := cat %>% as.character] %>%
            .[, Type     := "Hypo-"]
      rbind(A, B, C)
    }
    res <- res[, File := fname]
    res
  }
}

plotPsychEncodeEnrichment <- function(dt, file = "E073") {
  dt <- copy(dt)
  dt[, File := 
          gsub("../input/Enh_Prom/", "", File) %>% 
          gsub("_25_imputed12marks_dense.bed", "", .)]
  dt[, Category := factor(Category,
          labels = sort(unique(Category)) %>% gsub("^.*_", "", .),
          levels = sort(unique(Category))
          )]

  dt[File == file] %>%
  ggplot(., 
    aes(Category, log2(OR), ymin=log2(Lo), ymax=log2(Hi), 
      group=Type, shape=Type)) + 
  geom_point(aes(color=P < 0.05),
    position=position_dodge(width=.5), size=4) + 
  geom_errorbar(width=0.25, position=position_dodge(width=.5), size=0.25) +
  geom_line(aes(linetype=Type),
    position=position_dodge(width=.5), color="grey", size=0.25) +
  geom_hline(yintercept=0, size=0.5) + 
  ggtitle(glue("EpiRoadmap ({file})")) + 
  xlab("") + ylab("Odds ratio, log2") +
  scale_color_brewer(palette="Dark2") + 
  theme_bw(base_size=14) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1))
}

getEpiRoadmapEnrichment <- function(fits, qThresh = 0.05) {
  # AC marks
  require("bedr")
  stopifnot(check.binary("bedtools"))
  myfits <- copy(fits)
  dAC <- fread("../input/Enh_Prom/EpiMap_H3K27ac_Neu_sort.txt") %>%
          setnames(c("Chr", "Start", "End", "Name"))
  merged <- 
    bedr(engine = "bedtools", input = list(i = dAC), method = "merge", 
        check.zero.based = FALSE, check.chr = FALSE, 
        check.valid = FALSE, check.merge = FALSE, check.sort = FALSE, 
        verbose = TRUE)  %>%
    setDT %>%
    setnames(c("Chr", "Start", "End"))
  merged %>%  setkey(Chr, Start, End)
  myfits[, AC := FALSE]
  ids <- foverlaps(
    myfits[, list(Chr, Start = SNP, End = SNP, ID)], 
    merged, 
    by.x = c("Chr", "Start", "End"),
    by.y = c("Chr", "Start", "End")) %>%
    .[!is.na(Start), ID]
  myfits[ID %in% ids, AC := TRUE]

  # ME marks
  dME <- fread("../input/Enh_Prom/EpiMap_H3K4me3_Neu_sort.txt") %>%
          setnames(c("Chr", "Start", "End", "Name"))
  merged <- 
    bedr(engine = "bedtools", input = list(i = dME), method = "merge", 
        check.zero.based = FALSE, check.chr = FALSE, 
        check.valid = FALSE, check.merge = FALSE, check.sort = FALSE, 
        verbose = TRUE)  %>%
    setDT %>%
    setnames(c("Chr", "Start", "End"))
  merged %>%  setkey(Chr, Start, End)
  myfits[, ME := FALSE]
  ids <- foverlaps(
    myfits[, list(Chr, Start = SNP, End = SNP, ID)], 
    merged, 
    by.x = c("Chr", "Start", "End"),
    by.y = c("Chr", "Start", "End")) %>%
    .[!is.na(Start), ID]
  myfits[ID %in% ids, ME := TRUE]

  # Now mark enhancers and promoters
  myfits %>%
    # Promoter is AC only
    .[, Promoter := ifelse( AC & !ME, TRUE, FALSE)] %>%
    # Poised enhancer is ME only
    .[, PoisedEnhancer := ifelse( !AC & ME, TRUE, FALSE)] %>%
    # Active enhancer is when both ME and AC are present
    .[, ActiveEnhancer := ifelse( ME & AC, TRUE, FALSE)] %>%
    # Drop unused cols
    .[, AC := NULL] %>%
    .[, ME := NULL]

  # Compute enrichments
  categories <- c("Promoter", "PoisedEnhancer", "ActiveEnhancer")
  pd <- 
    foreach (cat = categories, .combine = rbind) %do% {
      dt <- myfits[, 
              list(ID, Sig = adj.P.Val < qThresh, 
                  Dir = sign(logFC), CAT = get(cat))]
      A <- dt[, table(CAT, Sig)] %>%
        myFisherTest %>%
        .[, Category := cat %>% as.character] %>%
        .[, Type     := "All"]

      B <- dt[, table(CAT, Sig & (Dir == 1))] %>%
        myFisherTest %>%
        .[, Category := cat %>% as.character] %>%
        .[, Type     := "Hyper-"]

      C <- dt[, table(CAT, Sig & (Dir == -1))] %>%
        myFisherTest %>%
        .[, Category := cat %>% as.character] %>%
        .[, Type     := "Hypo-"]
      rbind(A, B, C)
    }
  pd  
}

plotEpiRoadmapEnrichment <- function(dt,
  wrap = c("none", "vertical", "horizontal")
) {
  dt <- copy(dt)
  map <- list(
    Promoter = "Active Promoter",
    PoisedEnhancer = "Enhancer Poised",
    ActiveEnhancer = "Enhancer Active")
  dt[, Category := map[as.character(Category)] %>% unlist]

  p <- ggplot(dt, 
    aes(Category, 
        log2(OR), 
        ymin=log2(Lo), 
        ymax=log2(Hi), 
        group=Type, 
        shape=Type)) + 
  geom_point(aes(color=P < 0.05), size=4,
    position=position_dodge(width=.5)) +
  geom_errorbar(width=0.25, size=0.25,
    position=position_dodge(width=.5)) +
  geom_line(aes(linetype=Type), color="grey", size=0.25,
    position=position_dodge(width=.5)) +
  geom_hline(yintercept=0, size=0.5) + 
  ggtitle("PsychENCODE") + 
  xlab("") + ylab("Odds ratio, log2") +
  scale_color_brewer(palette="Dark2") + 
  theme_bw(base_size=14) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1))

  wrap <- match.arg(wrap)
  if (wrap == "vertical") {
    p <- p + facet_grid(Type~.)
  } else if (wrap == "horizontal") {
    p <- p + facet_grid(.~Type)
  }
  p  
}



getAppendixChipEnrichments <- function(fits, qThresh = 0.05) {
  d <- Sys.glob("../input/PD_Padlock_Appendix_ChIP/")
  fAC <- file.path(d, "ACvsinput_peaks.narrowPeak_IDR0.05_conservative_final.bed")
  fME <- file.path(d, "MEvsinput_peaks.narrowPeak_IDR0.05_conservative_final.bed")

  myfits <- copy(fits)

  # AC marks
  dAC <- fread(fAC) %>%
          setnames(c("Chr", "Start", "End", "Name", "Score", "Strand")) %>% 
          .[, list(Chr, Start, End)] %>%
          setkey(Chr, Start, End)
  myfits[, AC := FALSE]
  ids <- foverlaps(
    myfits[, list(Chr, Start = SNP, End = SNP, ID)], 
    dAC, 
    by.x = c("Chr", "Start", "End"),
    by.y = c("Chr", "Start", "End")) %>%
    .[!is.na(Start), ID]
  myfits[ID %in% ids, AC := TRUE]

  # ME marks
  dME <- fread(fME) %>%
          setnames(c("Chr", "Start", "End", "Name", "Score", "Strand")) %>% 
          .[, list(Chr, Start, End)] %>%
          setkey(Chr, Start, End)
  myfits[, ME := FALSE]
  ids <- foverlaps(
    myfits[, list(Chr, Start = SNP, End = SNP, ID)], 
    dME, 
    by.x = c("Chr", "Start", "End"),
    by.y = c("Chr", "Start", "End")) %>%
    .[!is.na(Start), ID]
  myfits[ID %in% ids, ME := TRUE]

  # Now mark enhancers and promoters
  myfits %>%
    # Promoter is AC only
    .[, Promoter := ifelse( AC & !ME, TRUE, FALSE)] %>%
    # Poised enhancer is ME only
    .[, PoisedEnhancer := ifelse( !AC & ME, TRUE, FALSE)] %>%
    # Active enhancer is when both ME and AC are present
    .[, ActiveEnhancer := ifelse( ME & AC, TRUE, FALSE)] %>%
    # Drop unused cols
    .[, AC := NULL] %>%
    .[, ME := NULL]

  # Compute enrichments
  categories <- c("Promoter", "PoisedEnhancer", "ActiveEnhancer")
  pd <- 
    foreach (cat = categories, .combine = rbind) %do% {
      dt <- myfits[, 
              list(ID, Sig = adj.P.Val < qThresh, 
                  Dir = sign(logFC), CAT = get(cat))]
      A <- dt[, table(CAT, Sig)] %>%
        myFisherTest %>%
        .[, Category := cat %>% as.character] %>%
        .[, Type     := "All"]

      B <- dt[, table(CAT, Sig & (Dir == 1))] %>%
        myFisherTest %>%
        .[, Category := cat %>% as.character] %>%
        .[, Type     := "Hyper-"]

      C <- dt[, table(CAT, Sig & (Dir == -1))] %>%
        myFisherTest %>%
        .[, Category := cat %>% as.character] %>%
        .[, Type     := "Hypo-"]
      rbind(A, B, C)
    }
  pd  
}

# This is actually the same as Epigenome Roadmap data
plotAppendixChipEnrichment <- function(dt, 
  wrap = c("none", "vertical", "horizontal")) {
  plotEpiRoadmapEnrichment(dt, wrap) +
  ggtitle("Appendix ChIP")
}