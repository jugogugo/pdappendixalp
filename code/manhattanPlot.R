if (basename(getwd()) == "code") {
	source("pformat.R")
} else {
	source("../code/pformat.R") 
}

manhattanPlot <- function(fit, qThresh = 0.05, labelGenes = TRUE) {

  fit <- copy(fit)
  # Mark the fit points which have to be labeled
  geneList <- c("GBA","NUCKS1","SLC41A1","SIPA1L2","TMEM163","CCNT2","STK39","CHMP2B","MCCC1","TMEM175","DGKQ","FAM200B","CD38","FAM47E","SNCA","HLA-DRB6","HLA-DQA1","KLHL7","NUPL2","GPNMB","MICU3","BAG3","DLG2","MIR4697","LRRK2","OGFOD2","GCH1","TMEM229B","VPS13C","ZNF646","KAT8","ARHGAP27","CRHR1","SPPL2C","MAPT","STH","KANSL1","SYT4","LSM7","DDRGK1","ITPKB","IL1R2","SCN3A","SATB1","NCKIPSD","CDC71","ALAS1","TLR9","DNAH1","BAP1","PHF7","NISCH","STAB1","ITIH3","ITIH4","ANK2","CAMK2D","ELOVL7","ZNF184","CTSB","SORBS3","PDLIM2","C8orf58","BIN3","SH3GL2","FAM171A1","GALC","COQ7","TOX3","ATP6V0A1","PSMC3IP","TUBG2")

  foo <- function(snps, pvals) {
    snps[ which.min(pvals)]
  }

  foo2 <- function(N, SNP) {
    i <- which.max(N)
    return(SNP[i])
  }


  getColorGroup <- function(Chr, Q, S) {
    Chr <- as.numeric(Chr)
    Q <- as.numeric(Q)
    S <- as.numeric(S)
    df <- data.frame(Chr, Q, S)
    apply(df, 1, function(I) {
      I <- as.list(I)
      if (I$Chr %% 2 == 0) {
        # Darker colors
        if (I$Q > qThresh) {
          "dark grey"
        } else {
          ifelse(I$S > 0, "dark positive", "dark negative")
        }
      } else {
        # Brighter colors
        if (I$Q > qThresh) {
          "light grey"
        } else {
          ifelse(I$S > 0, "light positive", "light negative")
        }
      }
    }) %>% unlist
  }

  if (labelGenes == TRUE) {
    effects <- 
      fit[ toupper(Gene) %in% geneList & adj.P.Val < qThresh, ] %>% 
        # count how many times each gene is negative of positive
        .[, list(.N, SNP = foo(SNP, P.Value), P = min(P.Value)), by=list(Gene, sign(logFC))] %>%
        # select the sign which has more representation
        .[, list(SNP=foo2(N, SNP)), Gene] %>%
        .[, Label := TRUE] %>% 
        merge(fit, ., by=c("Gene", "SNP"), all.x=TRUE)
    effects[is.na(Label), Label := FALSE]

    effects <- effects[Chr != "chrM"]    
  } else {
    effects <- fit[Chr != "chrM"] %>%
                .[, Label := FALSE] %>%
                .[, Gene := ""]
  }


  # Format the input data
  effects <- effects %>%
    .[, Diff := logFC ] %>%
    .[, SNP := strsplit(ID, split="_") %>% sapply(., '[[', 2) %>% as.numeric] %>%
    .[, P := P.Value] %>%
    .[, Q := adj.P.Val] %>%
    # Fix chromosome name 
    .[, Chr := gsub("chr", "", Chr) ] %>% 
    # Give chromosomes numbers
    .[, ChrNo :=  Chr %>% gsub("X", "23", .) %>% gsub("Y", "24", .) %>% gsub("M", "25", .) %>%as.numeric]  %>%
    setkey(ChrNo, SNP)

  # Prepare plot data
  effects <- effects[!is.na(Q)]
  pd <- effects %>% 
    # Subsample 
    .[ Q < 0.05 | ID %in% sample(ID, .N*0.3, prob=-log10(P))] %>%  
    # Compute chromosome size
    .[, list(ChrLen = max(SNP)), ChrNo] %>%
    # Calculate cumulative position of each chromosome
    .[, list(ChrNo, Tot = cumsum(ChrLen) - ChrLen)] %>%
    # Add this info to the initial dataset
    merge(effects, ., by=c("ChrNo")) %>%
    # Sort by ChrNo and position 
    .[ order(ChrNo, SNP), ] %>%
    # Compute position for each SNP on x axis
    .[, Loc := SNP + Tot] %>%
    # Get color group
    .[, Color := getColorGroup(ChrNo, Q, Diff) %>% as.factor ] %>%
    setkey(Loc)


  # Prepare the x axis
  pdAxis <- pd %>% 
    .[, list(center = ( max(Loc) + min(Loc) ) / 2 ), by=list(Chr, ChrNo)] %>%
    # omit some chromosomes
    .[, Chr := ifelse(ChrNo %in% seq(15, 21, 2), "", Chr)]


  # Prepare the core plot
  p1 <- pd %>% 
    ggplot(., aes(x=Loc, y=-1 * sign(Diff) * log10(P))) +
      # Show all points
      geom_point( aes(color=Color), alpha=0.5, size=0.5) + 
      # Plot bigger target points?
      geom_point(data = pd[Q < qThresh], aes(color=Color), size = 1) + 
      # custom X axis:
      scale_x_continuous( label = pdAxis$Chr, breaks= pdAxis$center ) +
      scale_color_manual( values = c("grey50", "#1f78b4", "#33a02c", 
                                     "grey", "#a6cee3","#b2df8a")) +
      guides(color = FALSE)
    

  # Add labels
  if (!is.null(pd$Label) & sum(pd$Label) > 0) {
    require("ggrepel")
    # Add highlighted points
    p1 <- p1 + 
      geom_point(data=pd[Label == TRUE], color="orange", size=1) +
      geom_label_repel(data=pd[Label == TRUE], aes(label=Gene), size=2.5)
  }

  # Label the axes
  p1 <- p1 + 
    xlab("Chromosome") + 
    ylab("SLP")


  # Enrichment barplot
  f <- 
  tryCatch({
    t <- fit[, table(
      Hypo=sign(logFC)<0, 
      Significant = adj.P.Val < qThresh)]
    f <- t %>% fisher.test    
    f$estimate <- format(f$estimate, digits=2, scientific=FALSE)
    f$p.value <- sprintf("p = %s", texP(f$p.value))
    f
  }, error = function(e) {
    return(list(
      estimate = 1,
      p.value = NA
      ))
  })

  p2 <- t %>% reshape2::melt() %>% as.data.table %>% 
    .[, list(Perc = value/sum(value), Hypo), by=Significant] %>%
    .[, Significant := factor(Significant, levels=c(FALSE, TRUE), labels=c("Background", "Significant"))] %>% 
    .[, Hypo := factor(Hypo, levels = c(FALSE, TRUE), labels=c("Hyper", "Hypo"))] %>% 
    .[Significant == "Significant"] %>% 
    ggplot(., aes(Hypo, Perc, fill=Hypo)) + 
    geom_bar(stat="identity") + 
    # scale_y_continuous(labels = scales::percent, limit=c(0, 0.75)) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual("", values=c("#33a02c", "#1f78b4")) + 
    xlab("-methylated") + 
    ylab("mC Proportion, %") + 
    ggtitle("Fisher test", latex2exp::TeX(f$p.value)) + 
    # ggtitle(glue::glue("Fisher test\np={texP(f$p.value)}")) + 
    guides(fill=FALSE)


  # Set theme
  p1 <- p1 +  
    theme_bw(base_size=8) +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

  p2 <- p2 + 
    theme_bw(base_size=8) + 
    theme( 
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  return(list(p1, p2))
}
