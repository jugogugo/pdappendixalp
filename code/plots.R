
plotHistograms <- function(p, observed, expected) {
	# P value histogram
	p1 <- ggplot(as.data.table(p), aes(p)) + 
		geom_histogram(breaks=seq(0, 1, 0.05), 
			fill="white", color="black") + 
		ggtitle("P value histogram") + 
		xlab("p-value") +
		theme_bw(base_size = 8)

	# Permutation 1
	p2 <- ggplot(expected, aes(P)) + 
		geom_histogram(bins=20, fill="white", color="black") + 
		geom_vline(data=observed, aes(xintercept=P), color="red") + 
		ggtitle(paste("Observed", format(observed$P, digits=2), "p =", 
			mean(expected$P >= observed$P)
			)) + 
		xlab("fraction of CGs with p < 0.05") + 
		theme_bw(base_size = 8)

	# Permutation 2
	p3 <- ggplot(expected, aes(PFDR)) + 
		geom_histogram(bins=20, fill="white", color="black") + 
		geom_vline(data=observed, aes(xintercept=PFDR), color="red") + 
		ggtitle(paste("Observed", format(observed$PFDR, digits=2), "p =", 
			mean(expected$PFDR >= observed$PFDR)
			)) + 
		xlab("number of FDR significant CGs") + 
		theme_bw(base_size = 8)

	# ggarrange(p1, 
	# 	ggarrange(p2, p3, ncol=2, nrow=1, labels=c("B", "C")),
	# 	labels = c("A"),
	#   	ncol = 1, nrow = 2)	

	ggarrange(p1, p2, p3, ncol=3, nrow=1, labels=c("A", "B", "C"))
}

plotManhattan <- function(chr, pos, p, coef, thresholds = NULL) {
	require(data.table)
	require(dplyr)
	require(ggplot2)


	pd <- data.table(Chr=gsub("chr", "", chr) %>% as.numeric, SNP=pos, P=p, S=coef)

	# significance thresholds
	if (is.null(thresholds)) {
		thresholds <- pd[p.adjust(P, "fdr") < 0.05, max(P)] %>% (function(x) {-log10(x)})
	}
	
	pd <- pd %>% 
		group_by(Chr) %>%
		summarize(chr_len = max(SNP)) %>%
		mutate(tot=cumsum(chr_len)-chr_len) %>%
		select(-chr_len) %>%
		left_join(pd, ., by=c("Chr"="Chr")) %>%
		arrange(Chr, SNP) %>%
		mutate( BPcum=SNP+tot) %>%
		setDT
	axisdf <- pd %>% 
		group_by(Chr) %>% 
		summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

	ggplot(pd, aes(x=BPcum, y=-log10(P) * sign(S), size=abs(S))) +
	    geom_point( aes(color=as.factor(Chr)), alpha=0.8, size=1.3) +
	    geom_hline(yintercept = abs(thresholds), color="red", size=0.3) + 
	    geom_hline(yintercept = -abs(thresholds), color="red", size=0.3) + 
	    geom_hline(yintercept = 0, color="black", size=0.3) + 
	    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
	    scale_x_continuous( label = axisdf$Chr, breaks= axisdf$center ) +
	    xlab("Chromosome") + 
	    ylab("Signed log P") + 
	    theme_bw(base_size=8) +
	    theme( 
	      legend.position="none",
	      panel.border = element_blank(),
	      panel.grid.major.x = element_blank(),
	      panel.grid.minor.x = element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1)
	    )
}

plotVolcano <- function(p, coef) {
	require(data.table)
	require(dplyr)
	require(ggplot2)

	data.table(P=p, S=coef) %>%
	ggplot(., aes(S, -log10(P), color=p.adjust(P, "fdr") <= 0.05)) + 
	geom_point(alpha=0.3) + 
	scale_color_brewer("FDR ", palette="Set1", labels=c("> 0.05", "<= 0.05"), direction = -1) + 
	xlab("Coefficient value") + ylab("-log10(P)") + 
	theme_bw(base_size=8) + 
	guides(size=FALSE) + 
	theme(legend.position = "topenv")
}