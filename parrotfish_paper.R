# Leïla - data parrotfish

## load packages
library(phyloseq)
library(biomformat)
library(taxize)
library(ggplot2)
library(vegan)
library(knitr)
library(dplyr)
library(data.table)
library(plyr)
library(DescTools)
## mixed model with random effect
library(lme4)
library(lmerTest)
library("DESeq2")
# install.packages("RColorBrewer")
library(RColorBrewer)
library("RColorBrewer")
display.brewer.all()
library(gridExtra)
library(gtable)
library(rcompanion)
library(Rcmdr)
library(devtools)
library(wesanderson)
library(nestedRanksTest)
library(FSA)
library(theseus)
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
## try http:// if https:// URLs are not supported
library("colorspace") 
library(microbiome)

otu_table = read.csv("otu_table.csv",sep=";",row.names=1)
otu_table = as.matrix(otu_table_trans)

##read in taxonomy
#separated by kingdom, phylum, class, order, family, genus, species
taxonomy = read.csv("transmission_taxo.csv", sep=";",row.names=1)
taxonomy= as.matrix(taxonomy_transmission)

## Import QIIME pipeline map-file
metadata="/Users.txt"
metadata<-import_qiime_sample_data(metadata)

## Import tree
tree_ <- read_tree("/Users.nwk")

OTU <- otu_table(otu_table_trans, taxa_are_rows=TRUE)
tax <- tax_table(taxonomy)

## Merge everything (mapping file, tree and OTUs)
merged<- merge_phyloseq((OTU), tax, (metadata), tree)

## Assign a real ID to each OTU instead of keeping the sequences
taxa_names(merged) <- paste("OTU", 1:ntaxa(merged), sep="_")

sort(sample_sums(merged))
all_r<- rarefy_even_depth(merged, sample.size=min(sample_sums(merged)), rngseed=678, replace=FALSE, trimOTUs=TRUE, verbose=TRUE)

# Create tables for VEGAN 
vegan_otu <- function(physeq){
  OTU <- otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
## apply function
vegan_matrix_transmission <- vegan_otu(all_r)
otu_tab <- otu_table(all_r)
tax_tab <- tax_table(all_r)

## Chao1
p_oBiom <- plot_richness(all_r, measures="Chao1", color="time", shape="what") +
  theme(axis.text.x=element_text(angle=70, hjust=1)) +
  xlab("") + ylab("Chao1") + ggtitle("Alpha Diveristy")
p_oBiom

alphadt<- data.table(p_oBiom$data)

chao1 <- aov((alphadt$value) ~ what, data=alphadt)

R1_chao1 <- resid(chao1); F1_chao1 <- fitted(chao1)
plot(F1_chao1, R1_chao1, xlab = "Fitted Values", ylab = "Normalized residuals"); abline(h = 0, lty = 2) 
summary(chao1)
hist(R1_chao1) # its pretty good distribution?
qqnorm(R1_chao1)
shapiro.test(R1_chao1)
LeveneTest(R1_chao1 ~ what, data=alphadt)
tukeychao1<-TukeyHSD(chao1)
tukeychao1

## Shannon
p_oBiom_Shannon <- plot_richness(all_r, measures="Shannon", shape="state", color="site") +
  theme(axis.text.x=element_text(angle=70, hjust=1)) +
  xlab("") + ylab("Shannon") + ggtitle("Alpha Diveristy")
p_oBiom_Shannon
alphadte = data.table(p_oBiom_Shannon$data)

shannon <- aov(sqrt(alphadte$value) ~state, data=alphadte)
##non-parametric Shannon
shannon_kruskal<- kruskal.test(alphadte$value ~state, data=alphadte)
dunnTest(alphadte$value ~ state, data=alphadte)

##Beta-diversity
##Community composition
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'euclidean', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
} 

ado<-adonis(vegan_matrix_transmission ~ what, data = alphadt)
pairwise.adonis(vegan_matrix_transmission, alphadtn$what, sim.method = 'bray', p.adjust.m = 'fdr', options(max.print=10000))

##Betadisper
qd_bc_<- phyloseq::distance(vegan_matrix_transmission, "bray")
bdisp<- betadisper(qd_bc, (alphadt$what), type=c("centroid"))
bdisp
aov.bdisp <-anova(bdisp)
permutest(bdisp)

labs <- paste("Dimension", 1:8, "(", 
              round(100*bdisp$eig / sum(bdisp$eig), 2), "%)")

plot(bdisp.mouth, cex=1, pch=15:17,
     main="Bray", cex.lab=1.25,
     xlab=labs[1], ylab=labs[2],
     hull=FALSE, ellipse=TRUE, conf=0.68, lwd=2)

box_test <- cbind(alphadt, dist=bdisp$distances)
options(max.print=1000000)

ggplot(box_test, aes(x=what, y=dist))+
  geom_boxplot() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##DESeq2 (basic - could be adjusted according to data)
##choose your data subset - unrarefied sOTU table!!
choix_taxglom <- tax_glom(physeq=xxx, taxrank=rank_names(xxx)[6], NArm=F)

treatmentdes = c("control","coral")
timedes="Tx"

oBiom_corals<- subset_samples(choix_taxglom, what %in% treatmentdes)

diagdds = phyloseq_to_deseq2(oBiom_corals, ~ what)

diagdds$Vegetation <- relevel(diagdds$Vegetation, ref = "no")
diagdds = DESeq(diagdds, test="Wald", fitType="mean", sfType="poscount")
res = results(diagdds, cooksCutoff = FALSE, alpha = 0.01)

sigtab = cbind(as(res, "data.frame"), as(tax_table(oBiom_corals)[rownames(res), ], "matrix"))
sigtab[order(sigtab$log2FoldChange),] %>% head(10)


#### Relative Abundance 
##subset 
oBiom <- subset_samples(all_r, what=="xx")
ra_fam<- tax_glom(physeq=xx, taxrank=rank_names(xx)[5], NArm=F)

physeq_sample_merge = merge_samples(ra_fam, "what")

transii  = transform_sample_counts(physeq_sample_merge, function(x) x / sum(x))
otu_table(trans, options(max.print=10000))
plot_bar(trans, fill="Family") + ggtitle("Relative Abundance %") + 
  scale_fill_brewer(palette="Paired") + theme_bw() +scale_colour_manual(values=cbPalette)


