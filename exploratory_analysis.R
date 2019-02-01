# source("https://bioconductor.org/biocLite.R")
# biocLite("ggtree")
# biocLite("phyloseq")
# biocLite("gridExtra")
# biocLite("dada2")
# biocLite("msa")
# biocLite("ggplot2")
# biocLite("plyr")
# biocLite("tidyverse")
# biocLite("DESeq2")
# biocLite("ggtree")
library(phyloseq)
library(RColorBrewer)
library(gridExtra)
library(dada2)
library(msa)
library(ggplot2)
library(plyr)
# library(phangorn)
library(tidyverse)
library(DESeq2)
library(ggtree)

setwd("M:/Projects/2018_MBCP/DADA2_outputs/")
experimentName=c("Name")
mainFactor="PMA_killed"
barcodePosition=7
plotFile=paste0("./",experimentName,".plots.pdf")
load(paste0("./",experimentName,".light.RData"))

################################################################################
# Function: plot abundance #####################################################
################################################################################

plot_abundance = function(physeq,
                          meta = sampledata.df$sampleNumber,
                          title = "",
                          FacetType = "wrap",
                          FacetBy = "Phylum",
                          Color = "Order",
                          legend.position = "none"){
  p1f = physeq ## Default - no subsetting.
  # Arbitrary subset, e.g., by Phylum, for plotting
  # p1f= subset_taxa(physeq, Phylum %in% c("p__Ascomycota"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  
  p <- ggplot(data = mphyseq,
         mapping = aes_string(x = meta,
                              y = "Abundance",
                              color = Color,
                              fill = Color)) +
    geom_point(size = 1,
               alpha = 0.3,
               position = position_jitter(width=0.3)) +
    scale_y_log10() +
    theme(legend.position = legend.position,
          axis.text.x = element_text(angle = 90)) +
    ggtitle(title) +
    theme(strip.text.y = element_text(angle = 0))
    if (FacetType == "wrap") {
      p + facet_wrap(facet = FacetBy)
    } else if (FacetType == "grid") {
      p + facet_grid(facet = FacetBy)
    } else { 
      print("Unknown FacetType! Use grid or wrap.")
    }
}
## END ABUNDANCE FUNCTION ######################################################
################################################################################

################################################################################
# PLOT MDS FUNCTION ############################################################
################################################################################
plot_my_MDS <- function(plotType,
                        myColor=mainFactor,
                        dist_methods=dist_methods,
                        ordMethod="PCoA",
                        wrapBy="distance") {
  
  print(paste0("Making ", ordMethod, " plots"))
  plist <- NULL
  for ( i in dist_methods ) {
    iDist <- distance(ps, method=i)
    myOrd <- tryCatch( 
      {
        message(paste0("Attempting ordination using ",i,"..."))
        ordinate(ps, ordMethod, distance=iDist)
      },
      error=function(cond)
      {
        message("Error!!! Returning NULL for this ordination.")
        message(cond)
        # Choose a return value in case of error
        return(NULL)
      },
      warning=function(cond)
      {
        message("Warning!!! This should still produce an ordination object though.")
        message(cond)
        # Choose a return value in case of warning
        return(NULL)
      },
      finally={ message(paste0("End ordination method ", i))  }
    ) # END TRYCATCH
    
    if ( is.null(myOrd) ) {
      p <- NULL
      plist[[i]] <- NULL
      print("Nothing here")
    } else {
      p <- NULL
      p <- plot_ordination(ps, myOrd, color=myColor, type=plotType)
      p <- p + ggtitle(paste(ordMethod, " using distance method ", i, sep=""))
      plist[[i]] = p
    }
  }
  
  message("Finished dist methods loop!")
  message("Making data frame. Dimensions:")
  df <- NULL
  df = plyr::ldply(plist, function(x) x$data)
  print(dim(df))
  message("Naming first column name as distance...")
  names(df)[1] <- "distance"
  message("Removing all NA entries from data frame for wrapping variable...")
  df <- df[!is.na(df[,wrapBy]),]
  message("Dimensions are now:")
  print(dim(df))
  message("Starting ggplot...")
  summaryPlot = ggplot(df,  aes(x=Axis.1, y=Axis.2, color=!!ensym(myColor))) # Quasiquotation instead of aes_string!
  message("Adding points...")
  summaryPlot = summaryPlot + geom_point(size=1, alpha=0.5)
  message("Adding facet wrap...")
  summaryPlot = summaryPlot + facet_wrap(vars(!!ensym(wrapBy)), scales="free")
  message("Adding title...")
  summaryPlot = summaryPlot + ggtitle(paste0(ordMethod, " using various distance metrics"), experimentName)
  message("Done")
  print(summaryPlot)
}
# END MDS PLOTTING FUNCTION ####################################################
################################################################################

################################################################################
# Subset to each V region function #############################################
################################################################################
split.to.regions <- function(ps=ps,region="V_region") {
  regions <- levels(ps@sam_data[[region]])
  
  # Initialize lists to store new Phyloseq objects
  ps_v_region_list <- NULL
  ps_v_region_list_transformed <- NULL
  ps_v_region_list_genus_glom <- NULL
  ps_v_region_list_genus_glom_transformed <- NULL
  ps_v_region_list_genus_glom_NArm <- NULL
  ps_v_region_list_genus_glom_NArm_transformed <- NULL
  
  ps_v_region_list <- vector("list", length(regions))
  ps_v_region_list_transformed <- vector("list", length(regions))
  ps_v_region_list_genus_glom <- vector("list", length(regions))
  ps_v_region_list_genus_glom_transformed <- vector("list", length(regions))
  ps_v_region_list_genus_glom_NArm <- vector("list", length(regions))
  ps_v_region_list_genus_glom_NArm_transformed <- vector("list", length(regions))
  
  for (i in 1:length(regions)) {
    print(regions[i])
    print(region)
    print(i)
    deparse(substitute(region))
    message("Subsetting phyloseq object into regions...")
    
    toPrune <- ps@sam_data[[region]]==regions[i]
    ps_v_region_list[[i]] <- prune_samples(toPrune, ps)
    names(ps_v_region_list[i]) <- regions[i]
    message("Computing relative abundance...")
    ps_v_region_list_transformed[[i]] <- transform_sample_counts(ps_v_region_list[[i]],
                                                                 function(x) 100*( x/sum(x)))
    names(ps_v_region_list_transformed[i]) <- regions[i]
    message("Agglomerating at genus level...")
    ps_v_region_list_genus_glom[[i]] <- tax_glom(ps_v_region_list[[i]],
                                                 taxrank="Genus",
                                                 NArm=FALSE)
    names(ps_v_region_list_genus_glom[i]) <- regions[i]
    message("Computing relative abundance of genus agglomerated...")
    ps_v_region_list_genus_glom_transformed[[i]] <- transform_sample_counts(ps_v_region_list_genus_glom[[i]],
                                                                            function(x) 100*( x/sum(x)))
    names(ps_v_region_list_genus_glom_transformed[i]) <- regions[i]
    message("Agglomerating at genus level, keeping NAs...")
    ps_v_region_list_genus_glom_NArm[[i]] <- tax_glom(ps_v_region_list[[i]],
                                                      taxrank="Genus",
                                                      NArm=TRUE)
    names(ps_v_region_list_genus_glom_NArm[i]) <- regions[i]
    message("Computing relative abundance of genus agglomerated, NAs kept...")
    ps_v_region_list_genus_glom_NArm_transformed[[i]] <- transform_sample_counts(ps_v_region_list_genus_glom_NArm[[i]],
                                                                                 function(x) 100*( x/sum(x)))
    names(ps_v_region_list_genus_glom_NArm_transformed[i]) <- regions[i]
  }
  
  myList <- list(regions=ps_v_region_list,
              regions_transformed=ps_v_region_list_transformed,
              regions_genus_glom=ps_v_region_list_genus_glom,
              regions_genus_glom_transformed=ps_v_region_list_genus_glom_transformed,
              regions_genus_glom_NArm=ps_v_region_list_genus_glom_NArm,
              regions_genus_glom_NArm_transformed=ps_v_region_list_genus_glom_NArm_transformed
              )
  return(myList)
}
## END REGION SPLITTING FUNCTION ###############################################
################################################################################

##################################
#### Create sample data table ####
##################################

# Not all experiments have the same samples (e.g., one V region may have no data)
# Start with sample names from seqtab.nochim
# Split on "." as delimiter
sampleList <- strsplit(row.names(seqtab.nochim), split="\\.")
# Sample number as taken from the name (i.e., use the barcode)
sampleNumber <- sapply(sampleList,"[[",barcodePosition) #### THIS DETERMINES WHICH PART OF THE NAME IS USED
# This creates the metadata dataframe
# First column of the table is essentially the barcode
sampledata.df$sampleNumber <- sampleNumber
# Make a lookup table to add metadata to correct samples
lookup <- data.frame(sampleNumber=unique(sampleNumber))
lookup$PMA <- c("No PMA","No PMA","No PMA","No PMA","No PMA","No PMA","PMA","PMA","PMA","PMA","PMA","PMA")
lookup$Killed <- c("Killed","Killed","Killed","Unkilled","Unkilled","Unkilled","Unkilled","Unkilled","Unkilled","Killed","Killed","Killed")
# Add in the metadata to the data frame
sampledata.df <- (merge(lookup, sampledata.df, by = 'sampleNumber'))
sampledata.df <- sampledata.df[order(sampledata.df$number),]
row.names(sampledata.df) <- row.names(seqtab.nochim)
sampledata.df$PMA_killed <- paste(sampledata.df$PMA,sampledata.df$Killed)

### Perform a sanity check on this!!!! Critical to have correct sample naming.
# row.names(seqtab.nochim) <- row.names(sampledata.df)
sampleOrder=row.names(sampledata.df)


###################################################
# Open a PDF file to record all plots as we go... #
###################################################

pdf(file=plotFile, width = 8.5, height = 11)

####################################
#### Run Phyloseq analysis, 16S ####
####################################

# Create phyloseq object of raw data
ps <- phyloseq(tax_table(taxa_silva),
               sample_data(sampledata.df),
               otu_table(seqtab.nochim,taxa_are_rows=FALSE),
               phy_tree(fitGTR$tree))

# Remove Unknown and Uncharacterized
ps0 <- subset_taxa(ps,  !is.na(Phylum) &! Phylum %in% c("", "uncharacterized"))

# Phyla before and after removal of unknown phyla
table(tax_table(ps)[, "Phylum"], exclude=NULL)
table(tax_table(ps0)[, "Phylum"], exclude=NULL)

# Compute prevalence of each feature, store as data frame
prevdf = apply(X = otu_table(ps0), MARGIN = ifelse(taxa_are_rows(ps0), yes =1, no =2), FUN = function(x) {sum(x>0)})
# Add taxonomy and total read counts to data frame
prevdf = data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(ps0),tax_table(ps0))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

#### If desired, manually remove certain taxa
# filterPhyla = c("Woesearchaeota", "Microgenomates", "Euryarchaeota", "Diapherotrites")
# ps1=subset_taxa(ps0, !Phylum %in% filterPhyla)
ps1=ps0 # Default, no filtering will be applied...

# Filter by prevalence
prevdf1 = subset(prevdf, Phylum  %in% get_taxa_unique(ps1,"Phylum"))
# Plot prevalance by abundance for each phylum
# This can help guide removal of taxa or setting a prevalence threshold
ggplot(prevdf1,
       aes(TotalAbundance, Prevalence/nsamples(ps1), color=Phylum)) + 
  geom_hline(yintercept=0.01, alpha=0.5, linetype=2) +
  geom_point(size=2, alpha=0.7) +
  scale_x_log10() +
  xlab("Total Abundance") +
  ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~ Phylum) +
  theme(legend.position="none") +
  ggtitle("Prevalence by abundance",experimentName)

# Prevalence filtering...
prevalenceThreshold = 0.01*nsamples(ps1)
keepTaxa=rownames(prevdf1)[(prevdf1$Prevalence>=prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps1)
length(get_taxa_unique(ps1,taxonomic.rank="Genus"))
length(get_taxa_unique(ps2,taxonomic.rank="Genus"))
# Removal of taxa observed less than 4 times in 10% of samples
# For data produced using the Ion Torrent 16S kit, it might be worth
# considering subsetting to only one V region at a time - e.g., V4
# ps_filtered <- filter_taxa(ps, function(x) sum(x > 3) > (0.2*length(x)), TRUE)


#################
# Plot richness #
#################

plot_richness(ps,
              x="PMA_killed",
              shape="PMA",
              color="V_region",
              measures=c("Observed", "Shannon", "Simpson", "Chao1"),
              title=paste0("Raw species richness for ",experimentName)) + 
  geom_point(size=5)

plot_richness(ps0,
              x="PMA_killed",
              shape="PMA",
              color="V_region",
              measures=c("Observed", "Shannon", "Simpson", "Chao1"),
              title=paste0("Species richness for ",experimentName,", unknown phyla removed")) + 
  geom_point(size=5)

plot_richness(ps2,
              x="PMA_killed",
              shape="PMA",
              color="V_region",
              measures=c("Observed", "Shannon", "Simpson", "Chao1")) + 
   geom_point(size=5) +
   ggtitle("Species richness after removal of low-prevalance samples",experimentName)

p1=plot_richness(ps,
              x="PMA_killed",
              shape="Killed",
              color="PMA",
              measures=c("Observed", "Shannon", "Simpson", "Chao1"))

ggplot(p1$data, aes(x=PMA_killed, y=value)) + 
  geom_point(size=3) + 
  facet_grid(variable~V_region) +
  ggtitle("Raw species richness, separated by V region",experimentName) +
  p1$mapping +
  theme(axis.text.x = element_text(angle=-90))

##############
# Plot trees #
##############

# Agglomerate at genus level, remove NAs
ps_glom_genus=tax_glom(ps2,"Genus" , NArm=TRUE)
ps_glom_genus_ra = transform_sample_counts(ps_glom_genus, function(x){x/sum(x)})

h1=0.2
ps_glom_height=tip_glom(ps2, h = h1)

multiPlotTitleTextSize=12

#plot_tree(ps2,method="treeonly",
#                  ladderize="left") + 
#  ggtitle("No agglomeration", experimentName) +
#  theme(plot.title=element_text(size= multiPlotTitleTextSize))

#plot_tree(ps_glom_genus,method="treeonly",
#                  ladderize="left") +
#  ggtitle("By Genus", experimentName) +
#  theme(plot.title=element_text(size= multiPlotTitleTextSize))

#plot_tree(ps_glom_height,method="treeonly",
#                 ladderize="left") +
#  ggtitle("By Height", experimentName)
#  theme(plot.title=element_text(size=multiPlotTitleTextSize))

# Plot together on one row
# grid.arrange(nrow=1, p2tree, p3tree, p4tree)

## Better Trees ##
  
ps_glom_species_na_kept <- tax_glom(ps, taxrank="Species", NArm=FALSE)
ps_glom_species_na_kept_transformed <-transform_sample_counts(ps_glom_species_na_kept, function(x) 100*( x/sum(x)))
  
ps_glom_genus_na_kept <- tax_glom(ps, taxrank="Genus", NArm=FALSE)
ps_glom_genus_na_kept_transformed <-transform_sample_counts(ps_glom_genus_na_kept, function(x) 100*( x/sum(x)))
  
plot_tree(ps_glom_species_na_kept_transformed,
            method="sampledodge",
            size="Abundance",
            justify="yes please",
            ladderize="left",
            color="PMA_killed") +
    scale_size_continuous(range = c(1, 3)) +
    ggtitle("Species-agglomerated tree, relative abundance, NAs kept",experimentName)
  

plot_tree(ps_glom_genus_na_kept_transformed,
          color="PMA_killed",
          label.tips="Phylum",
          method="sampledodge",
          justify="yes please",
          ladderize="left")

plot_tree(ps_glom_genus_na_kept_transformed,
          color="PMA_killed",
          method="sampledodge",
          justify="yes please",
          ladderize="left",
          size="abundance")

plot_tree(ps_glom_genus_na_kept_transformed,
            method="sampledodge",
            size="Abundance",
            justify="yes please",
            ladderize="left",
            color="PMA_killed") +
    scale_size_continuous(range = c(1, 3)) +
    ggtitle("Genus-agglomerated tree, relative abundance, NAs kept",experimentName)

###################
# Plot Abundances #
###################

ps2ra = transform_sample_counts(ps2, function(x){x/sum(x)})

# Non-agglomerated
plot_abundance(ps1,
               meta="PMA_killed",
               title="Asolute Abundance, PS1 (unfiltered)",
               FacetBy="Class")
plot_abundance(ps2,
               meta="PMA_killed",
               title="Absolute Abundance, PS2 (filtered by prevalence threshold)",
               FacetBy="Class")
plot_abundance(ps2ra,
               meta="PMA_killed",
               title="Relative Abundance, PS2 (filtered by prevalence threshold)",
               FacetBy="Class")
plot_abundance(ps2ra,
               meta="PMA_killed",
               title="Relative Abundance, PS2 (filtered by prevalence threshold), by V region and Class",
               FacetType="grid",
               FacetBy="Class~V_region")

# Agglomerated
plot_abundance(ps_glom_genus,
               meta="PMA_killed",
               title="Absolute abundance by treatment, Agglomerated at genus level",
               FacetBy="Class")
plot_abundance(ps_glom_genus_ra,
               meta="PMA_killed",
               title="Relative abundance by treatment, Agglomerated at genus level",
               FacetBy="Class")
plot_abundance(ps_glom_genus,
               meta="V_region",
               title="Absolute abundance by V region, Agglomerated at genus level",
               FacetBy="Class")
plot_abundance(ps_glom_genus_ra,
               meta="V_region",
               title="Relative abundance by V region, Agglomerated at genus level",
               FacetBy="Class")

plot_abundance(ps_glom_genus,
               meta="PMA_killed",
               title="Absolute abundance by V region, separated by Class, Agglomerated at genus level",
               FacetType="grid",
               FacetBy="Class~V_region")
plot_abundance(ps_glom_genus_ra,
               meta="PMA_killed",
               title="Relative abundance by V region, separated by Class, Agglomerated at genus level",
               FacetType="grid",
               FacetBy="Class~V_region")

 dev.off()

################################################################
#### GO BACK TO ORIGINAL PS OBJECT FOR STATISTICAL ANALYSIS ####
################################################################

##################
### ORDINATION ###
##################
 
#######################################################
#### Multiple distance metrics combined in single plots
ps.bak <- ps
# Remove samples with zero abundance
ps <- prune_samples(sample_sums(ps)>=1, ps)

dist_methods <- unlist(distanceMethodList)
print(dist_methods)
dist_methods = dist_methods[-which(dist_methods=="ANY")]
dist_methods = dist_methods[-which(dist_methods=="manhattan")]
print(dist_methods)
# plist <- vector("list", length(dist_methods)) #### LENGTH MUST AGREE WITH WHAT YOU DO BELOW



pdf(file="Ordination.pdf", width=12, height=12)

ordMethods=c("PCoA","MDS","NMDS","DCA","RDA","CCA","DPCoA")
ordMethods=c("MDS","PCoA")

for (ordMethod in ordMethods) {
    
    print(ordMethod)
    print(myColor)
    
    # Samples
    plot_my_MDS( plotType="samples",
                 myColor=mainFactor,
                 dist_methods=dist_methods,
                 ordMethod=ordMethod)
    
    plot_my_MDS( plotType="samples",
                 myColor=mainFactor,
                 dist_methods=dist_methods,
                 ordMethod=ordMethod,
                 wrapBy="V_region")
    
    plot_my_MDS( plotType="samples",
                 myColor="V_region",
                 dist_methods=dist_methods,
                 ordMethod=ordMethod)
    
    # Taxa
    plot_my_MDS( plotType="taxa",
                 myColor="Phylum",
                 dist_methods=dist_methods,
                 ordMethod=ordMethod)
    
    plot_my_MDS( plotType="taxa",
                 myColor="Order",
                 dist_methods=dist_methods,
                 ordMethod=ordMethod)
    
    # Split
    plot_my_MDS( plotType="split",
                 myColor=mainFactor,
                 dist_methods=dist_methods,
                 ordMethod=ordMethod)
    
    plot_my_MDS( plotType="split",
                 myColor="Phylum",
                 dist_methods=dist_methods,
                 ordMethod=ordMethod)
    
    plot_my_MDS( plotType="split",
                 myColor="Order",
                 dist_methods=dist_methods,
                 ordMethod=ordMethod)

}

dev.off()
ps <- ps.bak

################################################################################


ps_v_region_objects <- split.to.regions(ps)

plot_bar(ps_glom_genus_na_kept_transformed, fill="Family", x="PMA_killed")



top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
#ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
#ps.top20 <- transform_sample_counts(ps, function(x) 100*(x/sum(x)))
ps.top20 <- prune_taxa(top20, ps.top20)

getPalette = colorRampPalette(brewer.pal(11, "Spectral"))
getPalette = colorRampPalette(brewer.pal(9, "Set1")) # not bad - but may need to rearrange
getPalette = colorRampPalette(brewer.pal(8, "Set2")) # too muted
getPalette = colorRampPalette(brewer.pal(8, "Set3"))
getPalette = colorRampPalette(brewer.pal(8, "Accent")) # Ugly
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

plot_bar(ps.top20, x="Total_treatment", fill="Genus") + 
  labs(x="Silver Concentration") + 
  facet_grid(~PMA_killed,scales="free") +
  scale_fill_manual(values=getPalette(22))

top100 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:100]
ps.top100 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top100 <- prune_taxa(top100, ps.top100)
plot_bar(ps.top100, x="Sample", fill="Genus") + labs(x="Silver Concentration") + facet_grid(~PMA_killed,scales="free")



#######################################
# Close the PDF file to write to disk #
#######################################
dev.off()



### Differential abundance, adapted for data from dada2
ps_pruned <- prune_samples(sample_sums(ps) > 500, ps)
dds = phyloseq_to_deseq2(ps_pruned, ~ PMA_killed)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)
dds = DESeq(dds, fitType="local")
res = results(dds)

res60 <- results(dds, contrast=c("PMA_killed","60","0"), cooksCutoff = FALSE)
res145 <- results(dds, contrast=c("PMA_killed","145","0"), cooksCutoff = FALSE)
res347 <- results(dds, contrast=c("PMA_killed","347","0"), cooksCutoff = FALSE)
res833 <- results(dds, contrast=c("PMA_killed","833","0"), cooksCutoff = FALSE)
res2000 <- results(dds, contrast=c("PMA_killed","2000","0"), cooksCutoff = FALSE)

alpha = 0.01

sigtab60 = res60[which(res60$padj < alpha), ]
sigtab145 = res145[which(res145$padj < alpha), ]
sigtab347 = res347[which(res347$padj < alpha), ]
sigtab833 = res833[which(res833$padj < alpha), ]
sigtab2000 = res2000[which(res2000$padj < alpha), ]

sigtab60 = cbind(as(sigtab60, "data.frame"), as(tax_table(ps)[rownames(sigtab60), ], "matrix"))
sigtab145 = cbind(as(sigtab145, "data.frame"), as(tax_table(ps)[rownames(sigtab145), ], "matrix"))
sigtab347 = cbind(as(sigtab347, "data.frame"), as(tax_table(ps)[rownames(sigtab347), ], "matrix"))
sigtab833 = cbind(as(sigtab833, "data.frame"), as(tax_table(ps)[rownames(sigtab833), ], "matrix"))
sigtab2000 = cbind(as(sigtab2000, "data.frame"), as(tax_table(ps)[rownames(sigtab2000), ], "matrix"))




# Plot results
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

sigtableList <- NULL
sigtableList <- list(sigtab60,sigtab145,sigtab347,sigtab833,sigtab2000)

for (sigtable in sigtableList) {
  
  # Phylum order
  x = tapply(sigtable$log2FoldChange, sigtable$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sigtable$Phylum = factor(as.character(sigtable$Phylum), levels=names(x))
  # Genus order
  x = tapply(sigtable$log2FoldChange, sigtable$Genus, function(x) max(x))
  x = sort(x, TRUE)
  sigtable$Genus = factor(as.character(sigtable$Genus), levels=names(x))
  
}

sigtableListNA.RM <-lapply(sigtableList, subset, !is.na(Genus))

combined_sigtableListNA.RM <- do.call(rbind, sigtableListNA.RM)


#ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
#  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

## Genus plots, no N/A
ggplot(subset(sigtab833, !is.na(Genus)), aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))



# Make phyloseq object from the differentially abundant genera
dataSamplesPruned <- prune_samples(sample_sums(ps)>=1000, ps)
dataSamplesPrunedRarefied <- rarefy_even_depth(dataSamplesPruned)
dataSamplesPrunedRarefiedSpeciesGlom <- tax_glom(dataSamplesPrunedRarefied, taxrank="Species", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
dataSamplesPrunedRarefiedSpeciesGlom <- tax_glom(dataSamplesPrunedRarefied, taxrank="Species", NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))
dataSamplesPrunedRarefiedGenusGlom <- tax_glom(dataSamplesPrunedRarefied, taxrank="Genus", NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))

plot_heatmap(dataSamplesPrunedRarefiedSpeciesGlom, "RDA", "unifrac", "PMA_killed", "Genus",    sample.order=row.names(sampledata.df))

dataSamplesSpeciesGlomRarefied <- tax_glom(dataSamplesPrunedRarefied, taxrank="Species", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
# genera_affected <- subset_taxa(ps, Genus %in% levels(subset(sigtab, !is.na(Genus))$Genus))
genera_affected_sp_glom <- subset_taxa(dataSamplesPrunedRarefiedSpeciesGlom, Genus %in% levels(subset(combined_sigtableListNA.RM, !is.na(Genus))$Genus))
genera_affected_sp_glom_rarefied <- subset_taxa(dataSamplesSpeciesGlomRarefied, Genus %in% levels(subset(combined_sigtableListNA.RM, !is.na(Genus))$Genus))

plot_heatmap(genera_affected_sp_glom, "NMDS", "jaccard", "PMA_killed", "Genus", sample.order=sampleOrder) ## Nice

plot_heatmap(genera_affected_sp_glom, "RDA", "none", "PMA_killed", "Genus", sample.order=sampleOrder)
plot_heatmap(genera_affected_sp_glom_rarefied, "MDS", "unifrac", "PMA_killed", "Genus",    sample.order=sampleOrder)
plot_heatmap(dataSamplesPrunedRarefiedSpeciesGlom, "RDA", "unifrac", "PMA_killed", "Genus",    sample.order=sampleOrder)


# ORDERING TAXA
sigtab2000$taxID <- row.names(sigtab2000)
sigtab833$taxID <- row.names(sigtab833)
sigtab347$taxID <- row.names(sigtab347)
sigtab145$taxID <- row.names(sigtab145)
sigtab60$taxID <- row.names(sigtab60)
combined_sigtableListNA.RM$taxID <- row.names(combined_sigtableListNA.RM)


orderedGenera833 <- as.character(arrange(sigtab833, log2FoldChange)$taxID[!is.na(arrange(sigtab833, log2FoldChange)$taxID)])
orderedGenera60 <- as.character(arrange(sigtab60, log2FoldChange)$taxID[!is.na(arrange(sigtab60, log2FoldChange)$taxID)])
orderedGenera145 <- as.character(arrange(sigtab145, log2FoldChange)$taxID[!is.na(arrange(sigtab145, log2FoldChange)$taxID)])
orderedGenera347 <- as.character(arrange(sigtab347, log2FoldChange)$taxID[!is.na(arrange(sigtab347, log2FoldChange)$taxID)])
orderedGenera2000 <- as.character(arrange(sigtab2000, log2FoldChange)$taxID[!is.na(arrange(sigtab2000, log2FoldChange)$taxID)])

orderedGeneraAll <- as.character(arrange(combined_sigtableListNA.RM, log2FoldChange)$taxID[!is.na(arrange(combined_sigtableListNA.RM, log2FoldChange)$taxID)])
orderedGeneraByPval <- as.character(arrange(combined_sigtableListNA.RM, padj)$taxID[!is.na(arrange(combined_sigtableListNA.RM, log2FoldChange)$taxID)])








### ADD THIS IN PLACE FOR THE HEATMAP TO ELIMINATE UNASSIGNED GENERA
# subset_taxa(dataSamplesPrunedRarefied,  !is.na(Genus)

plot_heatmap(subset_taxa(dataSamplesPrunedRarefied,  !is.na(Genus)), "none", "none", "PMA_killed", 
             taxa.label="Genus", sample.order=sampleOrder, taxa.order=orderedGeneraAll, max.label=5000) +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10)) +
  labs(x="Silver Concentration")











nitrifying_bacteria <- c("Nitrosomonas", "Nitrosococcus", "Nitrobacter", "Nitrococcus", "Nitrospina", "Nitrospira", "Nitrosophaera")
ps_nitrifying <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
#ps_nitrifying <- subset_taxa(ps, Genus %in% nitrifying_bacteria)
ps_nitrifying <- subset_taxa(ps_nitrifying, Genus %in% nitrifying_bacteria)
ps_nitrifying_glom <- tax_glom(ps_nitrifying,taxrank="Genus")

plot_heatmap(ps_nitrifying_glom,  "PMA_killed", 
             taxa.label="Genus", sample.order=sampleOrder,  max.label=5000) +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10)) +
  labs(x="Silver Concentration")
plot_bar(ps_nitrifying_glom, fill="Genus")
plot_bar(ps_nitrifying_glom, x="PMA_killed", fill="Genus")

taxdf <- as.data.frame(tax_table(ps))
taxdf$taxid <- row.names(taxdf)
ps_transformed <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
nitro_tax <- subset(taxdf, grepl("Nitro", taxdf$Kingdom) | grepl("Nitro", taxdf$Phylum) | grepl("Nitro", taxdf$Class) | grepl("Nitro", taxdf$Order) | grepl("Nitro", taxdf$Genus) | grepl("Nitro", taxdf$Species) )
nitro_tax_ids <- row.names(nitro_tax)
ps_nitro=prune_taxa(nitro_tax_ids, ps_transformed)
ps_nitro_glom <- tax_glom(ps_nitro, taxrank="Genus")
plot_bar(ps_nitro_glom, x="PMA_killed", fill="Genus")
plot_heatmap(ps_nitro,  "PMA_killed", 
             taxa.label="Genus", sample.order=sampleOrder,  max.label=5000) +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10)) +
  labs(x="Silver Concentration")




