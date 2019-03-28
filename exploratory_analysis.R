### Author: Matt Meier
### This script will take an RData file with the output from DADA2
### and produce some useful plots using phyloseq as well as DESeq2
### Experiment-specific modification is required!!!

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
# library(phangorn)
#library(ggtree)
library(phyloseq)
library(RColorBrewer)
library(gridExtra)
library(dada2)
library(msa)
library(ggplot2)
library(plyr)
library(DESeq2)

#######################################################
## CONFIGURATION PARAMETERS REQUIRED TO RUN PROPERLY ##
## USER SUPPLIES THE FOLLOWING: #######################
setwd("/path/to/files")
experimentName=c("Name")
mainFactor="Factor"
myLevels=c("Level1",
           "Level2",
           "Level3",
           "Level4")
barcodePosition=7
regionPosition=12
#######################################################
#######################################################
## Non-user supplied values
timeOfScript=format(Sys.time(), "%Y-%m-%d.%H.%M")
plotFile=paste0("./",experimentName,".",timeOfScript,".plots.pdf")
load(paste0("./",experimentName,".light.RData"))

## For optimizing plot aesthetics
update_geom_defaults("point", list(size=1))
theme_set(theme_grey(base_size=10))

################################################################################
# Color palettes for plots... ##################################################
################################################################################

# getPalette = colorRampPalette(brewer.pal(11, "Spectral"))
# getPalette = colorRampPalette(brewer.pal(9, "Set1")) # not bad - but may need to rearrange
# getPalette = colorRampPalette(brewer.pal(8, "Set2")) # too muted
# getPalette = colorRampPalette(brewer.pal(8, "Set3"))
# getPalette = colorRampPalette(brewer.pal(8, "Accent")) # Ugly
# getPalette = colorRampPalette(brewer.pal(12, "Paired"))

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
run_my_MDS <- function(plotType,
                       dist_methods="None",
                       ordMethod="PCoA",
                       wrapBy="distance") {
  
  print(paste0("Making ", ordMethod, " plots"))
  plist <- NULL
  
  for ( i in dist_methods ) {
    message(paste0("Distance method: ", i))
    myOrd <- tryCatch( 
      {
        if (dist_methods=="None") {
          message(paste0("Attempting ordination using ",ordMethod,"..."))
          suppressWarnings(ordinate(ps.ord, ordMethod))
        } else {
          message(paste0("Attempting ordination using distance ",i," and ordination with ",ordMethod,"..."))
          suppressWarnings(ordinate(ps.ord, ordMethod, distance=iDist))
        }
      },
      error=function(cond)
      {
        message("Error: Returning NULL for this ordination.")
        message(cond)
        # Return NULL with error
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
      p <- plot_ordination(ps.ord, myOrd, type=plotType) # The colors here will be irrelevant
      # This is just used to make a dataframe for plotting
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
  return(df)
}

# Plot the MDS using the dataframe returned above; Could add shape variable to this?
plot_my_MDS <- function(df,
                        myColor=mainFactor,
                        myShape=NULL,
                        ordMethod="",
                        wrapBy="distance") {
  message("Starting ggplot...")
  colnames(df)[2:3] <- c("Axis.1","Axis.2")
  message(paste0("Coloring by ",myColor))
  print(myColor)
  # myColor <- syms(myColor)
  print(myColor)
  myColor2 <- ensym(myColor)
  print(ensym(myColor))
  # print(!!myColor2) # Won't work
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
    print(paste0("Processing region ",regions[i]," from ", region, " variable"))
    print(paste0("Loop number ",i, " of ",length(regions)))
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
sampleRegion <- factor(sapply(sampleList,"[[",regionPosition)) ####THIS DETERMINES WHERE V REGION IS LISTED
# This creates the metadata dataframe
# First column of the table is essentially the barcode
sampledata.df <- data.frame(sampleNumber=sampleNumber)
sampledata.df$V_region <- sampleRegion
# Make a lookup table to add metadata to correct samples
lookup <- data.frame(sampleNumber=unique(sampleNumber))
lookup$PMA <- c("No PMA","No PMA","No PMA","No PMA","No PMA","No PMA",
                "PMA","PMA","PMA","PMA","PMA","PMA")
lookup$Killed <- c("Killed","Killed","Killed",
                   "Unkilled","Unkilled","Unkilled","Unkilled","Unkilled","Unkilled",
                   "Killed","Killed","Killed")
# Add in the metadata to the data frame
sampledata.df <- dplyr::inner_join(lookup, sampledata.df, by = 'sampleNumber')
row.names(sampledata.df) <- row.names(seqtab.nochim)
sampledata.df$PMA_killed <- paste(sampledata.df$PMA,sampledata.df$Killed)

write.table(sampledata.df, file=paste0("metadata.",timeOfScript,".",experimentName,".txt"), quote=F, sep='\t', col.names=NA)

### Perform a sanity check on this!!!! Critical to have correct sample naming.
# row.names(seqtab.nochim) <- row.names(sampledata.df)
# sampleOrder=row.names(sampledata.df)

###################################################
# Open a PDF file to record all plots as we go... #
###################################################

pdf(file=plotFile, width = 8.5, height = 8.5)

####################################
#### Run Phyloseq analysis, 16S ####
####################################

# Create phyloseq object of raw data
ps <- phyloseq(tax_table(taxa_silva),
               sample_data(sampledata.df),
               otu_table(seqtab.nochim,taxa_are_rows=FALSE),
               phy_tree(fitGTR$tree))
# Make main factor an actual factor
# It is wise to customize this for each project to make sure it makes sense
ps@sam_data[[mainFactor]] <- factor(ps@sam_data[[mainFactor]], levels=myLevels)
sampleOrder <- sample_names(ps@sam_data[order(ps@sam_data[[mainFactor]]),])


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
ps=ps0

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

# Agglomerate at genus level, remove NAs
ps_glom_genus=tax_glom(ps2,"Genus" , NArm=TRUE)
ps_glom_genus_ra = transform_sample_counts(ps_glom_genus, function(x){x/sum(x)})

# Number of replicates, for transforming sample counts
reps=3

h1=0.2
ps_glom_height=tip_glom(ps2, h = h1)
ps_glom_species_na_kept <- tax_glom(ps2, taxrank="Species", NArm=FALSE)
ps_glom_species_na_kept_transformed <-transform_sample_counts(ps_glom_species_na_kept, function(x) 100/reps*( x/sum(x)))

ps_glom_genus_na_kept <- tax_glom(ps2, taxrank="Genus", NArm=FALSE)
ps_glom_genus_na_kept_transformed <-transform_sample_counts(ps_glom_genus_na_kept, function(x) 100/reps*( x/sum(x)))

# Make objects encompassing only the top 20 or top 100 OTUs
top20 <- names(sort(taxa_sums(ps2), decreasing=TRUE))[1:20]
ps.top20 <- prune_taxa(top20, ps2)
ps.top20 <- transform_sample_counts(ps.top20, function(x) 100/reps*( x/sum(x)) )

top100 <- names(sort(taxa_sums(ps2), decreasing=TRUE))[1:100]
ps.top100 <- prune_taxa(top100, ps2)
ps.top100 <- transform_sample_counts(ps.top100, function(x) 100/reps*( x/sum(x)) )

#################
# Plot richness #
#################

plot_richness(ps,
              x=mainFactor,
              shape="PMA",
              color="V_region",
              measures=c("Observed", "Shannon", "Simpson", "Chao1"),
              title=paste0("Raw species richness for ",experimentName)) + 
  geom_point(size=5)

plot_richness(ps0,
              x=mainFactor,
              shape="PMA",
              color="V_region",
              measures=c("Observed", "Shannon", "Simpson", "Chao1"),
              title=paste0("Species richness for ",experimentName,", unknown phyla removed")) + 
  geom_point(size=5)

plot_richness(ps2,
              x=mainFactor,
              shape="PMA",
              color="V_region",
              measures=c("Observed", "Shannon", "Simpson", "Chao1")) + 
  geom_point(size=5) +
  ggtitle("Species richness after removal of low-prevalance samples",experimentName)

p1=plot_richness(ps,
                 x=mainFactor,
                 shape="Killed",
                 color="PMA",
                 measures=c("Observed", "Shannon", "Simpson", "Chao1"))

ggplot(p1$data, aes(x=!!ensym(mainFactor), y=value)) + 
  geom_point(size=3) + 
  facet_grid(variable~V_region) +
  ggtitle("Raw species richness, separated by V region",experimentName) +
  p1$mapping +
  theme(axis.text.x = element_text(angle=-90))

##############
# Plot trees #
##############

plot_tree(ps_glom_species_na_kept_transformed,
          method="sampledodge",
          size="Abundance",
          justify="yes please",
          ladderize="left",
          color=mainFactor) +
  scale_size_continuous(range = c(1, 3)) +
  ggtitle("Species-agglomerated tree, relative abundance, NAs kept",experimentName)

plot_tree(ps_glom_genus_na_kept_transformed,
          color=mainFactor,
          label.tips="Phylum",
          method="sampledodge",
          justify="yes please",
          ladderize="left")

plot_tree(ps_glom_genus_na_kept_transformed,
          color=mainFactor,
          method="sampledodge",
          justify="yes please",
          ladderize="left",
          size="abundance")

plot_tree(ps_glom_genus_na_kept_transformed,
          method="sampledodge",
          size="Abundance",
          justify="yes please",
          ladderize="left",
          color=mainFactor) +
  scale_size_continuous(range = c(1, 3)) +
  ggtitle("Genus-agglomerated tree, relative abundance, NAs kept",experimentName)

#####################
# Other tree plots...
# multiPlotTitleTextSize=12

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

###################
# Plot Abundances #
###################

ps2ra = transform_sample_counts(ps2, function(x){x/sum(x)})

# Non-agglomerated
plot_abundance(ps1,
               meta=mainFactor,
               title="Asolute Abundance, PS1 (unfiltered)",
               FacetBy="Class")
plot_abundance(ps2,
               meta=mainFactor,
               title="Absolute Abundance, PS2 (filtered by prevalence threshold)",
               FacetBy="Class")
plot_abundance(ps2ra,
               meta=mainFactor,
               title="Relative Abundance, PS2 (filtered by prevalence threshold)",
               FacetBy="Class")
plot_abundance(ps2ra,
               meta=mainFactor,
               title="Relative Abundance, PS2 (filtered by prevalence threshold), by V region and Class",
               FacetType="grid",
               FacetBy="Class~V_region")

# Agglomerated
# Facet by class
plot_abundance(ps_glom_genus,
               meta=mainFactor,
               title="Absolute abundance by treatment, Agglomerated at genus level",
               FacetBy="Class")
plot_abundance(ps_glom_genus_ra,
               meta=mainFactor,
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

# Facet by Class and V region
plot_abundance(ps_glom_genus,
               meta=mainFactor,
               title="Absolute abundance by V region, separated by Class, Agglomerated at genus level",
               FacetType="grid",
               FacetBy="Class~V_region")
plot_abundance(ps_glom_genus_ra,
               meta=mainFactor,
               title="Relative abundance by V region, separated by Class, Agglomerated at genus level",
               FacetType="grid",
               FacetBy="Class~V_region")


################################################################
#### GO BACK TO ORIGINAL PS OBJECT FOR STATISTICAL ANALYSIS ####
################################################################

##################
### ORDINATION ###
##################

#######################################################
#### Multiple distance metrics combined in single plots

# pdf(file=paste0(experimentName,"ordination.pdf"), width = 8.5, height = 8.5)

# Remove samples with abundance less than N
N=10
ps.ord <- prune_samples(sample_sums(ps)>=N, ps)
ps.ord <- transform_sample_counts(ps.ord, function(x){x/sum(x)})

dist_methods <- unlist(distanceMethodList)
print(dist_methods)
dist_methods = dist_methods[-which(dist_methods=="ANY")]
dist_methods = dist_methods[-which(dist_methods=="manhattan")]
print(dist_methods)

ordMethods.dist_required=c("PCoA","MDS","NMDS")
ordMethods.no_dist=c("DCA","RDA","CCA","DPCoA")
# ordMethods=c("PCoA","MDS","NMDS","DCA","RDA","CCA","DPCoA")
# ordMethods=c("MDS","PCoA")
# ordMethod="PCoA"

for (ordMethod in ordMethods.dist_required) {
  
  print(ordMethod)
  
  # Run distance and ordination
  message("Making required data frames for plotting...")
  message("Samples...")
  myMDS.samples <- run_my_MDS( plotType="samples",
                               dist_methods=dist_methods,
                               ordMethod=ordMethod)
  message("Taxa...")
  myMDS.taxa <- run_my_MDS( plotType="taxa",
                            dist_methods=dist_methods,
                            ordMethod=ordMethod)
  message("Split...")
  myMDS.split <- run_my_MDS( plotType="split",
                             dist_methods=dist_methods,
                             ordMethod=ordMethod)
  # Samples
  message("Plotting by sample")
  plot_my_MDS( df=myMDS.samples,
               myColor=mainFactor,
               ordMethod=ordMethod)
  
  plot_my_MDS( df=myMDS.samples,
               myColor=mainFactor,
               wrapBy="V_region",
               ordMethod=ordMethod_no_dist)
  
  plot_my_MDS( df=myMDS.samples,
               myColor="V_region",
               ordMethod=ordMethod_no_dist)
  
  # Taxa
  message("Plotting by taxa")
  plot_my_MDS( df=myMDS.taxa,
               myColor="Phylum",
               ordMethod=ordMethod_no_dist)
  
  plot_my_MDS( df=myMDS.taxa,
               myColor="Order",
               ordMethod=ordMethod_no_dist)
  
  # Split
  message("Plotting split plots")
  plot_my_MDS( df=myMDS.split,
               myColor=mainFactor,
               ordMethod=ordMethod_no_dist)
  
  plot_my_MDS( df=myMDS.split,
               myColor="Phylum",
               ordMethod=ordMethod_no_dist)
  
  plot_my_MDS( df=myMDS.split,
               myColor="Order",
               ordMethod=ordMethod_no_dist)
  
}

for (ordMethod_no_dist in ordMethods.no_dist) {
  
  print(ordMethod_no_dist)
  
  # Run distance and ordination
  message("Making required data frames for plotting...")
  message("Samples...")
  myMDS.samples.no_dist <- run_my_MDS( plotType="samples",
                                       ordMethod=ordMethod_no_dist)
  message("Taxa...")
  myMDS.taxa.no_dist <- run_my_MDS( plotType="taxa",
                                    ordMethod=ordMethod_no_dist)
  message("Split...")
  myMDS.split.no_dist <- run_my_MDS( plotType="split",
                                     ordMethod=ordMethod_no_dist)
  # Samples
  message("Plotting by sample")
  plot_my_MDS( df=myMDS.samples.no_dist,
               myColor=mainFactor,
               ordMethod=ordMethod_no_dist)
  
  plot_my_MDS( df=myMDS.samples.no_dist,
               myColor=mainFactor,
               wrapBy="V_region",
               ordMethod=ordMethod_no_dist)
  
  plot_my_MDS( df=myMDS.samples.no_dist,
               myColor="V_region",
               ordMethod=ordMethod_no_dist)
  
  # Taxa
  message("Plotting by taxa")
  plot_my_MDS( df=myMDS.taxa.no_dist,
               myColor="Phylum",
               ordMethod=ordMethod_no_dist)
  
  plot_my_MDS( df=myMDS.taxa.no_dist,
               myColor="Order",
               ordMethod=ordMethod_no_dist)
  
  # Split
  message("Plotting split plots")
  plot_my_MDS( df=myMDS.split.no_dist,
               myColor=mainFactor,
               ordMethod=ordMethod_no_dist)
  
  plot_my_MDS( df=myMDS.split.no_dist,
               myColor="Phylum",
               ordMethod=ordMethod_no_dist)
  
  plot_my_MDS( df=myMDS.split.no_dist,
               myColor="Order",
               ordMethod=ordMethod_no_dist)
  
}
# End ordination
################################################################################


################################################################################
## Plot bar graphs
################################################################################

# Plot bar graphs, faceted by V region and factor of interest...
plot_bar(ps_glom_genus_na_kept_transformed,
         fill="Class",
         x="Phylum") +
  geom_bar( stat="identity", position="stack") +
  ggtitle(paste0("Relative abundance of Phyla, by ",mainFactor), experimentName) +
  facet_grid(expr(V_region ~ !!ensym(mainFactor)), scales="free_x")

plot_bar(ps_glom_genus_na_kept_transformed,
         fill="Order",
         x="Phylum") +
  geom_bar( stat="identity", position="stack") +
  ggtitle(paste0("Relative abundance of Phyla, by ",mainFactor), experimentName) +
  facet_grid(expr(V_region ~ !!ensym(mainFactor)), scales="free_x")

plot_bar(ps_glom_genus_na_kept_transformed,
         fill="Family",
         x="Phylum") +
  geom_bar( stat="identity", position="stack") +
  ggtitle(paste0("Relative abundance of Phyla, by ",mainFactor), experimentName) +
  facet_grid(expr(V_region ~ !!ensym(mainFactor)), scales="free_x")

plot_bar(ps_glom_genus_na_kept_transformed,
         fill="Genus",
         x="Phylum") +
  geom_bar( stat="identity", position="stack") +
  ggtitle(paste0("Relative abundance of Phyla, by ",mainFactor), experimentName) +
  facet_grid(expr(V_region ~ !!ensym(mainFactor)), scales="free_x") +
  theme(legend.position = "none")


plot_bar(ps_glom_genus_na_kept_transformed,
         fill="Order",
         x="Class") +
  geom_bar( stat="identity", position="stack") +
  ggtitle(paste0("Relative abundance of Class, by ",mainFactor), experimentName) +
  facet_grid(expr(V_region ~ !!ensym(mainFactor)), scales="free_x")

plot_bar(ps_glom_genus_na_kept_transformed,
         fill="Family",
         x="Class") +
  geom_bar( stat="identity", position="stack") +
  ggtitle(paste0("Relative abundance of Class, by ",mainFactor), experimentName) +
  facet_grid(expr(V_region ~ !!ensym(mainFactor)), scales="free_x")

plot_bar(ps_glom_genus_na_kept_transformed,
         fill="Genus",
         x="Class") +
  geom_bar( stat="identity", position="stack") +
  ggtitle(paste0("Relative abundance of Class, by ",mainFactor), experimentName) +
  facet_grid(expr(V_region ~ !!ensym(mainFactor)), scales="free_x")


plot_bar(ps_glom_genus_na_kept_transformed,
         fill="Family",
         x="Order") +
  geom_bar( stat="identity", position="stack") +
  ggtitle(paste0("Relative abundance of Order, by ",mainFactor), experimentName) +
  facet_grid(expr(V_region ~ !!ensym(mainFactor)), scales="free_x")

plot_bar(ps_glom_genus_na_kept_transformed,
         fill="Genus",
         x="Order") +
  geom_bar( stat="identity", position="stack") +
  ggtitle(paste0("Relative abundance of Order, by ",mainFactor), experimentName) +
  facet_grid(expr(V_region ~ !!ensym(mainFactor)), scales="free_x")

plot_bar(ps_glom_genus_na_kept_transformed,
         fill="Genus",
         x="Family") +
  geom_bar( stat="identity", position="stack") +
  ggtitle(paste0("Relative abundance of Family, by ",mainFactor), experimentName) +
  facet_grid(expr(V_region ~ !!ensym(mainFactor)), scales="free_x")


plot_bar(ps.top20,
         x="Family",
         fill="Genus") +
  labs(x="PMA Treatment") +
  geom_bar( stat="identity", position="stack") +
  ggtitle(paste0("Relative abundance of top 20 OTUs, by ",mainFactor), experimentName) +
  facet_grid(expr(V_region ~ !!ensym(mainFactor)), scales="free_x")

plot_bar(ps.top100,
         x="Family",
         fill="Genus") +
  labs(x="PMA Treatment") +
  geom_bar( stat="identity", position="stack") +
  ggtitle(paste0("Relative abundance of top 100 OTUs, by ",mainFactor), experimentName) +
  facet_grid(expr(V_region ~ !!ensym(mainFactor)), scales="free_x")

######################################################
#### Known pathogenic bacteria #######################
######################################################

# Work off these objects:
# ps_glom_genus
# ps_glom_genus_ra

pathogenGenusList <- scan("M:/Projects/general_purpose_files/pathogenic_genus.txt", what="", sep="\n")
pathogenGenusListNoBacillus <- scan("M:/Projects/general_purpose_files/pathogenic_genus_no_bacillus.txt", what="", sep="\n")
pathogenSpeciesList <- scan("M:/Projects/general_purpose_files/pathogenic_species.txt", what="", sep="\n")

pathogenicGenera <- subset_taxa(ps_glom_genus_ra,
                                Genus %in% pathogenGenusList)

pathogenicGeneraNoBacillus <- subset_taxa(ps_glom_genus_ra,
                                          Genus %in% pathogenGenusListNoBacillus)

pathogenicSpecies <- subset_taxa(ps_glom_genus_ra,
                                 Species %in% pathogenSpeciesList)

plot_bar(pathogenicGenera,
         x="Family",
         fill="Genus") +
  labs(x="PMA Treatment") +
  geom_bar( stat="identity", position="stack") +
  ggtitle(paste0("Relative abundance of genera with known pathogenic strains, by ",mainFactor), experimentName) +
  facet_grid(expr(V_region ~ !!ensym(mainFactor)), scales="free_x")

plot_bar(pathogenicGeneraNoBacillus,
         x="Family",
         fill="Genus") +
  labs(x="PMA Treatment") +
  geom_bar( stat="identity", position="stack") +
  ggtitle(paste0("Relative abundance of genera with known pathogenic strains (except for Bacillus), by ",mainFactor), experimentName) +
  facet_grid(expr(V_region ~ !!ensym(mainFactor)), scales="free_x")

plot_bar(pathogenicSpecies,
         x="Genus",
         fill="Species") +
  labs(x="PMA Treatment") +
  geom_bar( stat="identity", position="stack") +
  ggtitle(paste0("Relative abundance of genera with known pathogenic strains, by ",mainFactor), experimentName) +
  facet_grid(expr(V_region ~ !!ensym(mainFactor)), scales="free_x")


# End pathogenic bacteria analysis
######################################################




# # THIS WILL PRODUCE A LOT OF PLOTS, UNCOMMENT IF NEEDED FOR EXPLORATION
# ## Split phyloseq object into list of lists, one for each V region
ps_v_region_objects <- split.to.regions(ps) # If not already done...

# ## Plot bar graphs for each region, coloring by various taxonomic levels
# regions <- levels(ps@sam_data[[region]])
#
# # Plots for each region, plotted as grouped by main factor
# # Agglomerated at genus, transformed
# for (region in regions) {
#   print(region)
#   plot_bar(ps_v_region_objects$regions_genus_glom_transformed[[1]],
#            fill="Phylum",
#            x=mainFactor) + 
#     ggtitle(paste0("Relative abundance, ",region,", by ",mainFactor),
#             experimentName)
#   plot_bar(ps_v_region_objects$regions_genus_glom_transformed[[1]],
#            fill="Class",
#            x=mainFactor) + 
#     ggtitle(paste0("Relative abundance, ",region,", by ",mainFactor),
#             experimentName)
#   plot_bar(ps_v_region_objects$regions_genus_glom_transformed[[1]],
#            fill="Order",
#            x=mainFactor) + 
#     ggtitle(paste0("Relative abundance, ",region,", by ",mainFactor),
#             experimentName)
#   plot_bar(ps_v_region_objects$regions_genus_glom_transformed[[1]],
#            fill="Family",
#            x=mainFactor) + 
#     ggtitle(paste0("Relative abundance, ",region,", by ",mainFactor),
#             experimentName)
#   plot_bar(ps_v_region_objects$regions_genus_glom_transformed[[1]],
#            fill="Genus",
#            x=mainFactor) + 
#     ggtitle(paste0("Relative abundance, ",region,", by ",mainFactor),
#             experimentName)
# }
# 
# # Plots for each region, plotted as grouped by main factor
# # Agglomerated at genus, NAs removed, transformed
# for (region in regions) {
#   print(region)
#   plot_bar(ps_v_region_objects$regions_genus_glom_NArm_transformed[[1]],
#            fill="Phylum",
#            x=mainFactor) + 
#     ggtitle(paste0("Relative abundance, ",region,", by ",mainFactor),
#             experimentName)
#   plot_bar(ps_v_region_objects$regions_genus_glom_NArm_transformed[[1]],
#            fill="Class",
#            x=mainFactor) + 
#     ggtitle(paste0("Relative abundance, ",region,", by ",mainFactor),
#             experimentName)
#   plot_bar(ps_v_region_objects$regions_genus_glom_NArm_transformed[[1]],
#            fill="Order",
#            x=mainFactor) + 
#     ggtitle(paste0("Relative abundance, ",region,", by ",mainFactor),
#             experimentName)
#   plot_bar(ps_v_region_objects$regions_genus_glom_NArm_transformed[[1]],
#            fill="Family",
#            x=mainFactor) + 
#     ggtitle(paste0("Relative abundance, ",region,", by ",mainFactor),
#             experimentName)
#   plot_bar(ps_v_region_objects$regions_genus_glom_NArm_transformed[[1]],
#            fill="Genus",
#            x=mainFactor) + 
#     ggtitle(paste0("Relative abundance, ",region,", by ",mainFactor),
#             experimentName)
# }
# End barplots
################################################################################


################################################################################
## Heatmaps
################################################################################


# Heatmap

for (i in 1:length(ps_v_region_objects$regions_genus_glom_NArm_transformed)) {
  psHeatMap <- ps_v_region_objects$regions_genus_glom_NArm_transformed[[i]]
  customOrder <- sample_names(psHeatMap@sam_data[order(psHeatMap@sam_data[[mainFactor]]),])
  plot_heatmap(psHeatMap,
               "NMDS",
               "bray",
               mainFactor,
               "Genus",
               sample.order=customOrder) +
    ggtitle(paste0("Variable region ", region[i]), experimentName)
}

# plot_heatmap(ps, "NMDS", "jaccard", mainFactor, "Genus", sample.order=sampleOrder) ## Nice


# End heatmaps
################################################################################


################################################################################
################################################################################
### Differential abundance, adapted for data from dada2
### Need to customize experimental design for each experiment!!!!
################################################################################
################################################################################

# Only consider samples which have greater than N reads
N=100
ps_pruned <- prune_samples(sample_sums(ps) > N, ps)

# Convert phyloseq object to DESeq2 object
dds = phyloseq_to_deseq2(ps_pruned, as.formula(paste0("~",mainFactor)))

# Change factor levels ### CUSTOMIZE!!
# This shouldn't be necessary with new code added at start.
# colData(dds)[[mainFactor]] <- factor(colData(dds)[[mainFactor]], levels=myLevels)

# Geometric mean function
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Calculate geometric means
geoMeans = apply(counts(dds), 1, gm_mean)

# Run DESeq2
dds = estimateSizeFactors(dds, geoMeans = geoMeans)
dds = DESeq(dds, fitType="local")
res = results(dds)
DESeq2::plotMA(dds)
# colData(dds)

# Create list of DESeq results for each contrast (factor)...
resList <- list()

for (i in 2:length(myLevels)) {
  message("Contrast index:")
  print(i)
  message("Contrast factor:")
  print(myLevels[i])
  message("List index:")
  print(i-1)
  message("Computing results...")
  ### This may need customization!
  ### Defaults to use mainFactor as the variable for contrast,
  ### Then the first level of that factor as control
  ### And each successive level as the contrasts
  myContrast <- c(mainFactor,
                  myLevels[1],
                  myLevels[i])
  myResults <- results(dds,
                       contrast=myContrast,
                       cooksCutoff = FALSE)
  message("Done computing results. Adding to list.")
  resList[[i-1]] <- myResults
  message("Done adding to list.")
}

# Filter results table using adjusted p-value of alpha...
alpha = 0.01
sigtabList <- list()

for (i in 1:length(resList)) {
  print(i)
  # Filter by alpha (first removing NAs)
  sigTab <- resList[[i]][!is.na(resList[[i]]$padj) & resList[[i]]$padj < alpha, ]
  # Add taxonomy
  if (nrow(sigTab) == 0) {
    next
  } else {
    sigTab <- cbind(as(sigTab, "data.frame"), as(tax_table(ps)[rownames(sigTab), ], "matrix"))
    # Add factor name
    sigTab <- cbind(sigTab, contrast=myLevels[[i+1]])
    sigtabList[[i]] <- sigTab
  }
}
sigtabList <- sigtabList[!sapply(sigtabList, is.null)] 

##############
# Plot results
theme_set(theme_bw())
#scale_fill_discrete <- function(palname = "Set1", ...) {
#  scale_fill_brewer(palette = palname, ...)
#}

# Sort results tables by LFC
for (sigtable in sigtabList) {
  
  # Phylum order
  x = tapply(sigtable$log2FoldChange, sigtable$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sigtable$Phylum = factor(as.character(sigtable$Phylum), levels=names(x))
  # Genus order
  x = tapply(sigtable$log2FoldChange, sigtable$Genus, function(x) max(x))
  x = sort(x, TRUE)
  sigtable$Genus = factor(as.character(sigtable$Genus), levels=names(x))
  
}

# Write dataframe of all results
significantResults <- do.call(rbind, sigtabList)

# Create list of results tables where NA genera are removed
sigtableListNA.RM <-lapply(sigtabList, subset, !is.na(Genus))

# Create list of combined results
combined_sigtableListNA.RM <- do.call(rbind, sigtableListNA.RM)

## Genus plots, no N/A
ggplot(subset(combined_sigtableListNA.RM,
              !is.na(Genus)),
       aes(x=Genus,
           y=log2FoldChange,
           color=Phylum)) +
  geom_point(size=3) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  facet_grid(~contrast) +
  ggtitle(paste0("Taxa with significantly different relative abundance, grouped by ", mainFactor), experimentName)


# Make phyloseq object from the DIFFERENTIALLY ABUNDANT genera
# This could be limited to include only certain V regions...
# INCLUDE NAs
genera_affected_NA_included <- subset_taxa(ps, Genus %in% levels(subset(significantResults)$Genus))
genera_affected_NA_included_genus_glom <- tax_glom(genera_affected_NA_included, NArm=F)
# REMOVE NAs
genera_affected_NArm <- subset_taxa(ps, Genus %in% levels(subset(significantResults, !is.na(Genus))$Genus))
genera_affected_NArm_genus_glom <- tax_glom(genera_affected)

# Heatmap of signifcantly different genera
plot_heatmap(genera_affected_NA_included,
             "RDA",
             "unifrac",
             mainFactor,
             "Genus",
             sample.order=sampleOrder) +
  ggtitle("Differentially abundant genera, NAs included",experimentName)

plot_heatmap(genera_affected_NA_included_genus_glom,
             "RDA",
             "unifrac",
             mainFactor,
             "Genus",
             sample.order=sampleOrder) +
  ggtitle("Differentially abundant genera, NAs included, agglomerated",experimentName)

plot_heatmap(genera_affected_NArm,
             "RDA",
             "unifrac",
             mainFactor,
             "Genus",
             sample.order=sampleOrder) +
  ggtitle("Differentially abundant genera, NArm",experimentName)

plot_heatmap(genera_affected_NArm_genus_glom,
             "RDA",
             "unifrac",
             mainFactor,
             "Genus",
             sample.order=sampleOrder) +
  ggtitle("Differentially abundant genera, NArm, agglomerated",experimentName)



#######################################
### Write results table from DESeq2
#######################################
write.table(significantResults, file=paste0(experimentName,".",timeOfScript,".DESeq_output.txt"), quote=F, sep='\t', col.names=NA)
#######################################


# ORDERING TAXA
significantResults$taxID <- row.names(significantResults)

significantReults.orderedGByenera <- as.character(arrange(significantResults, log2FoldChange)$taxID[!is.na(arrange(significantResults, log2FoldChange)$taxID)])

significantReults.orderedGeneraByPval <- as.character(arrange(significantResults, padj)$taxID[!is.na(arrange(significantResults, log2FoldChange)$taxID)])


### ADD THIS IN PLACE FOR THE HEATMAP TO ELIMINATE UNASSIGNED GENERA
# subset_taxa(dataSamplesPrunedRarefied,  !is.na(Genus)

plot_heatmap(subset_taxa(ps,  !is.na(Genus)),
             "none",
             "none",
             mainFactor, 
             taxa.label="Genus",
             sample.order=sampleOrder,
             taxa.order=significantReults.orderedGByenera,
             max.label=5000) +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10)) +
  ggtitle(paste0("Statistically significant results, ordered by genus name, ", mainFactor, ", ", experimentName))


plot_heatmap(subset_taxa(ps,  !is.na(Genus)),
             "none",
             "none",
             mainFactor, 
             taxa.label="Genus",
             sample.order=sampleOrder,
             taxa.order=significantReults.orderedGeneraByPval,
             max.label=5000) +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10)) +
  ggtitle(paste0("Statistically significant results, ordered by adjusted P value, ", mainFactor, ", ", experimentName))


#######################################
# Close the PDF file to write to disk #
#######################################
dev.off()
#######################################
save.image(paste0("./",experimentName,".",timeOfScript,".RData"))

