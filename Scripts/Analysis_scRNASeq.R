#######################



# neurotransmitter receptors
total_ephys <- readRDS("total_ionchannels_gapjunctions_gproteincoupledreceptors.rds")
neurotransmitter_classes <- c("Gap junction proteins",
                              "5-hydroxytryptamine receptors, ionotropic ",
                              "5-hydroxytryptamine receptors, G protein-coupled",
                              "Glutamate ionotropic receptor delta type subunits",
                              "Glutamate ionotropic receptor AMPA type subunits",
                              "Glutamate ionotropic receptor NMDA type subunits",
                              "Glutamate ionotropic receptor kainate type subunits" ,
                              "Glutamate metabotropic receptors",
                              "Glycine receptors",
                              "Cholinergic receptors nicotinic subunits",
                              "Cholinergic receptors muscarinic",
                              "Adrenoceptors",
                              "Dopamine receptors",
                              "Gamma-aminobutyric acid type A receptor subunits",
                              "Gamma-aminobutyric acid type B receptor subunits",
                              "Histamine receptors",
                              "Opioid receptors")
neurotransmitter_genes <- total_ephys[total_ephys$Group.name %in% neurotransmitter_classes,]
neurotransmitter_genes <- neurotransmitter_genes[!neurotransmitter_genes$Group.name %in% ("Gap junction proteins"),]
neurotransmitter_genes$Group.name
neurotransmitter_genes$label <- str_remove_all(neurotransmitter_genes$Group.name, " " )
neurotransmitter_genes$label <- str_replace_all(neurotransmitter_genes$label, "-", "_")
neurotransmitter_genes$label <- str_replace_all(neurotransmitter_genes$label, ",", "" )
neurotransmitter_genes$label <- paste0("NT_", neurotransmitter_genes$label            )

yu <- readRDS("yu.rds")
neftel <- readRDS("neftel.rds")
pdx <- readRDS("pdx.rds")



long_short <- readRDS("long_short.rds")
neurotransmitter_genes_col <- merge(neurotransmitter_genes, long_short, by = "label" )
neurotransmitter_genes_col$Ensembl.gene.ID
neurotransmitter_genes_col %>%
  select(ntgroup, Approved.symbol, Approved.name, Ensembl.gene.ID) %>%
  write.csv("neurotransmitter_SITable.csv")



#colors
garofanocolors = c("#E32419", "#45A048", "#005EDB", "#69E3DC")
rColors <- c("#F2756D", "#19BCC1")
suvacolors <- c("#af5c7d", "#e2e055", "#6c86f0","#7da83b" )

############ Functions #############

getSuvaClass <- function(seuratObject) {
  cellStateAssignments <- read_xlsx("1-s2.0-S0092867419306877-mmc2.xlsx")
  colnames(cellStateAssignments) <- cellStateAssignments[4,]
  cellStateAssignments <- cellStateAssignments[-c(1:4),]
  colnames(cellStateAssignments)
  seuratObject <- AddModuleScore(seuratObject, features = list(cellStateAssignments$AC), name = "AC_score", search = T)
  print("added ac..")
  seuratObject <- AddModuleScore(seuratObject, features = list(cellStateAssignments$MES1), name = "MES1_score", search = T)
  print("added mes1..")
  seuratObject <- AddModuleScore(seuratObject, features = list(cellStateAssignments$MES2), name = "MES2_score", search = T)
  print("added mes2..")
  seuratObject <- AddModuleScore(seuratObject, features = list(cellStateAssignments$OPC), name = "OPC_score", search = T)
  print("added opc..")
  seuratObject <- AddModuleScore(seuratObject, features = list(cellStateAssignments$NPC1), name = "NPC1_score", search = T)

  seuratObject <- AddModuleScore(seuratObject, features = list(cellStateAssignments$NPC2), name = "NPC2_score", search = T)
  print("added npcs")
  metintern <- seuratObject@meta.data
  metintern <- metintern[,(ncol(metintern)-5):(ncol(metintern))]
  seuratObject@meta.data$CellState <- colnames(metintern)[(apply(metintern, 1, which.max))]
  seuratObject@meta.data$CellState <- factor(seuratObject@meta.data$CellState)

  print("done.")
  seuratObject@meta.data$CellState_short <- "default"
  seuratObject@meta.data$CellState_short[which(seuratObject@meta.data$CellState == "AC_score1")] <- "AClike"
  seuratObject@meta.data$CellState_short[which(seuratObject@meta.data$CellState == "MES1_score1")] <- "MESlike"
  seuratObject@meta.data$CellState_short[which(seuratObject@meta.data$CellState == "MES2_score1")] <- "MESlike"
  seuratObject@meta.data$CellState_short[which(seuratObject@meta.data$CellState == "OPC_score1")] <- "OPClike"
  seuratObject@meta.data$CellState_short[which(seuratObject@meta.data$CellState == "NPC1_score1")] <- "NPClike"
  seuratObject@meta.data$CellState_short[which(seuratObject@meta.data$CellState == "NPC2_score1")] <- "NPClike"
  seuratObject@meta.data$CellState_short <- factor(seuratObject@meta.data$CellState_short)
  return(seuratObject)
}
addInvasivityScore <- function(seu) {
  correlationTable <- readRDS("final_invasivity_score_correlationsTable.rds")
  seu <- addScoreCorrelationAsSplit(seu, correlationTable, "InvasivityScore")
  return(seu)
}


addScoreCorrelationAsSplit <- function(dataset, markerTable, title) {

  upregulatedMarkers <- markerTable %>%
    filter(correlation == "Correlated")

  downregulatedMarkers <- markerTable %>%
    filter(correlation == "Anticorrelated")
  upregulatedMarkers <- upregulatedMarkers$gene
  downregulatedMarkers <- downregulatedMarkers$gene

  name_up <- paste(title, "_Upregulated", sep = "")
  name_down <- paste(title, "_Downregulated", sep = "")
  name_up1 <- paste(name_up, "1", sep = "")
  name_down1 <- paste(name_down, "1", sep = "")
  dataset <-   AddModuleScore(dataset, features = list(upregulatedMarkers), name = name_up, search = T)
  print("Added Score for Upregulated Genes.")
  dataset <-   AddModuleScore(dataset, features = list(downregulatedMarkers), name = name_down, search = T)
  print("Added Score for Downregulated Genes.")
  getElement(dataset@meta.data, name_up1)
  upreg_score <- getElement(dataset@meta.data, name_up1)
  downreg_score <- getElement(dataset@meta.data, name_down1)

  combined_score <- data.frame(title = downreg_score - upreg_score)
  rownames(combined_score) <- rownames(dataset@meta.data)
  colnames(combined_score) <- title
  print("Calculated Combined Score. Adding to Seurat Object...")

  dataset <- AddMetaData(dataset,combined_score)

  return(dataset)
}

############ Functions #############

pdx <- getSuvaClass(pdx)
neftel <- addInvasivityScore(neftel)
yu <- addInvasivityScore(yu)


########## Prepare dataset lists for analyses
analysis_set <- list(
  list(name = "PDX Dataset",
    dataset = pdx,
    groupby = "SampleGroup",
    colors = rColors),

  list(name = "PDX Dataset",
       dataset = pdx,
       groupby = "Pathway_based",
       colors = garofanocolors),

  list(name = "PDX Dataset",
       dataset = pdx,
       groupby = "CellState_short",
       colors = suvacolors),



  list(name = "Yu Dataset",
    dataset = yu,
    groupby = "SampleGroup",
    colors = rColors),
  list(name = "Yu Dataset",
    dataset = yu,
    groupby = "Neftel_et_al.",
    colors = suvacolors
    ),
  list(name = "Yu Dataset",
    dataset = yu,
    groupby = "Pathway_based",
    colors = garofanocolors
    ),

  list(name = "Neftel Dataset",
    dataset = neftel,
    groupby = "CellState_short",
    colors = suvacolors
    )


)
names(analysis_set) <- c("PDX_conuncon", "PDX_garo", "PDX_suva", "Yu_rimcore", "Yu_suva", "Yu_garo", "Neftel")

########################################



##### Neurotransmitter receptors per group#############

splitted <- split(neurotransmitter_genes_col, neurotransmitter_genes_col$ntgroup)
addNeurotransmitter_Groups <- function(seu) {

  for (x in 1:length(splitted)) {
    genes <- unique(splitted[[x]]$Approved.symbol)
    print(name)
    print(genes)
    name <- names(splitted)[x]
    seu <- AddModuleScore(seu, features = list(genes), name = name)

  }
  return(seu)
}

groups <- colnames(datasets[[1]]$dataset@meta.data)[str_starts(colnames(datasets[[1]]$dataset@meta.data), "NT_")]
analysis_set$PDX_suva$dataset <- getSuvaClass(analysis_set$PDX_suva$dataset)
plots <- pblapply(analysis_set[[7]], gene = groups, FUN = function(set, gene) {
  # violin plot, set values from dataset

  p <- DotPlot(set$dataset, group.by = set$groupby, features = gene, cols = suvacolors) +
    ggtitle(paste(set$name, "-", "Neurotransmitter Receptors")) +
    theme(text = element_text(size = 6),
          axis.title =  element_text(size = 6),
          axis.text =  element_text(size = 6),
    ) +
    ylab("") +
    coord_flip()

})
ggarrange(plotlist = plots, nrow = 1)
ggsave(paste0("plots/Fig_", "neurotransmitters1_groups", ".png"), height = 11, width = 8.5)


##Synaptogenesis - Synapse assembly Score  vs Invasivity Score ############

data <- read.csv("synapse_assembly_goterm.txt", sep = "\t", header = F)


datasets <- list(analysis_set$Neftel,
                 analysis_set$PDX_garo,
                 analysis_set$Yu_rimcore
)

datasets <- pblapply(datasets, function(set) {
  set$dataset <- AddModuleScore(set$dataset, features = list(unique(toupper(data$V1))), name = "synapse_assembly")
  return(set)
})


getCorrelation_plot <- function(set) {
  dats <- FetchData(set$dataset, vars = c("InvasivityScore", "synapse_assembly1", set$groupby ))
  p1 <- ggplot(dats, aes(x = InvasivityScore, y = synapse_assembly1)) +
    geom_point( shape = 21, ) +
    stat_cor() +
    geom_smooth()  +
    scale_color_manual(values = set$colors) +
    theme_bw() +
    ggtitle(paste(set$name, "- Correlation Synaptogenic vs Invasivity Score"))  +
    ylab("Synaptogenic Score")
  return(p1)
}
correlation_plots <- pblapply(datasets, getCorrelation_plot)
ggarrange(plotlist = correlation_plots[c(1,3,2)])
ggsave(paste0("plots/Fig_", "synapseAssembly_vs_InvasivityCorrelation", ".png"), height = 8, width = 10)


sample_ranking <- datas %>%
  group_by(Sample) %>%
  summarize(mean = mean(synapse_assembly1),
            mean2 = mean(scaleInv + scaleSynp)
            ) %>%
  arrange(mean)
datas$Sample <- factor(datas$Sample, levels = sample_ranking$Sample)
datas <- datas %>%
  filter(!is.na(Sample))

p1 <- ggplot(datas, aes(x = Sample, y= synapse_assembly1, fill = Sample )) +
  geom_violin() +
  theme_bw() +
  RotatedAxis()+
  theme(legend.position = "none") +
  stat_summary() +
  ylab("Synaptogenesis Score")
p2 <- ggplot(datas, aes(x = Sample, y= InvasivityScore, fill = Sample )) +
  geom_violin() +
  theme_bw() +
  RotatedAxis()+
  theme(legend.position = "none") +
  stat_summary() +
  ylab("Invasivity Score")

p1 / p2
ggsave("plots/Synaptogenesis_vsinvasivity_perPatient2.png")



########## Cholinergic Aceytcholin receptors only######

plots <- pblapply(analysis_set[-c(3:4)], function(set) {
  p1 <- DotPlot(set$dataset, features = c("CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5"),
                group.by = set$groupby,
                  cols = suvacolors
                ) +
    ggtitle(paste(set$name, "- Acetylcholin")) +
    coord_flip() +
    xlab("") +
    ylab("")
  p1
  return(p1)
})
ggarrange(plotlist = plots, ncol = 2, nrow = 3)
ggsave("plots/CHRMs.png")









