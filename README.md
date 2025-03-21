# ASCVD_project

Gut bacteria derived odd chain fatty acid modulates cholesterol homeostasis and alleviates atherosclerosis


# SOURCE CODE

Refer to the Methods section for details.






# ANALYSIS PIPELINE

# I.	Megatnomic analysis

1.	Quality control of metagenomic reads

sickle se -t sanger -q 25 -f $entry -o $entry.t.fq

2.	Remove human reads

bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path=/lustre/home/mwcai/data/remove_human_reads_from_MG_MT qtrim=rl trimq=10 untrim -Xmx23g in=$entry outu=$entry.clean.fq outm=$entry.human.fq

3.	Fastq to fasta format

seqtk seq -a $entry > $entry.fa

4.	Taxonomic and abundance of the high-quality reads

metawrap kraken2 -t 32 -o kraken2_SRR11461968 SRR11461968.fastq.t.fq.clean.fq.fa

5.	Assessment of PA synthesis gene abundance in the two published cohorts

coverm contig --single $entry -r ref_gene.fa --min-read-percent-identity 95 --min-read-aligned-percent 50 -o $entry.coverm.ref_gene -m rpkm





# II.	Identification of PA synthesis genes in Unified Human Gastrointestinal Genome (UHGG) collection

blastp -db gene.fasta.faa -query NBT_total.faa -outfmt 6 -out NBT_total.faa_1E3_50.blastp -evalue 1E-3 -num_threads 64 -qcov_hsp_perc 50





# III.	Metatranscriptomic analysis of mouse liver

1.	derep analysis

#Load necessary library

library(DESeq2)

#Read CSV file

data <- read.csv("gene_expression.csv", row.names = 1)

#Create a data frame to store group information

group_info <- data.frame(
  sample = colnames(data),
  condition = factor(rep(c("BU", "CHD", "PBS"), each = 4))
)

#Check for NA values in the data and remove those rows

data <- na.omit(data)

#Convert FPKM values to pseudocounts

#Ensure all values are within the integer range

data_counts <- round(data * 1e2)
data_counts[data_counts < 0] <- 0  # Convert negative values to 0

#Create DESeq dataset

dds <- DESeqDataSetFromMatrix(countData = data_counts, colData = group_info, design = ~ condition)

#Run DESeq analysis

dds <- DESeq(dds)

#Define comparison function

perform_comparison <- function(dds, group1, group2) {
  res <- results(dds, contrast = c("condition", group1, group2))
  res$gene <- rownames(res)
  res <- as.data.frame(res)
  res <- res[, c("gene", "log2FoldChange", "pvalue")]
  res$log2FoldChange <- 2^res$log2FoldChange  # Convert to FC
  res <- res[order(res$pvalue), ]
  res$cluster <- paste(group1, "vs", group2)
  colnames(res) <- c("Name", "FC", "Pvalue", "cluster")
  return(res)
}

#Perform comparisons between different groups

res_BU_vs_CHD <- perform_comparison(dds, "BU", "CHD")
res_BU_vs_PBS <- perform_comparison(dds, "BU", "PBS")
res_PBS_vs_CHD <- perform_comparison(dds, "PBS", "CHD")

#Merge results

result <- rbind(res_BU_vs_CHD, res_BU_vs_PBS, res_PBS_vs_CHD)

#Output results to a CSV file

output_file <- "differential_expression_results.csv"
write.csv(result, file = output_file, row.names = FALSE)

#Print output file path

cat("Results saved to:", output_file, "\n")

2.	multi_volcano_plot

#Load packages

library(tidyverse)
library(ggrepel)   # For labeling

#Read CSV file

data <- read.csv("differential_expression_results.csv")

#Define a custom function for later use

mutiVolcano = function(df,         # Data for plotting
                       P = 0.05,   # P-value cutoff
                       FC = 1.5,   # Fold change cutoff
                       GroupName = c("Sig", "Not Sig"),      # Group labels
                       pointColor = c("#CC3333", "#0099CC"), # Colors for points
                       barFill = "#efefef",  # Color for bars
                       pointSize = 0.9,      # Size of points
                       labeltype = "1",      # Option for labeling differentially expressed genes. Options are "1" and "2"
                       labelNum = 5,         # Number of points to label when labeltype is "1"
                       labelName = NULL,     # Names of points to label when labeltype is "2"
                       tileLabel = "Label",  # Option for labeling comparison pairs. Options are "Label" and "Num". "Label" shows group names, "Num" shows numbers to avoid clutter
                       tileColor = NULL      # Colors for comparison pairs
){
 
  #Group data based on P-value cutoff

  dfSig = df %>% 
    mutate(log2FC = log2(FC)) %>%
    filter(FC > {{FC}} | FC < (1/{{FC}})) %>%
    mutate(Group = ifelse(Pvalue < 0.05, GroupName[[1]], GroupName[[2]])) %>%
    mutate(Group = factor(Group, levels = GroupName)) %>%
    mutate(Cluster = factor(Cluster, levels = unique(Cluster)))   # Cluster order follows the order in the file
  
  #Prepare data for bar plot

  dfBar = dfSig %>%
    group_by(Cluster) %>%
    summarise(min = min(log2FC, na.rm = T),
              max = max(log2FC, na.rm = T)
    )
  
  #Prepare data for scatter plot

  dfJitter = dfSig %>%
    mutate(jitter = jitter(as.numeric(Cluster), factor = 2))
  
  #Prepare data for labeling differentially expressed genes

  if(labeltype == "1"){
    # Label type 1: Top N points with smallest P-values
    dfLabel = dfJitter %>%
      group_by(Cluster) %>%
      slice_min(Pvalue, n = labelNum, with_ties = F) %>%
      ungroup()
  }else if(labeltype == "2"){
    # Label type 2: Specify points to label
    dfLabel = dfJitter %>%
      filter(Name %in% labelName)
  }else{
    dfLabel = dfJitter %>% slice()
  }
  
  #Create the plot

  p = ggplot() +
    # Draw bar plot
    geom_col(data = dfBar, aes(x = Cluster, y = max), fill = barFill) +
    geom_col(data = dfBar, aes(x = Cluster, y = min), fill = barFill) +
    # Draw scatter plot
    geom_point(data = dfJitter,
               aes(x = jitter, y = log2FC, color = Group),
               size = pointSize,
               show.legend = NA
    ) +
    #Draw middle label tiles
  ggplot2::geom_tile(data = dfSig,
                       ggplot2::aes(x = Cluster, y = 0, fill = Cluster), 
                       color = "black",
                       height = log2(FC) * 1.5,
                       show.legend = NA
    ) + 
    #Label differentially expressed genes
  ggrepel::geom_text_repel(
      data = dfLabel,
      aes(x = jitter,                   # geom_text_repel labeling function
          y = log2FC,          
          label = Name),        
      min.segment.length = 0.1,
      max.overlaps = 10000,             # Maximum overlap, adjust to avoid overlapping labels
      size = 3,                         # Font size
      box.padding = unit(0.5, 'lines'), # Margin around labels
      point.padding = unit(0.1, 'lines'), 
      segment.color = 'black',          # Color of label lines
      show.legend = F)
  
  if(tileLabel == "Label"){
    p =
      p +
      geom_text(data = dfSig, aes(x = Cluster, y = 0, label = Cluster)) +
      ggplot2::scale_fill_manual(values = tileColor,
                                 guide = NULL # Hide legend
      )
  }else if(tileLabel == "Num"){
    # If comparison pair names are too long, use numeric labels
    p =
      p +
      geom_text(data = dfSig, aes(x = Cluster, y = 0, label = as.numeric(Cluster)), show.legend = NA) +
      ggplot2::scale_fill_manual(values = tileColor,
                                 labels = c(paste0(1:length(unique(dfSig$Cluster)), ": ", unique(dfSig$Cluster))))
  }
  
  #Customize theme

  p = p + ggplot2::scale_color_manual(values = pointColor) +
    theme_classic() +
    ggplot2::scale_y_continuous(n.breaks = 5) + 
    ggplot2::theme(
      legend.position = "right", 
      legend.title = ggplot2::element_blank(), 
      legend.background = ggplot2::element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank()
    ) + 
    ggplot2::xlab("Clusters") + ggplot2::ylab("log2FC") + 
    guides(color = guide_legend(override.aes = list(size = 3)))
  
  return(p)
}

#Read data

df = read.csv("differential_expression_results.csv") %>%  # Read local CSV file
  as_tibble() %>%
  set_names(c("Name", "FC", "Pvalue", "Cluster"))

#Call the function to plot

mutiVolcano(
  df = df,    # Data for plotting
  P = 0.05,   # P-value cutoff
  FC = 1.5,   # Fold change cutoff
  GroupName = c("Sig", "Not Sig"),      # Group labels
  pointColor = c("#CC3333", "#0099CC"), # Colors for points
  barFill = "#efefef",   # Color for bars
  pointSize = 0.9,       # Size of points
  labeltype = "1",       # Option for labeling differentially expressed genes
  labelNum = 5,          # Number of points to label when labeltype is "1"
  labelName = c("Slc1a2", "Glul"),   # Names of points to label when labeltype is "2"
  tileLabel = "Label",           # Option for labeling comparison pairs
  tileColor = RColorBrewer::brewer.pal(length(unique(df$Cluster)), "Set3")   # Colors for comparison pairs
)

#Save the plot

ggsave("mutiVolcano.pdf", width = 8, height = 6)
