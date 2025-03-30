## gene presence/absence using heatmap in R

#required packages
install.packages("tidyverse")
install.packages("pheatmap")
library(tidyverse)
library(pheatmap)
install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)

# Read the CSV and fix column names
gene_data <- read.csv("gene_matrix.csv", row.names = 1, check.names = FALSE)
head(gene_data)


#annotating strain type
strain_types <- data.frame(
  Condition = c("Control", "Control", "Ethanol-Adapted", "Ethanol-Adapted", "Ethanol-Adapted")
)
rownames(strain_types) <- colnames(gene_matrix)
ann_colors <- list(
  Condition = c("Control" = "#b1cbee", "Ethanol-Adapted" = "#5a77f0")
)

# Set up your matrix

# GENE PRESENCE/ABSENCE HEATMAP IN R
mat <- as.matrix(gene_data)

# Define Sunrise palette (colors for presence/absence)
col_fun <- c("0" = "white", "1" = "#dd4c24")  # white = Absent, red = Present

# Define the condition colors using your Sunrise palette
strain_conditions <- c("Control", "Control", "Ethanol-Adapted", "Ethanol-Adapted", "Ethanol-Adapted")
names(strain_conditions) <- colnames(mat)

ha_bottom <- HeatmapAnnotation(
  Condition = strain_conditions,
  col = list(Condition = c("Control" = "#b1cbee", "Ethanol-Adapted" = "#5a77f0")),
  annotation_name_side = "left"
)

# Convert matrix to numeric (must be numeric for ComplexHeatmap)
mat_numeric <- apply(mat, 2, as.numeric)
rownames(mat_numeric) <- rownames(mat)

# Create the heatmap
Heatmap(mat_numeric,
        name = "Presence",
        col = col_fun,
        rect_gp = gpar(col = "black", lwd = 1),  # black border around cells
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        bottom_annotation = ha_bottom,
        heatmap_legend_param = list(
          at = c(0, 1),
          labels = c("Absent", "Present"),
          border = "black",  # legend box outline
          legend_height = unit(3, "cm")
        )
        # no title
)

# dN/dS VISUALIZATION
# Started with running 1 gene through this process before doing more, as suggested by ChatGPT
mamba install bioconda::hyphy -c conda-forge
sed 's/|.*//' accC_codon_alignment.fasta > accC_codon_cleaned.fasta
sed -E 's/\"//g' accC_newick | sed -E 's/([,(])[^,():|"]+\|([^,():|"]+)/\1\2/g' > accC_tree_fixed.nwk
# first try -> led to error that only 3 tree tips were being read from .nwk
hyphy FEL --alignment accC_codon_cleaned.fasta --tree accC_tree_fixed.nwk --branches All
cat accC_tree_fixed.nwk | tr '(),:' '\n'| grep -v '^$' | grep -v '[0-9]' | sort
grep "^>" accC_codon_cleaned.fasta | sed 's/^>//' | sort
# second try, this was supposed to force hyphy to use all tree tips no matter their similarity
# got same error, result is .json files that are empty
hyphy FEL --alignment accC_codon_cleaned.fasta --tree accC_tree_fixed.nwk --branches All --full-model Yes
# it is a problem that with accC there are 2 duplicate sequences

## Ligand bar graph and statistical testing
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Load your data
df <- read.csv("FabG_Docking_Data.csv")

# Summarize data (Mean and SD)
df_summary <- df %>%
  group_by(Protein_Version, Ligand) %>%
  summarize(
    Mean_Affinity = ifelse(all(is.na(Binding_Affinity)), NA, mean(Binding_Affinity, na.rm = TRUE)),
    SD_Affinity = ifelse(all(is.na(Binding_Affinity)), NA, sd(Binding_Affinity, na.rm = TRUE))
  )

# Perform two-sample t-tests (Glucose and Ribitol separately)
ttest_glucose <- t.test(Binding_Affinity ~ Protein_Version, 
                        data = subset(df, Ligand == "Glucose"),
                        var.equal = FALSE)

ttest_ribitol <- t.test(Binding_Affinity ~ Protein_Version, 
                        data = subset(df, Ligand == "Ribitol"),
                        var.equal = FALSE)

# Create dataframe for statistical results (dynamically use p-values)
stat_results <- data.frame(
  Ligand = c("Glucose", "Ribitol"),
  p_value = c(ttest_glucose$p.value, ttest_ribitol$p.value)
)

# Add significance labels clearly based on p-values
stat_results <- stat_results %>%
  mutate(Significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

stat_results$Ligand <- factor(stat_results$Ligand, levels = c("Glucose", "Ribitol"))

# Merge summary with statistical results to determine annotation positions
df_signif <- df_summary %>%
  group_by(Ligand) %>%
  summarize(max_affinity = max(Mean_Affinity + SD_Affinity, na.rm = TRUE)) %>%
  left_join(stat_results, by = "Ligand") %>%
  mutate(y_position = ifelse(max_affinity + 2 > 9, 9, max_affinity + 2)) # keeps labels clearly within graph

# Prepare data for "N/D" labels explicitly (retain Protein_Version)
df_nd <- df_summary %>%
  filter(is.na(Mean_Affinity)) %>%
  select(Ligand, Protein_Version) %>%
  mutate(y_nd = 9)

# Generate the bar graph clearly with significance annotations (no LOD line)
ggplot(df_summary, aes(x = Ligand, y = Mean_Affinity, fill = Protein_Version)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.8, na.rm = TRUE) +
  geom_errorbar(aes(ymin = Mean_Affinity - SD_Affinity, ymax = Mean_Affinity + SD_Affinity),
                position = position_dodge(0.9), width = 0.25, linewidth = 0.7, na.rm = TRUE) +
  theme_minimal() +
  labs(x = "Ligand",
       y = "Mean Binding Affinity (ΔG, kcal/mol)",
       fill = "FabG Version") +
  scale_fill_manual(values = c("#6BAED6", "#FB6A4A")) +
  coord_cartesian(ylim = c(-10, 10)) +
  geom_text(data = df_nd,
            aes(x = Ligand, y = y_nd, label = "N/D", group = Protein_Version),
            position = position_dodge(width = 0.9),
            color = "black",
            fontface = "bold",
            size = 4) +
  geom_text(data = df_signif,
            aes(x = Ligand, y = y_position, label = Significance),
            color = "black",
            size = 5,
            fontface = "bold")