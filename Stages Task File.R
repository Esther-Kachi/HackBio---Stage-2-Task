
#Stage 2 Task

#importing dataframes

#HBR vs UHR
fileOnline <- 'https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_top_deg_normalized_counts.csv'
#read in the file
read.csv(fileOnline, header = T)
#assign
HBR.vs.UHR <- read.csv(fileOnline, header = T)

#Deg_CHR22
fileOnline <- 'https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_deg_chr22_with_significance.csv'
#read in the file
read.csv(fileOnline, header = T)
#assign
Deg_CHR22 <- read.csv(fileOnline, header = T)

#Breast_CancerDS
fileOnline <- 'https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv'
#read in the file
read.csv(fileOnline, header = T)
#assign
Breast_CancerDS <- read.csv(fileOnline, header = T)


#Stage 2 Task 

#Part 1 -Gene Expression
#(a) Heatmap
pheatmap::pheatmap(mat = HBR.vs.UHR[ , 2:7],
                   border_color = 'black',
                   legend = T,
                   labels_row = HBR.vs.UHR$X,
                   fontsize_row = 5,
                   cluster_rows = T,
                   cluster_cols = T,
                   color = colorRampPalette(c('white', 'lightblue', 'navyblue'))(100)
                 
)

#(b) Volcano plot
#
ggplot(Deg_CHR22,
       aes(x = log2FoldChange, y = X.log10PAdj, color = significance)) +
  geom_point(size = 1.8, alpha = 1) +
  scale_color_manual(
    values = c("up" = "green",
      "down" = "orange",
      "ns" = "grey")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  labs(
    x = "log2 Fold Change",
    y = "X.log10PAdj",
    color = "Significance"
  ) +
  theme_minimal(base_size = 12)


#Part 2 - Breast Cancer Data Exploration
#(c) Scatter plot - radius vs texture
#colour grade
col_vec <- c("B" = "blue", # Benign
             "M" = "red"  # Malignant
             )

plot(x = Breast_CancerDS$radius_mean,
     y = Breast_CancerDS$texture_mean,
     xlab = 'radius_mean',
     ylab = 'texture_mean',
     xlim = c(7, 28),
     ylim = c(9, 40),
     las = 1,
     main = 'Radius vs Texture',
     col = col_vec[Breast_CancerDS$diagnosis],
     pch = 19,
     cex = 0.5)


#(d) Correlation Heatmap
key_features <- Breast_CancerDS[, c("radius_mean", "texture_mean", 
                                    "perimeter_mean","area_mean", 
                                    "smoothness_mean", "compactness_mean")]

# Compute correlation matrix
corr_matrix <- cor(key_features, use = "complete.obs")

# Plot with pheatmap
pheatmap::pheatmap( mat = corr_matrix,
                    display_numbers = TRUE, # annotate correlation values
                    number_format = "%.1f", # round to 1 decimal place
                    fontsize_number = 11,
                    border_color = 'black',
                    cluster_rows = FALSE,
                    cluster_cols = FALSE,
                    color = colorRampPalette(c("white", "lightblue", "navyblue"))(100),
                    main = "Correlation Heatmap of Key Features")



#(e) Scatter plot - Smoothness vs Compactness
plot(x = Breast_CancerDS$smoothness_mean,
     y = Breast_CancerDS$compactness_mean,
     xlab = 'smoothness_mean',
     ylab = 'compactness_mean',
     xlim = c(0.05, 0.17),
     ylim = c(0.01, 0.37),
     las = 1,
     main = 'Smoothness vs Compactness',
     col = as.factor(Breast_CancerDS$diagnosis),
     pch = 19,
     cex = 0.7)
grid(col = "lightgray", lty = "dotted", lwd = 1)


#(f) Density plot
# Split the data
area_Mal <- Breast_CancerDS$area_mean[Breast_CancerDS$diagnosis == "M"]
area_Be <- Breast_CancerDS$area_mean[Breast_CancerDS$diagnosis == "B"]

# Compute densities
Dens_M <- density(area_Mal)
Dens_B <- density(area_Be)

# Plot the first density
plot(Dens_M,
     col = "red",
     lwd = 2,
     xlab = "Area Mean",
     ylab = "Density",
     main = "Kernel Density Estimates of Area Mean",
     ylim = c(0, 0.003)
)

# Add the second density
lines(Dens_B, col = "blue", lwd = 2)

# Add legend
legend("topright",
       legend = c("Malignant (M)", "Benign (B)"),
       col = c("red", "blue"),
       lwd = 1,
       title = "Diagnosis"
)


#part 3
#Task 0
install.packages('readxl')
install.packages('ggplot2')
install.packages('pheatmap') 
install.packages('igraph')

#transparent colour
#use this function to create transparent colors
transparent_color <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

hb_pal <- c("#4e79a7", 
            "#8cd17d", 
            "#e15759", 
            "#fabfd2", 
            "#a0cbe8", 
            "#59a14f", 
            "#b07aa1", 
            "#ff9d9a", 
            "#f28e2b", 
            "#f1ce63",
            "#79706e",
            "#d4a6c8",
            "#e9e9e9",
            "#ffbe7d",
            "#bab0ac",
            "#9d7660",
            "#d37295",
            "#86bcb6",
            "#362a39",
            "#cd9942")

#test the color pallete
plot(1:length(hb_pal), 1:length(hb_pal), col = hb_pal, pch = 19)

#Importing XL
library(readxl)
read_excel('hb_stage_2.xlsx')

#Task 1 - Cell-type ratio distribution
sheet_a <- read_excel('hb_stage_2.xlsx', sheet = 'a')

boxplot(new_ratio ~ cell_type,
        data = sheet_a,
        notch = F,
        ylim = c(0, 0.5),
        las = 2,
        main = 'Cell-type ratio distribution',
        ylab = 'Ratio',
        xlab = 'Cell type',
        outpch = 10,      # solid dots
        outcex = 0.6,     # larger
        col = hb_pal)

#Conceptual check - Taking the log2 of half‑life stabilizes variance and makes differences between cell types easier to compare, especially when the raw values span several orders of magnitude. Boxplots summarize the distribution within each cell type, allowing us to visually assess whether certain immune compartments tend to have faster or slower mRNA decay. This helps identify cell types with distinct RNA stability profiles.



#Task 2 - Half-life vs alpha-life scatter
sheet_b <- read_excel('hb_stage_2.xlsx', sheet = 'b')

# (You can adjust these values if your paper uses specific cutoffs)
thresh_alpha <- median(sheet_b$alpha, na.rm = TRUE)
thresh_hl <- median(sheet_b$half_life, na.rm = TRUE)

sheet_b <- sheet_b %>%
  mutate(
    log_alpha = log2(alpha),
    log_hl = log2(half_life),
    regime = case_when(
      alpha > thresh_alpha & half_life > thresh_hl ~ "High/Stable",
      alpha > thresh_alpha & half_life <= thresh_hl ~ "High/Unstable",
      alpha <= thresh_alpha & half_life > thresh_hl ~ "Low/Stable",
      alpha <= thresh_alpha & half_life <= thresh_hl ~ "Low/Unstable"
    )
  )

#Create the Scatter Plot
Plot2b <- ggplot(sheet_b, 
                 aes(x = log_alpha, y = log_hl, color = regime)) +
  geom_point(alpha = 0.4, size = 1) +
  # Add Cutoff lines
  geom_vline(xintercept = log2(thresh_alpha), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = log2(thresh_hl), linetype = "dashed", color = "grey40") +
  # Label Exemplar Genes (Camp, Ccr2)
  geom_label_repel(data = filter(sheet_b, cell %in% c("Camp", "Ccr2")),
                   aes(label = cell),
                   color = "black", fontface = "bold",
                   vjust = -1) +
  # Custom Colors and Styling
  scale_color_manual(values = c("#e15759", "#4e79a7", "#f28e2b", "#76b7b2")) +
  theme_classic() +
  labs(
    title = "Half-life vs Alpha (log2 scale)",
    x = "log2(alpha)",
    y = "log2(half_life)",
    color = "Kinetic Regime"
  ) +
  theme_minimal(base_size = 14)

print(Plot2b)

#Conceptual Check - Plotting log2(half‑life) against log2(alpha) reveals how transcriptional activation (alpha) relates to mRNA stability (half‑life). The log transformation makes the relationship more linear and interpretable. Vertical and horizontal cutoffs divide the plot into four kinetic regimes, highlighting genes that are rapidly induced, stable, transient, or weakly expressed. Labeling exemplar genes like Camp and Ccr2 illustrates contrasting regulatory behaviors and anchors the interpretation in known biology.



#Task 3 - Heatmap across cell types and time
sheet_c <- read_excel('hb_stage_2.xlsx', sheet = 'c')
#inspect dataset
str(sheet_c)
#extract numeric matrix
mat <- as.matrix(sheet_c[ , -1])   # remove genes
rownames(mat) <- sheet_c$genes

#Build column annotations
#seperate the celltypes and time
celltypes <- str_extract(colnames(mat), "^[A-Za-z]+")
celltypes <- str_replace(celltypes, "n$", "")
print(celltypes)

times <- str_extract(colnames(mat), "[0-9]+h")
print(times)

annotation_col <- data.frame(CellType = celltypes,
                             Time = times
)
rownames(annotation_col) <- colnames(mat)

#Plot Heatmap
pheatmap::pheatmap( mat,
                    annotation_col = annotation_col,
                    fontsize = 11,
                    cluster_rows = TRUE,
                    cluster_cols = FALSE,
                    scale = "row",
                    color = colorRampPalette(c("white", "lightblue", "navyblue"))(100),
                    show_rownames = FALSE,
                    show_colnames = FALSE
)

#Conceptualcheck - A heatmap allows us to visualize temporal expression patterns across multiple immune cell types simultaneously. Clustering only the rows (genes) groups together genes with similar dynamic behavior, revealing coordinated transcriptional programs. Columns (timepoints) are not clustered because time has a natural biological order that should be preserved. Column annotations for cell type and time help contextualize the expression patterns and highlight compartment‑specific responses.


#Task 4 - Pathway enrichment heatmap
sheet_d <- read_excel('hb_stage_2.xlsx', sheet = 'd_1')

#using pathway names as rownames
rownames <- sheet_d$pathway
Mat_d <- as.matrix(sheet_d[ , -1])   # remove pathway column

#centering color scale
range_val <- max(abs(Mat_d))

pheatmap::pheatmap(Mat_d,
  fontsize = 10,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("red", "white", "blue"))(200),
  breaks = seq(-range_val, range_val, length.out = 200),
  main = "Pathway Enrichment Across Timepoints",
  annotation_legend = TRUE,
  legend = TRUE,
  labels_row = rownames(Mat_d),
  annotation_names_row = TRUE
)

#Conceptual Check - Pathway‑level heatmaps summarize complex gene expression changes into interpretable biological processes. Using pathway names as rownames preserves biological meaning. Clustering is disabled because pathway order is often curated or grouped by functional categories, and reordering them would obscure interpretation. A diverging color scale centered at zero clearly distinguishes activated pathways (positive enrichment) from suppressed pathways (negative enrichment), making temporal trends easy to interpret.


#Task 5 - Bubble plot of kinetic regimes
sheet_e <- read_excel('hb_stage_2.xlsx', sheet = 'e')

#Bubble plot
p2e <- ggplot(sheet_e, aes(x = half_life, y = alpha, size = count, color = stage)) +
  geom_point(alpha = 0.7) +
  
  # Apply project palette for stages
  scale_color_manual(values = c("6h" = "#e15759", "72h" = "#8cd17d")) +
  # Customize size range for readability
  scale_size_continuous(range = c(2, 10)) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "Functional Kinetic Regimes",
    x = "Half Life",
    y = "Alpha Life",
    size = "Gene Count",
    color = "Stage"
  )

# Display plot
print(p2e)

#Conceptual check - This plot compares transcriptional activation (alpha) with mRNA stability (half‑life) while encoding an additional variable (e.g., timepoint or expression level) through bubble size and color. Using a bubble plot allows three dimensions of information to be visualized simultaneously. This helps identify genes that are both strongly induced and long‑lived, or transient responders with short half‑lives. The color and size cues highlight how gene behavior changes across time, revealing dynamic regulatory patterns that would be missed in a simple scatterplot.



#Task 6 - Stacked Proportions
sheet_f <- read_excel('hb_stage_2.xlsx', sheet = 'f')

#Subset the data
#Filter for start and end timepoint
sheet_f <- sheet_f %>% 
  filter(stage %in% c("s00h", "s72h"))

#Create Stacked Barplot
p2f <- ggplot(sheet_f, aes(x = stage, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  # Fixed y-axis limit as per requirements
  scale_y_continuous(limits = c(0, 0.3), expand = c(0, 0)) +
  # Apply colors from the project palette (B and Plasma)
  scale_fill_manual(values = c("B" = "#4e79a7", "Plasma" = "#b07aa1")) +
  theme(
    plot.title = element_text(face = "bold", size = 12)
  ) +
  labs(
    title = "Figure 2f: Stacked Proportion",
    x = "Timepoint",
    y = "Proportion of Total Cells",
    fill = "Cell Type"
  )

#Display plot
print(p2f)

#Conceptual Check - A stacked barplot shows how two related cell populations jointly contribute to a fixed total proportion at each timepoint. Because the y‑axis represents a proportion, stacking emphasizes composition rather than absolute counts. This makes it easy to see how B cells and Plasma cells share the “proportion budget” over time and how their relative contributions shift during the immune response. Side‑by‑side bars would compare magnitudes, but stacked bars better highlight relative balance within each timepoint.


sheet_g <- read_excel('hb_stage_2.xlsx', sheet = 'g')

#Convert to Adjacency Matrix
mat_g <- as.matrix(sheet_g[,-1])
rownames(mat_g) <- sheet_g[[1]]
colnames(mat_g) <- sheet_g[[1]]

#Build Directed Graph
net <- graph_from_adjacency_matrix(mat_g, mode = "directed", weighted = TRUE, diag = FALSE)

#Remove Zero-weight Edges
net <- delete.edges(net, which(E(net)$weight <= 0))

#Define Layout and Styling
set.seed(42) 
layout_fr <- layout_with_fr(net)

#Map edge width and arrow size to the weight
E(net)$width <- E(net)$weight * 5
E(net)$arrow.size <- E(net)$weight * 1.5

#Map node colors from the palette
node_colors <- hb_pal[1:vcount(net)]

#Plot the Network
plot(net, 
     layout = layout_fr,
     vertex.color = node_colors,
     vertex.label.color = "black",
     vertex.label.cex = 0.8,
     vertex.label.dist = 2,
     vertex.frame.color = "white",
     edge.color = transparent_color("grey", 30), # Now this function will work!
     edge.curved = 0.2,
     main = "Figure 2g: Cell-Cell Interaction Network")


#Conceptual check - The interaction network represents communication between immune cell types using ligand–receptor signaling. The graph is directed because signaling has a biological direction: one cell type sends a ligand, and another receives it through a receptor. Edge weights encode the strength of communication, often based on ligand–receptor expression or interaction scores. Heavier edges indicate stronger or more influential signaling routes. This visualization highlights which cell types dominate communication and how information flows through the immune system during activation.



#Task 8. Final assembly (mandatory)
# 1. Load All Necessary Libraries
library(readxl)
library(ggplot2)
library(pheatmap)
library(igraph)
library(dplyr)
library(gridExtra)
library(grid)

# 2. Global Setup
file_path <- "C:\\Users\\user\\Downloads\\hb_stage_2.xlsx"

# --- RECREATE PANELS AS OBJECTS ---

# Panel 2a: Boxplot
p2a <- ggplot(sheet_a, aes(x = reorder(cell_type, -new_ratio, FUN = median), y = new_ratio, fill = cell_type)) +
  geom_boxplot(outlier.size = 0.5) + scale_fill_manual(values = hb_pal) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "none") +
  labs(title = "a) Cell-type Ratio distribution", x = "", y = "Ratio")

# Panel 2b: Scatter
p2b <- ggplot(sheet_b, aes(x = log2(alpha), y = log2(half_life))) +
  geom_point(alpha = 0.2, color = "grey50") + theme_classic() +
  geom_text(data = filter(sheet_b, cell %in% c("Camp", "Ccr2")), aes(label=cell), color="red", fontface="bold") +
  labs(title = "b) Kinetic Landscape (log2 scale)", x = "log2 (Alpha)", y = "log2 (Half-life)")

# Panel 2c: Gene Heatmap (Convert to Grob)
mat_c <- as.matrix(sheet_c[,-1]); rownames(mat_c) <- sheet_c[[1]]
p2c <- pheatmap::pheatmap(mat_c, cluster_cols = FALSE, scale = "row", silent = TRUE,
                color = colorRampPalette(c("white", "lightblue", "navyblue"))(100),
                main = "c) Gene Expression")$gtable

# Panel 2d: Pathway Heatmap (Convert to Grob)
mat_d <- as.matrix(sheet_d[,-1]); rownames(mat_d) <- sheet_d[[1]]
p2d <- pheatmap::pheatmap(mat_d, cluster_cols = FALSE, cluster_rows = FALSE, silent = TRUE,
                color = colorRampPalette(c("red", "white", "blue"))(200),
                main = "d) Pathway Enrichment Activity")$gtable

# Panel 2e: Bubble Plot
p2e <- ggplot(sheet_e, aes(x = half_life, y = alpha, size = count, color = stage)) +
  geom_point(alpha = 0.6) + scale_x_log10() + scale_y_log10() + theme_bw() +
  theme(legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size=7)) +
  labs(title = "e) Functional Kinetic Regimes", x = "Half-life", y = "Alpha")

# Panel 2f: Stacked Bar
sheet_f <- read_excel(file_path, sheet = "f") %>% filter(stage %in% c("s00h", "s72h"))
p2f <- ggplot(sheet_f, aes(x = stage, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.5) + theme_classic() +
  labs(title = "f) Lineage Shift", x = "", y = "Prop.")

# --- FINAL ASSEMBLY ---

# Create the layout grid
# We arrange them in a 3-row layout
final_plot <- grid.arrange(
  p2a, p2b, 
  p2c, p2d, 
  p2e, p2f,
  ncol = 2,
  widths = c(1.2, 0.8),
  top = textGrob("Figure 2: Global Immune Response Kinetics", gp = gpar(fontsize=16, fontface="bold"))
)

# Export as high-resolution PNG
ggsave("Figure2_Publication_Ready.png", final_plot, width = 12, height = 16, dpi = 300)

