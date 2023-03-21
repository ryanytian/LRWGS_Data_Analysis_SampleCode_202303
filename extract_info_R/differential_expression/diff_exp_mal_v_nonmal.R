library(Seurat)
library(limma)
library(EnhancedVolcano)
library(DESeq2)

DATA_DIR <- "/Users/ryanyutian/Desktop/NucSeq_data"       
SAVE_DIR <- "/Users/ryanyutian/Desktop/TRI_Brain_diff_exp/min_pct_50"
SAVE_DIR_DESEQ <- "/Users/ryanyutian/Desktop/TRI_Brain_diff_exp/deseq2"

dat <- readRDS(file.path(DATA_DIR,'GBM_Cohort_NoDoublets_Annotated_July152021_Seurat.rds'))

comparable_sample_id <- c("A_R_GBM607_180329", "A_RR_GBM809_180413", "B_P_GBM593.1_180413", 
                                  "B_P_GBM593.2_181205", "B_R_GBM898_181205", "C_P_GBM577.1_180718",
                                  "C_P_GBM577.2_181205", "C_R_GBM625_181205", "F_P_GBM620_200206", 
                                  "F_R_GBM691_181026", "G_P_GBM454_181102", "G_R_GBM833_181102", 
                                  "H_P_GBM460_181205", "H_R_GBM492_181026", "I_P_GBM440_190712", 
                                  "I_R_GBM532_181026", "J_P_GBM401_181129", "J_R_GBM498_190703", 
                                  "J_RR_GBM551_181026", "K_P_GBM529_190712", "K_R_GBM832_200218",
                                  "L_P_GBM618_190703", "L_R_SMTB152_190613",
                                  "M_P_GBM672_190620", "M_R_GBM828_190703", "N_P_BT2013110_190703",
                                  "N_R_GBM745_200218", "X3_R_SMTB265_200115", "O_P_SMTB665_190703",
                                  "O_R_GBM1070_190430", "X20_R_SMTB302_190923", "X23_R_SMTB814_190923")

# comparable_sample_id_archive <- c("L_P_GBM618_190703", "K_R_GBM832_200218")

### Seurat built-in differential expression

for (sample_id in comparable_sample_id) {
  
  temp_sample.cell_use <- colnames(dat)[which(dat@meta.data[, 1] == sample_id)]
  temp_sample <- subset(dat, cells = temp_sample.cell_use)
  
  Idents(object = temp_sample) <- temp_sample@meta.data$Final_Annotation
  print(levels(temp_sample))
  
  temp_sample.de.markers.minpct <- FindMarkers(temp_sample, ident.1 = "Malignant", ident.2 = NULL, min.pct = 0.5)
         
  write.csv(x=temp_sample.de.markers.minpct, 
            file=paste(SAVE_DIR, "/", substr(sample_id, 1, nchar(sample_id)-7), "_diff_exp.csv", sep=""),
            row.names = TRUE)
  
  print(EnhancedVolcano(temp_sample.de.markers.minpct,
                  title = sample_id,
                  lab = rownames(temp_sample.de.markers.minpct),
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  xlab = bquote(~Log[2]~ 'fold change'),
                  pointSize = 4.0,
                  labSize = 4.0,
                  colAlpha = 1,
                  legendPosition = 'right',
                  legendLabSize = 12,
                  legendIconSize = 4.0,
                  drawConnectors = TRUE,
                  widthConnectors = 0.75,
                  maxoverlapsConnectors = 20))
  
}

plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="/Users/ryanyutian/Desktop/TRI_Brain_diff_exp/min_pct_50/png")



### DESeq2

for (sample_id in comparable_sample_id) {
  
  temp_sample.cell_use <- colnames(dat)[which(dat@meta.data[, 1] == sample_id)]
  temp_sample <- subset(dat, cells = temp_sample.cell_use)
  
  Idents(object = temp_sample) <- temp_sample@meta.data$Final_Annotation
  print(levels(temp_sample))
  
  temp_sample@assays$RNA@counts <- as.matrix(temp_sample@assays$RNA@counts)+1
  temp_sample.de.markers.deseq <- FindMarkers(temp_sample, ident.1 = "Malignant", ident.2 = NULL,  test.use = "DESeq2")
  
  write.csv(x=temp_sample.de.markers.deseq, 
            file=paste(SAVE_DIR_DESEQ, "/", substr(sample_id, 1, nchar(sample_id)-7), "_diff_exp.csv", sep=""),
            row.names = TRUE)
  
  print(EnhancedVolcano(temp_sample.de.markers.deseq,
                        title = sample_id,
                        lab = rownames(temp_sample.de.markers.deseq),
                        x = 'avg_log2FC',
                        y = 'p_val_adj',
                        xlab = bquote(~Log[2]~ 'fold change'),
                        pointSize = 4.0,
                        labSize = 4.0,
                        colAlpha = 1,
                        legendPosition = 'right',
                        legendLabSize = 12,
                        legendIconSize = 4.0,
                        drawConnectors = TRUE,
                        widthConnectors = 0.75,
                        maxoverlapsConnectors = 20))
  
}

plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="/Users/ryanyutian/Desktop/TRI_Brain_diff_exp/deseq2/png")








### Archive code for one example

A_R_GBM607_cell.use <- colnames(dat)[which(dat@meta.data[, 1] == "A_R_GBM607_180329")]
A_R_GBM607 <- subset(dat, cells = A_R_GBM607_cell.use)

Idents(object = A_R_GBM607) <- A_R_GBM607@meta.data$Final_Annotation
levels(A_R_GBM607)

A_R_GBM607.de.markers <- FindMarkers(A_R_GBM607, ident.1 = "Malignant", ident.2 = NULL)

EnhancedVolcano(A_R_GBM607.de.markers,
                title = 'N061011 versus N61311',
                lab = rownames(A_R_GBM607.de.markers),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'fold change'),
                pointSize = 4.0,
                labSize = 4.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                maxoverlapsConnectors = 20)

A_R_GBM607.de.markers.minpct <- FindMarkers(A_R_GBM607, ident.1 = "Malignant", ident.2 = NULL, min.pct = 0.5)

EnhancedVolcano(A_R_GBM607.de.markers.minpct,
                lab = rownames(A_R_GBM607.de.markers.minpct),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'fold change'),
                pointSize = 4.0,
                labSize = 4.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                maxoverlapsConnectors = 20)

A_R_GBM607.malmarker = subset(A_R_GBM607.de.markers.minpct, -log10(p_val_adj)>1 & avg_log2FC>1)
A_R_GBM607.nmmarker = subset(A_R_GBM607.de.markers.minpct, -log10(p_val_adj)>1 & avg_log2FC<(-1))


A_R_GBM607.de.markers.DESeq <- FindMarkers(A_R_GBM607, ident.1 = "Malignant", ident.2 = NULL, test.use = "DESeq2")

EnhancedVolcano(A_R_GBM607.de.markers.DESeq,
                lab = rownames(A_R_GBM607.de.markers.DESeq),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'fold change'),
                pointSize = 4.0,
                labSize = 4.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)