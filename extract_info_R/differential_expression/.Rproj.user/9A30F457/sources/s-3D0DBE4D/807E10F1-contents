library(Seurat)
library(limma)
library(EnhancedVolcano)
library(DESeq2)

DATA_DIR <- "/Users/ryanyutian/Desktop/NucSeq_data"        

dat_mal <- readRDS(file.path(DATA_DIR,'20210715_MalignantCells_seurat.rds'))

Idents(object = dat_mal) <- dat_mal@meta.data[, 1]
levels(dat_mal)

A_R_vs_RR.de.markers <- FindMarkers(dat_mal, ident.1 = "A_R_GBM607_180329", ident.2 = "A_RR_GBM809_180413")

EnhancedVolcano(A_R_vs_RR.de.markers,
                lab = rownames(A_R_vs_RR.de.markers),
                x = 'avg_log2FC',
                y = 'p_val',
                xlab = bquote(~Log[2]~ 'fold change'),
                pointSize = 4.0,
                labSize = 4.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)

A_R_vs_RR.de.markers.minpct <- FindMarkers(dat_mal, ident.1 = "A_R_GBM607_180329", ident.2 = "A_RR_GBM809_180413", min.pct = 0.5)

EnhancedVolcano(A_R_vs_RR.de.markers.minpct,
                lab = rownames(A_R_vs_RR.de.markers.minpct),
                x = 'avg_log2FC',
                y = 'p_val',
                xlab = bquote(~Log[2]~ 'fold change'),
                pointSize = 4.0,
                labSize = 4.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)


A_R_vs_RR.de.markers.DESeq <-  FindMarkers(dat_mal, ident.1 = "A_R_GBM607_180329", ident.2 = "A_RR_GBM809_180413", test.use = "DESeq2")

EnhancedVolcano(A_R_vs_RR.de.markers.DESeq,
                lab = rownames(A_R_vs_RR.de.markers.DESeq),
                x = 'avg_log2FC',
                y = 'p_val',
                xlab = bquote(~Log[2]~ 'fold change'),
                pointSize = 4.0,
                labSize = 4.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)
