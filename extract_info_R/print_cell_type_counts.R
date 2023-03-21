library(Seurat)
library(limma)
library(EnhancedVolcano)
library(DESeq2)

DATA_DIR <- "/Users/ryanyutian/Desktop/NucSeq_data"       

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

for (sample_id in comparable_sample_id) {
  
  print(sample_id)
  temp_sample.cell_use <- colnames(dat)[which(dat@meta.data[, 1] == sample_id)]
  temp_sample <- subset(dat, cells = temp_sample.cell_use)
  
  Idents(object = temp_sample) <- temp_sample@meta.data$Final_Annotation
  print(table(Idents(object = temp_sample)))

}