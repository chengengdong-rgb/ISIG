# Load required packages
library(randomForest)
library(pROC)
library(ROCit)
library(survminer)
library(survival)
library(plyr)  
library(caret)
library(tidyr)
library(GSVA)
library(stringr)
library(ggplot2)
library(gridExtra)

# ROC curve visualization (best model)
set.seed(3252)
best_test_idx <- c(3,5,6,10,13)
train_best <- clustotal01[-best_test_idx, ]
test_best <- clustotal01[best_test_idx, ]

# Train final model
rf_final <- randomForest(as.factor(V14) ~ ., data = train_best, importance = TRUE)

# Predict and plot ROC
pred_best <- predict(rf_final, newdata = test_best[,1:13], type = "prob")
roc_best <- roc(test_best$V14, pred_best[,1])
roc_ext_best <- roc(clustotal_ext$V14, predict(rf_final, clustotal_ext[,1:13], type = "prob")[,1])

# Save ROC plot
pdf("AUC_curve.pdf", width = 6, height = 6)
plot(roc_best, col = "#2F4F4F", main = "Validation set", print.auc = FALSE, lwd = 3)
lines(roc_ext_best, col = "brown")
legend("bottomright", 
       legend = c(paste0('Validation AUC:', round(auc(roc_best),2)),
                  paste0('Testing AUC:', round(auc(roc_ext_best),2))),
       col = c('#2F4F4F','brown'), lwd = 4, cex = 0.9)
dev.off()

# ROCit style plots
pdf("ROCit_Validation.pdf", width = 6, height = 6)
rocit_val <- rocit(score = pred_best[,1], class = test_best$V14, negref = '1')
plot(rocit_val, legend = F, YIndex = F)
title('Validation set')
text(x = 0.6, y = 0.4, labels = paste0("AUC: ", round(rocit_val$AUC, 2)))
dev.off()

pdf("ROCit_Testing.pdf", width = 6, height = 6)
rocit_test <- rocit(score = predict(rf_final, clustotal_ext[,1:13], type = "prob")[,1],
                    class = clustotal_ext$V14, negref = '1')
plot(rocit_test, legend = F, YIndex = F)
title('Testing set')
text(x = 0.6, y = 0.4, labels = paste0("AUC: ", round(rocit_test$AUC, 2)))
dev.off()

# Cross-validation optimization
set.seed(448)
cv_test_idx <- c(5,7,9,10,11)
cv_train <- clustotal01[-cv_test_idx, ]

# 5-fold cross validation
ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE, 
                     summaryFunction = twoClassSummary)
rf_cv <- train(x = as.matrix(cv_train[,1:13]), 
               y = make.names(as.factor(cv_train$V14)),
               method = "rf", trControl = ctrl,
               preProcess = c("center", "scale"),
               tuneLength = 10, metric = "ROC")

# Apply CV model to TCGA data (same visualization logic as above)
# ... (reuse plot_survival and plot_stage_dist functions)

# SSGSEA score calculation for bulk RNA-seq
marker_genes <- read.csv("marker_genes.csv")
marker_list <- split(marker_genes$markergene, marker_genes$celltype)

# Read bulk RNA-seq data
bulk_data <- read.table("bulk_RNAseq.txt", header = TRUE, row.names = 1)

# Calculate SSGSEA scores
ssgsea_scores <- gsva(as.matrix(bulk_data), marker_list, 
                      method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)
ssgsea_df <- as.data.frame(t(ssgsea_scores))
ssgsea_df$sample <- rownames(ssgsea_df)

# Merge with clinical data
clinical_bulk <- read_excel("bulk_clinical.xlsx", skip = 2)
merged_bulk <- merge(ssgsea_df, clinical_bulk, by.x = "sample", by.y = "title")

# Calculate total score and grouping
merged_bulk$total_score <- rowMeans(merged_bulk[,2:14])
merged_bulk$score_group <- ifelse(merged_bulk$total_score > median(merged_bulk$total_score), 
                                  "high", "low")
merged_bulk$response <- ifelse(merged_bulk$BR %in% c("PR", "CR"), "R", "NR")


# Response rate visualization
pdf("Response_rate.pdf", width = 6, height = 5)
ggplot(merged_bulk, aes(x = score_group, fill = response)) +
  geom_bar(width = 0.5, position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#AEC7E8", "#FF9E4A")) +
  labs(x = "Score Group", y = "Percentage", fill = "Response") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
dev.off()

# Survival analysis for OS and PFS
pdf("Bulk_Survival.pdf", width = 10, height = 5)
par(mfrow = c(1,2))

# OS survival
os_surv <- Surv(merged_bulk$OS, merged_bulk$dead)
os_km <- survfit(os_surv ~ score_group, data = merged_bulk)
os_plot <- ggsurvplot(os_km, pval = TRUE, conf.int = F,
                      title = "OS", palette = c("#E7B800", "#2E9FDF"))
print(os_plot)

# PFS survival
pfs_surv <- Surv(merged_bulk$PFS, merged_bulk$dead)
pfs_km <- survfit(pfs_surv ~ score_group, data = merged_bulk)
pfs_plot <- ggsurvplot(pfs_km, pval = TRUE, conf.int = F,
                      title = "PFS", palette = c("#E7B800", "#2E9FDF"))
print(pfs_plot)
dev.off()

# Fisher's exact test
fisher_test <- fisher.test(xtabs(~score_group + response, data = merged_bulk))
print(fisher_test)