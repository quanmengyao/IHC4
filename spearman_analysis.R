library(readxl)
library(dplyr)

data <- read_excel("spearman.xlsx")
colnames(data)
ihc4 <- data$IHC4_score
features <- data[,3:ncol(data)]
results <- data.frame()
##########Spearman相关性分析############
for (i in 1:ncol(features)) {
  
  test <- cor.test(
    ihc4,
    features[[i]],
    method = "spearman"
  )
  
  results <- rbind(
    results,
    data.frame(
      Feature = colnames(features)[i],
      Spearman_r = test$estimate,
      p_value = test$p.value
    )
  )
}
results$FDR <- p.adjust(results$p_value, method = "fdr")
results <- results %>%
  arrange(FDR)
significant_features <- results %>%
  filter(FDR < 0.05)
write.csv(results,
          "Spearman_IHC4_radiomics_results.csv",
          row.names = FALSE)

write.csv(significant_features,
          "Spearman_significant_features.csv",
          row.names = FALSE)

##########特征相关性绘图############
library(stringr)

plot_data <- significant_features %>%
  arrange(Spearman_r) %>%
  mutate(
    Feature = str_replace_all(Feature, "\\.", "_"),
    Feature = str_extract(Feature, "[^_]+$"),
    Feature = factor(Feature, levels = Feature)
  )

feature_colors <- c(
  "#7BA6C9", "#C98C8C", "#A8C686", "#C7B37A",
  "#9B8CC2", "#6FB1B6", "#D1A56F", "#8FA39A",
  "#C27D9A", "#A7A7A7", "#6E9ECF", "#B7C36B"
)

feature_colors <- feature_colors[1:nrow(plot_data)]

p <- ggplot(plot_data, aes(x = Feature, y = Spearman_r, fill = Feature)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = feature_colors) +
  theme_classic() +
  labs(
    x = "",
    y = "correlation coefficient",
    title = ""
  ) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 12)
  )

p

ggsave(
  "Spearman_correlation_IHC4.pdf",
  plot = p,
  width = 7,
  height = 4,
  dpi = 600,
  device = cairo_pdf
)

##########相关性散点图############
library(ggplot2)
library(dplyr)
library(stringr)

# 确保输出文件夹存在
out_dir <- "Spearman_scatterplots"
if (!dir.exists(out_dir)) dir.create(out_dir)

# 使用显著特征表
plot_features <- significant_features

for (i in 1:nrow(plot_features)) {
  
  feature_name <- plot_features$Feature[i]
  p_value <- plot_features$p_value[i]
  
  # 简化特征名：保留最后一段 + 提取wavelet滤波名
  feature_clean <- str_replace_all(feature_name, "\\.", "_")
  feature_short <- str_extract(feature_clean, "[^_]+$")
  
  if (str_detect(feature_clean, "^wavelet_")) {
    wavelet_part <- str_extract(feature_clean, "(?<=wavelet_)[A-Z]{3}")
    y_label <- paste0(feature_short, " (", wavelet_part, ")")
  } else if (str_detect(feature_clean, "^original_")) {
    y_label <- paste0(feature_short, " (Original)")
  } else {
    y_label <- feature_short
  }
  
  # P值格式
  p_label <- ifelse(
    p_value < 0.001,
    "P < 0.001",
    paste0("P = ", format(round(p_value, 3), nsmall = 3))
  )
  
  p <- ggplot(data, aes(x = IHC4_score, y = .data[[feature_name]])) +
    geom_point(
      color = "#2C7FB8",
      size = 1.0,
      alpha = 0.75
    ) +
    geom_smooth(
      method = "lm",
      color = "gray",
      se = FALSE,
      linewidth = 0.9
    ) +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = p_label,
      hjust = -0.1, vjust = 1.2,
      size = 4
    ) +
    theme_classic() +
    labs(
      x = "IHC4 score",
      y = y_label,
      title = ""
    ) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12)
    )
  
  # 文件名安全化
  file_name <- feature_clean %>%
    str_replace_all("[^A-Za-z0-9_]+", "_") %>%
    paste0(".pdf")
  
  ggsave(
    filename = file.path(out_dir, file_name),
    plot = p,
    width = 4,
    height = 3.5,
    dpi = 600
  )
}

