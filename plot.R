#!/usr/bin/env Rscript

# ======================== #
#      MR 分析后续处理      #
# ======================== #

# 0. 设置工作路径
setwd("/Users/wilian/git/heihaheihaha")

# 1. 加载所需的 R 包
required_packages <- c("dplyr", "ggplot2", "readr", "stringr", "ggpubr", "cowplot")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# 2. 读取所有 MR 结果文件
# 假设所有结果文件命名为 "*_MRresults.csv" 并存储在 "data" 目录下
result_files <- list.files(path = "data", pattern = "_MRresults\\.csv$", full.names = TRUE)

# 检查是否找到任何结果文件
if (length(result_files) == 0) {
  stop("未在 'data' 目录中找到任何符合 '*_MRresults.csv' 模式的文件。请检查文件路径和命名。")
}

# 3. 读取并合并所有结果
all_results <- lapply(result_files, function(file) {
  message("读取文件: ", file)
  read_csv(file, show_col_types = FALSE)
})

# 合并所有结果为一个数据框
merged_results <- bind_rows(all_results)

# 检查必要的列是否存在
required_columns <- c("outcome", "exposure", "method", "nsnp", "b", "se", "pval", 
                      "id.outcome", "id.exposure", "or", "or_lci95", "or_uci95")
missing_columns <- setdiff(required_columns, colnames(merged_results))
if (length(missing_columns) > 0) {
  stop("以下必要的列在合并的结果中缺失: ", paste(missing_columns, collapse = ", "))
}

# 4. 多重检验校正
# 这里使用 Benjamini-Hochberg 方法调整 p 值以控制 FDR
merged_results <- merged_results %>%
  mutate(pval_adj = p.adjust(pval, method = "BH"))

# 5. 筛选显著结果
# 设定阈值，例如调整后的 p 值 < 0.05
significance_threshold <- 0.05
significant_results <- merged_results %>%
  filter(pval_adj < significance_threshold)

# 打印显著结果的数量
message("共有 ", nrow(significant_results), " 个显著结果（调整后的 p 值 < ", significance_threshold, "）。")

# 6. 保存显著结果
significant_outfile <- file.path("data", "significant_MR_results.csv")
write_csv(significant_results, significant_outfile)
message("显著结果已保存至: ", significant_outfile)

# 7. 可视化

## 7.1 检查是否有显著结果
if (nrow(significant_results) == 0) {
  message("没有显著结果可供绘图。")
} else {
  
  # 7.2 绘制森林图（Forest Plot）
  # 选择常用的 MR 方法，例如 IVW
  forest_data <- significant_results %>%
    filter(str_detect(method, "IVW")) %>%
    mutate(exposure_outcome = paste(exposure, "→", outcome))
  
  if (nrow(forest_data) > 0) {
    forest_plot <- ggplot(forest_data, aes(x = or, y = reorder(exposure_outcome, or))) +
      geom_point() +
      geom_errorbarh(aes(xmin = or_lci95, xmax = or_uci95), height = 0.2) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
      xlab("Odds Ratio (OR)") +
      ylab("Exposure → Outcome") +
      ggtitle("Forest Plot of Significant MR Results (IVW Method)") +
      theme_minimal()
    
    # 保存森林图
    ggsave(filename = file.path("data", "forest_plot_IVW.png"), plot = forest_plot, width = 8, height = 6)
    message("森林图已保存至: data/forest_plot_IVW.png")
  } else {
    message("没有符合条件的 IVW 方法的显著结果用于绘制森林图。")
  }
  
  # 7.3 绘制火山图（Volcano Plot）
  volcano_data <- merged_results %>%
    mutate(significant = pval_adj < significance_threshold)
  
  volcano_plot <- ggplot(volcano_data, aes(x = -log10(pval), y = -log10(pval_adj), color = significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("grey", "red")) +
    geom_hline(yintercept = -log10(significance_threshold), linetype = "dashed", color = "blue") +
    xlab("-log10(p-value)") +
    ylab("-log10(adjusted p-value)") +
    ggtitle("Volcano Plot of MR Results") +
    theme_minimal()
  
  # 保存火山图
  ggsave(filename = file.path("data", "volcano_plot_MR.png"), plot = volcano_plot, width = 8, height = 6)
  message("火山图已保存至: data/volcano_plot_MR.png")
  
  # 7.4 绘制 QQ 图
  # 计算预期 p 值
  observed_pvals <- merged_results$pval
  expected_pvals <- -log10(ppoints(length(observed_pvals)))
  
  qq_data <- data.frame(
    expected = expected_pvals,
    observed = -log10(sort(observed_pvals))
  )
  
  qq_plot <- ggplot(qq_data, aes(x = expected, y = observed)) +
    geom_point(size = 1, alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    xlab("Expected -log10(p)") +
    ylab("Observed -log10(p)") +
    ggtitle("QQ Plot of MR Results") +
    theme_minimal()
  
  # 保存 QQ 图
  ggsave(filename = file.path("data", "qq_plot_MR.png"), plot = qq_plot, width = 6, height = 6)
  message("QQ 图已保存至: data/qq_plot_MR.png")
  
  # 7.5 绘制条形图展示不同方法的显著结果数量
  method_counts <- significant_results %>%
    group_by(method) %>%
    summarise(count = n()) %>%
    arrange(desc(count))
  
  bar_plot <- ggplot(method_counts, aes(x = reorder(method, count), y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    xlab("MR 方法") +
    ylab("显著结果数量") +
    ggtitle("不同 MR 方法的显著结果数量") +
    theme_minimal() +
    coord_flip()
  
  # 保存条形图
  ggsave(filename = file.path("data", "significant_results_barplot.png"), plot = bar_plot, width = 8, height = 6)
  message("条形图已保存至: data/significant_results_barplot.png")
  
  # 7.6 生成综合图表
  # 如果存在多个图表，可以将它们组合在一起
  # 这里示例组合森林图和条形图
  if (exists("forest_plot") & exists("bar_plot")) {
    combined_plot <- plot_grid(forest_plot, bar_plot, labels = c("A", "B"), ncol = 2)
    ggsave(filename = file.path("data", "combined_plot.png"), plot = combined_plot, width = 16, height = 8)
    message("综合图表已保存至: data/combined_plot.png")
  }
}

# 8. 输出完成信息
message("MR 分析的后续处理和可视化已完成。所有结果和图表已保存在 'data' 目录下。")
