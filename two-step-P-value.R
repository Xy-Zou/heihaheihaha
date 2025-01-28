#!/usr/bin/env Rscript

# ========================== #
#  多重检验分两步校正示例脚本  #
# ========================== #

# ========== 0. 设置工作路径并加载包 ==========
setwd("/Users/wilian/git/heihaheihaha")

# 加载所需 R 包
required_packages <- c("dplyr", "ggplot2", "readr", "stringr", "purrr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# ========== 1. 读取 MR 结果文件 ==========
# 假设之前的 MR 脚本输出了若干个 "_MRresults.csv" 文件到 data 目录
result_files <- list.files(
  path        = "data",
  pattern     = "_MRresults\\.csv$",
  full.names  = TRUE
)

if (length(result_files) == 0) {
  stop("未找到任何 '_MRresults.csv' 结尾的文件，请检查 data 目录。")
}

# 逐一读取
message("读取以下 MR 结果文件：")
print(result_files)

all_results <- lapply(result_files, function(file) {
  read_csv(file, show_col_types = FALSE)
})

# 合并为一个大表
merged_results <- bind_rows(all_results)

# ========== 2. 第 1 阶段校正：在每个 outcome 内部做多重比较校正 ==========
# 这里使用 Benjamini-Hochberg (BH) 方法，可根据需求改成其他方法
# 假设 merged_results 中包含以下列：
# - outcome, exposure, SNP, pval, (可选: b, se, or, or_lci95, or_uci95 等)

# 先检查是否有 outcome / pval 列
if (!all(c("outcome", "pval") %in% colnames(merged_results))) {
  stop("merged_results 中缺失 'outcome' 或 'pval' 列，请检查数据格式。")
}

# 在每个 outcome 内部对 p 值进行校正
merged_results_intra <- merged_results %>%
  group_by(outcome) %>%
  mutate(pval_adj_intra = p.adjust(pval, method = "BH")) %>%
  ungroup()

# 保存第 1 阶段校正后的完整结果
write_csv(merged_results_intra, file = "data/MRresults_intra_corrected.csv")
message("已输出第 1 阶段校正结果：data/MRresults_intra_corrected.csv")

# ========== 3. 筛选在任意 outcome 下显著的 SNP （第 1 阶段结果）==========
# 例如，我们这里用 pval_adj_intra < 0.05 作为阈值
significance_threshold <- 0.05
significant_intra <- merged_results_intra %>%
  filter(pval_adj_intra < significance_threshold)

# 输出第 1 阶段显著结果
write_csv(significant_intra, "data/MRresults_intra_significant.csv")
message("已输出第 1 阶段显著结果：data/MRresults_intra_significant.csv")

# 如果没有任何显著结果，可视具体需求选择是否继续第二阶段
if (nrow(significant_intra) == 0) {
  message("第 1 阶段校正后无显著结果，不进行第二阶段校正。脚本结束。")
  quit(save = "no")
}

# ========== 4. 准备第 2 阶段校正的“重点 SNP”集合 ==========
# 这里以“在任意 outcome 中显著的 SNP 都纳入第二阶段”作为示范
top_snps <- unique(significant_intra$SNP)

# 在 merged_results_intra 中筛选只包含这些 SNP 的所有行（包括所有 outcome）
merged_results_top <- merged_results_intra %>%
  filter(SNP %in% top_snps)

# ========== 5. 第 2 阶段校正：在这些“重点 SNP”上做更严格校正 ==========
# 这里为了示例，演示“跨 outcome 的统一校正”，
# 即把所有与这些 SNP 相关的 p 值合并后再做一次 BH 校正。
# 也可以根据需求做别的组合方式，如只针对相同 outcome 进行二次校正等。
merged_results_top <- merged_results_top %>%
  mutate(pval_adj_second = p.adjust(pval, method = "BH"))

# 保存第二阶段整体结果
write_csv(merged_results_top, "data/MRresults_top_snps_secondwave.csv")
message("已输出纳入第二阶段校正的所有结果：data/MRresults_top_snps_secondwave.csv")

# ========== 6. 筛选第 2 阶段显著结果并输出 ==========
second_significant <- merged_results_top %>%
  filter(pval_adj_second < significance_threshold)

write_csv(second_significant, "data/MRresults_top_snps_secondwave_significant.csv")
message("已输出第 2 阶段显著结果：data/MRresults_top_snps_secondwave_significant.csv")

# ========== 7. 提示结束 ==========
message("分两步校正流程已完成！")
message("第一阶段：每个 outcome 内部校正；第二阶段：对通过第一阶段筛选的 SNP 再做更严格校正。")
