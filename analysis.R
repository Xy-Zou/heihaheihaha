#! /usr/bin/env Rscript

# ========== 0. 工作路径与环境 ========== 
# 请根据实际情况修改 setwd() 路径
setwd("/Users/wilian/git/heihaheihaha")

# 如果没有安装 R.utils，会导致无法直接 fread gz 文件，所以做个检测
if (!requireNamespace("R.utils", quietly = TRUE)) {
  install.packages("R.utils", repos = "http://cran.us.r-project.org")
}
library(R.utils)

# ========== 1. 加载所需 R 包 ========== 
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
})

# ========== 2. 读取免疫细胞 GWAS ID 列表 ========== 
# 该文件位于 data/731_immune_cell/ 目录下
immune_list_file <- "data/731_immune_cell/ICgwasid.csv"
datalist <- read.table(
  file         = immune_list_file,
  header       = TRUE,
  sep          = ",",
  quote        = "",
  comment.char = ""
)
imc_ids <- as.vector(datalist$id)

# ========== 3. 定义三个结局文件并循环分析 ========== 
# 这里的结局文件都位于 data/ 目录下
outcome_files <- c(
  "summary_stats_release_finngen_R12_FE.gz",
  "summary_stats_release_finngen_R12_FE_MODE.gz",
  "summary_stats_release_finngen_R12_FE_STRICT.gz"
)

# ========== 3.1. 循环处理每个结局文件 ========== 
for (out_file in outcome_files) {
  
  message("========== 开始处理结局文件: ", out_file, " ==========")
  
  # 定义结果存储的 data.frame
  final_results <- data.frame()
  
  # ========== 3.2. 循环处理每一个免疫细胞 GWAS 暴露 ========== 
  for (exposure_id in imc_ids) {
    
    # ========== 3.2.1. 读取暴露数据 ========== 
    # 暴露文件放在 data/731_immune_cell 中，形如 ebi-a-GCSTxxxxxx.csv
    exposure_file_path <- file.path("data", "731_immune_cell", paste0(exposure_id, ".csv"))
    
    # 如果文件不存在，跳过该循环
    if (!file.exists(exposure_file_path)) {
      message("找不到暴露文件: ", exposure_file_path, "，跳过...")
      next
    }
    
    exposure_dat <- read.table(
      file   = exposure_file_path,
      header = TRUE,
      sep    = ",",
      quote  = "",
      comment.char = ""
    )
    
    # 如果暴露数据本身需要重命名列，可以在这里进行
    # exposure_dat <- exposure_dat %>%
    #   rename(
    #     SNP = "snp_col_name",
    #     beta.exposure = "beta_col_name",
    #     se.exposure   = "se_col_name",
    #     effect_allele.exposure = "alt_col_name",
    #     other_allele.exposure  = "ref_col_name"
    #   )
    
    # ========== 3.2.2. 从大文件中读取匹配 SNP 的结局数据 ========== 
    outcome_dat <- read_outcome_data(
      snps               = exposure_dat$SNP,
      filename           = file.path("data", out_file),
      sep                = "\t",
      snp_col            = "rsids",
      beta_col           = "beta",
      se_col             = "sebeta",
      effect_allele_col  = "alt",
      other_allele_col   = "ref",
      eaf_col            = "af_alt",
      pval_col           = "pval"
    )
    
    # 如果没有匹配到数据，跳过
    if (nrow(outcome_dat) == 0) {
      message("  暴露: ", exposure_id, " 与结局文件: ", out_file, " 没有匹配到 SNP。跳过...")
      next
    }
    
    # ========== 3.2.3. 和谐化数据（Harmonise Data） ========== 
    harm_dat <- harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat  = outcome_dat,
      action       = 2  # 如果存在正负链不一致的等位基因，则进行一定校正
    )
    
    # ========== 3.2.4. 进行 MR 分析（以 IVW 为例） ========== 
    mr_result <- mr(
      harm_dat,
      method_list = c("mr_ivw")
    )
    
    # ========== 3.2.5. 转化为 OR（若暴露-结局是二分类） ========== 
    mr_result_or <- generate_odds_ratios(mr_result)
    
    # 合并到总结果中
    final_results <- bind_rows(final_results, mr_result_or)
  }
  
  # ========== 4. 输出结果 ========== 
  # 将 .gz 后缀替换为 _MRresults.csv，例如：summary_stats_release_finngen_R12_FE.gz -> summary_stats_release_finngen_R12_FE_MRresults.csv
  out_filename  <- gsub(".gz", "_MRresults.csv", out_file)
  out_full_path <- file.path("data", out_filename)
  
  write.csv(
    final_results,
    file      = out_full_path,
    row.names = FALSE
  )
  
  message("输出结果: ", out_full_path)
  message("========== 完成结局文件: ", out_file, " ==========")
}

message("所有文件处理完毕！")
