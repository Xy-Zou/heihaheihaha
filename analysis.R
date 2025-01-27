#! /usr/bin/env Rscript

# ========== 0. 工作路径与环境 ========== 
setwd("/Users/wilian/git/heihaheihaha")

# ========== 1. 加载所需 R 包 ========== 
if (!requireNamespace("R.utils", quietly = TRUE)) {
  install.packages("R.utils", repos = "http://cran.us.r-project.org")
}
library(R.utils)

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
})

# ========== 2. 读取免疫细胞 GWAS ID 列表 ========== 
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
outcome_files <- c(
  "summary_stats_release_finngen_R12_FE.gz",
  "summary_stats_release_finngen_R12_FE_MODE.gz",
  "summary_stats_release_finngen_R12_FE_STRICT.gz"
)

for (out_file in outcome_files) {
  
  message("========== 开始处理结局文件: ", out_file, " ==========")
  
  # 给这一次循环的 outcome 一个更“短而有意义”的名字
  # 例如你可以只保留去掉后缀 .gz 后的名称
  outcome_label <- gsub("\\.gz$", "", out_file)
  
  # 定义结果存储的 data.frame
  final_results <- data.frame()
  
  # ---------- 3.1. 循环处理每一个免疫细胞 GWAS 暴露 ----------
  for (exposure_id in imc_ids) {
    
    # ========== 3.1.1 读取暴露数据 ========== 
    exposure_file_path <- file.path("data", "731_immune_cell", paste0(exposure_id, ".csv"))
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
    
    # 指定暴露名称，让后续不会出现随机 ID
    # 在 TwoSampleMR 数据结构中，通常使用这两个字段
    exposure_dat$id.exposure <- exposure_id
    exposure_dat$exposure    <- exposure_id
    
    # ========== 3.1.2. 从大文件中读取匹配 SNP 的结局数据 ========== 
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
    
    if (nrow(outcome_dat) == 0) {
      message("  暴露: ", exposure_id, " 与结局文件: ", out_file, " 没有匹配到 SNP。跳过...")
      next
    }
    
    # 指定结局的名称，避免出现“FhDktu”或默认的 "outcome" 等随机标识
    outcome_dat$id.outcome <- outcome_label
    outcome_dat$outcome    <- outcome_label
    
    # ========== 3.1.3. 和谐化数据 ========== 
    harm_dat <- harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat  = outcome_dat,
      action       = 2
    )
    
    # ========== 3.1.4. MR 分析（仅 IVW） ========== 
    mr_result <- mr(
      harm_dat,
      method_list = c("mr_ivw")
    )
    
    # 转化为 OR
    mr_result_or <- generate_odds_ratios(mr_result)
    
    final_results <- bind_rows(final_results, mr_result_or)
  }
  
  # ========== 4. 输出结果 ========== 
  out_filename  <- gsub(".gz", "_MRresults.csv", out_file)
  out_full_path <- file.path("data", out_filename)
  
  write.csv(final_results, file = out_full_path, row.names = FALSE)
  
  message("输出结果: ", out_full_path)
  message("========== 完成结局文件: ", out_file, " ==========")
}

message("所有文件处理完毕！")
