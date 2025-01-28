#! /usr/bin/env Rscript

# ========== 0. 工作路径与环境 ========== 
setwd("~/git/heihaheihaha")

# ========== 1. 加载所需 R 包 ========== 
if (!requireNamespace("R.utils", quietly = TRUE)) {
  install.packages("R.utils", repos = "http://cran.us.r-project.org")
}
library(R.utils)

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
})

if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
  install.packages("MRPRESSO", repos = "http://cran.us.r-project.org")
}
library(MRPRESSO)

library(parallel)

# ========== 2. 读取 immune cell GWAS ID 列表 ========== 
immune_list_file <- "data/731_immune_cell/ICgwasid.csv"
datalist <- read.table(
  file         = immune_list_file,
  header       = TRUE,
  sep          = ",",
  quote        = "",
  comment.char = ""
)
imc_ids <- as.vector(datalist$id)

# ========== 2.1 先读取所有 exposure 文件中的 SNP 并去重 ========== 
all_exposure_snps <- c()

for (exposure_id in imc_ids) {
  exposure_file_path <- file.path("data", "731_immune_cell", paste0(exposure_id, ".csv"))
  if (!file.exists(exposure_file_path)) {
    message("[警告] 找不到暴露文件: ", exposure_file_path, "，跳过...")
    next
  }
  exposure_dat <- read.csv(exposure_file_path, header = TRUE, sep = ",", quote = "", comment.char = "")
  
  # 提取 SNP 列
  all_exposure_snps <- c(all_exposure_snps, exposure_dat$SNP)
}

all_exposure_snps <- unique(all_exposure_snps)
message("汇总得到的 SNP 总数: ", length(all_exposure_snps))

# ========== 3. 分析新的结局 ieu-b-10.vcf.gz ========== 

# 指定结局标签
outcome_label <- "ieu-b-10"

# 指向下载好的结局文件 (标准化 IEU VCF 格式)
outcome_file_path <- file.path("data", "ieu-b-10.vcf.gz")

# ========== 3.1 读取所有需要的 SNP ========== 
#   snp_col            = "ID"   (VCF ID 列包含 rsID)
#   beta_col           = "BETA"
#   se_col             = "SE"
#   effect_allele_col  = "ALT"
#   other_allele_col   = "REF"
#   eaf_col            = "AF"
#   pval_col           = "P"
outcome_dat_all <- read_outcome_data(
  snps               = all_exposure_snps,
  filename           = outcome_file_path,
  sep                = "\t",
  snp_col            = "ID",
  beta_col           = "BETA",
  se_col             = "SE",
  effect_allele_col  = "ALT",
  other_allele_col   = "REF",
  eaf_col            = "AF",
  pval_col           = "P"
)

if (nrow(outcome_dat_all) == 0) {
  stop("该结局文件中没有匹配到任何需要的 SNP，请检查 SNP 列名或文件内容。")
}

# 给这批数据打上 outcome 标识
outcome_dat_all$id.outcome <- outcome_label
outcome_dat_all$outcome    <- outcome_label

message("\n========== 开始处理结局: ", outcome_label, " ==========")

# ========== 3.2 并行处理每一个 exposure ========== 
# 这里把 mc.cores 限制为 4
final_results_list <- mclapply(
  X         = imc_ids,
  FUN       = function(exposure_id) {
    
    exposure_file_path <- file.path("data", "731_immune_cell", paste0(exposure_id, ".csv"))
    if (!file.exists(exposure_file_path)) {
      message("[警告] 找不到暴露文件: ", exposure_file_path, "，跳过...")
      return(NULL)
    }
    
    # 读取 exposure 数据
    exposure_dat <- read.csv(
      file         = exposure_file_path,
      header       = TRUE,
      sep          = ",",
      quote        = "",
      comment.char = ""
    )
    
    # 指定暴露名称
    exposure_dat$id.exposure <- exposure_id
    exposure_dat$exposure    <- exposure_id
    
    # 从 outcome_dat_all 中筛选本 exposure 所需 SNP
    snps_needed    <- exposure_dat$SNP
    outcome_subset <- subset(outcome_dat_all, SNP %in% snps_needed)
    
    if (nrow(outcome_subset) == 0) {
      return(NULL)
    }
    
    # 进行 harmonise
    harm_dat <- harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat  = outcome_subset,
      action       = 2
    )
    
    # 如果有效 SNP 太少，也可能导致后面某些方法报错，可以根据需要做过滤
    if (nrow(harm_dat) < 3) {
      # MR-Egger 和其他一些方法需要至少 3~4 个及以上的独立 SNP
      return(NULL)
    }
    
    # ========== 3.2.1: 使用 TwoSampleMR 包内置的方法 ========== 
    mr_result <- mr(
      harm_dat,
      method_list = c(
        "mr_ivw",              # IVW 
        "mr_egger_regression", # MR-Egger
        "mr_weighted_median"   # Weighted median
        # 可酌情加入其他方法, 如 "mr_weighted_mode"
      )
    )
    mr_result_or <- generate_odds_ratios(mr_result)
    
    # ========== 3.2.2: 使用 MR-PRESSO ========== 
    presso_res <- NULL
    try({
      presso_out <- mr_presso(
        BetaOutcome    = "beta.outcome", 
        BetaExposure   = "beta.exposure",
        SdOutcome      = "se.outcome", 
        SdExposure     = "se.exposure", 
        OUTLIERtest    = TRUE,
        DISTORTIONtest = TRUE,
        data           = harm_dat,
        NbDistribution = 1000, 
        SignifThreshold = 0.05
      )
      
      # 提取 MR-PRESSO 的主要结果
      presso_main <- as.data.frame(presso_out$`Main MR results`)
      
      if (nrow(presso_main) > 0) {
        temp_df <- data.frame(
          outcome       = unique(harm_dat$outcome),
          exposure      = unique(harm_dat$exposure),
          method        = paste0("MR-PRESSO (", rownames(presso_main), ")"),
          nsnp          = length(unique(harm_dat$SNP)),
          b             = as.numeric(presso_main[ , "Causal.Estimate"]),
          se            = as.numeric(presso_main[ , "Sd"]),
          pval          = as.numeric(presso_main[ , "Pvalue"]),
          id.outcome    = unique(harm_dat$id.outcome),
          id.exposure   = unique(harm_dat$id.exposure)
        )
        
        # 转换为 OR
        temp_df$or        = exp(temp_df$b)
        temp_df$or_lci95  = exp(temp_df$b - 1.96 * temp_df$se)
        temp_df$or_uci95  = exp(temp_df$b + 1.96 * temp_df$se)
        
        presso_res <- temp_df
      }
    }, silent = TRUE)
    
    # 合并结果
    if (is.null(presso_res)) {
      final_res <- mr_result_or
    } else {
      final_res <- bind_rows(mr_result_or, presso_res)
    }
    
    return(final_res)
  },
  mc.cores = 4  # 调整为你希望的并行核心数
)

final_results <- bind_rows(final_results_list)

# ========== 3.3 写出结果 ========== 
out_filename  <- paste0(outcome_label, "_MRresults.csv")
out_full_path <- file.path("data", out_filename)

write.csv(final_results, file = out_full_path, row.names = FALSE)
message("输出结果: ", out_full_path)
message("========== 完成结局: ", outcome_label, " ==========")

message("\n所有文件处理完毕！")
