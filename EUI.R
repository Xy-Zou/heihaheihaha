#! /usr/bin/env Rscript

# ========== 0. 工作路径与环境 ========== 
setwd("~/git/heihaheihaha")  # 请根据需要修改工作路径

# ========== 1. 加载所需 R 包 ========== 
# 若某些包未安装，会先尝试安装
if (!requireNamespace("R.utils", quietly = TRUE)) {
  install.packages("R.utils", repos = "http://cran.us.r-project.org")
}
library(R.utils)

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(dplyr)
  library(parallel)
})

# 如果需要使用 MR-PRESSO，请确保安装/加载
if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
  install.packages("MRPRESSO", repos = "http://cran.us.r-project.org")
}
library(MRPRESSO)

# 如需从 IEU OpenGWAS 拉取数据，需要以下几个包
if (!requireNamespace("ieugwasr", quietly = TRUE)) {
  install.packages("ieugwasr", repos = "http://cran.us.r-project.org")
}
if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
  install.packages("VariantAnnotation", repos = "http://cran.us.r-project.org")
}
if (!requireNamespace("gwasvcf", quietly = TRUE)) {
  install.packages("gwasvcf", repos = "http://cran.us.r-project.org")
}
library(ieugwasr)
library(VariantAnnotation)
library(gwasvcf)

# 如果你有自己的 JWT（或不需要），可根据需要设置
# 这里示例展示如何设置 OPENGWAS_JWT
Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.ey...")

# ========== 2. 读取 immune cell GWAS ID 列表 ========== 
immune_list_file <- "data/731_immune_cell/ICgwasid.csv"  # 里面保存了各个 immune cell 暴露的文件 ID
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
  exposure_dat <- read.csv(
    file         = exposure_file_path,
    header       = TRUE,
    sep          = ",",
    quote        = "",
    comment.char = ""
  )
  
  # 提取 SNP 列
  all_exposure_snps <- c(all_exposure_snps, exposure_dat$SNP)
}

all_exposure_snps <- unique(all_exposure_snps)
message("汇总得到的 SNP 总数: ", length(all_exposure_snps))

# ========== 3. 定义 outcome = "ieu-b-10" 并一次性提取结局数据 ========== 
outcome_id <- "ieu-b-10"

# 使用 TwoSampleMR 中的 extract_outcome_data() 函数从 IEU OpenGWAS 提取
# 可以根据需要调整 proxies / rsq / align_alleles / palindromes / maf_threshold 等参数
outcome_dat_all <- extract_outcome_data(
  snps          = all_exposure_snps,
  outcomes      = outcome_id,
  proxies       = 1,    # 是否寻找代理 SNP，0=不找，1=默认找
  rsq           = 0.8,  # 与原 SNP 的 LD 要求
  align_alleles = 1,
  palindromes   = 1,
  maf_threshold = 0.01
)

if (nrow(outcome_dat_all) == 0) {
  stop("在结局 'ieu-b-10' 中没有找到匹配的 SNP，脚本终止...")
}

# 为这批数据打上 outcome 的标识
outcome_dat_all$id.outcome <- outcome_id
outcome_dat_all$outcome    <- outcome_id

# ========== 4. 并行处理每一个 exposure，执行 MR 分析 ========== 
final_results_list <- mclapply(
  X         = imc_ids,
  FUN       = function(exposure_id) {
    
    exposure_file_path <- file.path("data", "731_immune_cell", paste0(exposure_id, ".csv"))
    if (!file.exists(exposure_file_path)) {
      message("[警告] 找不到暴露文件: ", exposure_file_path, "，跳过...")
      return(NULL)
    }
    
    # 读取 exposure
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
    
    # 从 outcome_dat_all 中筛选本暴露所需 SNP
    snps_needed    <- exposure_dat$SNP
    outcome_subset <- subset(outcome_dat_all, SNP %in% snps_needed)
    
    if (nrow(outcome_subset) == 0) {
      # 说明该暴露对应的 SNP 在结局里均无匹配
      return(NULL)
    }
    
    # 进行 harmonise（对齐等位基因、方向等）
    harm_dat <- harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat  = outcome_subset,
      action       = 2
    )
    
    # 如果有效 SNP 太少，也可能导致后面方法报错，这里设个简单判断
    if (nrow(harm_dat) < 3) {
      # 如果 <3 个独立 SNP，Egger/IVW 等很多方法都不稳，这里直接返回 NULL
      return(NULL)
    }
    
    # ========== 4.1: 使用 TwoSampleMR 包内置的方法（IVW / Egger / Weighted median） ========== 
    mr_result <- mr(
      harm_dat,
      method_list = c(
        "mr_ivw",
        "mr_egger_regression",
        "mr_weighted_median"
        # 若想加入 Weighted mode, 可添加 "mr_weighted_mode"
      )
    )
    mr_result_or <- generate_odds_ratios(mr_result)
    
    # ========== 4.2: 使用 MR-PRESSO 包进行稳健分析 ========== 
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
      
      presso_main <- as.data.frame(presso_out$`Main MR results`)
      
      # 如果结果不是空，可以将其整理成与 TwoSampleMR 类似的输出结构
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
        
        # 将因果估计转换为 OR
        temp_df$or       <- exp(temp_df$b)
        temp_df$or_lci95 <- exp(temp_df$b - 1.96 * temp_df$se)
        temp_df$or_uci95 <- exp(temp_df$b + 1.96 * temp_df$se)
        
        presso_res <- temp_df
      }
    }, silent = TRUE)
    
    # 合并各方法结果
    if (is.null(presso_res)) {
      final_res <- mr_result_or
    } else {
      final_res <- dplyr::bind_rows(mr_result_or, presso_res)
    }
    
    return(final_res)
  },
  mc.cores = 4  # 根据机器配置，可以调整并行核心数
)

# 整合所有结果
final_results <- dplyr::bind_rows(final_results_list)

# ========== 5. 写出结果 ========== 
out_filename  <- paste0("MRresults_", outcome_id, ".csv")  # 例如: MRresults_ieu-b-10.csv
out_full_path <- file.path("data", out_filename)

write.csv(final_results, file = out_full_path, row.names = FALSE)
message("输出结果: ", out_full_path)
message("所有文件处理完毕！")