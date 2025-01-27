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

# 并行包
library(parallel)   # 提供 mclapply
# 或者你也可以这样：
# library(doParallel)
# registerDoParallel(cores = detectCores())

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

# ========== 3. 定义三个结局文件并循环分析 ========== 
outcome_files <- c(
  "summary_stats_release_finngen_R12_FE.gz",
  "summary_stats_release_finngen_R12_FE_MODE.gz",
  "summary_stats_release_finngen_R12_FE_STRICT.gz"
)

for (out_file in outcome_files) {
  
  message("\n========== 开始处理结局文件: ", out_file, " ==========")
  
  # 确定一个简短的 outcome label
  outcome_label <- gsub("\\.gz$", "", out_file)
  
  # ========== 3.1 一次性读取该 outcome 文件中所有需要的 SNP ========== 
  # 注意: 这里只调用一次 read_outcome_data()， snps 参数传入 all_exposure_snps
  outcome_dat_all <- read_outcome_data(
    snps               = all_exposure_snps,
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
  
  # 如果没有任何 SNP 被匹配到，那么这个 outcome 文件就可以直接跳过
  if (nrow(outcome_dat_all) == 0) {
    message("该结局文件中没有匹配到任何需要的 SNP，跳过...")
    next
  }
  
  # 给这批数据打上明确的 outcome 标识
  outcome_dat_all$id.outcome <- outcome_label
  outcome_dat_all$outcome    <- outcome_label
  
  # ========== 3.2 并行处理每一个 exposure ========== 
  # 注意，这里我们用 mclapply 并行，默认为所有可用核心
  # mclapply 在 macOS / Linux 下可以用多核并行；Windows 下则只能用 mc.cores=1 或切换其他并行方式
  final_results_list <- mclapply(
    X         = imc_ids,
    FUN       = function(exposure_id) {
      
      exposure_file_path <- file.path("data", "731_immune_cell", paste0(exposure_id, ".csv"))
      if (!file.exists(exposure_file_path)) {
        # 直接返回NULL, 后面会滤掉
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
        # 没有匹配到的话返回空
        return(NULL)
      }
      
      # 进行 harmonise
      harm_dat <- harmonise_data(
        exposure_dat = exposure_dat,
        outcome_dat  = outcome_subset,
        action       = 2
      )
      
      # 只做 IVW
      mr_result <- mr(harm_dat, method_list = c("mr_ivw"))
      
      # 转换为 OR
      mr_result_or <- generate_odds_ratios(mr_result)
      return(mr_result_or)
    },
    mc.cores = detectCores()  # 或者手动指定一个合适的核心数
  )
  
  # 把 list 合并为一个 data.frame
  final_results <- bind_rows(final_results_list)
  
  # ========== 3.3 写出结果 ========== 
  out_filename  <- gsub(".gz", "_MRresults.csv", out_file)
  out_full_path <- file.path("data", out_filename)
  
  write.csv(final_results, file = out_full_path, row.names = FALSE)
  message("输出结果: ", out_full_path)
  message("========== 完成结局文件: ", out_file, " ==========")
}

message("\n所有文件处理完毕！")
