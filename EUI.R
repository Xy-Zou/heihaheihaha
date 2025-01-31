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
library(stringr)

# ========== 2. 读取 immune cell GWAS ID 列表 ========== 
immune_list_file <- "data/731_immune_cell/ICgwas-cd8-id.csv"
datalist <- read.table(
  file         = immune_list_file,
  header       = TRUE,
  sep          = ",",
  quote        = "",
  comment.char = ""
)
imc_ids <- as.vector(datalist$id)

# ========== 2.1 汇总所有 exposure 文件中的 SNP 并去重 ========== 
all_exposure_snps <- c()
for (exposure_id in imc_ids) {
  exposure_file_path <- file.path("data", "731_immune_cell", paste0(exposure_id, ".csv"))
  if (!file.exists(exposure_file_path)) {
    message("[警告] 找不到暴露文件: ", exposure_file_path, "，跳过...")
    next
  }
  exposure_dat <- read.csv(exposure_file_path, header = TRUE, sep = ",", quote = "", comment.char = "")
  all_exposure_snps <- c(all_exposure_snps, exposure_dat$SNP)
}
all_exposure_snps <- unique(all_exposure_snps)
message("汇总得到的 SNP 总数: ", length(all_exposure_snps))

# ========== 3. 本次只分析单个结局文件: ieu-b-10.vcf.gz ========== 
out_file <- "ieu-b-10.vcf.gz"
message("\n========== 开始处理结局文件: ", out_file, " ==========")

## 3.0 先统计多少行是##开头，用于skip
all_lines <- readLines(file.path("data", out_file), n = 200)  # 读前200行来统计
skip_n <- sum(grepl("^##", all_lines))
message("检测到 ", skip_n, " 行 ## 注释行，将在读表时跳过...")

# ========== 3.1 读取原始 VCF 并做初步解析 ========== 
raw_vcf <- read.table(
  file   = file.path("data", out_file),
  header = TRUE,
  comment.char = "",  # 不用默认#，否则#CHROM那行也跳过
  skip   = skip_n,    # 跳过##注释行
  sep    = "\t",
  stringsAsFactors = FALSE
)

# 假设读入后列为 [#CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, ieu-b-10]
# 重命名前10列方便操作
colnames(raw_vcf)[1:10] <- c("CHROM","POS","SNP","REF","ALT","QUAL","FILTER","INFO","FORMAT","GWASdata")

# 拆分 FORMAT (ES:SE:P:AF等) 和 GWASdata (例如 "0.05:0.01:0.001:0.30")
# 下面的 n=4 & colnames 需与你文件实际对应
splitted <- str_split_fixed(raw_vcf$GWASdata, ":", n=4)
colnames(splitted) <- c("BETA","SE","P","AF")

# 拼成一个解析好的 data.frame
raw_vcf_parsed <- cbind(raw_vcf[, c("CHROM","POS","SNP","REF","ALT")], splitted)
raw_vcf_parsed$BETA <- as.numeric(raw_vcf_parsed$BETA)
raw_vcf_parsed$SE   <- as.numeric(raw_vcf_parsed$SE)
raw_vcf_parsed$P    <- as.numeric(raw_vcf_parsed$P)
raw_vcf_parsed$AF   <- as.numeric(raw_vcf_parsed$AF)

# ========== 3.2 写出到临时文件，再用 read_outcome_data() 读取 ========== 
temp_file <- tempfile(pattern = "outcome_tmp_", fileext = ".txt")

# 将解析好的表写到临时文件
write.table(
  x = raw_vcf_parsed,
  file = temp_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# 通过 filename= 来读取，并指定列对应
outcome_dat_all <- read_outcome_data(
  snps                = all_exposure_snps,
  filename            = temp_file,
  sep                 = "\t",
  snp_col             = "SNP",
  beta_col            = "BETA",
  se_col              = "SE",
  effect_allele_col   = "ALT",   # 如果确认ALT就是效应等位基因
  other_allele_col    = "REF",
  eaf_col             = "AF",
  pval_col            = "P"
)

# 清理临时文件
file.remove(temp_file)

if (nrow(outcome_dat_all) == 0) {
  message("该结局文件中没有匹配到任何需要的 SNP，脚本退出...")
  quit("no")
}

# 给数据打上 outcome 标识
outcome_dat_all$id.outcome <- "ieu-b-10"
outcome_dat_all$outcome    <- "ieu-b-10"

# ========== 3.3 并行处理每一个 exposure ========== 
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
    
    exposure_dat$id.exposure <- exposure_id
    exposure_dat$exposure    <- exposure_id
    
    # 筛选本 exposure 需要的 SNP
    snps_needed    <- exposure_dat$SNP
    outcome_subset <- subset(outcome_dat_all, SNP %in% snps_needed)
    if (nrow(outcome_subset) == 0) {
      return(NULL)
    }
    
    # harmonise
    harm_dat <- harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat  = outcome_subset,
      action       = 2
    )
    
    if (nrow(harm_dat) < 3) {
      return(NULL)
    }
    
    # ========== (1) 使用 TwoSampleMR 包进行常规 MR ========== 
    mr_result <- mr(
      harm_dat,
      method_list = c(
        "mr_ivw",
        "mr_egger_regression",
        "mr_weighted_median"
      )
    )
    mr_result_or <- generate_odds_ratios(mr_result)
    
    # ========== (2) 使用 MR-PRESSO ========== 
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
        # 转成OR
        temp_df$or       <- exp(temp_df$b)
        temp_df$or_lci95 <- exp(temp_df$b - 1.96 * temp_df$se)
        temp_df$or_uci95 <- exp(temp_df$b + 1.96 * temp_df$se)
        
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
  mc.cores = 4
)
final_results <- bind_rows(final_results_list)

# ========== 3.4 写出结果 ========== 
out_filename  <- "ieu-b-10_MRresults.csv"
out_full_path <- file.path("data", out_filename)
write.csv(final_results, file = out_full_path, row.names = FALSE)
message("输出结果: ", out_full_path)
message("========== 完成结局文件: ", out_file, " ==========")

message("\n所有分析处理完毕！")

