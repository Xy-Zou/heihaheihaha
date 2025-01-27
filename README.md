# heihaheihaha
```bash
mamba env create -f environment.yml
mamba activate mr_env
```
Test R packages
```bash
Rscript -e "library(TwoSampleMR); library(dplyr)"
```
##直接分析-----
```r
# 2.Iterate over all immune cells data to get MR analysis results
#2.1
setwd()
# 2.2 load R package
library(TwoSampleMR)

# 2.3 读入731种免疫细胞的list
datalist <- read.table("ICgwasid.csv",header = TRUE,sep = ",",quote = "", comment.char = "")
imc <- as.vector(datalist$id)

# 2.4 读入结局数据
outcomefile <- 'summary_stats_release_finngen_R12_FE.gz'

# 2.5 分析
result <- data.frame()
for (i in imc ) {
  exposure_rt<-read.table(file =paste0(i,".csv"),header = TRUE,sep = ",")
  #下面几个名称你看着帮我改好了（狗头）
  outcome_rt <- read_outcome_data(
    snps = exposure_rt$SNP,
    filename = outcomefile,
    sep = "\t",
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    eaf_col = "af_alt",
    pval_col = "pval")
  harm_rt <- harmonise_data(
    exposure_dat =  exposure_rt, 
    outcome_dat = outcome_rt,action=2)
  mr_result<- mr(harm_rt,method_list=c("mr_ivw"))
  result_or=generate_odds_ratios(mr_result) 
  result <- rbind(result, result_or)
  
}

# 2.6 Output the results of Mendelian randomization
write.csv(result,"MRresults.csv",sep = "")
```
