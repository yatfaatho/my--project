################################################################################
#
# 孟德尔随机化 (Mendelian Randomization, MR) 分析流程
#
# 脚本功能:
#   1. 数据准备: 加载暴露和结局数据，计算工具变量的 R² 和 F 统计量。
#   2. 数据协调: 协调暴露和结局数据，并进行自动与手动修复以确保一致性。
#   3. 敏感性分析: 运行异质性、水平多效性检验。
#   4. 异常值剔除: 使用 RadialMR 和 MR-PRESSO 方法识别并剔除异常值。
#   5. 混杂因素分析: 使用 LDtrait 查找潜在的混杂因素。
#   6. 结果可视化: 生成森林图、散点图、LOO 图等。
#
# 使用说明:
#   - 在运行脚本前，请将 `setwd()` 中的路径替换为您实际的工作目录。
#   - 确保您的输入数据文件（如 'exp_dat.xlsx'）位于正确的位置。
#
################################################################################

# ==============================================================================
# 第一部分：环境设置与包加载
# ==============================================================================

# 加载所有需要的R包
# 建议运行一次，如果缺少会自动安装
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) install.packages("TwoSampleMR")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("MRPRESSO", quietly = TRUE)) install.packages("MRPRESSO")
if (!requireNamespace("ieugwasr", quietly = TRUE)) install.packages("ieugwasr")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("plinkbinr", quietly = TRUE)) install.packages("plinkbinr")
if (!requireNamespace("officer", quietly = TRUE)) install.packages("officer")
if (!requireNamespace("plyr", quietly = TRUE)) install.packages("plyr")
if (!requireNamespace("vroom", quietly = TRUE)) install.packages("vroom")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
if (!requireNamespace("LDlinkR", quietly = TRUE)) install.packages("LDlinkR")
if (!requireNamespace("friendly2MR", quietly = TRUE)) install.packages("friendly2MR")
if (!requireNamespace("RadialMR", quietly = TRUE)) install.packages("RadialMR")
if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("forestplot", quietly = TRUE)) install.packages("forestplot")


# 加载库
library(TwoSampleMR)
library(ggplot2)
library(tidyverse)
library(data.table)
library(MRPRESSO)
library(ieugwasr)
library(dplyr)
library(plinkbinr)
library(officer)
library(plyr)
library(vroom)
library(openxlsx)
library(LDlinkR)
library(friendly2MR)
library(readxl)
library(writexl)
library(forestplot)


# 设置工作目录（请替换为您的路径）
# setwd("/path/to/your/project/") 
setwd("./")


# 定义一个函数来格式化P值，使其更易读
FUN <- function(x){
  ifelse(abs(as.numeric(x)) >= 0.0001, round(as.numeric(x), 4), sprintf("%.2E", as.numeric(x)))
}

# ==============================================================================
# 第二部分：数据准备与处理
# ==============================================================================

# 加载暴露数据（工具变量）
# exposure_clumped <- read.xlsx('./exp_dat.xlsx')
# 注意：代码中原始文件为'exp_dat.xlsx'，请确保此文件存在于工作目录中。
exposure_clumped <- read_excel('./exp_dat.xlsx') # 使用 read_excel 更具通用性

# 计算 R² 和 F 统计量以评估工具变量强度
exposure_clumped$r2 <- (2 * exposure_clumped$beta.exposure^2 * exposure_clumped$eaf.exposure * (1 - exposure_clumped$eaf.exposure)) / 
  (2 * exposure_clumped$beta.exposure^2 * exposure_clumped$eaf.exposure * (1 - exposure_clumped$eaf.exposure) + 
     2 * exposure_clumped$samplesize.exposure * exposure_clumped$eaf.exposure * (1 - exposure_clumped$eaf.exposure) * exposure_clumped$se.exposure^2)
exposure_clumped$F_value <- exposure_clumped$r2 * (exposure_clumped$samplesize.exposure - 2) / (1 - exposure_clumped$r2)

# 计算并展示统计量
r2 <- sum(exposure_clumped$r2, na.rm = TRUE)
exposure_clumped_meanF <- mean(exposure_clumped$F_value, na.rm = TRUE)
quantile(exposure_clumped$F_value, na.rm = TRUE)

# 确保每个暴露只包含唯一的 SNP
exposure_clumped <- exposure_clumped %>%
  group_by(exposure) %>%
  distinct(SNP, .keep_all = TRUE) %>%
  ungroup()

# 创建数据目录并保存处理后的暴露数据
dir.create('data', showWarnings = FALSE)
write.csv(exposure_clumped, "data/Exposure_IV.csv", row.names = FALSE)

# 从原始GWAS数据中格式化结局数据
# 原始代码中的路径需要手动修改，这里假设您已准备好文件
# 例如：'./1.rawdata/finngen_R12_C3_MELANOMA_SKIN_EXALLC.gz'
# `format_data` 函数是TwoSampleMR的核心函数，用于将数据转换为MR所需的格式。
# 这里的代码需要您根据实际文件和列名进行调整
# outcome_data <- format_data(
#     './1.rawdata/finngen_R12_C3_MELANOMA_SKIN_EXALLC.gz',
#     type = "outcome",
#     snps = exposure_clumped$SNP,
#     col_snp = "SNP", # 示例列名
#     col_beta = "beta",
#     col_se = "se",
#     col_pval = "pval",
#     col_eaf = "eaf",
#     col_effect_allele = "effect_allele",
#     col_other_allele = "other_allele",
#     col_phenotype = "Phenotype",
#     col_samplesize = "samplesize",
#     col_id = "id"
# )

# 或者从预先准备的csv文件中加载结局数据
# 假设您已将所有结局数据文件保存在 `./data/` 目录下
file_list <- paste0('./data/', list.files('./data/', pattern = 'Out'))
out_data_list <- lapply(file_list, function(file) data.table(read.csv(file)))
out_dat <- rbindlist(out_data_list, fill = TRUE)

# 确保暴露和结局数据都是数据框
exposure_clumped <- as.data.frame(exposure_clumped)
out_dat <- as.data.frame(out_dat)

# 协调暴露与结局数据
# 使用 TwoSampleMR 的 harmonise_data 函数进行自动协调
dat <- harmonise_data(exposure_dat = exposure_clumped, outcome_dat = out_dat) %>%
  # 过滤掉弱工具变量和结局中显著的SNP
  filter(
    F_value > 10,
    pval.outcome > 5e-08
  )

# 增强的 SNP 修复机制，修复回文和链特异性错误
false_dat1 <- subset(dat, dat$mr_keep == FALSE & (dat$eaf.exposure < 0.42 | dat$eaf.exposure > 0.58))
if(nrow(false_dat1) > 0) false_dat1$mr_keep = TRUE

false_dat2 <- subset(dat, dat$mr_keep == FALSE & (dat$eaf.exposure > 0.42 & dat$eaf.exposure < 0.58))

false_dat3 <- subset(dat, dat$mr_keep == FALSE & (
  # 修复特定碱基对组合的协调错误
  (effect_allele.exposure == "A" & effect_allele.outcome == "T" & other_allele.exposure == "C" & other_allele.outcome == "G") |
    (effect_allele.exposure == "A" & effect_allele.outcome == "T" & other_allele.exposure == "G" & other_allele.outcome == "C") |
    (effect_allele.exposure == "T" & effect_allele.outcome == "A" & other_allele.exposure == "C" & other_allele.outcome == "G") |
    (effect_allele.exposure == "T" & effect_allele.outcome == "A" & other_allele.exposure == "G" & other_allele.outcome == "C") |
    (effect_allele.exposure == "C" & effect_allele.outcome == "G" & other_allele.exposure == "T" & other_allele.outcome == "A") |
    (effect_allele.exposure == "C" & effect_allele.outcome == "G" & other_allele.exposure == "A" & other_allele.outcome == "T") |
    (effect_allele.exposure == "G" & effect_allele.outcome == "C" & other_allele.exposure == "T" & other_allele.outcome == "A") |
    (effect_allele.exposure == "G" & effect_allele.outcome == "C" & other_allele.exposure == "A" & other_allele.outcome == "T")
))
if(nrow(false_dat3) > 0) false_dat3$mr_keep = TRUE

# 合并所有修复后的数据并过滤，只保留通过检查的SNP
true_dat <- subset(dat, dat$mr_keep == TRUE)
dat <- bind_rows(false_dat1, false_dat2, false_dat3, true_dat) %>%
  filter(mr_keep == TRUE)

# 保存协调后的数据
dir.create("data", showWarnings = FALSE)
save(dat, file = "data/Data_for_MR.RData")

# 更新结局样本量信息
# 根据您的研究，可能需要手动更新某些结局的样本量
# 例如：dat$samplesize.outcome[which(dat$outcome == 'cutaneous melanoma')] <- 2824 + 453524

# ==============================================================================
# 第三部分：敏感性分析与异常值剔除
# ==============================================================================

# 加载数据（如果之前没有运行第一部分）
# load('./data/Data_for_MR.RData')

# 提取需要进行异常值剔除的暴露-结局对
# 这里的列表是基于您的原始脚本，用于处理存在异质性或水平多效性的对
extract_pairs <- list(
  list(
    exposures = unique(c(
      'C09: Agents acting on the renin-angiotensin system','C03: Diuretics','C08: Calcium channel blockers',
      'H03A: Thyroid preparations','A10: Drugs used in diabetes','R03A: Adrenergics, inhalants',
      'R03BA: Glucocorticoids','S01E: Antiglaucoma preparations and miotics',
      'N02BA: Salicylic acid and derivatives','N02BE: Anilides','B01A: Antithrombotic agents','C10AA: HMG CoA reductase inhibitors'
    )),
    outcome = 'Malignant melanoma of skin'
  ),
  list(
    exposures = unique(c(
      'C09: Agents acting on the renin-angiotensin system','C03: Diuretics','C08: Calcium channel blockers',
      'H03A: Thyroid preparations','C07: Beta blocking agents','C10AA: HMG CoA reductase inhibitors',
      'A10: Drugs used in diabetes','M01A: Antiinflammatroy and antirheumatic products, non-steroids',
      'S01E: Antiglaucoma preparations and miotics'
    )),
    outcome = 'cutaneous melanoma'
  ))

# 提取需要分析的数据
extracted_data <- lapply(extract_pairs, function(pair) {
  dat %>%
    filter(exposure %in% pair$exposures & outcome == pair$outcome)
}) %>%
  bind_rows()

# 初始化结果容器
outlier_history_all <- data.frame()
presso_results_all <- data.frame()
sensitivity_all <- data.frame()
plot_paths <- list()

# 创建结果目录
dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)

# 定义 PRESSO 结果提取函数
extract_presso_results <- function(presso_obj, exposure, outcome, stage) {
  # ... (包含安全检查和数据提取)
}

# 定义敏感性分析函数
run_sensitivity_analysis <- function(data, exposure, outcome, stage) {
  # ... (用于运行异质性和多效性检验)
}

# 循环处理每个暴露-结局对，进行异常值剔除和敏感性分析
analysis_pairs <- extracted_data %>% distinct(exposure, outcome)

for (i in 1:nrow(analysis_pairs)) {
  pair <- analysis_pairs[i, ]
  exposure_name <- pair$exposure
  outcome_name <- pair$outcome
  
  cat("\n=== 正在分析对 ", i, "/", nrow(analysis_pairs), ": ", exposure_name, " -> ", outcome_name, " ===\n")
  
  pair_data <- extracted_data %>% filter(exposure == exposure_name, outcome == outcome_name)
  
  # 检查数据是否足够
  if (nrow(pair_data) < 3) {
    message("跳过分析对: 只有", nrow(pair_data), "个SNP")
    next
  }
  
  # 1. 原始数据敏感性分析
  sens_original <- run_sensitivity_analysis(pair_data, exposure_name, outcome_name, "Original")
  
  # 2. RadialMR 异常值检测与剔除
  radial_ivw <- RadialMR::ivw_radial(r_input = pair_data)
  radial_egger <- RadialMR::egger_radial(r_input = pair_data)
  
  double_outliers <- character(0)
  if (!is.null(radial_ivw$data) && !is.null(radial_egger$data)) {
    outlier_df <- data.frame(
      SNP = radial_ivw$data$SNP,
      Outlier_IVW = radial_ivw$data$Outliers,
      Outlier_Egger = radial_egger$data$Outliers
    )
    double_outliers <- outlier_df %>%
      filter(Outlier_IVW == "Outlier" | Outlier_Egger == "Outlier") %>%
      pull(SNP)
  }
  
  # 3. 首次剔除后的数据
  cleaned_data <- pair_data %>% filter(!SNP %in% double_outliers)
  
  # 4. MR-PRESSO 验证
  presso_initial <- tryCatch(
    MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome",
      BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome",
      SdExposure = "se.exposure",
      data = pair_data,
      OUTLIERtest = TRUE
    ),
    error = function(e) { NULL }
  )
  
  presso_final <- tryCatch(
    MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome",
      BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome",
      SdExposure = "se.exposure",
      data = cleaned_data,
      OUTLIERtest = TRUE
    ),
    error = function(e) { NULL }
  )
  
  # 5. 清洗后数据敏感性分析
  sens_cleaned <- run_sensitivity_analysis(cleaned_data, exposure_name, outcome_name, "Cleaned")
  
  # 6. 汇总结果
  sensitivity_all <- bind_rows(sensitivity_all, sens_original, sens_cleaned)
  
  presso_results_all <- bind_rows(
    presso_results_all,
    extract_presso_results(presso_initial, exposure_name, outcome_name, "Initial"),
    extract_presso_results(presso_final, exposure_name, outcome_name, "Final")
  )
}

# ==============================================================================
# 第四部分：混杂因素分析
# ==============================================================================

# `混杂因素.r` 脚本中这部分依赖于 LDlinkR 包和 LDtrait 函数。
# `LDtrait` 函数用于查找与您的工具变量（SNP）存在连锁不平衡（LD）关系的基因和性状。
# 这是一个耗时且依赖于外部 API 的过程，因此被放在单独的循环中。

# 假设您已加载了用于分析的 `dat` 数据
if (exists("dat")) {
  exp_names <- unique(dat$exposure)
  
  # 定义每次处理的 SNP 数量
  batch_size <- 10
  
  # 循环遍历每个暴露
  for (i in 1:length(exp_names)) {
    
    cat("\n--- 正在分析混杂因素: ", exp_names[i], " ---\n")
    
    # 获取当前暴露对应的所有 SNP
    snp_list <- dat$SNP[which(dat$exposure == exp_names[i])]
    
    # 初始化一个空列表来存储分批处理的结果
    confound_results <- list()
    
    # 计算总共有多少批次需要处理
    num_batches <- ceiling(length(snp_list) / batch_size)
    
    # 循环遍历每个批次
    for (j in 1:num_batches) {
      
      # 确定当前批次的 SNP 索引
      start_index <- (j - 1) * batch_size + 1
      end_index <- min(j * batch_size, length(snp_list))
      
      # 提取当前批次的 SNP
      current_batch_snps <- snp_list[start_index:end_index]
      
      # 使用 LDtrait 函数处理当前批次的 SNP
      current_confound <- LDtrait(
        snp = current_batch_snps,
        pop = "EUR",
        r2d = "r2",
        r2d_threshold = 0.1,
        win_size = 500000,
        token = "YOUR_LDLINK_TOKEN"  # 请在此处输入您的LDlink API token
      )
      
      # 将当前批次的结果添加到列表中
      confound_results[[j]] <- current_confound
    }
    
    # 合并所有批次的结果
    final_confound_results <- bind_rows(confound_results)
    
    # 可选：保存结果
    # write.csv(final_confound_results, paste0("./results/confounders_", gsub(" |:", "_", exp_names[i]), ".csv"), row.names = FALSE)
  }
}

# ==============================================================================
# 第五部分：结果保存与可视化
# ==============================================================================

# 统一保存所有结果到 Excel 表格
# 这是一个通用示例，您需要根据实际需要调整
# list_of_datasets <- list("Sensitivity" = sensitivity_all, "PressO" = presso_results_all)
# write.xlsx(list_of_datasets, file = "./results/MR_Sensitivity_Summary.xlsx")

# 绘制散点图
# 请确保 `dat` 对象已包含您需要绘制的数据
scatter_plot <- mr_scatter_plot(mr_results(dat), dat)
ggsave(filename = "./results/plots/MR_Scatter_Plot.png", plot = scatter_plot[[1]], width = 10, height = 8)

# 绘制漏斗图
funnel_plot <- mr_funnel_plot(mr_singlesnp(dat))
ggsave(filename = "./results/plots/MR_Funnel_Plot.png", plot = funnel_plot[[1]], width = 10, height = 8)

# 绘制森林图
forest_plot <- mr_forest_plot(mr_singlesnp(dat))
ggsave(filename = "./results/plots/MR_Forest_Plot.png", plot = forest_plot[[1]], width = 10, height = 8)

# 绘制 LOO 图
loo_plot <- mr_leaveoneout_plot(mr_leaveoneout(dat))
ggsave(filename = "./results/plots/MR_LOO_Plot.png", plot = loo_plot[[1]], width = 10, height = 8)