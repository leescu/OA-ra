# +================================================+ ####
# +====Section 0. Packages and setwd===============+ ####
# +================================================+ #### 
setwd("K:/UK biobank/Cohort data")

library(dplyr)
library(lubridate)
library(readr)
library(survival)
# +================================================+ ####
# +====Section 1. Oral disease to RA ==============+ ####
# +================================================+ ####
# >>>>> section 1.1 Combine ####
# ---- 1.1 读取与合并 ----
time_df          <- read_csv("Original data/Date_data.csv")          # 包含 center_date, death_date, censor_date
cov_df           <- read_csv("Original data/Covariate_data.csv")     # 协变量
oral_disease_df  <- read_csv("Original data/Oral_disease_data.csv")  # 含 K02_date/K04_date/K05_date/K12_date
oral_health_df   <- read_csv("Original data/Oral_health_data.csv")   # 口腔健康自评
RA_df            <- read_csv("Original data/Arthritis_data.csv")    # 含 rheumatoid_arthritis_date/define

# 协变量（含 Diabetes，便于基线排除用）
cov_df <- cov_df %>%
  select(eid, Age, Sex, Ethnicity, Education, TDI_quantile, Smoke, Alcohol, PA, BMI,
         Diabetes, Hypertension)

# 口腔健康（保留以备敏感性分析用）
oral_health_df <- oral_health_df %>%
  mutate(oral_health = ifelse(oral == "None of the above", 0, 1)) %>%
  select(eid, oral_health)

# RA 信息
RA_df <- RA_df %>%
  mutate(rheumatoid_arthritis_date = ymd(rheumatoid_arthritis_date)) %>%
  select(eid,
         RA_date  = rheumatoid_arthritis_date,
         RA_define = rheumatoid_arthritis_define)

# 合并
analysis_df <- time_df %>%
  left_join(oral_disease_df, by = "eid") %>%
  left_join(oral_health_df,  by = "eid") %>%
  left_join(RA_df,           by = "eid") %>%
  left_join(cov_df,          by = "eid")

# 主分析：排除基线已 RA
analysis_df <- analysis_df %>%
  filter(is.na(RA_date) | RA_date >= center_date) %>%
  filter(is.na(Diabetes) | Diabetes != 1)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++ ####
# >>>>> section 1.2 Dataset ####
# 函数：生成 Main 和 Sensitivity td 数据集
make_td_datasets_oral_to_DM <- function(df, exposure_var) {
  exp_date_col <- paste0(exposure_var, "_date")
  
  # Main：标准样本
  d_main <- df
  
  # Sensitivity：去掉“暴露日期缺失但 oral_health=1”的个体（未暴露组更‘干净’）
  d_sens <- df %>%
    filter(!(is.na(.data[[exp_date_col]]) & oral_health == 1))
  
  build_td <- function(d) {
    td_list <- vector("list", length = nrow(d))
    for (i in seq_len(nrow(d))) {
      exp_date <- d[[exp_date_col]][i]
      start1   <- d$center_date[i]
      end_all  <- min(c(d$RA_date[i], d$death_date[i], d$censor_date[i]), na.rm = TRUE)
      if (is.infinite(end_all)) end_all <- d$censor_date[i]
      
      if (!is.na(exp_date) && exp_date > start1 && exp_date < end_all) {
        td_list[[i]] <- bind_rows(
          data.frame(eid=d$eid[i], center_date=start1, start=start1, stop=exp_date,
                     exposure=0, status=0),
          data.frame(eid=d$eid[i], center_date=start1, start=exp_date, stop=end_all,
                     exposure=1,
                     status=ifelse(!is.na(d$RA_date[i]) && d$RA_date[i] <= end_all, 1, 0))
        )
      } else {
        td_list[[i]] <- data.frame(
          eid=d$eid[i], center_date=start1, start=start1, stop=end_all,
          exposure = ifelse(!is.na(exp_date) && exp_date <= start1, 1, 0),
          status   = ifelse(!is.na(d$RA_date[i]) && d$RA_date[i] <= end_all, 1, 0)
        )
      }
    }
    bind_rows(td_list) %>%
      mutate(
        start_num = as.numeric(start - center_date),
        stop_num  = as.numeric(stop  - center_date)
      ) %>%
      select(-center_date)
  }
  
  list(
    Main        = build_td(d_main) %>% left_join(cov_df, by = "eid"),
    Sensitivity = build_td(d_sens) %>% left_join(cov_df, by = "eid")
  )
}

# 批量构建
oral_exposure_list <- c("K02", "K04", "K05", "K12")
td_datasets_oral_to_DM <- lapply(oral_exposure_list, function(x)
  make_td_datasets_oral_to_DM(analysis_df, x))
names(td_datasets_oral_to_DM) <- oral_exposure_list
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++ ####
# >>>>> section 1.3 COX ####
# 协变量分组
model0_cov <- c()
model1_cov <- c("Age", "Sex", "Ethnicity")
model2_cov <- c(model1_cov, "Smoke", "Alcohol", "PA", "Education", "TDI_quantile")
model3_cov <- c(model2_cov, "BMI", "Hypertension","Diabetes")

run_models <- function(data) {
  cov_sets <- list(Model0=model0_cov, Model1=model1_cov, Model2=model2_cov, Model3=model3_cov)
  bind_rows(lapply(names(cov_sets), function(m){
    covars <- cov_sets[[m]]
    fml <- if (length(covars)) {
      as.formula(paste("Surv(start_num, stop_num, status) ~ exposure +", paste(covars, collapse=" + ")))
    } else {
      Surv(start_num, stop_num, status) ~ exposure
    }
    fit <- coxph(fml, data=data)
    hr  <- exp(coef(fit)["exposure"])
    ci  <- exp(confint(fit)["exposure", ])
    p   <- summary(fit)$coefficients["exposure", "Pr(>|z|)"]
    data.frame(Model=m, HR=hr, CI_lower=ci[1], CI_upper=ci[2], p_value=p)
  }))
}

cox_results_oral_to_DM <- lapply(oral_exposure_list, function(expname){
  list(
    Main        = run_models(td_datasets_oral_to_DM[[expname]]$Main),
    Sensitivity = run_models(td_datasets_oral_to_DM[[expname]]$Sensitivity)
  )
})
names(cox_results_oral_to_DM) <- oral_exposure_list

# 汇总表
icd10_labels <- c(K02="Dental caries",
                  K04="Pulp & periapical diseases",
                  K05="Gingivitis & periodontal diseases",
                  K12="Stomatitis & related lesions")

cox_summary_oral_to_DM <- bind_rows(lapply(names(cox_results_oral_to_DM), function(expname){
  bind_rows(
    cox_results_oral_to_DM[[expname]]$Main        %>% mutate(Exposure=expname, Outcome="Rheumatoid arthritis", Analysis="Main"),
    cox_results_oral_to_DM[[expname]]$Sensitivity %>% mutate(Exposure=expname, Outcome="Rheumatoid arthritis", Analysis="Sensitivity")
  )
})) %>%
  mutate(Disease = icd10_labels[Exposure]) %>%
  select(Disease, Exposure, Outcome, Analysis, Model, HR, CI_lower, CI_upper, p_value)

write.csv(cox_summary_oral_to_DM,
          "K:/UK biobank/Cohort data/Cox_Oral_disease_to_RA.csv",
          row.names = FALSE)
# +================================================+ ####
# +====Section 2. RA to Oral disease ==============+ ####
# +================================================+ ####
# >>>>> section 2.1 Combind ####
time_df          <- read_csv("Original data/Date_data.csv")
cov_df           <- read_csv("Original data/Covariate_data.csv")
oral_disease_df  <- read_csv("Original data/Oral_disease_data.csv")
oral_health_df   <- read_csv("Original data/Oral_health_data.csv")
RA_df            <- read_csv("Original data/Arthritis_data.csv")

cov_df <- cov_df %>%
  select(eid, Age, Sex, Ethnicity, Education, TDI_quantile, Smoke, Alcohol, PA, BMI,
         Diabetes, Hypertension)

oral_health_df <- oral_health_df %>%
  mutate(oral_health = ifelse(oral == "None of the above", 0, 1)) %>%
  select(eid, oral_health)

RA_df <- RA_df %>%
  mutate(rheumatoid_arthritis_date = ymd(rheumatoid_arthritis_date)) %>%
  select(eid,
         RA_date  = rheumatoid_arthritis_date,
         RA_define = rheumatoid_arthritis_define)

analysis_df <- time_df %>%
  left_join(oral_disease_df, by = "eid") %>%
  left_join(oral_health_df,  by = "eid") %>%
  left_join(RA_df,           by = "eid") %>%
  left_join(cov_df,          by = "eid")

# 口腔ICD列表 & 对应日期列
oral_list      <- c("K02","K04","K05","K12")
oral_date_cols <- paste0(oral_list, "_date")

# 可选：确保日期列为 Date
# analysis_df <- analysis_df %>% mutate(across(all_of(oral_date_cols), ~ ymd(.)))

# 关键：排除基线前存在任一口腔病者
analysis_df <- analysis_df %>%
  rowwise() %>%
  mutate(any_oral_prebaseline = any(c_across(all_of(oral_date_cols)) < center_date, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!any_oral_prebaseline) %>%
  select(-any_oral_prebaseline)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++ ####
# >>>>> section 2.2 Dataset ####
# 函数：生成 RA 暴露 → Oral disease 结局的 td 数据集（Main + Sensitivity）
make_td_datasets_RA_to_oral <- function(df, outcome_var) {
  out_date_col <- paste0(outcome_var, "_date")
  
  # 保留在基线前未发生该结局者
  d_main <- df %>%
    filter(is.na(.data[[out_date_col]]) | .data[[out_date_col]] >= center_date)
  
  # 敏感性：去掉“无DM但 oral_health=1”的个体（选择性，与你原逻辑一致）
  d_sens <- d_main %>%
    filter(!(is.na(RA_date) & oral_health == 1))
  
  build_td <- function(d) {
    td_list <- vector("list", length = nrow(d))
    for (i in seq_len(nrow(d))) {
      RA_date <- d$RA_date[i]
      start1  <- d$center_date[i]
      end_all <- min(c(d[[out_date_col]][i], d$death_date[i], d$censor_date[i]), na.rm = TRUE)
      if (is.infinite(end_all)) end_all <- d$censor_date[i]
      
      if (!is.na(RA_date) && RA_date <= start1) {
        # 基线前/当日已 RA：全程暴露
        td_list[[i]] <- data.frame(
          eid=d$eid[i], center_date=start1, start=start1, stop=end_all,
          exposure=1,
          status=ifelse(!is.na(d[[out_date_col]][i]) && d[[out_date_col]][i] <= end_all, 1, 0)
        )
      } else if (!is.na(RA_date) && RA_date > start1 && RA_date < end_all) {
        # 随访中新发 RA：拆段
        td_list[[i]] <- bind_rows(
          data.frame(eid=d$eid[i], center_date=start1, start=start1, stop=RA_date,
                     exposure=0, status=0),
          data.frame(eid=d$eid[i], center_date=start1, start=RA_date, stop=end_all,
                     exposure=1,
                     status=ifelse(!is.na(d[[out_date_col]][i]) && d[[out_date_col]][i] <= end_all, 1, 0))
        )
      } else {
        # 从未 RA 或 RA 在终点日/之后
        td_list[[i]] <- data.frame(
          eid=d$eid[i], center_date=start1, start=start1, stop=end_all,
          exposure=0,
          status=ifelse(!is.na(d[[out_date_col]][i]) && d[[out_date_col]][i] <= end_all, 1, 0)
        )
      }
    }
    bind_rows(td_list) %>%
      mutate(
        start_num = as.numeric(start - center_date),
        stop_num  = as.numeric(stop  - center_date)
      ) %>%
      select(-center_date)
  }
  
  list(
    Main        = build_td(d_main) %>% left_join(cov_df, by = "eid"),
    Sensitivity = build_td(d_sens) %>% left_join(cov_df, by = "eid")
  )
}

td_datasets_RA_to_oral <- lapply(oral_list, function(outc)
  make_td_datasets_RA_to_oral(analysis_df, outc))
names(td_datasets_RA_to_oral) <- oral_list

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++ ####
# >>>>> section 2.3 COX ####
model0_cov <- c()
model1_cov <- c("Age", "Sex", "Ethnicity")
model2_cov <- c(model1_cov, "Smoke", "Alcohol", "PA", "Education", "TDI_quantile")
model3_cov <- c(model2_cov, "BMI", "Hypertension","Diabetes")

run_models <- function(data) {
  cov_sets <- list(Model0=model0_cov, Model1=model1_cov, Model2=model2_cov, Model3=model3_cov)
  bind_rows(lapply(names(cov_sets), function(m){
    covars <- cov_sets[[m]]
    fml <- if (length(covars)) {
      as.formula(paste("Surv(start_num, stop_num, status) ~ exposure +", paste(covars, collapse=" + ")))
    } else {
      Surv(start_num, stop_num, status) ~ exposure
    }
    fit <- coxph(fml, data=data)
    hr  <- exp(coef(fit)["exposure"])
    ci  <- exp(confint(fit)["exposure", ])
    p   <- summary(fit)$coefficients["exposure", "Pr(>|z|)"]
    data.frame(Model=m, HR=hr, CI_lower=ci[1], CI_upper=ci[2], p_value=p)
  }))
}

cox_results_RA_to_oral <- lapply(oral_list, function(outc){
  list(
    Main        = run_models(td_datasets_RA_to_oral[[outc]]$Main),
    Sensitivity = run_models(td_datasets_RA_to_oral[[outc]]$Sensitivity)
  )
})
names(cox_results_RA_to_oral) <- oral_list

icd10_labels <- c(K02="Dental caries",
                  K04="Pulp & periapical diseases",
                  K05="Gingivitis & periodontal diseases",
                  K12="Stomatitis & related lesions")

cox_summary_RA_to_oral <- bind_rows(lapply(names(cox_results_RA_to_oral), function(outc){
  bind_rows(
    cox_results_RA_to_oral[[outc]]$Main        %>% mutate(Exposure="Rheumatoid arthritis", Outcome=outc, Analysis="Main"),
    cox_results_RA_to_oral[[outc]]$Sensitivity %>% mutate(Exposure="Rheumatoid arthritis", Outcome=outc, Analysis="Sensitivity")
  )
})) %>%
  mutate(Disease = icd10_labels[Outcome]) %>%
  select(Disease, Exposure, Outcome, Analysis, Model, HR, CI_lower, CI_upper, p_value)

write.csv(cox_summary_RA_to_oral,
          "K:/UK biobank/Cohort data/Cox_RA_to_Oral_disease.csv",
          row.names = FALSE)
# +================================================+ ####
# +====Section 3. Oral health to RA ==============+ ####
# +================================================+ ####
# >>>>> section 1 COMBINE ####
RA_df   <- read_csv("Original data/Arthritis_data.csv")   
oral_df <- read_csv("Original data/Oral_data.csv")
time_df <- read_csv("Original data/Date_data.csv")

RA_df <- RA_df %>%
  mutate(rheumatoid_arthritis_date = ymd(rheumatoid_arthritis_date)) %>%
  select(eid, rheumatoid_arthritis_date, rheumatoid_arthritis_define)

time_df <- time_df %>%
  mutate(center_date = ymd(center_date),
         death_date  = ymd(death_date),
         censor_date = ymd(censor_date))

df <- time_df %>%
  inner_join(oral_df, by = "eid") %>%
  inner_join(RA_df, by = "eid")

# 排除基线前 RA
exclude_DM <- df %>%
  filter(!is.na(rheumatoid_arthritis_date) & rheumatoid_arthritis_date < center_date)
write_csv(exclude_DM, "excluded_RA_before_baseline.csv")

df <- df %>%
  filter(is.na(rheumatoid_arthritis_date) | rheumatoid_arthritis_date >= center_date) %>%
  mutate(
    end_date        = pmin(rheumatoid_arthritis_date, death_date, censor_date, na.rm = TRUE),
    end_date        = if_else(is.infinite(as.numeric(end_date)), censor_date, end_date),
    rheumatoid_arthritis = ifelse(!is.na(rheumatoid_arthritis_date) & rheumatoid_arthritis_date == end_date, 1, 0),
    follow_up_days  = as.numeric(end_date - center_date),
    follow_up_years = follow_up_days / 365.25
  )

# 合并协变量
covar_df <- read_csv("K:/UK biobank/Original data/Covariate_data.csv") %>%
  select(eid, Age, Sex, Ethnicity, Education, TDI_quantile, Smoke, Alcohol, PA, BMI, Hypertension,Diabetes)

analysis_df <- df %>% left_join(covar_df, by = "eid")

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ####
# >>>>> section 2 COX ####
oral_vars <- c("Mouth_ulcers", "Painful_gums", "Bleeding_gums",
               "Loose_teeth", "Toothache", "Dentures", "Periodontal_disease")

model0_cov <- c()
model1_cov <- c("Age", "Sex", "Ethnicity")
model2_cov <- c(model1_cov, "Education", "TDI_quantile", "Smoke", "Alcohol", "PA")
model3_cov <- c(model2_cov, "BMI", "Hypertension","Diabetes")

cov_sets <- list(Model0=model0_cov, Model1=model1_cov, Model2=model2_cov, Model3=model3_cov)

cox_results <- list()

for (var in oral_vars) {
  for (m in names(cov_sets)) {
    covars <- cov_sets[[m]]
    fml <- if (length(covars)) {
      as.formula(paste0("Surv(follow_up_years, rheumatoid_arthritis) ~ ", var, " + ", paste(covars, collapse=" + ")))
    } else {
      as.formula(paste0("Surv(follow_up_years, rheumatoid_arthritis) ~ ", var))
    }
    fit <- coxph(fml, data = analysis_df)
    summ <- summary(fit)
    hr   <- exp(coef(fit)[var])
    ci   <- exp(confint(fit)[var, ])
    p    <- summ$coefficients[var, "Pr(>|z|)"]
    cox_results[[paste(var, m, sep = "_")]] <- data.frame(
      Exposure = var, Model = m, HR = hr, CI_lower = ci[1], CI_upper = ci[2], p_value = p
    )
  }
}

Cox_Oral_health_to_DM <- bind_rows(cox_results)

write_csv(Cox_Oral_health_to_DM,
          "K:/UK biobank/Cohort data/Cox_RA_oral_health_to_RA.csv")





