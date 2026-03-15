##########单因素cox分析############
library(survival)
library(dplyr)
library(readxl)

data <- read_excel("uni_cox.xlsx")

data$patient_id <- as.character(data$patient_id)
str(data)
vars <- c(
  "radiogenomic_score",
  "age",
  "lymph_node",
  "tumor_size_cm",
  "tumor_grade",
  "menopause_status",
  "vascular_invasion",
  "wavelet.LLH_firstorder_Kurtosis",
  "wavelet.HLL_ngtdm_Contrast",
  "wavelet.LLL_firstorder_RootMeanSquared"
)
results <- data.frame()

for (v in vars) {
  
  formula <- as.formula(
    paste("Surv(time_months, PFS_event) ~", v)
  )
  
  model <- coxph(formula, data = data)
  s <- summary(model)
  
  beta <- s$coefficients[,"coef"]
  HR <- s$coefficients[,"exp(coef)"]
  CI_low <- s$conf.int[,"lower .95"]
  CI_high <- s$conf.int[,"upper .95"]
  wald <- s$coefficients[,"z"]
  p <- s$coefficients[,"Pr(>|z|)"]
  
  results <- rbind(
    results,
    data.frame(
      Feature = v,
      beta = round(beta,3),
      HR_CI = paste0(
        round(HR,3),
        " (",
        round(CI_low,3),
        "–",
        round(CI_high,3),
        ")"
      ),
      Wald = round(wald,3),
      P = p
    )
  )
}
results
results <- results %>%
  mutate(
    P = ifelse(P < 0.001, paste0(format(P, scientific = TRUE)," ***"),
               ifelse(P < 0.01, paste0(round(P,3)," **"),
                      ifelse(P < 0.05, paste0(round(P,3)," *"),
                             round(P,3))))
  )
colnames(results) <- c(
  "Feature name",
  "β",
  "HR (95% CI)",
  "Wald",
  "P"
)
write.csv(results,
          "Univariate_cox_results.csv",
          row.names = FALSE)

##########多因素cox分析############
data$`wavelet.LLH_firstorder_Kurtosis` <- scale(data$`wavelet.LLH_firstorder_Kurtosis`)
data$`wavelet.LLL_firstorder_RootMeanSquared` <- scale(data$`wavelet.LLL_firstorder_RootMeanSquared`)
multi_model <- coxph(
  Surv(time_months, PFS_event) ~
    radiogenomic_score +
    `lymph_node` +
    `wavelet.LLH_firstorder_Kurtosis` +
    `wavelet.LLL_firstorder_RootMeanSquared`,
  data = data
)
summary(multi_model)
s <- summary(multi_model)

multi_results <- data.frame(
  `Feature name` = rownames(s$coefficients),
  `β` = round(s$coefficients[, "coef"], 3),
  `HR(95%CI)` = paste0(
    round(s$coefficients[, "exp(coef)"], 3),
    " (",
    round(s$conf.int[, "lower .95"], 3),
    "–",
    round(s$conf.int[, "upper .95"], 3),
    ")"
  ),
  `Wald` = round(s$coefficients[, "z"], 3),
  `P` = s$coefficients[, "Pr(>|z|)"]
)
multi_results

multi_results$P <- ifelse(
  multi_results$P < 0.001,
  paste0(format(multi_results$P, scientific = TRUE), " ***"),
  ifelse(
    multi_results$P < 0.01,
    paste0(round(multi_results$P, 3), " **"),
    ifelse(
      multi_results$P < 0.05,
      paste0(round(multi_results$P, 3), " *"),
      round(multi_results$P, 3)
    )
  )
)

multi_results
write.csv(multi_results,
          "multivariate_cox_results.csv",
          row.names = FALSE)

##########森林图############
library(survival)
library(ggplot2)

s <- summary(multi_model)

forest_data <- data.frame(
  Variable = rownames(s$coefficients),
  HR = s$coefficients[,2],
  lower = s$conf.int[,3],
  upper = s$conf.int[,4]
)

ggplot(forest_data, aes(x = HR, y = Variable)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  theme_classic() +
  xlab("Hazard Ratio") +
  ylab("")
ggsave("Forest_plot.pdf", width = 5.5, height = 4.5)
##########nomogram############
install.packages("rms")
library(rms)
data <- read_excel("uni_cox.xlsx")
dd <- datadist(data)
options(datadist = "dd")

cox_model <- cph(
  Surv(time_months, PFS_event) ~ 
    radiogenomic_score +
    wavelet.LLH_firstorder_Kurtosis,
  data = data,
  x = TRUE,
  y = TRUE,
  surv = TRUE
)

surv <- Survival(cox_model)

nom <- nomogram(
  cox_model,
  fun=list(function(x) surv(360,x),
           function(x) surv(1095,x),
           function(x) surv(1825,x)),
  funlabel=c("1-year PFS","3-year PFS","5-year PFS"),
  fun.at=list(
    c(0.95,0.97,0.99),
    c(0.85,0.90,0.95),
    c(0.70,0.80,0.90)
  )
)
pdf("nomogram_plot.pdf", width = 10, height = 6)
plot(nom)
dev.off()

##########DCA############
install.packages("rmda")
library(rmda)

# 线性预测值
data$lp <- predict(cox_model, type = "lp")
# 相对风险
data$risk <- exp(data$lp)
# 固定时间点生存概率
surv <- Survival(cox_model)
# 5-year PFS probability
data$pred_5y_pfs <- surv(1825, data$lp)
# 5-year event probability
data$pred_5y_event <- 1 - data$pred_5y_pfs
dca_model <- decision_curve(
  PFS_event ~ pred_5y_event,
  data = data,
  family = binomial(link = "logit"),
  thresholds = seq(0,1, by = 0.01),
  confidence.intervals = 0.95
)
pdf("DCA_curve.pdf", width = 5.5, height = 4.5)

plot_decision_curve(
  dca_model,
  curve.names = "Radiogenomic model",
  col = "#377EB8",
  lwd = 2,
  cost.benefit.axis = FALSE,
  standardize = TRUE,
  confidence.intervals = FALSE
)

dev.off()


