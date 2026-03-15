library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
########################一、RFS##########################
setwd("/Users/celeste/PhD/Rcode/Training/IHC4_project/KM_PFS")
##############Training_RFS################
head(Training_RFS)
colnames(Training_RFS)
# 构建Surv对象，注意 event=1 表示发生事件
Training_RFS$time_months <- Training_RFS$time_months / 30.44
surv_object <- Surv(time = Training_RFS$time_months, event = Training_RFS$PFS_event)

# 按 score 分组进行生存分析
fit <- survfit(surv_object ~ score, data = Training_RFS)

# 绘制KM曲线
km_plot <- ggsurvplot(
  fit,
  data = Training_RFS,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  legend.labs = c("Low score", "High score"),
  legend.title = NULL,
  xlab = "Time since surgery (months)",
  ylab = "Progression-Free Survival Probability",
  palette = c("#377EB8", "#E41A1C"),  # 蓝 + 红
  risk.table.height = 0.25,
  risk.table.title = "Number at risk",
  risk.table.fontsize = 3,
  font.main = c(14, "bold"),
  font.x = c(12),
  font.y = c(12),
  font.tickslab = c(10),
  legend = "top",
  ggtheme = theme_classic()
)
km_plot$plot <- km_plot$plot + 
  ggtitle("Training cohort")+
  theme(plot.title = element_text(hjust = 0.5))
png("KM_Training_RFS.png", width = 5.5, height = 4.5, units = "in", res = 600)
print(km_plot)  # 包含主图和风险表
dev.off()

pdf("KM_Training_RFS.pdf", width = 5.5, height = 4.5)
print(km_plot)
dev.off()

##############Test_RFS################
Test_RFS$time_months <- Test_RFS$time_months / 30.44
surv_object <- Surv(time = Test_RFS$time_months, event = Test_RFS$PFS_event)

fit <- survfit(surv_object ~ score, data = Test_RFS)

km_plot <- ggsurvplot(
  fit,
  data = Test_RFS,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  legend.labs = c("Low score", "High score"),
  legend.title = NULL,
  xlab = "Time since surgery (months)",
  ylab = "Progression-Free Survival Probability",
  palette = c("#377EB8", "#E41A1C"),
  risk.table.height = 0.25,
  risk.table.title = "Number at risk",
  risk.table.fontsize = 3,
  font.main = c(14, "bold"),
  font.x = c(12),
  font.y = c(12),
  font.tickslab = c(10),
  legend = "top",
  ggtheme = theme_classic()
)
km_plot$plot <- km_plot$plot + 
  ggtitle("Internal validation cohort")+
  theme(plot.title = element_text(hjust = 0.5))
png("KM_Test_RFS.png", width = 5.5, height = 4.5, units = "in", res = 600)
print(km_plot)
dev.off()

pdf("KM_Test_RFS.pdf", width = 5.5, height = 4.5)
print(km_plot)
dev.off()

##############WCH_RFS################
WCH_RFS$time_months <- WCH_RFS$time_months / 30.44
# 生存对象
surv_object <- Surv(time = WCH_RFS$time_months, event = WCH_RFS$PFS_event)

# 生存拟合
fit <- survfit(surv_object ~ score, data = WCH_RFS)

# 绘图
km_plot <- ggsurvplot(
  fit,
  data = WCH_RFS,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  legend.labs = c("Low score", "High score"),
  legend.title = NULL,
  xlab = "Time since surgery (months)",
  ylab = "Progression-Free Survival Probability",
  palette = c("#377EB8", "#E41A1C"),
  risk.table.height = 0.25,
  risk.table.title = "Number at risk",
  risk.table.fontsize = 3,
  font.main = c(14, "bold"),
  font.x = c(12),
  font.y = c(12),
  font.tickslab = c(10),
  legend = "top",
  ggtheme = theme_classic()
)
km_plot$plot <- km_plot$plot + 
  ggtitle("WCH Cohort")+
  theme(plot.title = element_text(hjust = 0.5))
# 保存图像
png("KM_WCH_RFS.png", width = 5.5, height = 4.5, units = "in", res = 600)
print(km_plot)
dev.off()

pdf("KM_WCH_RFS.pdf", width = 5.5, height = 4.5)
print(km_plot)
dev.off()
##############SYSUCC_RFS################
SYSUCC_RFS$time_months <- SYSUCC_RFS$time_months / 30.44
# 生存对象
surv_object <- Surv(time = SYSUCC_RFS$time_months, event = SYSUCC_RFS$PFS_event)

fit <- survfit(surv_object ~ score, data = SYSUCC_RFS)

km_plot <- ggsurvplot(
  fit,
  data = SYSUCC_RFS,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  legend.labs = c("Low score", "High score"),
  legend.title = NULL,
  xlab = "Time since surgery (months)",
  ylab = "Progression-Free Survival Probability",
  palette = c("#377EB8", "#E41A1C"),
  risk.table.height = 0.25,
  risk.table.title = "Number at risk",
  risk.table.fontsize = 3,
  font.main = c(14, "bold"),
  font.x = c(12),
  font.y = c(12),
  font.tickslab = c(10),
  legend = "top",
  ggtheme = theme_classic()
)
km_plot$plot <- km_plot$plot + 
  ggtitle("SYSUCC Cohort")+
  theme(plot.title = element_text(hjust = 0.5))
png("KM_SYSUCC_RFS.png", width = 5.5, height = 4.5, units = "in", res = 600)
print(km_plot)
dev.off()

pdf("KM_SYSUCC_RFS.pdf", width = 5.5, height = 4.5)
print(km_plot)
dev.off()

########################二、late DR##########################
graphics.off()
setwd("/Users/celeste/PhD/Rcode/Training/IHC4_project/KM_DR")
##############Training_lateDR################
# 构建生存对象 & 拟合KM模型
Training_lateDR$time_months <- Training_lateDR$time_months / 30.44
surv <- Surv(time = Training_lateDR$time_months, event = Training_lateDR$late_DR_event)
fit <- survfit(surv ~ score, data = Training_lateDR)
summary(fit)
# 绘图
dr_plot <- ggsurvplot(
  fit,
  data = Training_lateDR,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  legend.labs = c("Low score", "High score"),
  legend.title = NULL,
  xlab = "Time since 5-year landmark (months)",
  ylab = "Late distant recurrence-free probability",
  palette = c("#377EB8", "#E41A1C"),
  risk.table.height = 0.25,
  risk.table.title = "Number at risk",
  font.main = c(14, "bold"),
  font.x = c(12),
  font.y = c(12),
  font.tickslab = c(10),
  legend = "top",
  ggtheme = theme_classic()
)
dr_plot$plot <- dr_plot$plot + 
  ggtitle("Training cohort")+
  theme(plot.title = element_text(hjust = 0.5))

png("KM_Training_lateDR.png", width = 5.5, height = 4.5, units = "in", res = 600)
print(dr_plot)
dev.off()

pdf("KM_Training_lateDR.pdf", width = 5.5, height = 4.5)
print(dr_plot)
dev.off()
##############Test_lateDR################
# 构建生存对象 & 拟合KM模型
Test_lateDR$time_months <- Test_lateDR$time_months / 30.44
surv <- Surv(time = Test_lateDR$time_months, event = Test_lateDR$late_DR_event)
fit <- survfit(surv ~ score, data = Test_lateDR)
summary(fit)

# 绘图
dr_plot <- ggsurvplot(
  fit,
  data = Test_lateDR,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  legend.labs = c("Low score", "High score"),
  legend.title = NULL,
  xlab = "Time since 5-year landmark (months)",
  ylab = "Late distant recurrence-free probability",
  palette = c("#377EB8", "#E41A1C"),
  risk.table.height = 0.25,
  risk.table.title = "Number at risk",
  font.main = c(14, "bold"),
  font.x = c(12),
  font.y = c(12),
  font.tickslab = c(10),
  legend = "top",
  ggtheme = theme_classic()
)

dr_plot$plot <- dr_plot$plot + 
  ggtitle("Internal validation cohort") +
  theme(plot.title = element_text(hjust = 0.5))

png("KM_Test_lateDR.png", width = 5.5, height = 4.5, units = "in", res = 600)
print(dr_plot)
dev.off()

pdf("KM_Test_lateDR.pdf", width = 5.5, height = 4.5)
print(dr_plot)
dev.off()


##############WCH_lateDR################
# 构建生存对象 & 拟合KM模型
WCH_lateDR$time_months <- WCH_lateDR$time_months / 30.44
surv <- Surv(time = WCH_lateDR$time_months, event = WCH_lateDR$late_DR_event)
fit <- survfit(surv ~ score, data = WCH_lateDR)
summary(fit)

# 绘图
dr_plot <- ggsurvplot(
  fit,
  data = WCH_lateDR,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  legend.labs = c("Low score", "High score"),
  legend.title = NULL,
  xlab = "Time since 5-year landmark (months)",
  ylab = "Late distant recurrence-free probability",
  palette = c("#377EB8", "#E41A1C"),
  risk.table.height = 0.25,
  risk.table.title = "Number at risk",
  font.main = c(14, "bold"),
  font.x = c(12),
  font.y = c(12),
  font.tickslab = c(10),
  legend = "top",
  ggtheme = theme_classic()
)

dr_plot$plot <- dr_plot$plot + 
  ggtitle("WCH cohort") +
  theme(plot.title = element_text(hjust = 0.5))

png("KM_WCH_lateDR.png", width = 5.5, height = 4.5, units = "in", res = 600)
print(dr_plot)
dev.off()

pdf("KM_WCH_lateDR.pdf", width = 5.5, height = 4.5)
print(dr_plot)
dev.off()


##############SYSUCC_lateDR################
# 构建生存对象 & 拟合KM模型
SYSUCC_lateDR$time_months <- SYSUCC_lateDR$time_months / 30.44
surv <- Surv(time = SYSUCC_lateDR$time_months, event = SYSUCC_lateDR$late_DR_event)
fit <- survfit(surv ~ score, data = SYSUCC_lateDR)
summary(fit)

# 绘图
dr_plot <- ggsurvplot(
  fit,
  data = SYSUCC_lateDR,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  legend.labs = c("Low score", "High score"),
  legend.title = NULL,
  xlab = "Time since 5-year landmark (months)",
  ylab = "Late distant recurrence-free probability",
  palette = c("#377EB8", "#E41A1C"),
  risk.table.height = 0.25,
  risk.table.title = "Number at risk",
  font.main = c(14, "bold"),
  font.x = c(12),
  font.y = c(12),
  font.tickslab = c(10),
  legend = "top",
  ggtheme = theme_classic()
)

dr_plot$plot <- dr_plot$plot + 
  ggtitle("SYSUCC cohort") +
  theme(plot.title = element_text(hjust = 0.5))

png("KM_SYSUCC_lateDR.png", width = 5.5, height = 4.5, units = "in", res = 600)
print(dr_plot)
dev.off()

pdf("KM_SYSUCC_lateDR.pdf", width = 5.5, height = 4.5)
print(dr_plot)
dev.off()

