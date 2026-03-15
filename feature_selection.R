library(glmnet)
library(readxl)

data <- read_excel("Training_radiomics.xlsx")
#######################LASSO_logstic##########################
# 结局变量
y <- as.matrix(data$label)

# radiomics 特征
x <- as.matrix(data[,4:ncol(data)])

# 标准化
x <- scale(x)

# LASSO
lasso_model <- glmnet(
  x,
  y,
  family = "binomial",
  alpha = 1
)

# 路径图
pdf("LASSO_path.pdf", width = 5.5, height = 4.5)
plot(lasso_model, xvar = "lambda", lwd = 2)
dev.off()
set.seed(123)
# 交叉验证
cv_model <- cv.glmnet(
  x,
  y,
  family = "binomial",
  alpha = 1,
  nfolds = 5
)

pdf("LASSO_CV.pdf", width = 5.5, height = 4.5)
plot(cv_model)
dev.off()

# lambda
lambda_min <- cv_model$lambda.min
lambda_1se <- cv_model$lambda.1se

lambda_min
lambda_1se

# 提取特征
coef_lasso <- coef(cv_model, s = "lambda.min")
coef_lasso
coef_table <- as.data.frame(as.matrix(coef_lasso))
coef_table$feature <- rownames(coef_table)

colnames(coef_table)[1] <- "coefficient"

selected_features <- coef_table[coef_table$coefficient != 0, ]
selected_features
nrow(selected_features)
write.csv(
  selected_features,
  "LASSO_selected_features.csv",
  row.names = FALSE
)

#######################LASSO_cox##########################
library(glmnet)#加载glmnet包
FUSCC_Training_LASSO$time_months <- FUSCC_Training_LASSO$time_months / 30.44
attach(FUSCC_Training_LASSO)
View(FUSCC_Training_LASSO)
colnames(FUSCC_Training_LASSO[,1:45])#查看前44列的列名（根据自己数据调整）
y <- Surv(
  FUSCC_Training_LASSO$time_months,
  FUSCC_Training_LASSO$PFS_event
) # 提取第2列作为结局（建议放在第一列）
x <- as.matrix(FUSCC_Training_LASSO[, 4:45])  # 第3至第44列为自变量
#后边的代码除了s值基本不需更改

lasso_model <- glmnet(x, y, family = "cox",
                      alpha = 1) # 表示采用L1正则化，即Lasso回归。
max(lasso_model$lambda)
print(lasso_model) 
#绘制LASSO图
pdf("LASSO_path.pdf", width = 6, height = 5)
plot(lasso_model,
     xvar = "lambda",
     lwd = 2)
dev.off()
#交叉验证并绘制可视化结果
set.seed(123)
cv_model <- cv.glmnet(x, y, family = "cox",alpha = 1,nfolds = 5)
pdf("LASSO_CV.pdf", width = 6, height = 5)
plot(cv_model)
dev.off()
#根据交叉验证结果，选择lambda值，lambda.min或lambda.1se。
lambda_min <- cv_model$lambda.min
lambda_min
lambda_1se <- cv_model$lambda.1se
lambda_1se
#s为Lambda大小，Lambda越大表示模型的正则化强度越大，选择的自变量也越少。
#这里选择的是刚刚得到的lambda_1se的值
coef_lasso <- coef(cv_model, s = "lambda.min")
coef_lasso
#结果显示后边带有数值的变量为筛选得到的变量
coef_table <- as.data.frame(as.matrix(coef_lasso))
coef_table$feature <- rownames(coef_table)
colnames(coef_table)[1] <- "coefficient"
selected_features <- coef_table[coef_table$coefficient != 0, ]
write.csv(selected_features,
          "LASSO_selected_features.csv",
          row.names = FALSE)
write.csv(coef_table,
          "LASSO_all_coefficients.csv",
          row.names = FALSE)

