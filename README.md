# =====================================================
# SFMS‑IHI 完整 R 代码示例
# 包含固定权重版 (FW) 和数据驱动版 (DDW) 计算
# 三层 HLM 估计 β 系数和 ICC
# 固定 SMEImax = 25.0 （跨队列可比性）
# =====================================================

# 安装并加载必要的包
if (!require("lme4")) install.packages("lme4")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
library(lme4)
library(dplyr)
library(ggplot2)

# 设置随机种子
set.seed(2026)

# =====================================================
# 固定常数：SMEImax (理论/经验上界)
# =====================================================
SMEImax <- 25.0
cat("固定归一化常数 SMEImax =", SMEImax, "\n")

# =====================================================
# 1. 模拟嵌套数据结构
# =====================================================
n_families <- 200          # 家庭数（Level 2）
n_members <- 5             # 每家庭人数
n_time <- 4                # 时间点
N <- n_families * n_members * n_time

# 层级标识
family_id <- rep(1:n_families, each = n_members * n_time)
region_id <- rep(sample(1:10, n_families, replace = TRUE), each = n_members * n_time)
time <- rep(1:n_time, times = n_families * n_members)

# =====================================================
# 2. 模拟暴露变量（与手册一致）
# =====================================================

# 2.1 社会暴露原始指数 SMEI_raw (范围理论上可达 -理论下限 到 +理论上限)
#     实际模拟使用截断正态，使得大多数 |SMEI_raw| < 25.0，少数可超过（用于演示截断）
SMEI_raw <- rnorm(N, mean = 0, sd = 8)   # sd=8 使得多数在 [-20,20] 区间
# 为保证与固定常数兼容，计算归一化 S_norm
S_norm <- SMEI_raw / SMEImax
# 将超出 [-1,1] 的值截断并标记
if(any(abs(S_norm) > 1)) {
  cat("警告：存在 |S_norm| > 1 的值，已截断至 ±1，并标记为 supra-maximal\n")
  S_norm_orig <- S_norm
  S_norm <- pmax(-1, pmin(1, S_norm))
  supra_flag <- abs(S_norm_orig) > 1
} else {
  supra_flag <- rep(FALSE, N)
}
# 转换为 S_risk (0-1，越高风险越大) 用于 HLM
S_risk <- (S_norm + 1) / 2

# 2.2 环境暴露 E (0-1)
E_diet <- 1 - rbeta(N, 2, 2)        # 模拟膳食不利度，可用 HEI 代替
E_density <- runif(N, 0.2, 0.8)
E <- (E_diet + E_density) / 2

# 2.3 遗传风险 G (0-1)
G <- rbeta(N, 3, 3)

# 2.4 抗生素暴露 A (0-1)
A <- rbeta(N, 1.5, 4)

# 2.5 饮食共享 DS (0-1)
DS <- rbeta(N, 2, 4)

# 2.6 垂直传播 V (仅对母乳婴儿非零，本模拟中设为大部分0)
V <- sample(c(rep(0, N-100), runif(100, 0, 0.5)))

# =====================================================
# 3. 模拟 GHI (Gut Health Index)
# 真实模型：GHI 受暴露变量及随机效应影响
# =====================================================
# 设定真实 β 系数（固定效应）
beta_E <- -0.15
beta_G <- -0.08
beta_A <- -0.30
beta_DS <- -0.10
beta_V <- -0.05
beta_S <- -0.20

# 随机效应方差 (真实值)
var_family_true <- 0.12
var_region_true <- 0.03
var_resid_true <- 0.15

# 生成随机效应
u_family <- rep(rnorm(n_families, 0, sqrt(var_family_true)), each = n_members * n_time)
v_region <- rep(rnorm(10, 0, sqrt(var_region_true)), each = N / 10)
epsilon <- rnorm(N, 0, sqrt(var_resid_true))

# 线性预测值
linear_pred <- 0.5 + 
               beta_E * E + beta_G * G + beta_A * A + 
               beta_DS * DS + beta_V * V + beta_S * S_risk +
               u_family + v_region + epsilon

# 转换为 GHI (0-1 之间，逻辑斯蒂变换)
GHI <- plogis(linear_pred)

# 构建数据框
df <- data.frame(GHI, E, G, A, DS, V, S_norm, S_risk, SMEI_raw,
                 family_id, region_id, time, supra_flag)

# =====================================================
# 4. 固定权重版 IHI (FW)
# =====================================================
alpha <- 0.2
beta_param <- 2.0

# 计算 BaseExposure (非社会维度平均)
BaseExposure <- (E + G + A + DS + V) / 5

# 计算 ExposureScore_FW
ExposureScore_FW <- BaseExposure + alpha * S_norm

# 计算 IHI_FW
IHI_FW <- GHI * exp(-beta_param * ExposureScore_FW) * 100

# 健康分级
grade_fw <- ifelse(IHI_FW >= 80, "Excellent",
            ifelse(IHI_FW >= 60, "Good",
            ifelse(IHI_FW >= 40, "Moderate",
            ifelse(IHI_FW >= 20, "Poor", "Critical"))))

df$IHI_FW <- IHI_FW
df$Grade_FW <- grade_fw

# =====================================================
# 5. 三层 HLM 估计 β 系数和 ICC (数据驱动权重版预备)
# =====================================================

# 拟合条件模型 (固定效应 + 随机截距)
hlm_model <- lmer(GHI ~ E + G + A + DS + V + S_risk + 
                  (1 | family_id) + (1 | region_id), data = df)

# 查看固定效应系数
summary(hlm_model)

# 提取固定效应系数 (未标准化)
beta_raw <- fixef(hlm_model)[c("E", "G", "A", "DS", "V", "S_risk")]

# 提取随机效应方差
var_comp <- as.data.frame(VarCorr(hlm_model))
var_family <- var_comp[var_comp$grp == "family_id", "vcov"]
var_region <- var_comp[var_comp$grp == "region_id", "vcov"]
var_resid <- attr(VarCorr(hlm_model), "sc")^2

# 计算 ICC
total_var <- var_family + var_region + var_resid
ICC_family <- var_family / total_var
ICC_region <- var_region / total_var
cat("Family ICC =", round(ICC_family, 3), "\n")
cat("Region ICC =", round(ICC_region, 3), "\n")

# =====================================================
# 6. 数据驱动权重 w_k 计算
# =====================================================
# 标准化系数 β*_k = β_k * sd(X_k) / sd(Y)
sd_x <- apply(df[, c("E", "G", "A", "DS", "V", "S_risk")], 2, sd)
sd_y <- sd(df$GHI)
beta_std <- beta_raw * sd_x / sd_y
w <- abs(beta_std) / sum(abs(beta_std))
names(w) <- c("E", "G", "A", "DS", "V", "S_risk")
print("Data-driven weights (based on fixed SMEImax normalisation):")
print(w)

# =====================================================
# 7. 数据驱动版 IHI (DDW)
# =====================================================
ExposureScore_DDW <- as.matrix(df[, c("E", "G", "A", "DS", "V", "S_risk")]) %*% w
IHI_DDW <- df$GHI * exp(-ExposureScore_DDW) * 100

df$IHI_DDW <- IHI_DDW

# =====================================================
# 8. 两种 IHI 的 Spearman 相关性
# =====================================================
cor_rank <- cor(df$IHI_FW, df$IHI_DDW, method = "spearman")
cat("Spearman correlation between IHI_FW and IHI_DDW:", round(cor_rank, 3), "\n")

# =====================================================
# 9. 可视化：两种 IHI 的箱线图对比
# =====================================================
df_long <- data.frame(
  IHI = c(df$IHI_FW, df$IHI_DDW),
  Version = rep(c("FW", "DDW"), each = N)
)
ggplot(df_long, aes(x = Version, y = IHI, fill = Version)) +
  geom_boxplot() +
  labs(title = "Comparison of IHI_FW and IHI_DDW",
       y = "SFMS-IHI Score (0-100)") +
  theme_minimal()

# =====================================================
# 10. 输出部分数据示例
# =====================================================
cat("\n前10行数据示例 (包含 SMEI_raw, S_norm, IHI_FW, IHI_DDW, Grade_FW):\n")
print(head(df[, c("GHI", "SMEI_raw", "S_norm", "IHI_FW", "IHI_DDW", "Grade_FW")], 10))

# 若有 supra-maximal 暴露，列出比例
if(any(supra_flag)) {
  cat("\n注意：", sum(supra_flag), "个观测的 |SMEI_raw| >", SMEImax, 
      "，S_norm 已被截断至 ±1\n")
}
