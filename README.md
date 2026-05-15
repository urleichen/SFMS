# SFMS
1.0
# =====================================================================
# SFMS‑IHI 完整 R 代码
# 包含：ENBI / GHI 计算 | 暴露维度 | 固定权重版 IHI (FW)
#       三层 HLM | 数据驱动版 IHI (DDW) | ICC | 权重 w | 可视化
# 基于：Corral Lopez et al. (Science 2026); SFMS 框架
# =====================================================================

# 加载必要的包
if (!require("SpiecEasi")) install.packages("SpiecEasi")
if (!require("lme4")) install.packages("lme4")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
library(SpiecEasi)
library(lme4)
library(dplyr)
library(ggplot2)

set.seed(2026)

# ==================== 1. 模拟宏基因组数据 ====================
# 生成物种相对丰度表 (100样本, 80物种)
n_samples <- 100
n_species <- 80

rdirichlet <- function(n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l * n, shape = alpha), ncol = l, byrow = TRUE)
  return(x / rowSums(x))
}

healthy_abund <- rdirichlet(40, rep(1, n_species))
disease_abund <- rdirichlet(60, rep(1, n_species))
abund <- rbind(healthy_abund, disease_abund)
rownames(abund) <- paste0("S", 1:n_samples)
colnames(abund) <- paste0("sp", 1:n_species)
group <- factor(c(rep("Healthy", 40), rep("Disease", 60)))

# 转换为整数计数 (SpiecEasi要求)
abund_counts <- round(abund * 1e6)

# ==================== 2. ENBI 和 GHI 计算 ====================
# 网络推断 (SPIEC-EASI)
se <- spiec.easi(abund_counts, method = 'mb', lambda.min.ratio = 1e-2, nlambda = 20)
cov_mat <- getOptCov(se)
cor_mat <- cov2cor(cov_mat)
W <- cor_mat
W[abs(W) < 0.1] <- 0
W <- sign(W)
diag(W) <- 0

# 计算每个样本的 ρ
calc_rho <- function(abund_vec, W) {
  n <- length(abund_vec)
  num <- den <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (W[i,j] != 0) {
        prod_ab <- abund_vec[i] * abund_vec[j]
        num <- num + W[i,j] * prod_ab
        den <- den + abs(W[i,j]) * prod_ab
      }
    }
  }
  if (den == 0) return(0) else return(num / den)
}
rho <- apply(abund, 1, calc_rho, W = W)
names(rho) <- rownames(abund)

median_rho_healthy <- median(rho[group == "Healthy"])
ENBI <- rho - median_rho_healthy
GHI <- 1 / (1 + exp(3 * ENBI))

# ==================== 3. 模拟暴露变量（个体层级） ====================
n_total <- n_samples
S_norm <- rnorm(n_total, 0, 0.3)                 # 净社会暴露 (-1~1)
S_risk <- (S_norm + 1) / 2                       # 映射到 [0,1] 风险方向
E_diet <- 1 - rbeta(n_total, 2, 2)
E_density <- runif(n_total, 0.2, 0.8)
E <- (E_diet + E_density) / 2
G <- rbeta(n_total, 3, 3)
A <- rbeta(n_total, 1.5, 4)
DS <- rbeta(n_total, 2, 4)
V <- sample(c(rep(0, n_total-20), runif(20, 0, 0.5)))   # 仅少数非零

# 生成家庭和区域嵌套结构 (用于 HLM，需足够样本)
n_families <- 200
n_members <- 5
n_time <- 4
N_nested <- n_families * n_members * n_time
family_id <- rep(1:n_families, each = n_members * n_time)
region_id <- rep(sample(1:10, n_families, replace = TRUE), each = n_members * n_time)

# 为嵌套数据重复生成暴露变量 (使每个个体有独立值)
E_nested <- runif(N_nested, 0, 1)
G_nested <- runif(N_nested, 0, 1)
A_nested <- runif(N_nested, 0, 1)
DS_nested <- runif(N_nested, 0, 1)
V_nested <- sample(c(rep(0, N_nested-100), runif(100, 0, 0.5)))
S_norm_nested <- rnorm(N_nested, 0, 0.3)
S_risk_nested <- (S_norm_nested + 1) / 2

# 模拟嵌套数据的 GHI (受暴露和随机效应影响)
beta_E <- -0.15; beta_G <- -0.08; beta_A <- -0.30
beta_DS <- -0.10; beta_V <- -0.05; beta_S <- -0.20
u_family <- rep(rnorm(n_families, 0, sqrt(0.12)), each = n_members * n_time)
v_region <- rep(rnorm(10, 0, sqrt(0.03)), each = N_nested / 10)
epsilon <- rnorm(N_nested, 0, sqrt(0.15))
linear_pred <- 0.5 + beta_E*E_nested + beta_G*G_nested + beta_A*A_nested +
               beta_DS*DS_nested + beta_V*V_nested + beta_S*S_risk_nested +
               u_family + v_region + epsilon
GHI_nested <- plogis(linear_pred)

df_nested <- data.frame(GHI = GHI_nested, E = E_nested, G = G_nested, A = A_nested,
                        DS = DS_nested, V = V_nested, S_risk = S_risk_nested,
                        family_id = family_id, region_id = region_id)

# ==================== 4. 固定权重版 IHI (FW) ====================
alpha_fw <- 0.2
beta_fw <- 2.0
BaseExposure <- (E_nested + G_nested + A_nested + DS_nested + V_nested) / 5
ExposureScore_FW <- BaseExposure + alpha_fw * S_norm_nested
IHI_FW <- GHI_nested * exp(-beta_fw * ExposureScore_FW) * 100
grade_fw <- ifelse(IHI_FW >= 80, "Excellent",
            ifelse(IHI_FW >= 60, "Good",
            ifelse(IHI_FW >= 40, "Moderate",
            ifelse(IHI_FW >= 20, "Poor", "Critical"))))
df_nested$IHI_FW <- IHI_FW
df_nested$Grade_FW <- grade_fw

# ==================== 5. 三层 HLM (估计 β 系数和 ICC) ====================
hlm_model <- lmer(GHI ~ E + G + A + DS + V + S_risk + (1 | family_id) + (1 | region_id),
                  data = df_nested)
summary(hlm_model)

# 提取固定效应系数 (未标准化)
beta_raw <- fixef(hlm_model)[c("E", "G", "A", "DS", "V", "S_risk")]

# 随机效应方差
var_comp <- as.data.frame(VarCorr(hlm_model))
var_family <- var_comp[var_comp$grp == "family_id", "vcov"]
var_region <- var_comp[var_comp$grp == "region_id", "vcov"]
var_resid <- attr(VarCorr(hlm_model), "sc")^2
total_var <- var_family + var_region + var_resid
ICC_family <- var_family / total_var
ICC_region <- var_region / total_var
cat("Family ICC =", round(ICC_family, 3), "\n")
cat("Region ICC =", round(ICC_region, 3), "\n")

# ==================== 6. 数据驱动版 IHI (DDW) 权重 ====================
sd_x <- apply(df_nested[, c("E", "G", "A", "DS", "V", "S_risk")], 2, sd)
sd_y <- sd(df_nested$GHI)
beta_std <- beta_raw * sd_x / sd_y
w <- abs(beta_std) / sum(abs(beta_std))
names(w) <- c("E", "G", "A", "DS", "V", "S_risk")
print("Data-driven weights:")
print(round(w, 3))

# ==================== 7. 计算 IHI_DDW ====================
ExposureScore_DDW <- as.matrix(df_nested[, c("E", "G", "A", "DS", "V", "S_risk")]) %*% w
IHI_DDW <- df_nested$GHI * exp(-ExposureScore_DDW) * 100
df_nested$IHI_DDW <- IHI_DDW

# ==================== 8. 两种 IHI 的相关性 ====================
cor_rank <- cor(df_nested$IHI_FW, df_nested$IHI_DDW, method = "spearman")
cat("Spearman correlation between IHI_FW and IHI_DDW:", round(cor_rank, 3), "\n")

# ==================== 9. 箱线图对比 ====================
df_plot <- data.frame(
  IHI = c(df_nested$IHI_FW, df_nested$IHI_DDW),
  Version = rep(c("FW", "DDW"), each = N_nested)
)
ggplot(df_plot, aes(x = Version, y = IHI, fill = Version)) +
  geom_boxplot() +
  labs(title = "Comparison of IHI_FW and IHI_DDW", y = "SFMS‑IHI (0-100)") +
  theme_minimal()

# ==================== 10. 输出结果示例 ====================
head(df_nested[, c("GHI", "IHI_FW", "IHI_DDW", "Grade_FW")], 10)
