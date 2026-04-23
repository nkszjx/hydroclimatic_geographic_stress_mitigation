
# The effect of urbanization on the impacts of geographic stresses on UGS coverage

rm(list=ls())  
library(parallel)
library(lmtest)
library(DescTools)
library(foreign)
library(Matrix)
library(lfe)  
library(magrittr)
library(margins)
library(naniar)
library(dplyr)
library(plotly)
library(zoo)
library(readxl)
library(mice)
library(rio)
library(orca)
library(DMwR2)
library(car)
library(AER)
library(lme4)
library(ggplot2)
library(brms)
library(mgcv)
library(future)
library(plm)
library(dplyr)
library(data.table)
library(fixest)
library(clubSandwich) 
library(ggpubr)
library(marginaleffects)
library(webshot)
library(htmlwidgets)
library(scales)
library(scatterplot3d)
library(extrafont)
library(showtext)
library(future)
library(future.apply)
library(stargazer)
library(data.table)
library(corrplot)
library(stringr)
library(patchwork)



################################################### 城市内外交互 ################################################## 
merged_dt <- fread("E:/2025/UrbanGreenSpaceWorks/Rcodes/dataset_2000_2020.csv")
dim(merged_dt) 

# 将纬度取绝对值
merged_dt$Latitude <- abs(merged_dt$Latitude)

#  VIF
lm_model <- lm(FVC ~  Precipitation + WaterDeficit + Seasonalityindex + AverageTemperature + TemperatureRange +  
                 WindSpeed + SoilMoisture +  SoilPH + 
                 Elevation + Latitude +   Longitude + 
                 GDP_per_capita_PPP + HDI + Population + 
                 ImperviousSurface + HumanSettlement + CityArea, 
               data = merged_dt)
vif_values <- vif(lm_model)
print(vif_values)

# 二次项
merged_dt$WaterDeficit_2 <- merged_dt$WaterDeficit * merged_dt$WaterDeficit
merged_dt$Seasonalityindex_2 <- merged_dt$Seasonalityindex * merged_dt$Seasonalityindex
merged_dt$Elevation_2 <- merged_dt$Elevation * merged_dt$Elevation
merged_dt$Latitude_2 <- merged_dt$Latitude * merged_dt$Latitude


################################################### elevation ################################################## 
cols_to_center <- c('WaterDeficit', 'Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_Elevation <- copy(merged_dt)
merged_dt_Elevation <- merged_dt_Elevation[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                           .SDcols = cols_to_center]
merged_dt_Elevation$urban <- factor(merged_dt_Elevation$urban, levels = c("out", "in"))

# linear regression
# Table 2 -- model 1
model_urban <- felm(FVC ~ 1 + WaterDeficit +
                      Seasonalityindex +
                      AverageTemperature +
                      TemperatureRange +
                      WindSpeed +
                      Elevation +
                      # Elevation_2 +
                      Latitude +
                      Precipitation +
                      SoilMoisture +
                      SoilPH + 
                      Longitude + 
                      
                      Elevation:urban +
                      # Elevation_2:urban +
                      
                      GDP_per_capita_PPP + HDI + Population +
                      ImperviousSurface + HumanSettlement + CityArea
                    |  year + CountryID, data = merged_dt_Elevation)

summary(model_urban)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_out.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

# nonlinear regression
# Table 2 -- model 3
model_urban <- felm(FVC ~ 1 + WaterDeficit +
                      Seasonalityindex +
                      AverageTemperature +
                      TemperatureRange +
                      WindSpeed +
                      Elevation +
                      Elevation_2 +
                      Latitude +
                      Precipitation +
                      SoilMoisture +
                      SoilPH + 
                      Longitude + 
                      
                      Elevation:urban +
                      Elevation_2:urban +
                      
                      GDP_per_capita_PPP + HDI + Population +
                      ImperviousSurface + HumanSettlement + CityArea
                    |  year + CountryID, data = merged_dt_Elevation)

summary(model_urban)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_out2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

######################################### 绘制效应图与边际效应图 #######################################
coeffs <- coef(model_urban)
beta1 <- coeffs["Elevation"]
beta2 <- coeffs["Elevation_2"]
beta1_city <- coeffs["Elevation:urbanin"]
beta2_city <- coeffs["Elevation_2:urbanin"]
beta1_total <- beta1 + beta1_city
beta2_total <- beta2 + beta2_city

P_min <- quantile(merged_dt_Elevation$Elevation, 0.001, na.rm = TRUE)
P_max <- quantile(merged_dt_Elevation$Elevation, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)

# 计算效应值
effect_non_urban <- beta1 * P_values + beta2 * (P_values^2)
effect_urban <- beta1_total * P_values + beta2_total * (P_values^2)

vcov_matrix <- vcovHC(model_urban, type = "HC1", cluster = "group", group = merged_dt_Elevation$CountryID)

cov_vars_non_urban <- c("Elevation", "Elevation_2")
vcov_non_urban <- vcov_matrix[cov_vars_non_urban, cov_vars_non_urban]
d_beta1 <- P_values
d_beta2 <- P_values^2
var_effect_non_urban <- (d_beta1^2) * vcov_non_urban["Elevation", "Elevation"] +
  (d_beta2^2) * vcov_non_urban["Elevation_2", "Elevation_2"] +
  2 * d_beta1 * d_beta2 * vcov_non_urban["Elevation", "Elevation_2"]
se_effect_non_urban <- sqrt(var_effect_non_urban)
ci_effect_non_urban_lower <- effect_non_urban - 1.96 * se_effect_non_urban
ci_effect_non_urban_upper <- effect_non_urban + 1.96 * se_effect_non_urban

cov_vars_urban <- c("Elevation", "Elevation_2", "Elevation:urbanin", "Elevation_2:urbanin")
vcov_urban <- vcov_matrix[cov_vars_urban, cov_vars_urban]
d_beta1_u <- P_values
d_beta2_u <- P_values^2
var_effect_urban <- (d_beta1_u^2) * vcov_urban["Elevation", "Elevation"] +
  (d_beta2_u^2) * vcov_urban["Elevation_2", "Elevation_2"] +
  (d_beta1_u^2) * vcov_urban["Elevation:urbanin", "Elevation:urbanin"] +
  (d_beta2_u^2) * vcov_urban["Elevation_2:urbanin", "Elevation_2:urbanin"] +
  2 * d_beta1_u^2 * vcov_urban["Elevation", "Elevation:urbanin"] +
  2 * d_beta2_u^2 * vcov_urban["Elevation_2", "Elevation_2:urbanin"] +
  2 * d_beta1_u * d_beta2_u * vcov_urban["Elevation", "Elevation_2"] +
  2 * d_beta1_u * d_beta2_u * vcov_urban["Elevation", "Elevation_2:urbanin"] +
  2 * d_beta1_u * d_beta2_u * vcov_urban["Elevation:urbanin", "Elevation_2"] +
  2 * d_beta1_u * d_beta2_u * vcov_urban["Elevation:urbanin", "Elevation_2:urbanin"]
se_effect_urban <- sqrt(var_effect_urban)
ci_effect_urban_lower <- effect_urban - 1.96 * se_effect_urban
ci_effect_urban_upper <- effect_urban + 1.96 * se_effect_urban

# 计算边际效应
marginal_non_urban <- beta1 + 2 * beta2 * P_values
marginal_urban <- beta1_total + 2 * beta2_total * P_values

d_marg_beta1 <- 1
d_marg_beta2 <- 2 * P_values
var_marginal_non_urban <- (d_marg_beta1^2) * vcov_non_urban["Elevation", "Elevation"] +
  (d_marg_beta2^2) * vcov_non_urban["Elevation_2", "Elevation_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_non_urban["Elevation", "Elevation_2"]
se_marginal_non_urban <- sqrt(var_marginal_non_urban)
ci_marginal_non_urban_lower <- marginal_non_urban - 1.96 * se_marginal_non_urban
ci_marginal_non_urban_upper <- marginal_non_urban + 1.96 * se_marginal_non_urban

d_marg_beta1_u <- 1
d_marg_beta2_u <- 2 * P_values
var_marginal_urban <- (d_marg_beta1_u^2) * vcov_urban["Elevation", "Elevation"] +
  (d_marg_beta2_u^2) * vcov_urban["Elevation_2", "Elevation_2"] +
  (d_marg_beta1_u^2) * vcov_urban["Elevation:urbanin", "Elevation:urbanin"] +
  (d_marg_beta2_u^2) * vcov_urban["Elevation_2:urbanin", "Elevation_2:urbanin"] +
  2 * d_marg_beta1_u^2 * vcov_urban["Elevation", "Elevation:urbanin"] +
  2 * d_marg_beta2_u^2 * vcov_urban["Elevation_2", "Elevation_2:urbanin"] +
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Elevation", "Elevation_2"] +
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Elevation", "Elevation_2:urbanin"] +
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Elevation:urbanin", "Elevation_2"] +
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Elevation:urbanin", "Elevation_2:urbanin"]
se_marginal_urban <- sqrt(var_marginal_urban)
ci_marginal_urban_lower <- marginal_urban - 1.96 * se_marginal_urban
ci_marginal_urban_upper <- marginal_urban + 1.96 * se_marginal_urban

effect_df <- rbind(
  data.frame(
    Elevation = P_values,
    Region = "Non-urban",
    Effect = effect_non_urban,
    SE = se_effect_non_urban,
    CI_Lower = ci_effect_non_urban_lower,
    CI_Upper = ci_effect_non_urban_upper
  ),
  data.frame(
    Elevation = P_values,
    Region = "Urban",
    Effect = effect_urban,
    SE = se_effect_urban,
    CI_Lower = ci_effect_urban_lower,
    CI_Upper = ci_effect_urban_upper
  )
)

marginal_df <- rbind(
  data.frame(
    Elevation = P_values,
    Region = "Non-urban",
    Marginal_Effect = marginal_non_urban,
    SE = se_marginal_non_urban,
    CI_Lower = ci_marginal_non_urban_lower,
    CI_Upper = ci_marginal_non_urban_upper
  ),
  data.frame(
    Elevation = P_values,
    Region = "Urban",
    Marginal_Effect = marginal_urban,
    SE = se_marginal_urban,
    CI_Lower = ci_marginal_urban_lower,
    CI_Upper = ci_marginal_urban_upper
  )
)

# 绘制效应值图
# Figure 3a
urban_colors <- c(
  "Urban" = "#FF7A3C", 
  "Non-urban" = "#00A850"
)

urban_colors_hist <- c(
  "in" = "#FF7A3C", 
  "out" = "#00A850"
)

effect_plot <- ggplot(effect_df, aes(x = Elevation, y = Effect, color = Region)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha = 0.2, color = NA
  ) +
  theme_classic() +
  labs(
    x = "Elevation (m)",
    y = "Effect on FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 4000), breaks = seq(0, 4000, by =1000), 
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.text.y = element_text(color = "black", size = 15),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 16),
    legend.position = c(0.3, 0.3), 
    legend.text = element_text(size = 14),
  ) +
  scale_color_manual(values = urban_colors) +
  scale_fill_manual(values = urban_colors) 

print(effect_plot)

# 绘制直方图
hist_plot <- ggplot(merged_dt_Elevation, aes(x = Elevation)) +
  geom_histogram(aes(y = ..density.., fill=urban),  binwidth =  15,  alpha = 0.5) +
  scale_x_continuous(limits = c(0, 4000), breaks = seq(0, 4000, by =1000)) +
  labs(x = "Elevation (m)",
       y = "",
       subtitle = paste0("")) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16),
    plot.margin = margin(t=-1.32,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text.y = element_text(color = "black", size = 15),
    axis.text.x = element_text(color = "black", size = 15),
    legend.position = "none") +
  scale_fill_manual(values = urban_colors_hist) 

print(hist_plot)

# 合并
effect_elevation <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(effect_elevation)


# 绘制边际效应图
# Figure S5a
marginal_plot <- ggplot(marginal_df, aes(x = Elevation, y = Marginal_Effect, color = Region)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha = 0.2, color = NA
  ) +
  labs(
    x = "Elevation (m)",
    y = "Marginal effect on FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  theme_classic() +
  scale_x_continuous(limits = c(0, 4000), breaks = seq(0, 4000, by =1000)) +
  scale_y_continuous(limits = c(-0.0004, 0.0002),
                     breaks = seq(-0.0004, 0.0002, by = 0.0002),
                     labels = scales::number_format(accuracy = 0.0001)) +
  
  theme(
    axis.text.y = element_text(color = "black", size = 15),  
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),  
    axis.title.y = element_text(size = 16), 
    legend.position = c(0.3,0.3), 
    legend.text = element_text(size = 14),  
    legend.background = element_blank()
  ) +
  scale_color_manual(values = urban_colors) +
  scale_fill_manual(values = urban_colors)

print(marginal_plot)

marginal_elevation <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(marginal_elevation)



################################################### latitude ################################################## 
cols_to_center <- c('WaterDeficit', 'Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation',
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_Latitude <- copy(merged_dt)
merged_dt_Latitude <- merged_dt_Latitude[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                         .SDcols = cols_to_center]
merged_dt_Latitude$urban <- factor(merged_dt_Latitude$urban, levels = c("out", "in"))

# linear regression
# Table 2 -- model 2
model_urban <- felm(FVC ~ 1 + WaterDeficit +
                      Seasonalityindex +
                      AverageTemperature +
                      TemperatureRange +
                      WindSpeed +
                      Elevation +
                      Latitude +
                      # Latitude_2 +
                      Precipitation +
                      SoilMoisture +
                      SoilPH + 
                      Longitude + 
                      
                      Latitude:urban +
                      # Latitude_2:urban +
                      
                      GDP_per_capita_PPP + HDI + Population +
                      ImperviousSurface + HumanSettlement + CityArea
                    |  year + CountryID, data = merged_dt_Latitude)

summary(model_urban)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_out.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

# nonlinear regression
# Table 2 -- model 4
model_urban <- felm(FVC ~ 1 + WaterDeficit +
                      Seasonalityindex +
                      AverageTemperature +
                      TemperatureRange +
                      WindSpeed +
                      Elevation +
                      Latitude +
                      Latitude_2 +
                      Precipitation +
                      SoilMoisture +
                      SoilPH + 
                      Longitude + 
                      
                      Latitude:urban +
                      Latitude_2:urban +
                      
                      GDP_per_capita_PPP + HDI + Population +
                      ImperviousSurface + HumanSettlement + CityArea
                    |  year + CountryID, data = merged_dt_Latitude)

summary(model_urban)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_out2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

######################################### 绘制效应图与边际效应图 #######################################
coeffs <- coef(model_urban)
beta1 <- coeffs["Latitude"]
beta2 <- coeffs["Latitude_2"]
beta1_city <- coeffs["Latitude:urbanin"]
beta2_city <- coeffs["Latitude_2:urbanin"]
beta1_total <- beta1 + beta1_city
beta2_total <- beta2 + beta2_city

P_min <- quantile(merged_dt_Latitude$Latitude, 0.001, na.rm = TRUE)
P_max <- quantile(merged_dt_Latitude$Latitude, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)

# 计算效应值
effect_non_urban <- beta1 * P_values + beta2 * (P_values^2)
effect_urban <- beta1_total * P_values + beta2_total * (P_values^2)

vcov_matrix <- vcovHC(model_urban, type = "HC1", cluster = "group", group = merged_dt_Latitude$CountryID)

cov_vars_non_urban <- c("Latitude", "Latitude_2")
vcov_non_urban <- vcov_matrix[cov_vars_non_urban, cov_vars_non_urban]
d_beta1 <- P_values
d_beta2 <- P_values^2
var_effect_non_urban <- (d_beta1^2) * vcov_non_urban["Latitude", "Latitude"] +
  (d_beta2^2) * vcov_non_urban["Latitude_2", "Latitude_2"] +
  2 * d_beta1 * d_beta2 * vcov_non_urban["Latitude", "Latitude_2"]
se_effect_non_urban <- sqrt(var_effect_non_urban)
ci_effect_non_urban_lower <- effect_non_urban - 1.96 * se_effect_non_urban
ci_effect_non_urban_upper <- effect_non_urban + 1.96 * se_effect_non_urban

cov_vars_urban <- c("Latitude", "Latitude_2", "Latitude:urbanin", "Latitude_2:urbanin")
vcov_urban <- vcov_matrix[cov_vars_urban, cov_vars_urban]
d_beta1_u <- P_values
d_beta2_u <- P_values^2
var_effect_urban <- (d_beta1_u^2) * vcov_urban["Latitude", "Latitude"] +
  (d_beta2_u^2) * vcov_urban["Latitude_2", "Latitude_2"] +
  (d_beta1_u^2) * vcov_urban["Latitude:urbanin", "Latitude:urbanin"] +
  (d_beta2_u^2) * vcov_urban["Latitude_2:urbanin", "Latitude_2:urbanin"] +
  2 * d_beta1_u^2 * vcov_urban["Latitude", "Latitude:urbanin"] +
  2 * d_beta2_u^2 * vcov_urban["Latitude_2", "Latitude_2:urbanin"] +
  2 * d_beta1_u * d_beta2_u * vcov_urban["Latitude", "Latitude_2"] +
  2 * d_beta1_u * d_beta2_u * vcov_urban["Latitude", "Latitude_2:urbanin"] +
  2 * d_beta1_u * d_beta2_u * vcov_urban["Latitude:urbanin", "Latitude_2"] +
  2 * d_beta1_u * d_beta2_u * vcov_urban["Latitude:urbanin", "Latitude_2:urbanin"]
se_effect_urban <- sqrt(var_effect_urban)
ci_effect_urban_lower <- effect_urban - 1.96 * se_effect_urban
ci_effect_urban_upper <- effect_urban + 1.96 * se_effect_urban

# 计算边际效应
marginal_non_urban <- beta1 + 2 * beta2 * P_values
marginal_urban <- beta1_total + 2 * beta2_total * P_values

d_marg_beta1 <- 1
d_marg_beta2 <- 2 * P_values
var_marginal_non_urban <- (d_marg_beta1^2) * vcov_non_urban["Latitude", "Latitude"] +
  (d_marg_beta2^2) * vcov_non_urban["Latitude_2", "Latitude_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_non_urban["Latitude", "Latitude_2"]
se_marginal_non_urban <- sqrt(var_marginal_non_urban)
ci_marginal_non_urban_lower <- marginal_non_urban - 1.96 * se_marginal_non_urban
ci_marginal_non_urban_upper <- marginal_non_urban + 1.96 * se_marginal_non_urban

d_marg_beta1_u <- 1
d_marg_beta2_u <- 2 * P_values
var_marginal_urban <- (d_marg_beta1_u^2) * vcov_urban["Latitude", "Latitude"] +
  (d_marg_beta2_u^2) * vcov_urban["Latitude_2", "Latitude_2"] +
  (d_marg_beta1_u^2) * vcov_urban["Latitude:urbanin", "Latitude:urbanin"] +
  (d_marg_beta2_u^2) * vcov_urban["Latitude_2:urbanin", "Latitude_2:urbanin"] +
  2 * d_marg_beta1_u^2 * vcov_urban["Latitude", "Latitude:urbanin"] +
  2 * d_marg_beta2_u^2 * vcov_urban["Latitude_2", "Latitude_2:urbanin"] +
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Latitude", "Latitude_2"] +
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Latitude", "Latitude_2:urbanin"] +
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Latitude:urbanin", "Latitude_2"] +
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Latitude:urbanin", "Latitude_2:urbanin"]
se_marginal_urban <- sqrt(var_marginal_urban)
ci_marginal_urban_lower <- marginal_urban - 1.96 * se_marginal_urban
ci_marginal_urban_upper <- marginal_urban + 1.96 * se_marginal_urban

effect_df <- rbind(
  data.frame(
    Latitude = P_values,
    Region = "Non-urban",
    Effect = effect_non_urban,
    SE = se_effect_non_urban,
    CI_Lower = ci_effect_non_urban_lower,
    CI_Upper = ci_effect_non_urban_upper
  ),
  data.frame(
    Latitude = P_values,
    Region = "Urban",
    Effect = effect_urban,
    SE = se_effect_urban,
    CI_Lower = ci_effect_urban_lower,
    CI_Upper = ci_effect_urban_upper
  )
)

marginal_df <- rbind(
  data.frame(
    Latitude = P_values,
    Region = "Non-urban",
    Marginal_Effect = marginal_non_urban,
    SE = se_marginal_non_urban,
    CI_Lower = ci_marginal_non_urban_lower,
    CI_Upper = ci_marginal_non_urban_upper
  ),
  data.frame(
    Latitude = P_values,
    Region = "Urban",
    Marginal_Effect = marginal_urban,
    SE = se_marginal_urban,
    CI_Lower = ci_marginal_urban_lower,
    CI_Upper = ci_marginal_urban_upper
  )
)

# 绘制效应值图
# Figure 3b
effect_plot <- ggplot(effect_df, aes(x = Latitude, y = Effect, color = Region)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha = 0.2, color = NA
  ) +
  theme_classic() +
  labs(
    x = "Latitude (°)",
    y = "Effect on FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, by =10), 
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.text.y = element_text(color = "black", size = 15),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 16),
    legend.position = "none", 
    legend.text = element_text(size = 14),
  ) +
  scale_color_manual(values = urban_colors) +
  scale_fill_manual(values = urban_colors) 

print(effect_plot)

# 绘制直方图
hist_plot <- ggplot(merged_dt_Latitude, aes(x = Latitude)) +
  geom_histogram(aes(y = ..density.., fill=urban),  binwidth =  0.05,  alpha = 0.5) +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, by =10)) +
  labs(x = "Abs(Latitude) (°)",
       y = "",
       subtitle = paste0("")) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16),
    plot.margin = margin(t=-1.32,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text.y = element_text(color = "black", size = 15),
    axis.text.x = element_text(color = "black", size = 15),
    legend.position = "none") +
  scale_fill_manual(values = urban_colors_hist) 

print(hist_plot)

# 合并
effect_latitude <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(effect_latitude)

# 绘制边际效应图
# Figure S5b
marginal_plot <- ggplot(marginal_df, aes(x = Latitude, y = Marginal_Effect, color = Region)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha = 0.2, color = NA
  ) +
  labs(
    x = "Abs(Latitude) (°)",
    y = "Marginal effect on FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  theme_classic() +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, by =10)) +
  # scale_y_continuous(limits = c(-0.0004, 0.0004),
  #                    breaks = seq(-0.0004, 0.0004, by = 0.0002),
  #                    labels = scales::number_format(accuracy = 0.0001)) +
  
  theme(
    axis.text.y = element_text(color = "black", size = 15), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),  
    axis.title.y = element_text(size = 16),  
    legend.position = "none", 
    legend.text = element_text(size = 14),  
  ) +
  scale_color_manual(values = urban_colors) +
  scale_fill_manual(values = urban_colors)

print(marginal_plot)

marginal_latitude <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(marginal_latitude)



################################################### 城市内异质性分析 ################################################## 

##################################################### continent ####################################################
################################################### elevation ################################################## 

urban_in_dt <- fread("E:/2025/UrbanGreenSpaceWorks/Rcodes/dataset_in_urban.csv")
dim(urban_in_dt)

# 将纬度取绝对值
urban_in_dt$Latitude <- abs(urban_in_dt$Latitude)

# 二次项 
urban_in_dt$WaterDeficit_2 <- urban_in_dt$WaterDeficit * urban_in_dt$WaterDeficit
urban_in_dt$Seasonalityindex_2 <- urban_in_dt$Seasonalityindex * urban_in_dt$Seasonalityindex
urban_in_dt$Elevation_2 <- urban_in_dt$Elevation * urban_in_dt$Elevation
urban_in_dt$Latitude_2 <- urban_in_dt$Latitude * urban_in_dt$Latitude

# 筛选条件
urban_in_dt0 <- copy(urban_in_dt)
urban_in_dt <- urban_in_dt0[FVC > 0 & AverageTemperature > 0]
urban_in_dt <- urban_in_dt[!(continent == "SouthAmerica" & TemperatureRange > 14.5)]
urban_in_dt <- urban_in_dt[!(continent == "SouthAmerica" & WaterDeficit > 2000)]
urban_in_dt <- urban_in_dt[!(continent == "Oceania" & Seasonalityindex > 0.04)]
urban_in_dt <- urban_in_dt[!(ClimateZone == "Cold" & Seasonalityindex > 0.045)]


################################################### elevation ################################################## 
cols_to_center <- c('WaterDeficit', 'Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea') 
urban_in_dt_Elevation <- copy(urban_in_dt)
urban_in_dt_Elevation <- urban_in_dt_Elevation[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                               .SDcols = cols_to_center]
urban_in_dt_Elevation$continent <- factor(urban_in_dt_Elevation$continent, levels = c( "Asia",
                                                                                       "NorthAmerica",
                                                                                       "Africa",
                                                                                       "Europe",
                                                                                       "Oceania",
                                                                                       "SouthAmerica"))
# linear regression
# Table S3 -- model 3
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         Seasonalityindex +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         # Elevation_2 +
                         Latitude +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         Elevation:continent +
                         # Elevation_2:continent +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea
                       |  year + CountryID, data = urban_in_dt_Elevation)

summary(model_urban_in)


######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_continent.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

# nonlinear regression
# Table S3 -- model 7
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         Seasonalityindex +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Elevation_2 +
                         Latitude +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         Elevation:continent +
                         Elevation_2:continent +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea
                       |  year + CountryID, data = urban_in_dt_Elevation)

summary(model_urban_in)


######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_continent2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

######################################### 绘制效应图与边际效应图 #######################################
beta1 <- coef(model_urban_in)["Elevation"]
beta2 <- coef(model_urban_in)["Elevation_2"]

continents <- c("NorthAmerica","Africa",  "Europe",  "Oceania", "SouthAmerica")

beta1_continent <- sapply(continents, function(cont) {
  coef(model_urban_in)[paste0("Elevation:continent", cont)]
})
beta2_continent <- sapply(continents, function(cont) {
  coef(model_urban_in)[paste0("Elevation_2:continent", cont)]
})

beta1_total <- beta1 + beta1_continent
beta2_total <- beta2 + beta2_continent

names(beta1_total) <- continents
names(beta2_total) <- continents

P_min <- quantile(urban_in_dt_Elevation$Elevation, 0.001, na.rm = TRUE)
P_max <- quantile(urban_in_dt_Elevation$Elevation, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)

# 计算效应值
effect_base <- beta1 * P_values + beta2 * (P_values^2)

effect_continent <- lapply(continents, function(cont) {
  beta1_total[[cont]] * P_values + beta2_total[[cont]] * (P_values^2)
})
names(effect_continent) <- continents

vcov_matrix <- vcovHC(model_urban_in, type = "HC1", cluster = "group",
                      group = urban_in_dt_Elevation$CountryID)

cov_vars_base <- c("Elevation", "Elevation_2")
vcov_base <- vcov_matrix[cov_vars_base, cov_vars_base]

d_beta1 <- P_values
d_beta2 <- P_values^2

var_effect_base <- (d_beta1^2) * vcov_base["Elevation", "Elevation"] +
  (d_beta2^2) * vcov_base["Elevation_2", "Elevation_2"] +
  2 * d_beta1 * d_beta2 * vcov_base["Elevation", "Elevation_2"]
se_effect_base <- sqrt(var_effect_base)
ci_effect_base_lower <- effect_base - 1.96 * se_effect_base
ci_effect_base_upper <- effect_base + 1.96 * se_effect_base

var_effect_continent <- list()
se_effect_continent <- list()
ci_effect_continent_lower <- list()
ci_effect_continent_upper <- list()

for (cont in continents) {
  cov_vars <- c("Elevation", "Elevation_2", 
                paste0("Elevation:continent", cont), 
                paste0("Elevation_2:continent", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  
  var_effect <- (d_beta1^2) * vcov_cont["Elevation", "Elevation"] +
    (d_beta2^2) * vcov_cont["Elevation_2", "Elevation_2"] +
    (d_beta1^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_beta2^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_beta1^2 * vcov_cont["Elevation", cov_vars[3]] +
    2 * d_beta2^2 * vcov_cont["Elevation_2", cov_vars[4]] +
    2 * d_beta1 * d_beta2 * vcov_cont["Elevation", "Elevation_2"] +
    2 * d_beta1 * d_beta2 * vcov_cont["Elevation", cov_vars[4]] +
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], "Elevation_2"] +
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], cov_vars[4]]
  
  var_effect_continent[[cont]] <- var_effect
  se_effect_continent[[cont]] <- sqrt(var_effect)
  ci_effect_continent_lower[[cont]] <- effect_continent[[cont]] - 1.96 * se_effect_continent[[cont]]
  ci_effect_continent_upper[[cont]] <- effect_continent[[cont]] + 1.96 * se_effect_continent[[cont]]
}

# 计算边际效应
marginal_base <- beta1 + 2 * beta2 * P_values

marginal_continent <- lapply(continents, function(cont) {
  beta1_total[[cont]] + 2 * beta2_total[[cont]] * P_values
})
names(marginal_continent) <- continents

d_marg_beta1 <- 1
d_marg_beta2 <- 2 * P_values

var_marginal_base <- (d_marg_beta1^2) * vcov_base["Elevation", "Elevation"] +
  (d_marg_beta2^2) * vcov_base["Elevation_2", "Elevation_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_base["Elevation", "Elevation_2"]
se_marginal_base <- sqrt(var_marginal_base)
ci_marginal_base_lower <- marginal_base - 1.96 * se_marginal_base
ci_marginal_base_upper <- marginal_base + 1.96 * se_marginal_base

var_marginal_continent <- list()
se_marginal_continent <- list()
ci_marginal_continent_lower <- list()
ci_marginal_continent_upper <- list()

for (cont in continents) {
  cov_vars <- c("Elevation", "Elevation_2",
                paste0("Elevation:continent", cont),
                paste0("Elevation_2:continent", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  
  d_marg_beta1_u <- 1
  d_marg_beta2_u <- 2 * P_values
  
  var_marginal <- (d_marg_beta1_u^2) * vcov_cont["Elevation", "Elevation"] +
    (d_marg_beta2_u^2) * vcov_cont["Elevation_2", "Elevation_2"] +
    (d_marg_beta1_u^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_marg_beta2_u^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_marg_beta1_u^2 * vcov_cont["Elevation", cov_vars[3]] +
    2 * d_marg_beta2_u^2 * vcov_cont["Elevation_2", cov_vars[4]] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["Elevation", "Elevation_2"] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["Elevation", cov_vars[4]] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], "Elevation_2"] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], cov_vars[4]]
  
  var_marginal_continent[[cont]] <- var_marginal
  se_marginal_continent[[cont]] <- sqrt(var_marginal)
  ci_marginal_continent_lower[[cont]] <- marginal_continent[[cont]] - 1.96 * se_marginal_continent[[cont]]
  ci_marginal_continent_upper[[cont]] <- marginal_continent[[cont]] + 1.96 * se_marginal_continent[[cont]]
}

effect_df <- data.frame(
  Elevation = P_values,
  Continent = "基准组",
  Effect = effect_base,
  SE = se_effect_base,
  CI_Lower = ci_effect_base_lower,
  CI_Upper = ci_effect_base_upper
)

for (cont in continents) {
  temp_df <- data.frame(
    Elevation = P_values,
    Continent = cont,
    Effect = effect_continent[[cont]],
    SE = se_effect_continent[[cont]],
    CI_Lower = ci_effect_continent_lower[[cont]],
    CI_Upper = ci_effect_continent_upper[[cont]]
  )
  effect_df <- rbind(effect_df, temp_df)
}

marginal_df <- data.frame(
  Elevation = P_values,
  Continent = "基准组",
  Marginal_Effect = marginal_base,
  SE = se_marginal_base,
  CI_Lower = ci_marginal_base_lower,
  CI_Upper = ci_marginal_base_upper
)

for (cont in continents) {
  temp_df <- data.frame(
    Elevation = P_values,
    Continent = cont,
    Marginal_Effect = marginal_continent[[cont]],
    SE = se_marginal_continent[[cont]],
    CI_Lower = ci_marginal_continent_lower[[cont]],
    CI_Upper = ci_marginal_continent_upper[[cont]]
  )
  marginal_df <- rbind(marginal_df, temp_df)
}

effect_df$Continent <- factor(effect_df$Continent, 
                              levels = c("基准组", continents),
                              labels = c("Asia", "North America", "Africa",  "Europe", "Oceania", "South America"))

marginal_df$Continent <- factor(marginal_df$Continent,
                                levels = c("基准组", continents),
                                labels = c("Asia", "North America", "Africa",  "Europe", "Oceania", "South America"))

continent_limits <- tibble(
  Continent = c("Asia", "Africa", "North America", "South America", "Europe", "Oceania"),
  min_WD = c(
    min(urban_in_dt_Elevation[urban_in_dt_Elevation$continent == "Asia", ]$Elevation, na.rm = TRUE),
    min(urban_in_dt_Elevation[urban_in_dt_Elevation$continent == "Africa", ]$Elevation, na.rm = TRUE),
    min(urban_in_dt_Elevation[urban_in_dt_Elevation$continent == "NorthAmerica", ]$Elevation, na.rm = TRUE),
    min(urban_in_dt_Elevation[urban_in_dt_Elevation$continent == "SouthAmerica", ]$Elevation, na.rm = TRUE),
    min(urban_in_dt_Elevation[urban_in_dt_Elevation$continent == "Europe", ]$Elevation, na.rm = TRUE),
    min(urban_in_dt_Elevation[urban_in_dt_Elevation$continent == "Oceania", ]$Elevation, na.rm = TRUE)
  ),
  max_WD = c(
    max(urban_in_dt_Elevation[urban_in_dt_Elevation$continent == "Asia", ]$Elevation, na.rm = TRUE),
    max(urban_in_dt_Elevation[urban_in_dt_Elevation$continent == "Africa", ]$Elevation, na.rm = TRUE),
    max(urban_in_dt_Elevation[urban_in_dt_Elevation$continent == "NorthAmerica", ]$Elevation, na.rm = TRUE),
    max(urban_in_dt_Elevation[urban_in_dt_Elevation$continent == "SouthAmerica", ]$Elevation, na.rm = TRUE),
    max(urban_in_dt_Elevation[urban_in_dt_Elevation$continent == "Europe", ]$Elevation, na.rm = TRUE),
    max(urban_in_dt_Elevation[urban_in_dt_Elevation$continent == "Oceania", ]$Elevation, na.rm = TRUE)
  )
)

print(continent_limits)

effect_df <- effect_df %>%
  left_join(continent_limits, by = "Continent") %>%
  filter(Elevation >= min_WD & Elevation <= max_WD) %>%
  select(-min_WD, -max_WD)

marginal_df <- marginal_df %>%
  left_join(continent_limits, by = "Continent") %>%
  filter(Elevation >= min_WD & Elevation <= max_WD) %>%
  select(-min_WD, -max_WD)

# 绘制效应值图
# Figure 3c
continent_colors <- c(
  "Asia" = "#3498db",
  "Africa" = "#95a5a6",
  "North America" = "#2ecc71",
  "South America" = "#9b59b6",
  "Europe" = "#e74c3c",
  "Oceania" = "#f39c12"
)

effect_df$Continent <- factor(effect_df$Continent, levels = c("Asia", "Africa", "North America", "South America", "Europe", "Oceania"))

effect_plot <- ggplot(effect_df, aes(x = Elevation, y = Effect, color = Continent)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Continent),
    alpha = 0.15, color = NA
  ) +
  theme_classic() +
  labs(
    x = "",
    y = "Effect on urban FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 4050), breaks = seq(0, 4050, by =1000), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.1, 0.06),
                     breaks = seq(-0.1, 0.06, by = 0.05),
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.text.y = element_text(color = "black", size = 13),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.35, 0.36),
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'), 
    legend.background = element_blank()
  ) +
  scale_color_manual(values = continent_colors) +
  scale_fill_manual(values = continent_colors)

print(effect_plot)

# 绘制直方图
hist_plot <- ggplot(urban_in_dt_Elevation, aes(x = Elevation)) +
  geom_histogram(aes(y = ..density.., fill=continent), binwidth = 2, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 4050), breaks = seq(0, 4050, by =1000)) +
  scale_y_continuous(limits = c(0, 0.09), breaks = seq(0, 0.09, by =0.04)) +
  labs(x = "",
       y = "") +
  theme_classic() +
  theme(
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text = element_text(color = "black", size = 13),
    axis.title = element_text(size = 14),
    legend.position = "none"
  ) + 
  scale_fill_manual(values = c("Africa" = "#95a5a6",
                               "Asia" = "#3498db",
                               "Europe" = "#e74c3c",
                               "NorthAmerica" = "#2ecc71",
                               "Oceania" = "#f39c12",
                               "SouthAmerica" = "#9b59b6"))

print(hist_plot)

# 合并
effect_continent_elevation <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_continent_elevation)


# 绘制边际效应图
# Figure S5c
marginal_df$Continent <- factor(marginal_df$Continent, levels = c("Asia", "Africa", "North America", "South America", "Europe", "Oceania"))

marginal_plot <- ggplot(marginal_df, aes(x = Elevation, y = Marginal_Effect, color = Continent)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Continent),
    alpha = 0.15, color = NA
  ) +
  labs(
    x = "",
    y = "Marginal effect on urban FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  theme_classic() +
  scale_x_continuous(limits = c(0, 4050), breaks = seq(0, 4050, by =1000), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.00008, 0.00016),
                     breaks = seq(-0.00008, 0.00016, by = 0.00008),
                     labels = scales::number_format(accuracy = 0.00001)) +
  theme(
    axis.text.y = element_text(color = "black", size = 12),  
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),  
    axis.title.y = element_text(size = 12),  
    legend.position = c(0.55,0.85), 
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'), 
    legend.background = element_blank()
  ) +
  scale_color_manual(values = continent_colors) +
  scale_fill_manual(values = continent_colors)

print(marginal_plot)

marginal_continent_elevation <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_continent_elevation)


##################################################### continent ####################################################
################################################### latitude ######################################
cols_to_center <- c('WaterDeficit', 'Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation',  
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
urban_in_dt_Latitude <- copy(urban_in_dt)
urban_in_dt_Latitude <- urban_in_dt_Latitude[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                             .SDcols = cols_to_center]
urban_in_dt_Latitude$continent <- factor(urban_in_dt_Latitude$continent, levels = c( "Asia",
                                                                                     "NorthAmerica",
                                                                                     "Africa",
                                                                                     "Europe",
                                                                                     "Oceania",
                                                                                     "SouthAmerica"))
# linear regression
# Table S3 -- model 4
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         Seasonalityindex +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Latitude +
                         # Latitude_2 +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         Latitude:continent +
                         # Latitude_2:continent +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea
                       |  year + CountryID, data = urban_in_dt_Latitude)

summary(model_urban_in)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_continent.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  


# nonlinear regression
# Table S3 -- model 8
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         Seasonalityindex +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Latitude +
                         Latitude_2 +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         Latitude:continent +
                         Latitude_2:continent +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea
                       |  year + CountryID, data = urban_in_dt_Latitude)

summary(model_urban_in)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_continent2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  


######################################### 绘制效应图与边际效应图 #######################################
beta1 <- coef(model_urban_in)["Latitude"]
beta2 <- coef(model_urban_in)["Latitude_2"]

continents <- c("NorthAmerica","Africa","Europe", "Oceania", "SouthAmerica")

beta1_continent <- sapply(continents, function(cont) {
  coef(model_urban_in)[paste0("Latitude:continent", cont)]
})
beta2_continent <- sapply(continents, function(cont) {
  coef(model_urban_in)[paste0("Latitude_2:continent", cont)]
})

beta1_total <- beta1 + beta1_continent
beta2_total <- beta2 + beta2_continent

names(beta1_total) <- continents
names(beta2_total) <- continents

P_min <- quantile(urban_in_dt_Latitude$Latitude, 0.001, na.rm = TRUE)
P_max <- quantile(urban_in_dt_Latitude$Latitude, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)

# 计算效应值
effect_base <- beta1 * P_values + beta2 * (P_values^2)

effect_continent <- lapply(continents, function(cont) {
  beta1_total[[cont]] * P_values + beta2_total[[cont]] * (P_values^2)
})
names(effect_continent) <- continents

vcov_matrix <- vcovHC(model_urban_in, type = "HC1", cluster = "group",
                      group = urban_in_dt_Latitude$CountryID)

cov_vars_base <- c("Latitude", "Latitude_2")
vcov_base <- vcov_matrix[cov_vars_base, cov_vars_base]

d_beta1 <- P_values
d_beta2 <- P_values^2

var_effect_base <- (d_beta1^2) * vcov_base["Latitude", "Latitude"] +
  (d_beta2^2) * vcov_base["Latitude_2", "Latitude_2"] +
  2 * d_beta1 * d_beta2 * vcov_base["Latitude", "Latitude_2"]
se_effect_base <- sqrt(var_effect_base)
ci_effect_base_lower <- effect_base - 1.96 * se_effect_base
ci_effect_base_upper <- effect_base + 1.96 * se_effect_base

var_effect_continent <- list()
se_effect_continent <- list()
ci_effect_continent_lower <- list()
ci_effect_continent_upper <- list()

for (cont in continents) {
  cov_vars <- c("Latitude", "Latitude_2", 
                paste0("Latitude:continent", cont), 
                paste0("Latitude_2:continent", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  
  var_effect <- (d_beta1^2) * vcov_cont["Latitude", "Latitude"] +
    (d_beta2^2) * vcov_cont["Latitude_2", "Latitude_2"] +
    (d_beta1^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_beta2^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_beta1^2 * vcov_cont["Latitude", cov_vars[3]] +
    2 * d_beta2^2 * vcov_cont["Latitude_2", cov_vars[4]] +
    2 * d_beta1 * d_beta2 * vcov_cont["Latitude", "Latitude_2"] +
    2 * d_beta1 * d_beta2 * vcov_cont["Latitude", cov_vars[4]] +
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], "Latitude_2"] +
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], cov_vars[4]]
  
  var_effect_continent[[cont]] <- var_effect
  se_effect_continent[[cont]] <- sqrt(var_effect)
  ci_effect_continent_lower[[cont]] <- effect_continent[[cont]] - 1.96 * se_effect_continent[[cont]]
  ci_effect_continent_upper[[cont]] <- effect_continent[[cont]] + 1.96 * se_effect_continent[[cont]]
}

# 计算边际效应
marginal_base <- beta1 + 2 * beta2 * P_values

marginal_continent <- lapply(continents, function(cont) {
  beta1_total[[cont]] + 2 * beta2_total[[cont]] * P_values
})
names(marginal_continent) <- continents

d_marg_beta1 <- 1
d_marg_beta2 <- 2 * P_values

var_marginal_base <- (d_marg_beta1^2) * vcov_base["Latitude", "Latitude"] +
  (d_marg_beta2^2) * vcov_base["Latitude_2", "Latitude_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_base["Latitude", "Latitude_2"]
se_marginal_base <- sqrt(var_marginal_base)
ci_marginal_base_lower <- marginal_base - 1.96 * se_marginal_base
ci_marginal_base_upper <- marginal_base + 1.96 * se_marginal_base

var_marginal_continent <- list()
se_marginal_continent <- list()
ci_marginal_continent_lower <- list()
ci_marginal_continent_upper <- list()

for (cont in continents) {
  cov_vars <- c("Latitude", "Latitude_2",
                paste0("Latitude:continent", cont),
                paste0("Latitude_2:continent", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  
  d_marg_beta1_u <- 1
  d_marg_beta2_u <- 2 * P_values
  
  var_marginal <- (d_marg_beta1_u^2) * vcov_cont["Latitude", "Latitude"] +
    (d_marg_beta2_u^2) * vcov_cont["Latitude_2", "Latitude_2"] +
    (d_marg_beta1_u^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_marg_beta2_u^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_marg_beta1_u^2 * vcov_cont["Latitude", cov_vars[3]] +
    2 * d_marg_beta2_u^2 * vcov_cont["Latitude_2", cov_vars[4]] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["Latitude", "Latitude_2"] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["Latitude", cov_vars[4]] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], "Latitude_2"] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], cov_vars[4]]
  
  var_marginal_continent[[cont]] <- var_marginal
  se_marginal_continent[[cont]] <- sqrt(var_marginal)
  ci_marginal_continent_lower[[cont]] <- marginal_continent[[cont]] - 1.96 * se_marginal_continent[[cont]]
  ci_marginal_continent_upper[[cont]] <- marginal_continent[[cont]] + 1.96 * se_marginal_continent[[cont]]
}

effect_df <- data.frame(
  Latitude = P_values,
  Continent = "基准组",
  Effect = effect_base,
  SE = se_effect_base,
  CI_Lower = ci_effect_base_lower,
  CI_Upper = ci_effect_base_upper
)

for (cont in continents) {
  temp_df <- data.frame(
    Latitude = P_values,
    Continent = cont,
    Effect = effect_continent[[cont]],
    SE = se_effect_continent[[cont]],
    CI_Lower = ci_effect_continent_lower[[cont]],
    CI_Upper = ci_effect_continent_upper[[cont]]
  )
  effect_df <- rbind(effect_df, temp_df)
}

marginal_df <- data.frame(
  Latitude = P_values,
  Continent = "基准组",
  Marginal_Effect = marginal_base,
  SE = se_marginal_base,
  CI_Lower = ci_marginal_base_lower,
  CI_Upper = ci_marginal_base_upper
)

for (cont in continents) {
  temp_df <- data.frame(
    Latitude = P_values,
    Continent = cont,
    Marginal_Effect = marginal_continent[[cont]],
    SE = se_marginal_continent[[cont]],
    CI_Lower = ci_marginal_continent_lower[[cont]],
    CI_Upper = ci_marginal_continent_upper[[cont]]
  )
  marginal_df <- rbind(marginal_df, temp_df)
}

effect_df$Continent <- factor(effect_df$Continent, 
                              levels = c("基准组", continents),
                              labels = c("Asia", "North America","Africa","Europe", "Oceania", "South America"))

marginal_df$Continent <- factor(marginal_df$Continent,
                                levels = c("基准组", continents),
                                labels = c("Asia", "North America","Africa","Europe", "Oceania", "South America"))

continent_limits <- tibble(
  Continent = c("Asia", "Africa", "North America", "South America", "Europe", "Oceania"),
  min_WD = c(
    min(urban_in_dt_Latitude[urban_in_dt_Latitude$continent == "Asia", ]$Latitude, na.rm = TRUE),
    min(urban_in_dt_Latitude[urban_in_dt_Latitude$continent == "Africa", ]$Latitude, na.rm = TRUE),
    min(urban_in_dt_Latitude[urban_in_dt_Latitude$continent == "NorthAmerica", ]$Latitude, na.rm = TRUE),
    min(urban_in_dt_Latitude[urban_in_dt_Latitude$continent == "SouthAmerica", ]$Latitude, na.rm = TRUE),
    min(urban_in_dt_Latitude[urban_in_dt_Latitude$continent == "Europe", ]$Latitude, na.rm = TRUE),
    min(urban_in_dt_Latitude[urban_in_dt_Latitude$continent == "Oceania", ]$Latitude, na.rm = TRUE)
  ),
  max_WD = c(
    max(urban_in_dt_Latitude[urban_in_dt_Latitude$continent == "Asia", ]$Latitude, na.rm = TRUE),
    max(urban_in_dt_Latitude[urban_in_dt_Latitude$continent == "Africa", ]$Latitude, na.rm = TRUE),
    max(urban_in_dt_Latitude[urban_in_dt_Latitude$continent == "NorthAmerica", ]$Latitude, na.rm = TRUE),
    max(urban_in_dt_Latitude[urban_in_dt_Latitude$continent == "SouthAmerica", ]$Latitude, na.rm = TRUE),
    max(urban_in_dt_Latitude[urban_in_dt_Latitude$continent == "Europe", ]$Latitude, na.rm = TRUE),
    max(urban_in_dt_Latitude[urban_in_dt_Latitude$continent == "Oceania", ]$Latitude, na.rm = TRUE)
  )
)

print(continent_limits)

effect_df <- effect_df %>%
  left_join(continent_limits, by = "Continent") %>%
  filter(Latitude >= min_WD & Latitude <= max_WD) %>%
  select(-min_WD, -max_WD)

marginal_df <- marginal_df %>%
  left_join(continent_limits, by = "Continent") %>%
  filter(Latitude >= min_WD & Latitude <= max_WD) %>%
  select(-min_WD, -max_WD)

# 绘制效应值图
# Figure 3d
effect_plot <- ggplot(effect_df, aes(x = Latitude, y = Effect, color = Continent)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Continent),
    alpha = 0.15, color = NA
  ) +
  theme_classic() +
  labs(
    x = "",
    y = "",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, by =10), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.2, 0.2),
                     breaks = seq(-0.2, 0.2, by = 0.1),
                     labels = scales::number_format(accuracy = 0.1)) +
  theme(
    axis.text.y = element_text(color = "black", size = 13),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),
  ) +
  scale_color_manual(values = continent_colors) +
  scale_fill_manual(values = continent_colors)

print(effect_plot)

# 绘制直方图
hist_plot <- ggplot(urban_in_dt_Latitude, aes(x = Latitude)) +
  geom_histogram(aes(y = ..density.., fill=continent), binwidth = 0.1, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, by =10)) +
  labs(x = "",
       y = "") +
  theme_classic() +
  theme(
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text = element_text(color = "black", size = 13),
    axis.title = element_text(size = 14),
    legend.position = "none"
  ) + 
  scale_fill_manual(values = c(
    "Asia" = "#3498db",  "Africa" = "#95a5a6",
    "Europe" = "#e74c3c",
    "NorthAmerica" = "#2ecc71",
    "Oceania" = "#f39c12",
    "SouthAmerica" = "#9b59b6"))

print(hist_plot)

# 合并
effect_continent_latitude <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_continent_latitude)

# 绘制边际效应图 
# Figure S5d
marginal_plot <- ggplot(marginal_df, aes(x = Latitude, y = Marginal_Effect, color = Continent)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Continent),
    alpha = 0.15, color = NA
  ) +
  labs(
    x = "",
    y = "",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  theme_classic() +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, by =10), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.01, 0.01),
                     breaks = seq(-0.01, 0.01, by = 0.01),
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.text.y = element_text(color = "black", size = 13),  
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),  
    axis.title.y = element_text(size = 14),  
    legend.position = "none", 
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),  
  ) +
  scale_color_manual(values = continent_colors) +
  scale_fill_manual(values = continent_colors)

print(marginal_plot)

marginal_continent_latitude <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))

print(marginal_continent_latitude)


################################################### 城市内异质性分析 ################################################## 
##################################################### climate ################################################## 
################################################### elevation ######################################
cols_to_center <- c('WaterDeficit', 'Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
urban_in_dt_Elevation <- copy(urban_in_dt)
urban_in_dt_Elevation <- urban_in_dt_Elevation[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                               .SDcols = cols_to_center]
urban_in_dt_Elevation$ClimateZone <- factor(urban_in_dt_Elevation$ClimateZone, levels = c("Temperate",
                                                                                          "Tropical",
                                                                                          "Arid",
                                                                                          "Cold"))
# linear regression
# Table S4 -- model 3
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         Seasonalityindex +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         # Elevation_2 +
                         Latitude +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         Elevation:ClimateZone +
                         # Elevation_2:ClimateZone +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea
                       |  year + CountryID, data = urban_in_dt_Elevation)

summary(model_urban_in)


######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_climate.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  


# nonlinear regression
# Table S4 -- model 7
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         Seasonalityindex +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Elevation_2 +
                         Latitude +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         Elevation:ClimateZone +
                         Elevation_2:ClimateZone +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea
                       |  year + CountryID, data = urban_in_dt_Elevation)

summary(model_urban_in)


######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_climate2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

######################################### 绘制效应图与边际效应图 #######################################
beta1 <- coef(model_urban_in)["Elevation"]
beta2 <- coef(model_urban_in)["Elevation_2"]

climates <- c("Tropical","Arid",  "Cold")

beta1_climate <- sapply(climates, function(cont) {
  coef(model_urban_in)[paste0("Elevation:ClimateZone", cont)]
})
beta2_climate <- sapply(climates, function(cont) {
  coef(model_urban_in)[paste0("Elevation_2:ClimateZone", cont)]
})

beta1_total <- beta1 + beta1_climate
beta2_total <- beta2 + beta2_climate

names(beta1_total) <- climates
names(beta2_total) <- climates

P_min <- quantile(urban_in_dt_Elevation$Elevation, 0.001, na.rm = TRUE)
P_max <- quantile(urban_in_dt_Elevation$Elevation, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)

# 计算效应值
effect_base <- beta1 * P_values + beta2 * (P_values^2)

effect_climate <- lapply(climates, function(cont) {
  beta1_total[[cont]] * P_values + beta2_total[[cont]] * (P_values^2)
})
names(effect_climate) <- climates

vcov_matrix <- vcovHC(model_urban_in, type = "HC1", cluster = "group",
                      group = urban_in_dt_Elevation$CountryID)

cov_vars_base <- c("Elevation", "Elevation_2")
vcov_base <- vcov_matrix[cov_vars_base, cov_vars_base]

d_beta1 <- P_values
d_beta2 <- P_values^2

var_effect_base <- (d_beta1^2) * vcov_base["Elevation", "Elevation"] +
  (d_beta2^2) * vcov_base["Elevation_2", "Elevation_2"] +
  2 * d_beta1 * d_beta2 * vcov_base["Elevation", "Elevation_2"]
se_effect_base <- sqrt(var_effect_base)
ci_effect_base_lower <- effect_base - 1.96 * se_effect_base
ci_effect_base_upper <- effect_base + 1.96 * se_effect_base

var_effect_climate <- list()
se_effect_climate <- list()
ci_effect_climate_lower <- list()
ci_effect_climate_upper <- list()

for (cont in climates) {
  cov_vars <- c("Elevation", "Elevation_2", 
                paste0("Elevation:ClimateZone", cont), 
                paste0("Elevation_2:ClimateZone", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  
  var_effect <- (d_beta1^2) * vcov_cont["Elevation", "Elevation"] +
    (d_beta2^2) * vcov_cont["Elevation_2", "Elevation_2"] +
    (d_beta1^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_beta2^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_beta1^2 * vcov_cont["Elevation", cov_vars[3]] +
    2 * d_beta2^2 * vcov_cont["Elevation_2", cov_vars[4]] +
    2 * d_beta1 * d_beta2 * vcov_cont["Elevation", "Elevation_2"] +
    2 * d_beta1 * d_beta2 * vcov_cont["Elevation", cov_vars[4]] +
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], "Elevation_2"] +
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], cov_vars[4]]
  
  var_effect_climate[[cont]] <- var_effect
  se_effect_climate[[cont]] <- sqrt(var_effect)
  ci_effect_climate_lower[[cont]] <- effect_climate[[cont]] - 1.96 * se_effect_climate[[cont]]
  ci_effect_climate_upper[[cont]] <- effect_climate[[cont]] + 1.96 * se_effect_climate[[cont]]
}

# 计算边际效应
marginal_base <- beta1 + 2 * beta2 * P_values

marginal_climate <- lapply(climates, function(cont) {
  beta1_total[[cont]] + 2 * beta2_total[[cont]] * P_values
})
names(marginal_climate) <- climates

d_marg_beta1 <- 1
d_marg_beta2 <- 2 * P_values

var_marginal_base <- (d_marg_beta1^2) * vcov_base["Elevation", "Elevation"] +
  (d_marg_beta2^2) * vcov_base["Elevation_2", "Elevation_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_base["Elevation", "Elevation_2"]
se_marginal_base <- sqrt(var_marginal_base)
ci_marginal_base_lower <- marginal_base - 1.96 * se_marginal_base
ci_marginal_base_upper <- marginal_base + 1.96 * se_marginal_base

var_marginal_climate <- list()
se_marginal_climate <- list()
ci_marginal_climate_lower <- list()
ci_marginal_climate_upper <- list()

for (cont in climates) {
  cov_vars <- c("Elevation", "Elevation_2",
                paste0("Elevation:ClimateZone", cont),
                paste0("Elevation_2:ClimateZone", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  
  d_marg_beta1_u <- 1
  d_marg_beta2_u <- 2 * P_values
  
  var_marginal <- (d_marg_beta1_u^2) * vcov_cont["Elevation", "Elevation"] +
    (d_marg_beta2_u^2) * vcov_cont["Elevation_2", "Elevation_2"] +
    (d_marg_beta1_u^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_marg_beta2_u^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_marg_beta1_u^2 * vcov_cont["Elevation", cov_vars[3]] +
    2 * d_marg_beta2_u^2 * vcov_cont["Elevation_2", cov_vars[4]] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["Elevation", "Elevation_2"] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["Elevation", cov_vars[4]] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], "Elevation_2"] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], cov_vars[4]]
  
  var_marginal_climate[[cont]] <- var_marginal
  se_marginal_climate[[cont]] <- sqrt(var_marginal)
  ci_marginal_climate_lower[[cont]] <- marginal_climate[[cont]] - 1.96 * se_marginal_climate[[cont]]
  ci_marginal_climate_upper[[cont]] <- marginal_climate[[cont]] + 1.96 * se_marginal_climate[[cont]]
}

effect_df <- data.frame(
  Elevation = P_values,
  Climate = "基准组",
  Effect = effect_base,
  SE = se_effect_base,
  CI_Lower = ci_effect_base_lower,
  CI_Upper = ci_effect_base_upper
)

for (cont in climates) {
  temp_df <- data.frame(
    Elevation = P_values,
    Climate = cont,
    Effect = effect_climate[[cont]],
    SE = se_effect_climate[[cont]],
    CI_Lower = ci_effect_climate_lower[[cont]],
    CI_Upper = ci_effect_climate_upper[[cont]]
  )
  effect_df <- rbind(effect_df, temp_df)
}

marginal_df <- data.frame(
  Elevation = P_values,
  Climate = "基准组",
  Marginal_Effect = marginal_base,
  SE = se_marginal_base,
  CI_Lower = ci_marginal_base_lower,
  CI_Upper = ci_marginal_base_upper
)

for (cont in climates) {
  temp_df <- data.frame(
    Elevation = P_values,
    Climate = cont,
    Marginal_Effect = marginal_climate[[cont]],
    SE = se_marginal_climate[[cont]],
    CI_Lower = ci_marginal_climate_lower[[cont]],
    CI_Upper = ci_marginal_climate_upper[[cont]]
  )
  marginal_df <- rbind(marginal_df, temp_df)
}

effect_df$Climate <- factor(effect_df$Climate, 
                            levels = c("基准组", climates),
                            labels = c("Temperate", "Tropical", "Arid",  "Cold"))

marginal_df$Climate <- factor(marginal_df$Climate,
                              levels = c("基准组", climates),
                              labels = c("Temperate", "Tropical", "Arid",  "Cold"))

climate_limits <- tibble(
  Climate = c("Temperate", "Tropical", "Arid",  "Cold"),
  min_WD = c(
    min(urban_in_dt_Elevation[urban_in_dt_Elevation$ClimateZone == "Temperate", ]$Elevation, na.rm = TRUE),
    min(urban_in_dt_Elevation[urban_in_dt_Elevation$ClimateZone == "Tropical", ]$Elevation, na.rm = TRUE),
    min(urban_in_dt_Elevation[urban_in_dt_Elevation$ClimateZone == "Arid", ]$Elevation, na.rm = TRUE),
    min(urban_in_dt_Elevation[urban_in_dt_Elevation$ClimateZone == "Cold", ]$Elevation, na.rm = TRUE)
  ),
  max_WD = c(
    max(urban_in_dt_Elevation[urban_in_dt_Elevation$ClimateZone == "Temperate", ]$Elevation, na.rm = TRUE),
    max(urban_in_dt_Elevation[urban_in_dt_Elevation$ClimateZone == "Tropical", ]$Elevation, na.rm = TRUE),
    max(urban_in_dt_Elevation[urban_in_dt_Elevation$ClimateZone == "Arid", ]$Elevation, na.rm = TRUE),
    max(urban_in_dt_Elevation[urban_in_dt_Elevation$ClimateZone == "Cold", ]$Elevation, na.rm = TRUE)
  )
)

print(climate_limits)

effect_df <- effect_df %>%
  left_join(climate_limits, by = "Climate") %>%
  filter(Elevation >= min_WD & Elevation <= max_WD) %>%
  select(-min_WD, -max_WD)

marginal_df <- marginal_df %>%
  left_join(climate_limits, by = "Climate") %>%
  filter(Elevation >= min_WD & Elevation <= max_WD) %>%
  select(-min_WD, -max_WD)

# 绘制效应值图
# Figure 3e
climate_colors <- c(
  "Temperate" = "#2ecc71",
  "Tropical" = "#e74c3c",
  "Arid" = "#95a5a6",
  "Cold" = "#3498db"
)

effect_df$Climate <- factor(effect_df$Climate, levels = c( "Tropical", "Temperate",  "Cold","Arid"))

effect_plot <- ggplot(effect_df, aes(x = Elevation, y = Effect, color = Climate)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Climate),
    alpha = 0.15, color = NA
  ) +
  theme_classic() +
  labs(
    x = "",
    y = "Effect on urban FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 4050), breaks = seq(0, 4050, by =1000), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.08, 0.04),
                     breaks = seq(-0.08, 0.04, by = 0.04),
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.text.y = element_text(color = "black", size = 13),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.35,0.35),
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),
    legend.background = element_blank()
  ) +
  scale_color_manual(values = climate_colors) +
  scale_fill_manual(values = climate_colors)

print(effect_plot)

# 绘制直方图
hist_plot <- ggplot(urban_in_dt_Elevation, aes(x = Elevation)) +
  geom_histogram(aes(y = ..density.., fill=ClimateZone), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 4050), breaks = seq(0, 4050, by =1000)) +
  scale_y_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, by =0.02)) +
  labs(x = "",
       y = "") +
  theme_classic() +
  theme(
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text = element_text(color = "black", size = 13),
    axis.title = element_text(size = 14),
    legend.position = "none"
  ) + 
  scale_fill_manual(values = climate_colors)

print(hist_plot)

# 合并
effect_climate_elevation <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_climate_elevation)

# 绘制边际效应图 
# Figure S5e
marginal_df$Climate <- factor(marginal_df$Climate, levels = c("Tropical", "Temperate",  "Cold","Arid"))

marginal_plot <- ggplot(marginal_df, aes(x = Elevation, y = Marginal_Effect, color = Climate)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Climate),
    alpha = 0.15, color = NA
  ) +
  labs(
    x = "",
    y = "Marginal effect on urban FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  theme_classic() +
  scale_x_continuous(limits = c(0, 4050), breaks = seq(0, 4050, by =1000), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.00006, 0.00006),
                     breaks = seq(-0.00006, 0.00006, by = 0.00003),
                     labels = scales::number_format(accuracy = 0.00001)) +
  theme(
    axis.text.y = element_text(color = "black", size = 12),  
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),  
    axis.title.y = element_text(size = 12),  
    legend.position = c(0.25,0.27), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),  
    legend.background = element_blank()
  ) +
  scale_color_manual(values = climate_colors) +
  scale_fill_manual(values = climate_colors)

print(marginal_plot)

marginal_climate_elevation <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_climate_elevation)

##################################################### climate ####################################################
################################################### latitude ######################################
cols_to_center <- c('WaterDeficit', 'Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
urban_in_dt_Latitude <- copy(urban_in_dt)
urban_in_dt_Latitude <- urban_in_dt_Latitude[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                             .SDcols = cols_to_center]
urban_in_dt_Latitude$ClimateZone <- factor(urban_in_dt_Latitude$ClimateZone, levels = c("Temperate",
                                                                                        "Tropical",
                                                                                        "Arid",
                                                                                        "Cold"))
# linear regression
# Table S4 -- model 4
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         Seasonalityindex +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Latitude +
                         # Latitude_2 +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         Latitude:ClimateZone +
                         # Latitude_2:ClimateZone +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea
                       |  year + CountryID, data = urban_in_dt_Latitude)

summary(model_urban_in)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_climate.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

# nonlinear regression
# Table S4 -- model 8
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         Seasonalityindex +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Latitude +
                         Latitude_2 +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         Latitude:ClimateZone +
                         Latitude_2:ClimateZone +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea
                       |  year + CountryID, data = urban_in_dt_Latitude)

summary(model_urban_in)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_climate2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

######################################### 绘制效应图与边际效应图 #######################################
beta1 <- coef(model_urban_in)["Latitude"]
beta2 <- coef(model_urban_in)["Latitude_2"]

climates <- c("Tropical","Arid",  "Cold")

beta1_climate <- sapply(climates, function(cont) {
  coef(model_urban_in)[paste0("Latitude:ClimateZone", cont)]
})
beta2_climate <- sapply(climates, function(cont) {
  coef(model_urban_in)[paste0("Latitude_2:ClimateZone", cont)]
})

beta1_total <- beta1 + beta1_climate
beta2_total <- beta2 + beta2_climate

names(beta1_total) <- climates
names(beta2_total) <- climates

P_min <- quantile(urban_in_dt_Latitude$Latitude, 0.001, na.rm = TRUE)
P_max <- quantile(urban_in_dt_Latitude$Latitude, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)

# 计算效应值
effect_base <- beta1 * P_values + beta2 * (P_values^2)

effect_climate <- lapply(climates, function(cont) {
  beta1_total[[cont]] * P_values + beta2_total[[cont]] * (P_values^2)
})
names(effect_climate) <- climates

vcov_matrix <- vcovHC(model_urban_in, type = "HC1", cluster = "group",
                      group = urban_in_dt_Latitude$CountryID)

cov_vars_base <- c("Latitude", "Latitude_2")
vcov_base <- vcov_matrix[cov_vars_base, cov_vars_base]

d_beta1 <- P_values
d_beta2 <- P_values^2

var_effect_base <- (d_beta1^2) * vcov_base["Latitude", "Latitude"] +
  (d_beta2^2) * vcov_base["Latitude_2", "Latitude_2"] +
  2 * d_beta1 * d_beta2 * vcov_base["Latitude", "Latitude_2"]
se_effect_base <- sqrt(var_effect_base)
ci_effect_base_lower <- effect_base - 1.96 * se_effect_base
ci_effect_base_upper <- effect_base + 1.96 * se_effect_base

var_effect_climate <- list()
se_effect_climate <- list()
ci_effect_climate_lower <- list()
ci_effect_climate_upper <- list()

for (cont in climates) {
  cov_vars <- c("Latitude", "Latitude_2", 
                paste0("Latitude:ClimateZone", cont), 
                paste0("Latitude_2:ClimateZone", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  
  var_effect <- (d_beta1^2) * vcov_cont["Latitude", "Latitude"] +
    (d_beta2^2) * vcov_cont["Latitude_2", "Latitude_2"] +
    (d_beta1^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_beta2^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_beta1^2 * vcov_cont["Latitude", cov_vars[3]] +
    2 * d_beta2^2 * vcov_cont["Latitude_2", cov_vars[4]] +
    2 * d_beta1 * d_beta2 * vcov_cont["Latitude", "Latitude_2"] +
    2 * d_beta1 * d_beta2 * vcov_cont["Latitude", cov_vars[4]] +
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], "Latitude_2"] +
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], cov_vars[4]]
  
  var_effect_climate[[cont]] <- var_effect
  se_effect_climate[[cont]] <- sqrt(var_effect)
  ci_effect_climate_lower[[cont]] <- effect_climate[[cont]] - 1.96 * se_effect_climate[[cont]]
  ci_effect_climate_upper[[cont]] <- effect_climate[[cont]] + 1.96 * se_effect_climate[[cont]]
}

# 计算边际效应
marginal_base <- beta1 + 2 * beta2 * P_values

marginal_climate <- lapply(climates, function(cont) {
  beta1_total[[cont]] + 2 * beta2_total[[cont]] * P_values
})
names(marginal_climate) <- climates

d_marg_beta1 <- 1
d_marg_beta2 <- 2 * P_values

var_marginal_base <- (d_marg_beta1^2) * vcov_base["Latitude", "Latitude"] +
  (d_marg_beta2^2) * vcov_base["Latitude_2", "Latitude_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_base["Latitude", "Latitude_2"]
se_marginal_base <- sqrt(var_marginal_base)
ci_marginal_base_lower <- marginal_base - 1.96 * se_marginal_base
ci_marginal_base_upper <- marginal_base + 1.96 * se_marginal_base

var_marginal_climate <- list()
se_marginal_climate <- list()
ci_marginal_climate_lower <- list()
ci_marginal_climate_upper <- list()

for (cont in climates) {
  cov_vars <- c("Latitude", "Latitude_2",
                paste0("Latitude:ClimateZone", cont),
                paste0("Latitude_2:ClimateZone", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  
  d_marg_beta1_u <- 1
  d_marg_beta2_u <- 2 * P_values
  
  var_marginal <- (d_marg_beta1_u^2) * vcov_cont["Latitude", "Latitude"] +
    (d_marg_beta2_u^2) * vcov_cont["Latitude_2", "Latitude_2"] +
    (d_marg_beta1_u^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_marg_beta2_u^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_marg_beta1_u^2 * vcov_cont["Latitude", cov_vars[3]] +
    2 * d_marg_beta2_u^2 * vcov_cont["Latitude_2", cov_vars[4]] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["Latitude", "Latitude_2"] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["Latitude", cov_vars[4]] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], "Latitude_2"] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], cov_vars[4]]
  
  var_marginal_climate[[cont]] <- var_marginal
  se_marginal_climate[[cont]] <- sqrt(var_marginal)
  ci_marginal_climate_lower[[cont]] <- marginal_climate[[cont]] - 1.96 * se_marginal_climate[[cont]]
  ci_marginal_climate_upper[[cont]] <- marginal_climate[[cont]] + 1.96 * se_marginal_climate[[cont]]
}

effect_df <- data.frame(
  Latitude = P_values,
  Climate = "基准组",
  Effect = effect_base,
  SE = se_effect_base,
  CI_Lower = ci_effect_base_lower,
  CI_Upper = ci_effect_base_upper
)

for (cont in climates) {
  temp_df <- data.frame(
    Latitude = P_values,
    Climate = cont,
    Effect = effect_climate[[cont]],
    SE = se_effect_climate[[cont]],
    CI_Lower = ci_effect_climate_lower[[cont]],
    CI_Upper = ci_effect_climate_upper[[cont]]
  )
  effect_df <- rbind(effect_df, temp_df)
}

marginal_df <- data.frame(
  Latitude = P_values,
  Climate = "基准组",
  Marginal_Effect = marginal_base,
  SE = se_marginal_base,
  CI_Lower = ci_marginal_base_lower,
  CI_Upper = ci_marginal_base_upper
)

for (cont in climates) {
  temp_df <- data.frame(
    Latitude = P_values,
    Climate = cont,
    Marginal_Effect = marginal_climate[[cont]],
    SE = se_marginal_climate[[cont]],
    CI_Lower = ci_marginal_climate_lower[[cont]],
    CI_Upper = ci_marginal_climate_upper[[cont]]
  )
  marginal_df <- rbind(marginal_df, temp_df)
}

effect_df$Climate <- factor(effect_df$Climate, 
                            levels = c("基准组", climates),
                            labels = c("Temperate", "Tropical", "Arid",  "Cold"))

marginal_df$Climate <- factor(marginal_df$Climate,
                              levels = c("基准组", climates),
                              labels = c("Temperate", "Tropical", "Arid",  "Cold"))

climate_limits <- tibble(
  Climate = c("Temperate", "Tropical", "Arid",  "Cold"),
  min_WD = c(
    min(urban_in_dt_Latitude[urban_in_dt_Latitude$ClimateZone == "Temperate", ]$Latitude, na.rm = TRUE),
    min(urban_in_dt_Latitude[urban_in_dt_Latitude$ClimateZone == "Tropical", ]$Latitude, na.rm = TRUE),
    min(urban_in_dt_Latitude[urban_in_dt_Latitude$ClimateZone == "Arid", ]$Latitude, na.rm = TRUE),
    min(urban_in_dt_Latitude[urban_in_dt_Latitude$ClimateZone == "Cold", ]$Latitude, na.rm = TRUE)
  ),
  max_WD = c(
    max(urban_in_dt_Latitude[urban_in_dt_Latitude$ClimateZone == "Temperate", ]$Latitude, na.rm = TRUE),
    max(urban_in_dt_Latitude[urban_in_dt_Latitude$ClimateZone == "Tropical", ]$Latitude, na.rm = TRUE),
    max(urban_in_dt_Latitude[urban_in_dt_Latitude$ClimateZone == "Arid", ]$Latitude, na.rm = TRUE),
    max(urban_in_dt_Latitude[urban_in_dt_Latitude$ClimateZone == "Cold", ]$Latitude, na.rm = TRUE)
  )
)

print(climate_limits)

effect_df <- effect_df %>%
  left_join(climate_limits, by = "Climate") %>%
  filter(Latitude >= min_WD & Latitude <= max_WD) %>%
  select(-min_WD, -max_WD)

marginal_df <- marginal_df %>%
  left_join(climate_limits, by = "Climate") %>%
  filter(Latitude >= min_WD & Latitude <= max_WD) %>%
  select(-min_WD, -max_WD)

# 绘制效应值图
# Figure 3f
effect_plot <- ggplot(effect_df, aes(x = Latitude, y = Effect, color = Climate)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Climate),
    alpha = 0.15, color = NA
  ) +
  theme_classic() +
  labs(
    x = "",
    y = "",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, by =10), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.1, 0.05),
                     breaks = seq(-0.1, 0.05, by = 0.05),
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.text.y = element_text(color = "black", size = 13),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),
  ) +
  scale_color_manual(values = climate_colors) +
  scale_fill_manual(values = climate_colors)

print(effect_plot)

# 绘制直方图
hist_plot <- ggplot(urban_in_dt_Latitude, aes(x = Latitude)) +
  geom_histogram(aes(y = ..density.., fill=ClimateZone), binwidth = 0.05, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, by =10)) +
  scale_y_continuous(limits = c(0, 0.25), breaks = seq(0, 0.25, by =0.1)) +
  labs(x = "",
       y = "") +
  theme_classic() +
  theme(
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text = element_text(color = "black", size = 13),
    axis.title = element_text(size = 14),
    legend.position = "none"
  ) + 
  scale_fill_manual(values = climate_colors)

print(hist_plot)

# 合并
effect_climate_latitude <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_climate_latitude)


# 绘制边际效应图 
# Figure S5f
marginal_plot <- ggplot(marginal_df, aes(x = Latitude, y = Marginal_Effect, color = Climate)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Climate),
    alpha = 0.15, color = NA
  ) +
  labs(
    x = "",
    y = "",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  theme_classic() +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, by =10), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.01, 0.005),
                     breaks = seq(-0.01, 0.005, by = 0.005),
                     labels = scales::number_format(accuracy = 0.001)) +
  theme(
    axis.text.y = element_text(color = "black", size = 13), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),  
    axis.title.y = element_text(size = 14),  
    legend.position = "none", 
    legend.direction = "vertical", 
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'), 
  ) +
  scale_color_manual(values = climate_colors) +
  scale_fill_manual(values = climate_colors)

print(marginal_plot)

marginal_climate_latitude <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_climate_latitude)

################################################### 城市内异质性分析 ################################################## 
##################################################### development level ################################################## 
################################################### elevation ######################################
cols_to_center <- c('WaterDeficit', 'Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
urban_in_dt_Elevation <- copy(urban_in_dt)
urban_in_dt_Elevation <- urban_in_dt_Elevation[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                               .SDcols = cols_to_center]
urban_in_dt_Elevation$DevelopmentLevel <- factor(urban_in_dt_Elevation$DevelopmentLevel, levels = c("developing", "developed"))

# linear regression
# Table S5 -- model 3
model_urban <- felm(FVC ~ 1 + WaterDeficit +
                      Seasonalityindex +
                      AverageTemperature +
                      TemperatureRange +
                      WindSpeed +
                      Elevation +
                      # Elevation_2 +
                      Latitude +
                      Precipitation +
                      SoilMoisture +
                      SoilPH + 
                      Longitude + 
                      
                      Elevation:DevelopmentLevel +
                      # Elevation_2:DevelopmentLevel +
                      
                      GDP_per_capita_PPP + HDI + Population + 
                      ImperviousSurface + HumanSettlement + CityArea
                    |  year + CountryID, data = urban_in_dt_Elevation)

summary(model_urban)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_level.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

# nonlinear regression
# Table S5 -- model 7
model_urban <- felm(FVC ~ 1 + WaterDeficit +
                      Seasonalityindex +
                      AverageTemperature +
                      TemperatureRange +
                      WindSpeed +
                      Elevation +
                      Elevation_2 +
                      Latitude +
                      Precipitation +
                      SoilMoisture +
                      SoilPH + 
                      Longitude + 
                      
                      Elevation:DevelopmentLevel +
                      Elevation_2:DevelopmentLevel +
                      
                      GDP_per_capita_PPP + HDI + Population + 
                      ImperviousSurface + HumanSettlement + CityArea
                    |  year + CountryID, data = urban_in_dt_Elevation)

summary(model_urban)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_level2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

######################################### 绘制效应图与边际效应图 #######################################
coeffs <- coef(model_urban)
beta1 <- coeffs["Elevation"]
beta2 <- coeffs["Elevation_2"]
beta1_city <- coeffs["Elevation:DevelopmentLeveldeveloped"]
beta2_city <- coeffs["Elevation_2:DevelopmentLeveldeveloped"]
beta1_total <- beta1 + beta1_city
beta2_total <- beta2 + beta2_city

P_min <- quantile(urban_in_dt_Elevation$Elevation, 0.001, na.rm = TRUE)
P_max <- quantile(urban_in_dt_Elevation$Elevation, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)

# 计算效应值
effect_non_urban <- beta1 * P_values + beta2 * (P_values^2)
effect_urban <- beta1_total * P_values + beta2_total * (P_values^2)

vcov_matrix <- vcovHC(model_urban, type = "HC1", cluster = "group", group = urban_in_dt_Elevation$CountryID)

cov_vars_non_urban <- c("Elevation", "Elevation_2")
vcov_non_urban <- vcov_matrix[cov_vars_non_urban, cov_vars_non_urban]
d_beta1 <- P_values
d_beta2 <- P_values^2
var_effect_non_urban <- (d_beta1^2) * vcov_non_urban["Elevation", "Elevation"] +
  (d_beta2^2) * vcov_non_urban["Elevation_2", "Elevation_2"] +
  2 * d_beta1 * d_beta2 * vcov_non_urban["Elevation", "Elevation_2"]
se_effect_non_urban <- sqrt(var_effect_non_urban)
ci_effect_non_urban_lower <- effect_non_urban - 1.96 * se_effect_non_urban
ci_effect_non_urban_upper <- effect_non_urban + 1.96 * se_effect_non_urban

cov_vars_urban <- c("Elevation", "Elevation_2", "Elevation:DevelopmentLeveldeveloped", "Elevation_2:DevelopmentLeveldeveloped")
vcov_urban <- vcov_matrix[cov_vars_urban, cov_vars_urban]
d_beta1_u <- P_values
d_beta2_u <- P_values^2
var_effect_urban <- (d_beta1_u^2) * vcov_urban["Elevation", "Elevation"] +
  (d_beta2_u^2) * vcov_urban["Elevation_2", "Elevation_2"] +
  (d_beta1_u^2) * vcov_urban["Elevation:DevelopmentLeveldeveloped", "Elevation:DevelopmentLeveldeveloped"] +
  (d_beta2_u^2) * vcov_urban["Elevation_2:DevelopmentLeveldeveloped", "Elevation_2:DevelopmentLeveldeveloped"] +
  2 * d_beta1_u^2 * vcov_urban["Elevation", "Elevation:DevelopmentLeveldeveloped"] +
  2 * d_beta2_u^2 * vcov_urban["Elevation_2", "Elevation_2:DevelopmentLeveldeveloped"] +
  2 * d_beta1_u * d_beta2_u * vcov_urban["Elevation", "Elevation_2"] +
  2 * d_beta1_u * d_beta2_u * vcov_urban["Elevation", "Elevation_2:DevelopmentLeveldeveloped"] +
  2 * d_beta1_u * d_beta2_u * vcov_urban["Elevation:DevelopmentLeveldeveloped", "Elevation_2"] +
  2 * d_beta1_u * d_beta2_u * vcov_urban["Elevation:DevelopmentLeveldeveloped", "Elevation_2:DevelopmentLeveldeveloped"]
se_effect_urban <- sqrt(var_effect_urban)
ci_effect_urban_lower <- effect_urban - 1.96 * se_effect_urban
ci_effect_urban_upper <- effect_urban + 1.96 * se_effect_urban

# 计算边际效应
marginal_non_urban <- beta1 + 2 * beta2 * P_values
marginal_urban <- beta1_total + 2 * beta2_total * P_values

d_marg_beta1 <- 1
d_marg_beta2 <- 2 * P_values
var_marginal_non_urban <- (d_marg_beta1^2) * vcov_non_urban["Elevation", "Elevation"] +
  (d_marg_beta2^2) * vcov_non_urban["Elevation_2", "Elevation_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_non_urban["Elevation", "Elevation_2"]
se_marginal_non_urban <- sqrt(var_marginal_non_urban)
ci_marginal_non_urban_lower <- marginal_non_urban - 1.96 * se_marginal_non_urban
ci_marginal_non_urban_upper <- marginal_non_urban + 1.96 * se_marginal_non_urban

d_marg_beta1_u <- 1
d_marg_beta2_u <- 2 * P_values
var_marginal_urban <- (d_marg_beta1_u^2) * vcov_urban["Elevation", "Elevation"] +
  (d_marg_beta2_u^2) * vcov_urban["Elevation_2", "Elevation_2"] +
  (d_marg_beta1_u^2) * vcov_urban["Elevation:DevelopmentLeveldeveloped", "Elevation:DevelopmentLeveldeveloped"] +
  (d_marg_beta2_u^2) * vcov_urban["Elevation_2:DevelopmentLeveldeveloped", "Elevation_2:DevelopmentLeveldeveloped"] +
  2 * d_marg_beta1_u^2 * vcov_urban["Elevation", "Elevation:DevelopmentLeveldeveloped"] +
  2 * d_marg_beta2_u^2 * vcov_urban["Elevation_2", "Elevation_2:DevelopmentLeveldeveloped"] +
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Elevation", "Elevation_2"] +
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Elevation", "Elevation_2:DevelopmentLeveldeveloped"] +
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Elevation:DevelopmentLeveldeveloped", "Elevation_2"] +
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Elevation:DevelopmentLeveldeveloped", "Elevation_2:DevelopmentLeveldeveloped"]
se_marginal_urban <- sqrt(var_marginal_urban)
ci_marginal_urban_lower <- marginal_urban - 1.96 * se_marginal_urban
ci_marginal_urban_upper <- marginal_urban + 1.96 * se_marginal_urban

effect_df <- rbind(
  data.frame(
    Elevation = P_values,
    Region = "Developing",
    Effect = effect_non_urban,
    SE = se_effect_non_urban,
    CI_Lower = ci_effect_non_urban_lower,
    CI_Upper = ci_effect_non_urban_upper
  ),
  data.frame(
    Elevation = P_values,
    Region = "Developed",
    Effect = effect_urban,
    SE = se_effect_urban,
    CI_Lower = ci_effect_urban_lower,
    CI_Upper = ci_effect_urban_upper
  )
)

marginal_df <- rbind(
  data.frame(
    Elevation = P_values,
    Region = "Developing",
    Marginal_Effect = marginal_non_urban,
    SE = se_marginal_non_urban,
    CI_Lower = ci_marginal_non_urban_lower,
    CI_Upper = ci_marginal_non_urban_upper
  ),
  data.frame(
    Elevation = P_values,
    Region = "Developed",
    Marginal_Effect = marginal_urban,
    SE = se_marginal_urban,
    CI_Lower = ci_marginal_urban_lower,
    CI_Upper = ci_marginal_urban_upper
  )
)

level_limits <- tibble(
  Region = c("Developing", "Developed"),
  min_WD = c(
    min(urban_in_dt_Elevation[urban_in_dt_Elevation$DevelopmentLevel == "developing", ]$Elevation, na.rm = TRUE),
    min(urban_in_dt_Elevation[urban_in_dt_Elevation$DevelopmentLevel == "developed", ]$Elevation, na.rm = TRUE)
  ),
  max_WD = c(
    max(urban_in_dt_Elevation[urban_in_dt_Elevation$DevelopmentLevel == "developing", ]$Elevation, na.rm = TRUE),
    max(urban_in_dt_Elevation[urban_in_dt_Elevation$DevelopmentLevel == "developed", ]$Elevation, na.rm = TRUE)
  )
)

print(level_limits)

effect_df <- effect_df %>%
  left_join(level_limits, by = "Region") %>%
  filter(Elevation >= min_WD & Elevation <= max_WD) %>%
  select(-min_WD, -max_WD)

marginal_df <- marginal_df %>%
  left_join(level_limits, by = "Region") %>%
  filter(Elevation >= min_WD & Elevation <= max_WD) %>%
  select(-min_WD, -max_WD)

# 绘制效应值图
# Figure 3g
level_colors <- c(
  "Developed" = "#0071bc", 
  "Developing" = "#d95218"
)

effect_plot <- ggplot(effect_df, aes(x = Elevation, y = Effect, color = Region)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha = 0.15, color = NA
  ) +
  theme_classic() +
  labs(
    x = "Water deficit",
    y = "Effect on urban FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 4050), breaks = seq(0, 4050, by =1000), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.1, 0.05),
                     breaks = seq(-0.1, 0.05, by = 0.05),
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.text.y = element_text(color = "black", size = 13),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.position = c(0.3,0.35), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),
    legend.background = element_blank()
  ) +
  scale_color_manual(values = level_colors) +
  scale_fill_manual(values = level_colors) 

print(effect_plot)

# 绘制直方图
hist_plot <- ggplot(urban_in_dt_Elevation, aes(x = Elevation)) +
  geom_histogram(aes(y = ..density.., fill = DevelopmentLevel), binwidth =  5,  alpha = 0.5) +
  scale_x_continuous(limits = c(0, 4050), breaks = seq(0, 4050, by =1000)) +
  scale_y_continuous(limits = c(0, 0.025), breaks = seq(0, 0.025, by =0.01)) +
  labs(x = "",
       y = "",
       subtitle = paste0("")) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 14),
    plot.margin = margin(t=-1.32,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text.y = element_text(color = "black", size = 13),
    axis.text.x = element_text(color = "black", size = 13),
    legend.position = "none") +
  scale_fill_manual(values = c(  "developed" = "#0071bc", 
                                 "developing" = "#d95218"))

print(hist_plot)

# 合并
effect_level_elevation <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(effect_level_elevation)


# 绘制边际效应图 
# Figure S5g
marginal_plot <- ggplot(marginal_df, aes(x = Elevation, y = Marginal_Effect, color = Region)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha = 0.15, color = NA
  ) +
  labs(
    x = "Elevation (m)",
    y = "Marginal effect on urban FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  theme_classic() +
  scale_x_continuous(limits = c(0, 4000), breaks = seq(0, 4000, by =1000)) +
  scale_y_continuous(limits = c(-0.0001, 0.0001),
                     breaks = seq(-0.0001, 0.0001, by = 0.0001),
                     labels = scales::number_format(accuracy = 0.0001)) +
  
  theme(
    axis.text.y = element_text(color = "black", size = 12),  
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),  
    axis.title.y = element_text(size = 12),  
    legend.position = c(0.25,0.3), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'), 
    legend.background = element_blank()
  ) +
  scale_color_manual(values = level_colors) +
  scale_fill_manual(values = level_colors)

print(marginal_plot)

marginal_level_elevation <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(marginal_level_elevation)


##################################################### level ####################################################
################################################### latitude ######################################
cols_to_center <- c('WaterDeficit', 'Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
urban_in_dt_Latitude <- copy(urban_in_dt)
urban_in_dt_Latitude <- urban_in_dt_Latitude[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                             .SDcols = cols_to_center]
urban_in_dt_Latitude$DevelopmentLevel <- factor(urban_in_dt_Latitude$DevelopmentLevel, levels = c("developing", "developed"))

# linear regression
# Table S5 -- model 4
model_urban <- felm(FVC ~ 1 + WaterDeficit +
                      Seasonalityindex +
                      AverageTemperature +
                      TemperatureRange +
                      WindSpeed +
                      Elevation +
                      Latitude +
                      # Latitude_2 +
                      Precipitation +
                      SoilMoisture +
                      SoilPH + 
                      Longitude + 
                      
                      Latitude:DevelopmentLevel +
                      # Latitude_2:DevelopmentLevel +
                      
                      GDP_per_capita_PPP + HDI + Population + 
                      ImperviousSurface + HumanSettlement + CityArea
                    |  year + CountryID, data = urban_in_dt_Latitude)

summary(model_urban)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_level.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

# nonlinear regression
# Table S5 -- model 8
model_urban <- felm(FVC ~ 1 + WaterDeficit +
                      Seasonalityindex +
                      AverageTemperature +
                      TemperatureRange +
                      WindSpeed +
                      Elevation +
                      Latitude +
                      Latitude_2 +
                      Precipitation +
                      SoilMoisture +
                      SoilPH + 
                      Longitude + 
                      
                      Latitude:DevelopmentLevel +
                      Latitude_2:DevelopmentLevel +
                      
                      GDP_per_capita_PPP + HDI + Population + 
                      ImperviousSurface + HumanSettlement + CityArea
                    |  year + CountryID, data = urban_in_dt_Latitude)

summary(model_urban)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_level2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  
######################################### 绘制效应图与边际效应图 #######################################
coeffs <- coef(model_urban)
beta1 <- coeffs["Latitude"]
beta2 <- coeffs["Latitude_2"]
beta1_city <- coeffs["Latitude:DevelopmentLeveldeveloped"]
beta2_city <- coeffs["Latitude_2:DevelopmentLeveldeveloped"]
beta1_total <- beta1 + beta1_city
beta2_total <- beta2 + beta2_city

P_min <- quantile(urban_in_dt_Latitude$Latitude, 0.001, na.rm = TRUE)
P_max <- quantile(urban_in_dt_Latitude$Latitude, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)

# 计算效应值
effect_non_urban <- beta1 * P_values + beta2 * (P_values^2)
effect_urban <- beta1_total * P_values + beta2_total * (P_values^2)

vcov_matrix <- vcovHC(model_urban, type = "HC1", cluster = "group", group = urban_in_dt_Latitude$CountryID)

cov_vars_non_urban <- c("Latitude", "Latitude_2")
vcov_non_urban <- vcov_matrix[cov_vars_non_urban, cov_vars_non_urban]
d_beta1 <- P_values
d_beta2 <- P_values^2
var_effect_non_urban <- (d_beta1^2) * vcov_non_urban["Latitude", "Latitude"] +
  (d_beta2^2) * vcov_non_urban["Latitude_2", "Latitude_2"] +
  2 * d_beta1 * d_beta2 * vcov_non_urban["Latitude", "Latitude_2"]
se_effect_non_urban <- sqrt(var_effect_non_urban)
ci_effect_non_urban_lower <- effect_non_urban - 1.96 * se_effect_non_urban
ci_effect_non_urban_upper <- effect_non_urban + 1.96 * se_effect_non_urban

cov_vars_urban <- c("Latitude", "Latitude_2", "Latitude:DevelopmentLeveldeveloped", "Latitude_2:DevelopmentLeveldeveloped")
vcov_urban <- vcov_matrix[cov_vars_urban, cov_vars_urban]
d_beta1_u <- P_values
d_beta2_u <- P_values^2
var_effect_urban <- (d_beta1_u^2) * vcov_urban["Latitude", "Latitude"] +
  (d_beta2_u^2) * vcov_urban["Latitude_2", "Latitude_2"] +
  (d_beta1_u^2) * vcov_urban["Latitude:DevelopmentLeveldeveloped", "Latitude:DevelopmentLeveldeveloped"] +
  (d_beta2_u^2) * vcov_urban["Latitude_2:DevelopmentLeveldeveloped", "Latitude_2:DevelopmentLeveldeveloped"] +
  2 * d_beta1_u^2 * vcov_urban["Latitude", "Latitude:DevelopmentLeveldeveloped"] +
  2 * d_beta2_u^2 * vcov_urban["Latitude_2", "Latitude_2:DevelopmentLeveldeveloped"] +
  2 * d_beta1_u * d_beta2_u * vcov_urban["Latitude", "Latitude_2"] +
  2 * d_beta1_u * d_beta2_u * vcov_urban["Latitude", "Latitude_2:DevelopmentLeveldeveloped"] +
  2 * d_beta1_u * d_beta2_u * vcov_urban["Latitude:DevelopmentLeveldeveloped", "Latitude_2"] +
  2 * d_beta1_u * d_beta2_u * vcov_urban["Latitude:DevelopmentLeveldeveloped", "Latitude_2:DevelopmentLeveldeveloped"]
se_effect_urban <- sqrt(var_effect_urban)
ci_effect_urban_lower <- effect_urban - 1.96 * se_effect_urban
ci_effect_urban_upper <- effect_urban + 1.96 * se_effect_urban

# 计算边际效应
marginal_non_urban <- beta1 + 2 * beta2 * P_values
marginal_urban <- beta1_total + 2 * beta2_total * P_values

d_marg_beta1 <- 1
d_marg_beta2 <- 2 * P_values
var_marginal_non_urban <- (d_marg_beta1^2) * vcov_non_urban["Latitude", "Latitude"] +
  (d_marg_beta2^2) * vcov_non_urban["Latitude_2", "Latitude_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_non_urban["Latitude", "Latitude_2"]
se_marginal_non_urban <- sqrt(var_marginal_non_urban)
ci_marginal_non_urban_lower <- marginal_non_urban - 1.96 * se_marginal_non_urban
ci_marginal_non_urban_upper <- marginal_non_urban + 1.96 * se_marginal_non_urban

d_marg_beta1_u <- 1
d_marg_beta2_u <- 2 * P_values
var_marginal_urban <- (d_marg_beta1_u^2) * vcov_urban["Latitude", "Latitude"] +
  (d_marg_beta2_u^2) * vcov_urban["Latitude_2", "Latitude_2"] +
  (d_marg_beta1_u^2) * vcov_urban["Latitude:DevelopmentLeveldeveloped", "Latitude:DevelopmentLeveldeveloped"] +
  (d_marg_beta2_u^2) * vcov_urban["Latitude_2:DevelopmentLeveldeveloped", "Latitude_2:DevelopmentLeveldeveloped"] +
  2 * d_marg_beta1_u^2 * vcov_urban["Latitude", "Latitude:DevelopmentLeveldeveloped"] +
  2 * d_marg_beta2_u^2 * vcov_urban["Latitude_2", "Latitude_2:DevelopmentLeveldeveloped"] +
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Latitude", "Latitude_2"] +
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Latitude", "Latitude_2:DevelopmentLeveldeveloped"] +
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Latitude:DevelopmentLeveldeveloped", "Latitude_2"] +
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Latitude:DevelopmentLeveldeveloped", "Latitude_2:DevelopmentLeveldeveloped"]
se_marginal_urban <- sqrt(var_marginal_urban)
ci_marginal_urban_lower <- marginal_urban - 1.96 * se_marginal_urban
ci_marginal_urban_upper <- marginal_urban + 1.96 * se_marginal_urban

effect_df <- rbind(
  data.frame(
    Latitude = P_values,
    Region = "Developing",
    Effect = effect_non_urban,
    SE = se_effect_non_urban,
    CI_Lower = ci_effect_non_urban_lower,
    CI_Upper = ci_effect_non_urban_upper
  ),
  data.frame(
    Latitude = P_values,
    Region = "Developed",
    Effect = effect_urban,
    SE = se_effect_urban,
    CI_Lower = ci_effect_urban_lower,
    CI_Upper = ci_effect_urban_upper
  )
)

marginal_df <- rbind(
  data.frame(
    Latitude = P_values,
    Region = "Developing",
    Marginal_Effect = marginal_non_urban,
    SE = se_marginal_non_urban,
    CI_Lower = ci_marginal_non_urban_lower,
    CI_Upper = ci_marginal_non_urban_upper
  ),
  data.frame(
    Latitude = P_values,
    Region = "Developed",
    Marginal_Effect = marginal_urban,
    SE = se_marginal_urban,
    CI_Lower = ci_marginal_urban_lower,
    CI_Upper = ci_marginal_urban_upper
  )
)

level_limits <- tibble(
  Region = c("Developing", "Developed"),
  min_WD = c(
    min(urban_in_dt_Latitude[urban_in_dt_Latitude$DevelopmentLevel == "developing", ]$Latitude, na.rm = TRUE),
    min(urban_in_dt_Latitude[urban_in_dt_Latitude$DevelopmentLevel == "developed", ]$Latitude, na.rm = TRUE)
  ),
  max_WD = c(
    max(urban_in_dt_Latitude[urban_in_dt_Latitude$DevelopmentLevel == "developing", ]$Latitude, na.rm = TRUE),
    max(urban_in_dt_Latitude[urban_in_dt_Latitude$DevelopmentLevel == "developed", ]$Latitude, na.rm = TRUE)
  )
)

print(level_limits)

effect_df <- effect_df %>%
  left_join(level_limits, by = "Region") %>%
  filter(Latitude >= min_WD & Latitude <= max_WD) %>%
  select(-min_WD, -max_WD)

marginal_df <- marginal_df %>%
  left_join(level_limits, by = "Region") %>%
  filter(Latitude >= min_WD & Latitude <= max_WD) %>%
  select(-min_WD, -max_WD)

# 绘制效应值图
# Figure 3h
effect_plot <- ggplot(effect_df, aes(x = Latitude, y = Effect, color = Region)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha = 0.15, color = NA
  ) +
  theme_classic() +
  labs(
    x = "",
    y = "",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, by =10), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.08, 0.04),
                     breaks = seq(-0.08, 0.04, by = 0.04),
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.text.y = element_text(color = "black", size = 13),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.position = "none", 
    legend.direction = "vertical", 
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),
  ) +
  scale_color_manual(values = level_colors) +
  scale_fill_manual(values = level_colors) 

print(effect_plot)

# 绘制直方图
hist_plot <- ggplot(urban_in_dt_Latitude, aes(x = Latitude)) +
  geom_histogram(aes(y = ..density.., fill = DevelopmentLevel), binwidth =  0.05,  alpha = 0.5) +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, by =10)) +
  scale_y_continuous(limits = c(0, 0.14), breaks = seq(0, 0.14, by =0.06)) +
  labs(x = "",
       y = "",
       subtitle = paste0("")) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 14),
    plot.margin = margin(t=-1.32,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text.y = element_text(color = "black", size = 13),
    axis.text.x = element_text(color = "black", size = 13),
    legend.position = "none") +
  scale_fill_manual(values = c(  "developed" = "#0071bc", 
                                 "developing" = "#d95218"))

print(hist_plot)

# 合并
effect_level_latitude <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(effect_level_latitude)

# 绘制边际效应图 
# Figure S5h
marginal_plot <- ggplot(marginal_df, aes(x = Latitude, y = Marginal_Effect, color = Region)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha = 0.15, color = NA
  ) +
  labs(
    x = "",
    y = "",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  theme_classic() +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, by =10)) +
  scale_y_continuous(limits = c(-0.005, 0.005),
                     breaks = seq(-0.005, 0.005, by = 0.005),
                     labels = scales::number_format(accuracy = 0.001)) +
  
  theme(
    axis.text.y = element_text(color = "black", size = 13),  
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),  
    axis.title.y = element_text(size = 14),  
    legend.position = "none", 
    legend.direction = "vertical", 
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'), 
  ) +
  scale_color_manual(values = level_colors) +
  scale_fill_manual(values = level_colors)

print(marginal_plot)

marginal_level_latitude <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(marginal_level_latitude)

################################################### 城市内异质性分析 ################################################## 
##################################################### expansion rate ################################################## 
################################################### elevation ######################################
cols_to_center <- c('WaterDeficit', 'Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
urban_in_dt_Elevation <- copy(urban_in_dt)
urban_in_dt_Elevation <- urban_in_dt_Elevation[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                               .SDcols = cols_to_center]
urban_in_dt_Elevation$expansionrate <- factor(urban_in_dt_Elevation$expansionrate, levels = c( "Slow",   "Moderate",    "Rapid",    "Intense"))

# linear regression
# Table S6 -- model 3
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         Seasonalityindex +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         # Elevation_2 +
                         Latitude +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         Elevation:expansionrate +
                         # Elevation_2:expansionrate +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea
                       | year + CountryID, data = urban_in_dt_Elevation)

summary(model_urban_in)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_rate.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

# nonlinear regression
# Table S6 -- model 7
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         Seasonalityindex +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Elevation_2 +
                         Latitude +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         Elevation:expansionrate +
                         Elevation_2:expansionrate +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea
                       | year + CountryID, data = urban_in_dt_Elevation)

summary(model_urban_in)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_rate2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  


######################################### 绘制效应图与边际效应图 #######################################
beta1 <- coef(model_urban_in)["Elevation"]
beta2 <- coef(model_urban_in)["Elevation_2"]

rates <- c( "Moderate",    "Rapid",    "Intense")

beta1_rate <- sapply(rates, function(cont) {
  coef(model_urban_in)[paste0("Elevation:expansionrate", cont)]
})
beta2_rate <- sapply(rates, function(cont) {
  coef(model_urban_in)[paste0("Elevation_2:expansionrate", cont)]
})

beta1_total <- beta1 + beta1_rate
beta2_total <- beta2 + beta2_rate

names(beta1_total) <- rates
names(beta2_total) <- rates

P_min <- quantile(urban_in_dt_Elevation$Elevation, 0.001, na.rm = TRUE)
P_max <- quantile(urban_in_dt_Elevation$Elevation, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)

# 计算效应值
effect_base <- beta1 * P_values + beta2 * (P_values^2)

effect_rate <- lapply(rates, function(cont) {
  beta1_total[[cont]] * P_values + beta2_total[[cont]] * (P_values^2)
})
names(effect_rate) <- rates

vcov_matrix <- vcovHC(model_urban_in, type = "HC1", cluster = "group", 
                      group = urban_in_dt_Elevation$CountryID)

cov_vars_base <- c("Elevation", "Elevation_2")
vcov_base <- vcov_matrix[cov_vars_base, cov_vars_base]

d_beta1 <- P_values
d_beta2 <- P_values^2

var_effect_base <- (d_beta1^2) * vcov_base["Elevation", "Elevation"] +
  (d_beta2^2) * vcov_base["Elevation_2", "Elevation_2"] +
  2 * d_beta1 * d_beta2 * vcov_base["Elevation", "Elevation_2"]
se_effect_base <- sqrt(var_effect_base)
ci_effect_base_lower <- effect_base - 1.96 * se_effect_base
ci_effect_base_upper <- effect_base + 1.96 * se_effect_base

var_effect_rate <- list()
se_effect_rate <- list()
ci_effect_rate_lower <- list()
ci_effect_rate_upper <- list()

for (cont in rates) {
  cov_vars <- c("Elevation", "Elevation_2", 
                paste0("Elevation:expansionrate", cont), 
                paste0("Elevation_2:expansionrate", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  
  var_effect <- (d_beta1^2) * vcov_cont["Elevation", "Elevation"] +
    (d_beta2^2) * vcov_cont["Elevation_2", "Elevation_2"] +
    (d_beta1^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_beta2^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_beta1^2 * vcov_cont["Elevation", cov_vars[3]] +
    2 * d_beta2^2 * vcov_cont["Elevation_2", cov_vars[4]] +
    2 * d_beta1 * d_beta2 * vcov_cont["Elevation", "Elevation_2"] +
    2 * d_beta1 * d_beta2 * vcov_cont["Elevation", cov_vars[4]] +
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], "Elevation_2"] +
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], cov_vars[4]]
  
  var_effect_rate[[cont]] <- var_effect
  se_effect_rate[[cont]] <- sqrt(var_effect)
  ci_effect_rate_lower[[cont]] <- effect_rate[[cont]] - 1.96 * se_effect_rate[[cont]]
  ci_effect_rate_upper[[cont]] <- effect_rate[[cont]] + 1.96 * se_effect_rate[[cont]]
}

# 计算边际效应
marginal_base <- beta1 + 2 * beta2 * P_values

marginal_rate <- lapply(rates, function(cont) {
  beta1_total[[cont]] + 2 * beta2_total[[cont]] * P_values
})
names(marginal_rate) <- rates

d_marg_beta1 <- 1
d_marg_beta2 <- 2 * P_values

var_marginal_base <- (d_marg_beta1^2) * vcov_base["Elevation", "Elevation"] +
  (d_marg_beta2^2) * vcov_base["Elevation_2", "Elevation_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_base["Elevation", "Elevation_2"]
se_marginal_base <- sqrt(var_marginal_base)
ci_marginal_base_lower <- marginal_base - 1.96 * se_marginal_base
ci_marginal_base_upper <- marginal_base + 1.96 * se_marginal_base

var_marginal_rate <- list()
se_marginal_rate <- list()
ci_marginal_rate_lower <- list()
ci_marginal_rate_upper <- list()

for (cont in rates) {
  cov_vars <- c("Elevation", "Elevation_2",
                paste0("Elevation:expansionrate", cont),
                paste0("Elevation_2:expansionrate", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  
  d_marg_beta1_u <- 1
  d_marg_beta2_u <- 2 * P_values
  
  var_marginal <- (d_marg_beta1_u^2) * vcov_cont["Elevation", "Elevation"] +
    (d_marg_beta2_u^2) * vcov_cont["Elevation_2", "Elevation_2"] +
    (d_marg_beta1_u^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_marg_beta2_u^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_marg_beta1_u^2 * vcov_cont["Elevation", cov_vars[3]] +
    2 * d_marg_beta2_u^2 * vcov_cont["Elevation_2", cov_vars[4]] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["Elevation", "Elevation_2"] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["Elevation", cov_vars[4]] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], "Elevation_2"] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], cov_vars[4]]
  
  var_marginal_rate[[cont]] <- var_marginal
  se_marginal_rate[[cont]] <- sqrt(var_marginal)
  ci_marginal_rate_lower[[cont]] <- marginal_rate[[cont]] - 1.96 * se_marginal_rate[[cont]]
  ci_marginal_rate_upper[[cont]] <- marginal_rate[[cont]] + 1.96 * se_marginal_rate[[cont]]
}

effect_df <- data.frame(
  Elevation = P_values,
  rate = "基准组",
  Effect = effect_base,
  SE = se_effect_base,
  CI_Lower = ci_effect_base_lower,
  CI_Upper = ci_effect_base_upper
)

for (cont in rates) {
  temp_df <- data.frame(
    Elevation = P_values,
    rate = cont,
    Effect = effect_rate[[cont]],
    SE = se_effect_rate[[cont]],
    CI_Lower = ci_effect_rate_lower[[cont]],
    CI_Upper = ci_effect_rate_upper[[cont]]
  )
  effect_df <- rbind(effect_df, temp_df)
}

marginal_df <- data.frame(
  Elevation = P_values,
  rate = "基准组",
  Marginal_Effect = marginal_base,
  SE = se_marginal_base,
  CI_Lower = ci_marginal_base_lower,
  CI_Upper = ci_marginal_base_upper
)

for (cont in rates) {
  temp_df <- data.frame(
    Elevation = P_values,
    rate = cont,
    Marginal_Effect = marginal_rate[[cont]],
    SE = se_marginal_rate[[cont]],
    CI_Lower = ci_marginal_rate_lower[[cont]],
    CI_Upper = ci_marginal_rate_upper[[cont]]
  )
  marginal_df <- rbind(marginal_df, temp_df)
}

effect_df$rate <- factor(effect_df$rate, 
                         levels = c("基准组", rates),
                         labels = c( "Low",   "Middle",    "High",    "Very high"))

marginal_df$rate <- factor(marginal_df$rate,
                           levels = c("基准组", rates),
                           labels = c( "Low",   "Middle",    "High",    "Very high"))

rate_limits <- tibble(
  rate = c("Low",   "Middle",    "High",    "Very high"),
  min_WD = c(
    min(urban_in_dt_Elevation[urban_in_dt_Elevation$expansionrate == "Slow", ]$Elevation, na.rm = TRUE),
    min(urban_in_dt_Elevation[urban_in_dt_Elevation$expansionrate == "Moderate", ]$Elevation, na.rm = TRUE),
    min(urban_in_dt_Elevation[urban_in_dt_Elevation$expansionrate == "Rapid", ]$Elevation, na.rm = TRUE),
    min(urban_in_dt_Elevation[urban_in_dt_Elevation$expansionrate == "Intense", ]$Elevation, na.rm = TRUE)
  ),
  max_WD = c(
    max(urban_in_dt_Elevation[urban_in_dt_Elevation$expansionrate == "Slow", ]$Elevation, na.rm = TRUE),
    max(urban_in_dt_Elevation[urban_in_dt_Elevation$expansionrate == "Moderate", ]$Elevation, na.rm = TRUE),
    max(urban_in_dt_Elevation[urban_in_dt_Elevation$expansionrate == "Rapid", ]$Elevation, na.rm = TRUE),
    max(urban_in_dt_Elevation[urban_in_dt_Elevation$expansionrate == "Intense", ]$Elevation, na.rm = TRUE)
  )
)

print(rate_limits)

effect_df <- effect_df %>%
  left_join(rate_limits, by = "rate") %>%
  filter(Elevation >= min_WD & Elevation <= max_WD) %>%
  select(-min_WD, -max_WD)

marginal_df <- marginal_df %>%
  left_join(rate_limits, by = "rate") %>%
  filter(Elevation >= min_WD & Elevation <= max_WD) %>%
  select(-min_WD, -max_WD)

# 绘制效应值图
# Figure 3i
rate_colors <- c(
  "Low" = "#27ae60",
  "Middle" = "#3498db",
  "High" = "#e67e22",
  "Very high" = "#c0392b"
)

effect_df$rate <- factor(effect_df$rate, levels = c( "Low",   "Middle",    "High",    "Very high"))

effect_plot <- ggplot(effect_df, aes(x = Elevation, y = Effect, color = rate)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = rate),
    alpha = 0.15, color = NA
  ) +
  theme_classic() +
  labs(
    x = "",
    y = "Effect on urban FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 4050), breaks = seq(0, 4050, by =1000), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.08, 0.04),
                     breaks = seq(-0.08, 0.04, by = 0.04),
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.text.y = element_text(color = "black", size = 13),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.3, 0.35),
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),
    legend.background = element_blank()
  ) +
  scale_color_manual(values = rate_colors) +
  scale_fill_manual(values = rate_colors)

print(effect_plot)

# 绘制直方图
rate_colors_hist <- c(
  "Slow" = "#27ae60",
  "Moderate" = "#3498db",
  "Rapid" = "#e67e22",
  "Intense" = "#c0392b"
)

hist_plot <- ggplot(urban_in_dt_Elevation, aes(x = Elevation)) +
  geom_histogram(aes(y = ..density.., fill=expansionrate), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 4050), breaks = seq(0, 4050, by =1000)) +
  scale_y_continuous(limits = c(0, 0.07), breaks = seq(0, 0.07, by =0.03)) +
  labs(x = "",
       y = "") +
  theme_classic() +
  theme(
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text = element_text(color = "black", size = 13),
    axis.title = element_text(size = 14),
    legend.position = "none"
  ) + 
  scale_fill_manual(values = rate_colors_hist)

print(hist_plot)

# 合并
effect_rate_elevation <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_rate_elevation)

# 绘制边际效应图
# Figure S5i
marginal_df$rate <- factor(marginal_df$rate, levels = c( "Low",   "Middle",    "High",    "Very high"))
marginal_plot <- ggplot(marginal_df, aes(x = Elevation, y = Marginal_Effect, color = rate)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = rate),
    alpha = 0.15, color = NA
  ) +
  labs(
    x = "",
    y = "Marginal effect on urban FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  theme_classic() +
  scale_x_continuous(limits = c(0, 4050), breaks = seq(0, 4050, by =1000), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.00008, 0.00004),
                     breaks = seq(-0.00008, 0.00004, by = 0.00004),
                     labels = scales::number_format(accuracy = 0.00001)) +
  theme(
    axis.text.y = element_text(color = "black", size = 12),  
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title.y = element_text(size = 12),  
    legend.position = c(0.3,0.3),    
    legend.direction = "vertical", 
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),  
    legend.background = element_blank()
  ) +
  scale_color_manual(values = rate_colors) +
  scale_fill_manual(values = rate_colors)

print(marginal_plot)

marginal_rate_elevation <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))

print(marginal_rate_elevation)

##################################################### expansion rate ####################################################
################################################### latitude ######################################
cols_to_center <- c('WaterDeficit', 'Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation',  
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
urban_in_dt_Latitude <- copy(urban_in_dt)
urban_in_dt_Latitude <- urban_in_dt_Latitude[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                             .SDcols = cols_to_center]
urban_in_dt_Latitude$expansionrate <- factor(urban_in_dt_Latitude$expansionrate, levels = c( "Slow",   "Moderate",    "Rapid",    "Intense"))

# linear regression
# Table S6 -- model 4
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         Seasonalityindex +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Latitude +
                         # Latitude_2 +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         Latitude:expansionrate +
                         # Latitude_2:expansionrate +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea
                       | year + CountryID, data = urban_in_dt_Latitude)

summary(model_urban_in)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_rate.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

# nonlinear regression
# Table S6 -- model 8
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         Seasonalityindex +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Latitude +
                         Latitude_2 +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         Latitude:expansionrate +
                         Latitude_2:expansionrate +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea
                       | year + CountryID, data = urban_in_dt_Latitude)

summary(model_urban_in)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_rate2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

######################################### 绘制效应图与边际效应图 #######################################
beta1 <- coef(model_urban_in)["Latitude"]
beta2 <- coef(model_urban_in)["Latitude_2"]

rates <- c( "Moderate",    "Rapid",    "Intense")

beta1_rate <- sapply(rates, function(cont) {
  coef(model_urban_in)[paste0("Latitude:expansionrate", cont)]
})
beta2_rate <- sapply(rates, function(cont) {
  coef(model_urban_in)[paste0("Latitude_2:expansionrate", cont)]
})

beta1_total <- beta1 + beta1_rate
beta2_total <- beta2 + beta2_rate

names(beta1_total) <- rates
names(beta2_total) <- rates

P_min <- quantile(urban_in_dt_Latitude$Latitude, 0.001, na.rm = TRUE)
P_max <- quantile(urban_in_dt_Latitude$Latitude, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)

# 计算效应值
effect_base <- beta1 * P_values + beta2 * (P_values^2)

effect_rate <- lapply(rates, function(cont) {
  beta1_total[[cont]] * P_values + beta2_total[[cont]] * (P_values^2)
})
names(effect_rate) <- rates

vcov_matrix <- vcovHC(model_urban_in, type = "HC1", cluster = "group", 
                      group = urban_in_dt_Latitude$CountryID)

cov_vars_base <- c("Latitude", "Latitude_2")
vcov_base <- vcov_matrix[cov_vars_base, cov_vars_base]

d_beta1 <- P_values
d_beta2 <- P_values^2

var_effect_base <- (d_beta1^2) * vcov_base["Latitude", "Latitude"] +
  (d_beta2^2) * vcov_base["Latitude_2", "Latitude_2"] +
  2 * d_beta1 * d_beta2 * vcov_base["Latitude", "Latitude_2"]
se_effect_base <- sqrt(var_effect_base)
ci_effect_base_lower <- effect_base - 1.96 * se_effect_base
ci_effect_base_upper <- effect_base + 1.96 * se_effect_base

var_effect_rate <- list()
se_effect_rate <- list()
ci_effect_rate_lower <- list()
ci_effect_rate_upper <- list()

for (cont in rates) {
  cov_vars <- c("Latitude", "Latitude_2", 
                paste0("Latitude:expansionrate", cont), 
                paste0("Latitude_2:expansionrate", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  
  var_effect <- (d_beta1^2) * vcov_cont["Latitude", "Latitude"] +
    (d_beta2^2) * vcov_cont["Latitude_2", "Latitude_2"] +
    (d_beta1^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_beta2^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_beta1^2 * vcov_cont["Latitude", cov_vars[3]] +
    2 * d_beta2^2 * vcov_cont["Latitude_2", cov_vars[4]] +
    2 * d_beta1 * d_beta2 * vcov_cont["Latitude", "Latitude_2"] +
    2 * d_beta1 * d_beta2 * vcov_cont["Latitude", cov_vars[4]] +
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], "Latitude_2"] +
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], cov_vars[4]]
  
  var_effect_rate[[cont]] <- var_effect
  se_effect_rate[[cont]] <- sqrt(var_effect)
  ci_effect_rate_lower[[cont]] <- effect_rate[[cont]] - 1.96 * se_effect_rate[[cont]]
  ci_effect_rate_upper[[cont]] <- effect_rate[[cont]] + 1.96 * se_effect_rate[[cont]]
}

# 计算边际效应
marginal_base <- beta1 + 2 * beta2 * P_values

marginal_rate <- lapply(rates, function(cont) {
  beta1_total[[cont]] + 2 * beta2_total[[cont]] * P_values
})
names(marginal_rate) <- rates

d_marg_beta1 <- 1
d_marg_beta2 <- 2 * P_values

var_marginal_base <- (d_marg_beta1^2) * vcov_base["Latitude", "Latitude"] +
  (d_marg_beta2^2) * vcov_base["Latitude_2", "Latitude_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_base["Latitude", "Latitude_2"]
se_marginal_base <- sqrt(var_marginal_base)
ci_marginal_base_lower <- marginal_base - 1.96 * se_marginal_base
ci_marginal_base_upper <- marginal_base + 1.96 * se_marginal_base

var_marginal_rate <- list()
se_marginal_rate <- list()
ci_marginal_rate_lower <- list()
ci_marginal_rate_upper <- list()

for (cont in rates) {
  cov_vars <- c("Latitude", "Latitude_2",
                paste0("Latitude:expansionrate", cont),
                paste0("Latitude_2:expansionrate", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  
  d_marg_beta1_u <- 1
  d_marg_beta2_u <- 2 * P_values
  
  var_marginal <- (d_marg_beta1_u^2) * vcov_cont["Latitude", "Latitude"] +
    (d_marg_beta2_u^2) * vcov_cont["Latitude_2", "Latitude_2"] +
    (d_marg_beta1_u^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_marg_beta2_u^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_marg_beta1_u^2 * vcov_cont["Latitude", cov_vars[3]] +
    2 * d_marg_beta2_u^2 * vcov_cont["Latitude_2", cov_vars[4]] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["Latitude", "Latitude_2"] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["Latitude", cov_vars[4]] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], "Latitude_2"] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], cov_vars[4]]
  
  var_marginal_rate[[cont]] <- var_marginal
  se_marginal_rate[[cont]] <- sqrt(var_marginal)
  ci_marginal_rate_lower[[cont]] <- marginal_rate[[cont]] - 1.96 * se_marginal_rate[[cont]]
  ci_marginal_rate_upper[[cont]] <- marginal_rate[[cont]] + 1.96 * se_marginal_rate[[cont]]
}

effect_df <- data.frame(
  Latitude = P_values,
  rate = "基准组",
  Effect = effect_base,
  SE = se_effect_base,
  CI_Lower = ci_effect_base_lower,
  CI_Upper = ci_effect_base_upper
)

for (cont in rates) {
  temp_df <- data.frame(
    Latitude = P_values,
    rate = cont,
    Effect = effect_rate[[cont]],
    SE = se_effect_rate[[cont]],
    CI_Lower = ci_effect_rate_lower[[cont]],
    CI_Upper = ci_effect_rate_upper[[cont]]
  )
  effect_df <- rbind(effect_df, temp_df)
}

marginal_df <- data.frame(
  Latitude = P_values,
  rate = "基准组",
  Marginal_Effect = marginal_base,
  SE = se_marginal_base,
  CI_Lower = ci_marginal_base_lower,
  CI_Upper = ci_marginal_base_upper
)

for (cont in rates) {
  temp_df <- data.frame(
    Latitude = P_values,
    rate = cont,
    Marginal_Effect = marginal_rate[[cont]],
    SE = se_marginal_rate[[cont]],
    CI_Lower = ci_marginal_rate_lower[[cont]],
    CI_Upper = ci_marginal_rate_upper[[cont]]
  )
  marginal_df <- rbind(marginal_df, temp_df)
}

effect_df$rate <- factor(effect_df$rate, 
                         levels = c("基准组", rates),
                         labels = c( "Low",   "Middle",    "High",    "Very high"))

marginal_df$rate <- factor(marginal_df$rate,
                           levels = c("基准组", rates),
                           labels = c( "Low",   "Middle",    "High",    "Very high"))

rate_limits <- tibble(
  rate = c("Low",   "Middle",    "High",    "Very high"),
  min_WD = c(
    min(urban_in_dt_Latitude[urban_in_dt_Latitude$expansionrate == "Slow", ]$Latitude, na.rm = TRUE),
    min(urban_in_dt_Latitude[urban_in_dt_Latitude$expansionrate == "Moderate", ]$Latitude, na.rm = TRUE),
    min(urban_in_dt_Latitude[urban_in_dt_Latitude$expansionrate == "Rapid", ]$Latitude, na.rm = TRUE),
    min(urban_in_dt_Latitude[urban_in_dt_Latitude$expansionrate == "Intense", ]$Latitude, na.rm = TRUE)
  ),
  max_WD = c(
    max(urban_in_dt_Latitude[urban_in_dt_Latitude$expansionrate == "Slow", ]$Latitude, na.rm = TRUE),
    max(urban_in_dt_Latitude[urban_in_dt_Latitude$expansionrate == "Moderate", ]$Latitude, na.rm = TRUE),
    max(urban_in_dt_Latitude[urban_in_dt_Latitude$expansionrate == "Rapid", ]$Latitude, na.rm = TRUE),
    max(urban_in_dt_Latitude[urban_in_dt_Latitude$expansionrate == "Intense", ]$Latitude, na.rm = TRUE)
  )
)

print(rate_limits)

effect_df <- effect_df %>%
  left_join(rate_limits, by = "rate") %>%
  filter(Latitude >= min_WD & Latitude <= max_WD) %>%
  select(-min_WD, -max_WD)

marginal_df <- marginal_df %>%
  left_join(rate_limits, by = "rate") %>%
  filter(Latitude >= min_WD & Latitude <= max_WD) %>%
  select(-min_WD, -max_WD)

# 绘制效应值图
# Figure 3j
effect_plot <- ggplot(effect_df, aes(x = Latitude, y = Effect, color = rate)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = rate),
    alpha = 0.15, color = NA
  ) +
  theme_classic() +
  labs(
    x = "",
    y = "",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, by =10), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.08, 0.04),
                     breaks = seq(-0.08, 0.04, by = 0.04),
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.text.y = element_text(color = "black", size = 13),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),
  ) +
  scale_color_manual(values = rate_colors) +
  scale_fill_manual(values = rate_colors)

print(effect_plot)

# 绘制直方图
hist_plot <- ggplot(urban_in_dt_Latitude, aes(x = Latitude)) +
  geom_histogram(aes(y = ..density.., fill=expansionrate), binwidth = 0.15, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, by =10)) +
  scale_y_continuous(limits = c(0, 0.25), breaks = seq(0, 0.25, by =0.1)) +
  labs(x = "",
       y = "") +
  theme_classic() +
  theme(
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text = element_text(color = "black", size = 13),
    axis.title = element_text(size = 14),
    legend.position = "none"
  ) + 
  scale_fill_manual(values = rate_colors_hist)

print(hist_plot)

# 合并
effect_rate_latitude <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_rate_latitude)

# 绘制边际效应图
# Figure S5j
marginal_plot <- ggplot(marginal_df, aes(x = Latitude, y = Marginal_Effect, color = rate)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = rate),
    alpha = 0.15, color = NA
  ) +
  labs(
    x = "",
    y = "",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  theme_classic() +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, by =10), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.005, 0.005),
                     breaks = seq(-0.005, 0.005, by = 0.005),
                     labels = scales::number_format(accuracy = 0.001)) +
  theme(
    axis.text.y = element_text(color = "black", size = 13),  
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title.y = element_text(size = 14),  
    legend.position = "none",   
    legend.direction = "vertical", 
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),  
  ) +
  scale_color_manual(values = rate_colors) +
  scale_fill_manual(values = rate_colors)

print(marginal_plot)

marginal_rate_latitude <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))

print(marginal_rate_latitude)



################################################### 城市内异质性分析 ################################################## 
##################################################### country ################################################## 
################################################### elevation ######################################
cols_to_center <- c('WaterDeficit', 'Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed',  'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
urban_in_dt_urbanization_Elevation <- copy(urban_in_dt)
urban_in_dt_urbanization_Elevation <- urban_in_dt_urbanization_Elevation[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                                         .SDcols = cols_to_center]
urban_in_dt_urbanization_Elevation$country <- factor(urban_in_dt_urbanization_Elevation$country, 
                                                     levels = c( "United States of America",
                                                                 "Europe",
                                                                 "Canada",
                                                                 "China",
                                                                 "India",
                                                                 "Brazil", "other"))

# linear regression
# Table S7 -- model 3
model_urban_in_country <- felm(FVC ~ 1 + WaterDeficit +
                                 Seasonalityindex +
                                 AverageTemperature +
                                 TemperatureRange +
                                 WindSpeed +
                                 Elevation +
                                 # Elevation_2 +
                                 Latitude +
                                 Precipitation +
                                 SoilMoisture +
                                 SoilPH + 
                                 Longitude + 
                                 
                                 Elevation:country +      
                                 # Elevation_2:country +
                                 
                                 GDP_per_capita_PPP + HDI + Population + 
                                 ImperviousSurface + HumanSettlement + CityArea
                               | year, 
                               data = urban_in_dt_urbanization_Elevation)

summary(model_urban_in_country)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_country.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in_country))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

# nonlinear regression
# Table S7 -- model 7
model_urban_in_country <- felm(FVC ~ 1 + WaterDeficit +
                                 Seasonalityindex +
                                 AverageTemperature +
                                 TemperatureRange +
                                 WindSpeed +
                                 Elevation +
                                 Elevation_2 +
                                 Latitude +
                                 Precipitation +
                                 SoilMoisture +
                                 SoilPH + 
                                 Longitude + 
                                 
                                 Elevation:country +      
                                 Elevation_2:country +
                                 
                                 GDP_per_capita_PPP + HDI + Population + 
                                 ImperviousSurface + HumanSettlement + CityArea
                               | year, 
                               data = urban_in_dt_urbanization_Elevation)

summary(model_urban_in_country)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_country2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in_country))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

######################################### 绘制效应图与边际效应图 #######################################
beta1 <- coef(model_urban_in_country)["Elevation"]
beta2 <- coef(model_urban_in_country)["Elevation_2"]

countries <- c("Europe", "Canada",  "China", "India", "Brazil", "other")

beta1_country <- sapply(countries, function(cty) {
  coef(model_urban_in_country)[paste0("Elevation:country", cty)]
})
beta2_country <- sapply(countries, function(cty) {
  coef(model_urban_in_country)[paste0("Elevation_2:country", cty)]
})

beta1_total <- beta1 + beta1_country
beta2_total <- beta2 + beta2_country

names(beta1_total) <- countries
names(beta2_total) <- countries

P_min <- quantile(urban_in_dt_urbanization_Elevation$Elevation, 0.001, na.rm = TRUE)
P_max <- quantile(urban_in_dt_urbanization_Elevation$Elevation, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)

# 计算效应值
effect_base <- beta1 * P_values + beta2 * (P_values^2)

effect_country <- lapply(countries, function(cty) {
  beta1_total[[cty]] * P_values + beta2_total[[cty]] * (P_values^2)
})
names(effect_country) <- countries

vcov_matrix <- vcovHC(model_urban_in_country, type = "HC1", cluster = "group",
                      group = urban_in_dt_urbanization_Elevation$CountryID)

cov_vars_base <- c("Elevation", "Elevation_2")
vcov_base <- vcov_matrix[cov_vars_base, cov_vars_base]

d_beta1 <- P_values
d_beta2 <- P_values^2

var_effect_base <- (d_beta1^2) * vcov_base["Elevation", "Elevation"] +
  (d_beta2^2) * vcov_base["Elevation_2", "Elevation_2"] +
  2 * d_beta1 * d_beta2 * vcov_base["Elevation", "Elevation_2"]
se_effect_base <- sqrt(var_effect_base)
ci_effect_base_lower <- effect_base - 1.96 * se_effect_base
ci_effect_base_upper <- effect_base + 1.96 * se_effect_base

var_effect_country <- list()
se_effect_country <- list()
ci_effect_country_lower <- list()
ci_effect_country_upper <- list()

for (cty in countries) {
  cov_vars <- c("Elevation", "Elevation_2", 
                paste0("Elevation:country", cty), 
                paste0("Elevation_2:country", cty))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  
  var_effect <- (d_beta1^2) * vcov_cont["Elevation", "Elevation"] +
    (d_beta2^2) * vcov_cont["Elevation_2", "Elevation_2"] +
    (d_beta1^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_beta2^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_beta1^2 * vcov_cont["Elevation", cov_vars[3]] +
    2 * d_beta2^2 * vcov_cont["Elevation_2", cov_vars[4]] +
    2 * d_beta1 * d_beta2 * vcov_cont["Elevation", "Elevation_2"] +
    2 * d_beta1 * d_beta2 * vcov_cont["Elevation", cov_vars[4]] +
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], "Elevation_2"] +
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], cov_vars[4]]
  
  var_effect_country[[cty]] <- var_effect
  se_effect_country[[cty]] <- sqrt(var_effect)
  ci_effect_country_lower[[cty]] <- effect_country[[cty]] - 1.96 * se_effect_country[[cty]]
  ci_effect_country_upper[[cty]] <- effect_country[[cty]] + 1.96 * se_effect_country[[cty]]
}

# 计算边际效应
marginal_base <- beta1 + 2 * beta2 * P_values

marginal_country <- lapply(countries, function(cty) {
  beta1_total[[cty]] + 2 * beta2_total[[cty]] * P_values
})
names(marginal_country) <- countries

d_marg_beta1 <- 1
d_marg_beta2 <- 2 * P_values

var_marginal_base <- (d_marg_beta1^2) * vcov_base["Elevation", "Elevation"] +
  (d_marg_beta2^2) * vcov_base["Elevation_2", "Elevation_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_base["Elevation", "Elevation_2"]
se_marginal_base <- sqrt(var_marginal_base)
ci_marginal_base_lower <- marginal_base - 1.96 * se_marginal_base
ci_marginal_base_upper <- marginal_base + 1.96 * se_marginal_base

var_marginal_country <- list()
se_marginal_country <- list()
ci_marginal_country_lower <- list()
ci_marginal_country_upper <- list()

for (cty in countries) {
  cov_vars <- c("Elevation", "Elevation_2",
                paste0("Elevation:country", cty),
                paste0("Elevation_2:country", cty))
  vcov_cty <- vcov_matrix[cov_vars, cov_vars]
  
  d_marg_beta1_u <- 1
  d_marg_beta2_u <- 2 * P_values
  
  var_marginal <- (d_marg_beta1_u^2) * vcov_cty["Elevation", "Elevation"] +
    (d_marg_beta2_u^2) * vcov_cty["Elevation_2", "Elevation_2"] +
    (d_marg_beta1_u^2) * vcov_cty[cov_vars[3], cov_vars[3]] +
    (d_marg_beta2_u^2) * vcov_cty[cov_vars[4], cov_vars[4]] +
    2 * d_marg_beta1_u^2 * vcov_cty["Elevation", cov_vars[3]] +
    2 * d_marg_beta2_u^2 * vcov_cty["Elevation_2", cov_vars[4]] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cty["Elevation", "Elevation_2"] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cty["Elevation", cov_vars[4]] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cty[cov_vars[3], "Elevation_2"] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cty[cov_vars[3], cov_vars[4]]
  
  var_marginal_country[[cty]] <- var_marginal
  se_marginal_country[[cty]] <- sqrt(var_marginal)
  ci_marginal_country_lower[[cty]] <- marginal_country[[cty]] - 1.96 * se_marginal_country[[cty]]
  ci_marginal_country_upper[[cty]] <- marginal_country[[cty]] + 1.96 * se_marginal_country[[cty]]
}

marginal_df <- data.frame(
  Elevation = P_values,
  Country = "基准组",
  Marginal_Effect = marginal_base,
  SE = se_marginal_base,
  CI_Lower = ci_marginal_base_lower,
  CI_Upper = ci_marginal_base_upper
)

for (cty in countries) {
  temp_df <- data.frame(
    Elevation = P_values,
    Country = cty,
    Marginal_Effect = marginal_country[[cty]],
    SE = se_marginal_country[[cty]],
    CI_Lower = ci_marginal_country_lower[[cty]],
    CI_Upper = ci_marginal_country_upper[[cty]]
  )
  marginal_df <- rbind(marginal_df, temp_df)
}

effect_df <- data.frame(
  Elevation = P_values,
  Country = "基准组",
  Effect = effect_base,
  SE = se_effect_base,
  CI_Lower = ci_effect_base_lower,
  CI_Upper = ci_effect_base_upper
)

for (cty in countries) {
  temp_df <- data.frame(
    Elevation = P_values,
    Country = cty,
    Effect = effect_country[[cty]],
    SE = se_effect_country[[cty]],
    CI_Lower = ci_effect_country_lower[[cty]],
    CI_Upper = ci_effect_country_upper[[cty]]
  )
  effect_df <- rbind(effect_df, temp_df)
}

effect_df$Country <- factor(effect_df$Country, 
                            levels = c("基准组", countries),
                            labels = c("United States", "Europe (excluding Russia)", "Canada",  "China", "India", "Brazil", "other"))

marginal_df$Country <- factor(marginal_df$Country, 
                              levels = c("基准组", countries),
                              labels = c("United States", "Europe (excluding Russia)", "Canada",  "China", "India", "Brazil", "other"))

country_limits <- tibble(
  Country = c("United States", "Europe (excluding Russia)", "Canada",  "China", "India", "Brazil", "other"),
  min_WD = c(
    min(urban_in_dt_urbanization_Elevation[country == "United States of America"]$Elevation, na.rm = TRUE),
    min(urban_in_dt_urbanization_Elevation[country == "Europe" ]$Elevation, na.rm = TRUE),
    min(urban_in_dt_urbanization_Elevation[country == "Canada" ]$Elevation, na.rm = TRUE),
    min(urban_in_dt_urbanization_Elevation[country == "China" ]$Elevation, na.rm = TRUE),
    min(urban_in_dt_urbanization_Elevation[country == "India" ]$Elevation, na.rm = TRUE),
    min(urban_in_dt_urbanization_Elevation[country == "Brazil" ]$Elevation, na.rm = TRUE),
    min(urban_in_dt_urbanization_Elevation[country == "other" ]$Elevation, na.rm = TRUE)
  ),
  max_WD = c(
    max(urban_in_dt_urbanization_Elevation[country == "United States of America" ]$Elevation, na.rm = TRUE),
    max(urban_in_dt_urbanization_Elevation[country == "Europe" ]$Elevation, na.rm = TRUE),
    max(urban_in_dt_urbanization_Elevation[country == "Canada"]$Elevation, na.rm = TRUE),
    max(urban_in_dt_urbanization_Elevation[country == "China"]$Elevation, na.rm = TRUE),
    max(urban_in_dt_urbanization_Elevation[country == "India"]$Elevation, na.rm = TRUE),
    max(urban_in_dt_urbanization_Elevation[country == "Brazil"]$Elevation, na.rm = TRUE),
    max(urban_in_dt_urbanization_Elevation[country == "other"]$Elevation, na.rm = TRUE)
  )
)

print(country_limits)

effect_df <- effect_df %>%
  left_join(country_limits, by = "Country") %>%
  filter(Elevation >= min_WD & Elevation <= max_WD) %>%
  select(-min_WD, -max_WD)

marginal_df <- marginal_df %>%
  left_join(country_limits, by = "Country") %>%
  filter(Elevation >= min_WD & Elevation <= max_WD) %>%
  select(-min_WD, -max_WD)

# 绘制效应值图
# Figure 3k
country_colors <- c(
  "United States" = "#2ecc71",
  "Europe (excluding Russia)" = "#e84393",
  "Canada" = "#3498db",
  "China" = "#CB0505",
  "India" = "#f39c12",
  "Brazil" = "#9b59b6",
  "other" = "#95a5a6"
)

effect_df$Country <- factor(effect_df$Country, levels = c("United States", "Europe (excluding Russia)", "Canada", "China", "India", "Brazil" ,"other"))

effect_plot <- ggplot(effect_df, aes(x = Elevation, y = Effect, color = Country)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Country),
    alpha = 0.15, color = NA
  ) +
  theme_classic() +
  labs(
    x = "",
    y = "Effect on urban FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 4050), breaks = seq(0, 4050, by =1000), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.3, 0.85),
                     breaks = seq(-0.3, 0.85, by = 0.3),
                     labels = scales::number_format(accuracy = 0.1)) +
  theme(
    axis.text.y = element_text(color = "black", size = 13),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.52,0.77),
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),
    legend.background = element_blank()
  ) +
  scale_color_manual(values = country_colors) +
  scale_fill_manual(values = country_colors)

print(effect_plot)

# 绘制直方图
hist_plot <- ggplot(urban_in_dt_urbanization_Elevation, aes(x = Elevation)) +
  geom_histogram(aes(y = ..density.., fill=country), binwidth = 3, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 4050), breaks = seq(0, 4050, by =1000)) +
  scale_y_continuous(limits = c(0, 0.07), breaks = seq(0, 0.07, by =0.03)) +
  labs(x = "Elevation (m)",
       y = "") +
  theme_classic() +
  theme(
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text = element_text(color = "black", size = 13),
    axis.title = element_text(size = 14),
    legend.position = "none"
  ) + 
  scale_fill_manual(values = c(
    "United States of America" = "#2ecc71",
    "Europe (excluding Russia)" = "#e84393",
    "Canada" = "#3498db",
    "China" = "#CB0505",
    "India" = "#f39c12",
    "Brazil" = "#9b59b6",
    "other" = "#95a5a6"
  ))

print(hist_plot)

# 合并
effect_country_elevation <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_country_elevation)


# 绘制边际效应图
# Figure S5k
marginal_df$Country <- factor(marginal_df$Country, levels = c("United States", "Europe (excluding Russia)", "Canada", "China", "India", "Brazil" ,"other"))

marginal_plot <- ggplot(marginal_df, aes(x = Elevation, y = Marginal_Effect, color = Country)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Country),
    alpha = 0.15, color = NA
  ) +
  labs(
    x = "",
    y = "Marginal effect on urban FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  theme_classic() +
  scale_x_continuous(limits = c(0, 4050), breaks = seq(0, 4050, by =1000), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.0004, 0.0008),
                     breaks = seq(-0.0004, 0.0008, by = 0.0004),
                     labels = scales::number_format(accuracy = 0.0001)) +
  theme(
    axis.text.y = element_text(color = "black", size = 12), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),  
    axis.title.y = element_text(size = 12),  
    legend.position = c(0.62,0.8),   
    legend.direction = "vertical", 
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),  
    legend.background = element_blank()
  ) +
  scale_color_manual(values = country_colors) +
  scale_fill_manual(values = country_colors)

print(marginal_plot)

marginal_country_elevation <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_country_elevation)

##################################################### country ####################################################
################################################### latitude ######################################
cols_to_center <- c('WaterDeficit', 'Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
urban_in_dt_urbanization_Latitude <- copy(urban_in_dt)
urban_in_dt_urbanization_Latitude <- urban_in_dt_urbanization_Latitude[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                                       .SDcols = cols_to_center]
urban_in_dt_urbanization_Latitude$country <- factor(urban_in_dt_urbanization_Latitude$country, 
                                                    levels = c( "United States of America",
                                                                "Europe",
                                                                "Canada",
                                                                "China",
                                                                "India",
                                                                "Brazil", "other"))

# linear regression
# Table S7 -- model 4
model_urban_in_country <- felm(FVC ~ 1 + WaterDeficit +
                                 Seasonalityindex +
                                 AverageTemperature +
                                 TemperatureRange +
                                 WindSpeed +
                                 Elevation +
                                 Latitude +
                                 # Latitude_2 +
                                 Precipitation +
                                 SoilMoisture +
                                 SoilPH + 
                                 Longitude +
                                 
                                 Latitude:country +      
                                 # Latitude_2:country +
                                 
                                 GDP_per_capita_PPP + HDI + Population + 
                                 ImperviousSurface + HumanSettlement + CityArea
                               | year, 
                               data = urban_in_dt_urbanization_Latitude)

summary(model_urban_in_country)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_country.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in_country))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con) 

# nonlinear regression
# Table S7 -- model 8
model_urban_in_country <- felm(FVC ~ 1 + WaterDeficit +
                                 Seasonalityindex +
                                 AverageTemperature +
                                 TemperatureRange +
                                 WindSpeed +
                                 Elevation +
                                 Latitude +
                                 Latitude_2 +
                                 Precipitation +
                                 SoilMoisture +
                                 SoilPH + 
                                 Longitude +
                                 
                                 Latitude:country +      
                                 Latitude_2:country +
                                 
                                 GDP_per_capita_PPP + HDI + Population + 
                                 ImperviousSurface + HumanSettlement + CityArea
                               | year, 
                               data = urban_in_dt_urbanization_Latitude)

summary(model_urban_in_country)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_country2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in_country))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con) 


######################################### 绘制效应图与边际效应图 #######################################
beta1 <- coef(model_urban_in_country)["Latitude"]
beta2 <- coef(model_urban_in_country)["Latitude_2"]

countries <- c("Europe", "Canada",  "China", "India", "Brazil", "other")

beta1_country <- sapply(countries, function(cty) {
  coef(model_urban_in_country)[paste0("Latitude:country", cty)]
})
beta2_country <- sapply(countries, function(cty) {
  coef(model_urban_in_country)[paste0("Latitude_2:country", cty)]
})

beta1_total <- beta1 + beta1_country
beta2_total <- beta2 + beta2_country

names(beta1_total) <- countries
names(beta2_total) <- countries

P_min <- quantile(urban_in_dt_urbanization_Latitude$Latitude, 0.001, na.rm = TRUE)
P_max <- quantile(urban_in_dt_urbanization_Latitude$Latitude, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)

# 计算效应值
effect_base <- beta1 * P_values + beta2 * (P_values^2)

effect_country <- lapply(countries, function(cty) {
  beta1_total[[cty]] * P_values + beta2_total[[cty]] * (P_values^2)
})
names(effect_country) <- countries

vcov_matrix <- vcovHC(model_urban_in_country, type = "HC1", cluster = "group",
                      group = urban_in_dt_urbanization_Latitude$CountryID)

cov_vars_base <- c("Latitude", "Latitude_2")
vcov_base <- vcov_matrix[cov_vars_base, cov_vars_base]

d_beta1 <- P_values
d_beta2 <- P_values^2

var_effect_base <- (d_beta1^2) * vcov_base["Latitude", "Latitude"] +
  (d_beta2^2) * vcov_base["Latitude_2", "Latitude_2"] +
  2 * d_beta1 * d_beta2 * vcov_base["Latitude", "Latitude_2"]
se_effect_base <- sqrt(var_effect_base)
ci_effect_base_lower <- effect_base - 1.96 * se_effect_base
ci_effect_base_upper <- effect_base + 1.96 * se_effect_base

var_effect_country <- list()
se_effect_country <- list()
ci_effect_country_lower <- list()
ci_effect_country_upper <- list()

for (cty in countries) {
  cov_vars <- c("Latitude", "Latitude_2", 
                paste0("Latitude:country", cty), 
                paste0("Latitude_2:country", cty))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  
  var_effect <- (d_beta1^2) * vcov_cont["Latitude", "Latitude"] +
    (d_beta2^2) * vcov_cont["Latitude_2", "Latitude_2"] +
    (d_beta1^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_beta2^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_beta1^2 * vcov_cont["Latitude", cov_vars[3]] +
    2 * d_beta2^2 * vcov_cont["Latitude_2", cov_vars[4]] +
    2 * d_beta1 * d_beta2 * vcov_cont["Latitude", "Latitude_2"] +
    2 * d_beta1 * d_beta2 * vcov_cont["Latitude", cov_vars[4]] +
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], "Latitude_2"] +
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], cov_vars[4]]
  
  var_effect_country[[cty]] <- var_effect
  se_effect_country[[cty]] <- sqrt(var_effect)
  ci_effect_country_lower[[cty]] <- effect_country[[cty]] - 1.96 * se_effect_country[[cty]]
  ci_effect_country_upper[[cty]] <- effect_country[[cty]] + 1.96 * se_effect_country[[cty]]
}

# 计算边际效应
marginal_base <- beta1 + 2 * beta2 * P_values

marginal_country <- lapply(countries, function(cty) {
  beta1_total[[cty]] + 2 * beta2_total[[cty]] * P_values
})
names(marginal_country) <- countries

d_marg_beta1 <- 1
d_marg_beta2 <- 2 * P_values

var_marginal_base <- (d_marg_beta1^2) * vcov_base["Latitude", "Latitude"] +
  (d_marg_beta2^2) * vcov_base["Latitude_2", "Latitude_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_base["Latitude", "Latitude_2"]
se_marginal_base <- sqrt(var_marginal_base)
ci_marginal_base_lower <- marginal_base - 1.96 * se_marginal_base
ci_marginal_base_upper <- marginal_base + 1.96 * se_marginal_base

var_marginal_country <- list()
se_marginal_country <- list()
ci_marginal_country_lower <- list()
ci_marginal_country_upper <- list()

for (cty in countries) {
  cov_vars <- c("Latitude", "Latitude_2",
                paste0("Latitude:country", cty),
                paste0("Latitude_2:country", cty))
  vcov_cty <- vcov_matrix[cov_vars, cov_vars]
  
  d_marg_beta1_u <- 1
  d_marg_beta2_u <- 2 * P_values
  
  var_marginal <- (d_marg_beta1_u^2) * vcov_cty["Latitude", "Latitude"] +
    (d_marg_beta2_u^2) * vcov_cty["Latitude_2", "Latitude_2"] +
    (d_marg_beta1_u^2) * vcov_cty[cov_vars[3], cov_vars[3]] +
    (d_marg_beta2_u^2) * vcov_cty[cov_vars[4], cov_vars[4]] +
    2 * d_marg_beta1_u^2 * vcov_cty["Latitude", cov_vars[3]] +
    2 * d_marg_beta2_u^2 * vcov_cty["Latitude_2", cov_vars[4]] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cty["Latitude", "Latitude_2"] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cty["Latitude", cov_vars[4]] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cty[cov_vars[3], "Latitude_2"] +
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cty[cov_vars[3], cov_vars[4]]
  
  var_marginal_country[[cty]] <- var_marginal
  se_marginal_country[[cty]] <- sqrt(var_marginal)
  ci_marginal_country_lower[[cty]] <- marginal_country[[cty]] - 1.96 * se_marginal_country[[cty]]
  ci_marginal_country_upper[[cty]] <- marginal_country[[cty]] + 1.96 * se_marginal_country[[cty]]
}

marginal_df <- data.frame(
  Latitude = P_values,
  Country = "基准组",
  Marginal_Effect = marginal_base,
  SE = se_marginal_base,
  CI_Lower = ci_marginal_base_lower,
  CI_Upper = ci_marginal_base_upper
)

for (cty in countries) {
  temp_df <- data.frame(
    Latitude = P_values,
    Country = cty,
    Marginal_Effect = marginal_country[[cty]],
    SE = se_marginal_country[[cty]],
    CI_Lower = ci_marginal_country_lower[[cty]],
    CI_Upper = ci_marginal_country_upper[[cty]]
  )
  marginal_df <- rbind(marginal_df, temp_df)
}

effect_df <- data.frame(
  Latitude = P_values,
  Country = "基准组",
  Effect = effect_base,
  SE = se_effect_base,
  CI_Lower = ci_effect_base_lower,
  CI_Upper = ci_effect_base_upper
)

for (cty in countries) {
  temp_df <- data.frame(
    Latitude = P_values,
    Country = cty,
    Effect = effect_country[[cty]],
    SE = se_effect_country[[cty]],
    CI_Lower = ci_effect_country_lower[[cty]],
    CI_Upper = ci_effect_country_upper[[cty]]
  )
  effect_df <- rbind(effect_df, temp_df)
}

effect_df$Country <- factor(effect_df$Country, 
                            levels = c("基准组", countries),
                            labels = c("United States", "Europe (excluding Russia)", "Canada",  "China", "India", "Brazil", "other"))

marginal_df$Country <- factor(marginal_df$Country, 
                              levels = c("基准组", countries),
                              labels = c("United States", "Europe (excluding Russia)", "Canada",  "China", "India", "Brazil", "other"))

country_limits <- tibble(
  Country = c("United States", "Europe (excluding Russia)", "Canada",  "China", "India", "Brazil", "other"),
  min_WD = c(
    min(urban_in_dt_urbanization_Latitude[country == "United States of America"]$Latitude, na.rm = TRUE),
    min(urban_in_dt_urbanization_Latitude[country == "Europe" ]$Latitude, na.rm = TRUE),
    min(urban_in_dt_urbanization_Latitude[country == "Canada" ]$Latitude, na.rm = TRUE),
    min(urban_in_dt_urbanization_Latitude[country == "China" ]$Latitude, na.rm = TRUE),
    min(urban_in_dt_urbanization_Latitude[country == "India" ]$Latitude, na.rm = TRUE),
    min(urban_in_dt_urbanization_Latitude[country == "Brazil" ]$Latitude, na.rm = TRUE),
    min(urban_in_dt_urbanization_Latitude[country == "other" ]$Latitude, na.rm = TRUE)
  ),
  max_WD = c(
    max(urban_in_dt_urbanization_Latitude[country == "United States of America" ]$Latitude, na.rm = TRUE),
    max(urban_in_dt_urbanization_Latitude[country == "Europe" ]$Latitude, na.rm = TRUE),
    max(urban_in_dt_urbanization_Latitude[country == "Canada"]$Latitude, na.rm = TRUE),
    max(urban_in_dt_urbanization_Latitude[country == "China"]$Latitude, na.rm = TRUE),
    max(urban_in_dt_urbanization_Latitude[country == "India"]$Latitude, na.rm = TRUE),
    max(urban_in_dt_urbanization_Latitude[country == "Brazil"]$Latitude, na.rm = TRUE),
    max(urban_in_dt_urbanization_Latitude[country == "other"]$Latitude, na.rm = TRUE)
  )
)

print(country_limits)

effect_df <- effect_df %>%
  left_join(country_limits, by = "Country") %>%
  filter(Latitude >= min_WD & Latitude <= max_WD) %>%
  select(-min_WD, -max_WD)

marginal_df <- marginal_df %>%
  left_join(country_limits, by = "Country") %>%
  filter(Latitude >= min_WD & Latitude <= max_WD) %>%
  select(-min_WD, -max_WD)

# 绘制效应值图
# Figure 3l
effect_plot <- ggplot(effect_df, aes(x = Latitude, y = Effect, color = Country)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Country),
    alpha = 0.15, color = NA
  ) +
  theme_classic() +
  labs(
    x = "",
    y = "",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, by =10), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.2, 0.1),
                     breaks = seq(-0.2, 0.1, by = 0.1),
                     labels = scales::number_format(accuracy = 0.1)) +
  theme(
    axis.text.y = element_text(color = "black", size = 13),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),
  ) +
  scale_color_manual(values = country_colors) +
  scale_fill_manual(values = country_colors)

print(effect_plot)

# 绘制直方图
hist_plot <- ggplot(urban_in_dt_urbanization_Latitude, aes(x = Latitude)) +
  geom_histogram(aes(y = ..density.., fill=country), binwidth = 0.1, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, by =10)) +
  labs(x = "Abs(Latitude) (°)") +
  theme_classic() +
  theme(
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text = element_text(color = "black", size = 13),
    axis.title = element_text(size = 14),
    legend.position = "none"
  ) + 
  scale_fill_manual(values = c(
    "United States of America" = "#2ecc71",
    "Europe (excluding Russia)" = "#e84393",
    "Canada" = "#3498db",
    "China" = "#CB0505",
    "India" = "#f39c12",
    "Brazil" = "#9b59b6",
    "other" = "#95a5a6"
  ))

print(hist_plot)

# 合并
effect_country_latitude <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_country_latitude)

# 绘制边际效应图
# Figure S5l
marginal_plot <- ggplot(marginal_df, aes(x = Latitude, y = Marginal_Effect, color = Country)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Country),
    alpha = 0.15, color = NA
  ) +
  labs(
    x = "",
    y = "",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  theme_classic() +
  scale_x_continuous(limits = c(0, 70), breaks = seq(0, 70, by =10), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.008, 0.004),
                     breaks = seq(-0.008, 0.004, by = 0.004),
                     labels = scales::number_format(accuracy = 0.001)) +
  theme(
    axis.text.y = element_text(color = "black", size = 13),  
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),  
    axis.title.y = element_text(size = 14), 
    legend.position = "none",  
    legend.direction = "vertical", 
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),  
  ) +
  scale_color_manual(values = country_colors) +
  scale_fill_manual(values = country_colors)

print(marginal_plot)

marginal_country_latitude <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))

print(marginal_country_latitude)


################################################合并图形###############################################################
# Figure 3: The effect of urbanization on the impacts of geographic stresses on UGS coverage
ii <- cowplot::plot_grid(effect_elevation, nrow = 1, ncol = 1, labels = c("a"),               
                         label_size = 18, scale=1, label_x = 0.02, label_y = 1.00)

jj <- cowplot::plot_grid(effect_latitude, nrow = 1, ncol = 1, labels = c("b"),
                         label_size = 18, scale=1, label_x = 0.02, label_y = 1.00)

kk <- cowplot::plot_grid(effect_continent_elevation, nrow = 1, ncol = 1, labels = c("c"),
                         label_size = 18, scale=1, label_x = -0.05, label_y = 1.15)

mm <- cowplot::plot_grid(effect_continent_latitude, nrow = 1, ncol = 1, labels = c("d"),
                         label_size = 18, scale=1, label_x = 0.02, label_y = 1.15)

oo <- cowplot::plot_grid(effect_climate_elevation, nrow = 1, ncol = 1, labels = c("e"),
                         label_size = 18, scale=1, label_x = -0.05, label_y = 1.15)

pp <- cowplot::plot_grid(effect_climate_latitude, nrow = 1, ncol = 1, labels = c("f"),
                         label_size = 18, scale=1, label_x = 0.02, label_y = 1.15)

oo2 <- cowplot::plot_grid(effect_level_elevation, nrow = 1, ncol = 1, labels = c("g"),
                          label_size = 18, scale=1, label_x = -0.05, label_y = 1.15)

pp2 <- cowplot::plot_grid(effect_level_latitude, nrow = 1, ncol = 1, labels = c("h"),
                          label_size = 18, scale=1, label_x = 0.02, label_y = 1.15)

oo3 <- cowplot::plot_grid(effect_rate_elevation, nrow = 1, ncol = 1, labels = c("i"),
                          label_size = 18, scale=1, label_x = -0.05, label_y = 1.15)

pp3 <- cowplot::plot_grid(effect_rate_latitude, nrow = 1, ncol = 1, labels = c("j"),
                          label_size = 18, scale=1, label_x = 0.02, label_y = 1.15)

oo4 <- cowplot::plot_grid(effect_country_elevation, nrow = 1, ncol = 1, labels = c("k"),
                          label_size = 18, scale=1, label_x = -0.05, label_y = 1.15)

pp4 <- cowplot::plot_grid(effect_country_latitude, nrow = 1, ncol = 1, labels = c("l"),
                          label_size = 18, scale=1, label_x = 0.02, label_y = 1.15)

left_col <- cowplot::plot_grid(
  ii, jj,
  nrow = 2,
  ncol = 1,
  align = "v",
  rel_heights = c(1,1)
)

right_grid <- cowplot::plot_grid(
  kk,   mm,
  oo,   pp,
  oo2,  pp2,
  oo3,  pp3,
  oo4,  pp4,
  nrow = 5,
  ncol = 2,
  align = "hv"
)

p_final <- cowplot::plot_grid(
  left_col, right_grid,
  nrow = 1,
  ncol = 2,
  align = "h",
  rel_widths = c(0.46, 0.54)
) +
  theme(plot.margin = margin(28, 5, 5, 5))

print(p_final)

ggsave("E:/2025/UrbanGreenSpaceWorks/Rcodes/Fig3.pdf", 
       plot = p_final, width = 320, height = 380, units = "mm",dpi=1000)



################################################合并图形###############################################################
# Figure S5: The marginal effect of urbanization on the impacts of geographic stresses on UGS coverage
ii <- cowplot::plot_grid(marginal_elevation, nrow = 1, ncol = 1, labels = c("a"),               
                         label_size = 18, scale=1, label_x = 0.02, label_y = 1.00)

jj <- cowplot::plot_grid(marginal_latitude, nrow = 1, ncol = 1, labels = c("b"),
                         label_size = 18, scale=1, label_x = 0.02, label_y = 1.00)

kk <- cowplot::plot_grid(marginal_continent_elevation, nrow = 1, ncol = 1, labels = c("c"),
                         label_size = 18, scale=1, label_x = -0.05, label_y = 1.15)

mm <- cowplot::plot_grid(marginal_continent_latitude, nrow = 1, ncol = 1, labels = c("d"),
                         label_size = 18, scale=1, label_x = 0.02, label_y = 1.15)

oo <- cowplot::plot_grid(marginal_climate_elevation, nrow = 1, ncol = 1, labels = c("e"),
                         label_size = 18, scale=1, label_x = -0.05, label_y = 1.15)

pp <- cowplot::plot_grid(marginal_climate_latitude, nrow = 1, ncol = 1, labels = c("f"),
                         label_size = 18, scale=1, label_x = 0.02, label_y = 1.15)

oo2 <- cowplot::plot_grid(marginal_level_elevation, nrow = 1, ncol = 1, labels = c("g"),
                          label_size = 18, scale=1, label_x = -0.05, label_y = 1.15)

pp2 <- cowplot::plot_grid(marginal_level_latitude, nrow = 1, ncol = 1, labels = c("h"),
                          label_size = 18, scale=1, label_x = 0.02, label_y = 1.15)

oo3 <- cowplot::plot_grid(marginal_rate_elevation, nrow = 1, ncol = 1, labels = c("i"),
                          label_size = 18, scale=1, label_x = -0.05, label_y = 1.15)

pp3 <- cowplot::plot_grid(marginal_rate_latitude, nrow = 1, ncol = 1, labels = c("j"),
                          label_size = 18, scale=1, label_x = 0.02, label_y = 1.15)

oo4 <- cowplot::plot_grid(marginal_country_elevation, nrow = 1, ncol = 1, labels = c("k"),
                          label_size = 18, scale=1, label_x = -0.05, label_y = 1.15)

pp4 <- cowplot::plot_grid(marginal_country_latitude, nrow = 1, ncol = 1, labels = c("l"),
                          label_size = 18, scale=1, label_x = 0.02, label_y = 1.15)

left_col <- cowplot::plot_grid(
  ii, jj,
  nrow = 2,
  ncol = 1,
  align = "v",
  rel_heights = c(1,1)
)

right_grid <- cowplot::plot_grid(
  kk,   mm,
  oo,   pp,
  oo2,  pp2,
  oo3,  pp3,
  oo4,  pp4,
  nrow = 5,
  ncol = 2,
  align = "hv"
)

p_final <- cowplot::plot_grid(
  left_col, right_grid,
  nrow = 1,
  ncol = 2,
  align = "h",
  rel_widths = c(0.46, 0.54)
) +
  theme(plot.margin = margin(28, 5, 5, 5))

print(p_final)

ggsave("E:/2025/UrbanGreenSpaceWorks/Rcodes/FigS5.pdf", 
       plot = p_final, width = 320, height = 380, units = "mm",dpi=1000)




