
# The effect of urbanization on the impacts of hydroclimatic stresses on UGS coverage.

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



################################################### water deficit ################################################## 
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_WaterDeficit <- copy(merged_dt)
merged_dt_WaterDeficit <- merged_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                 .SDcols = cols_to_center]
merged_dt_WaterDeficit$urban <- factor(merged_dt_WaterDeficit$urban, levels = c("out", "in"))

# linear regression
# Table 1 -- model 1
model_urban <- felm(FVC ~ 1 + WaterDeficit +
                      # WaterDeficit_2 +
                      Seasonalityindex +
                      AverageTemperature +
                      TemperatureRange +
                      WindSpeed +
                      Elevation +
                      Latitude +
                      Precipitation +
                      SoilMoisture +
                      SoilPH + 
                      Longitude + 
                      
                      WaterDeficit:urban +
                      # WaterDeficit_2:urban +
                      
                      GDP_per_capita_PPP + HDI + Population + 
                      ImperviousSurface + HumanSettlement + CityArea
                    |  year + CountryID, data = merged_dt_WaterDeficit)  
summary(model_urban)


######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_out.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

# nonlinear regression
# Table 1 -- model 3
model_urban <- felm(FVC ~ 1 + WaterDeficit +
                      WaterDeficit_2 +
                      Seasonalityindex +
                      AverageTemperature +
                      TemperatureRange +
                      WindSpeed +
                      Elevation +
                      Latitude +
                      Precipitation +
                      SoilMoisture +
                      SoilPH + 
                      Longitude + 
                      
                      WaterDeficit:urban +
                      WaterDeficit_2:urban +
                      
                      GDP_per_capita_PPP + HDI + Population + 
                      ImperviousSurface + HumanSettlement + CityArea
                    |  year + CountryID, data = merged_dt_WaterDeficit)  
summary(model_urban)


######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_out2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

######################################### 绘制效应图与边际效应图 #######################################
coeffs <- coef(model_urban)
beta1 <- coeffs["WaterDeficit"]              
beta2 <- coeffs["WaterDeficit_2"]             
beta1_city <- coeffs["WaterDeficit:urbanin"]  
beta2_city <- coeffs["WaterDeficit_2:urbanin"]
beta1_total <- beta1 + beta1_city
beta2_total <- beta2 + beta2_city
P_min <- quantile(merged_dt_WaterDeficit$WaterDeficit, 0.001, na.rm = TRUE)
P_max <- quantile(merged_dt_WaterDeficit$WaterDeficit, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200) 

# 计算效应值
effect_non_urban <- beta1 * P_values + beta2 * (P_values^2)  
effect_urban <- beta1_total * P_values + beta2_total * (P_values^2)  
vcov_matrix <- vcovHC(model_urban, type = "HC1", cluster = "group", group = merged_dt_WaterDeficit$CountryID)
cov_vars_non_urban <- c("WaterDeficit", "WaterDeficit_2")
vcov_non_urban <- vcov_matrix[cov_vars_non_urban, cov_vars_non_urban]
d_beta1 <- P_values  
d_beta2 <- P_values^2 
var_effect_non_urban <- (d_beta1^2) * vcov_non_urban["WaterDeficit", "WaterDeficit"] +
  (d_beta2^2) * vcov_non_urban["WaterDeficit_2", "WaterDeficit_2"] +
  2 * d_beta1 * d_beta2 * vcov_non_urban["WaterDeficit", "WaterDeficit_2"]
se_effect_non_urban <- sqrt(var_effect_non_urban)
ci_effect_non_urban_lower <- effect_non_urban - 1.96 * se_effect_non_urban
ci_effect_non_urban_upper <- effect_non_urban + 1.96 * se_effect_non_urban
cov_vars_urban <- c("WaterDeficit", "WaterDeficit_2", "WaterDeficit:urbanin", "WaterDeficit_2:urbanin")
vcov_urban <- vcov_matrix[cov_vars_urban, cov_vars_urban]
d_beta1_u <- P_values  
d_beta2_u <- P_values^2  
var_effect_urban <- (d_beta1_u^2) * vcov_urban["WaterDeficit", "WaterDeficit"] +
  (d_beta2_u^2) * vcov_urban["WaterDeficit_2", "WaterDeficit_2"] +
  (d_beta1_u^2) * vcov_urban["WaterDeficit:urbanin", "WaterDeficit:urbanin"] +
  (d_beta2_u^2) * vcov_urban["WaterDeficit_2:urbanin", "WaterDeficit_2:urbanin"] +
  2 * d_beta1_u^2 * vcov_urban["WaterDeficit", "WaterDeficit:urbanin"] +  
  2 * d_beta2_u^2 * vcov_urban["WaterDeficit_2", "WaterDeficit_2:urbanin"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["WaterDeficit", "WaterDeficit_2"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["WaterDeficit", "WaterDeficit_2:urbanin"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["WaterDeficit:urbanin", "WaterDeficit_2"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["WaterDeficit:urbanin", "WaterDeficit_2:urbanin"]  
se_effect_urban <- sqrt(var_effect_urban)
ci_effect_urban_lower <- effect_urban - 1.96 * se_effect_urban
ci_effect_urban_upper <- effect_urban + 1.96 * se_effect_urban

# 计算边际效应 
marginal_non_urban <- beta1 + 2 * beta2 * P_values  
marginal_urban <- beta1_total + 2 * beta2_total * P_values  
d_marg_beta1 <- 1  
d_marg_beta2 <- 2 * P_values  
var_marginal_non_urban <- (d_marg_beta1^2) * vcov_non_urban["WaterDeficit", "WaterDeficit"] +
  (d_marg_beta2^2) * vcov_non_urban["WaterDeficit_2", "WaterDeficit_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_non_urban["WaterDeficit", "WaterDeficit_2"]
se_marginal_non_urban <- sqrt(var_marginal_non_urban)
ci_marginal_non_urban_lower <- marginal_non_urban - 1.96 * se_marginal_non_urban
ci_marginal_non_urban_upper <- marginal_non_urban + 1.96 * se_marginal_non_urban
d_marg_beta1_u <- 1  
d_marg_beta2_u <- 2 * P_values  
var_marginal_urban <- (d_marg_beta1_u^2) * vcov_urban["WaterDeficit", "WaterDeficit"] +
  (d_marg_beta2_u^2) * vcov_urban["WaterDeficit_2", "WaterDeficit_2"] +
  (d_marg_beta1_u^2) * vcov_urban["WaterDeficit:urbanin", "WaterDeficit:urbanin"] +
  (d_marg_beta2_u^2) * vcov_urban["WaterDeficit_2:urbanin", "WaterDeficit_2:urbanin"] +
  2 * d_marg_beta1_u^2 * vcov_urban["WaterDeficit", "WaterDeficit:urbanin"] +  
  2 * d_marg_beta2_u^2 * vcov_urban["WaterDeficit_2", "WaterDeficit_2:urbanin"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["WaterDeficit", "WaterDeficit_2"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["WaterDeficit", "WaterDeficit_2:urbanin"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["WaterDeficit:urbanin", "WaterDeficit_2"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["WaterDeficit:urbanin", "WaterDeficit_2:urbanin"]  
se_marginal_urban <- sqrt(var_marginal_urban)
ci_marginal_urban_lower <- marginal_urban - 1.96 * se_marginal_urban
ci_marginal_urban_upper <- marginal_urban + 1.96 * se_marginal_urban

# 效应值结果
effect_df <- rbind(
  data.frame(
    WaterDeficit = P_values,
    Region = "Non-urban",
    Effect = effect_non_urban,
    SE = se_effect_non_urban,
    CI_Lower = ci_effect_non_urban_lower,
    CI_Upper = ci_effect_non_urban_upper
  ),
  data.frame(
    WaterDeficit = P_values,
    Region = "Urban",
    Effect = effect_urban,
    SE = se_effect_urban,
    CI_Lower = ci_effect_urban_lower,
    CI_Upper = ci_effect_urban_upper
  )
)

# 边际效应结果
marginal_df <- rbind(
  data.frame(
    WaterDeficit = P_values,
    Region = "Non-urban",
    Marginal_Effect = marginal_non_urban,
    SE = se_marginal_non_urban,
    CI_Lower = ci_marginal_non_urban_lower,
    CI_Upper = ci_marginal_non_urban_upper
  ),
  data.frame(
    WaterDeficit = P_values,
    Region = "Urban",
    Marginal_Effect = marginal_urban,
    SE = se_marginal_urban,
    CI_Lower = ci_marginal_urban_lower,
    CI_Upper = ci_marginal_urban_upper
  )
)


# 绘制效应值图
# Figure 2a
urban_colors <- c(
  "Urban" = "#FF7A3C", 
  "Non-urban" = "#00A850"
)
urban_colors_hist <- c(
  "in" = "#FF7A3C", 
  "out" = "#00A850"
)
effect_plot <- ggplot(effect_df, aes(x = WaterDeficit, y = Effect, color = Region)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha = 0.2, color = NA
  ) +
  theme_classic() + 
  labs(
    x = "Water deficit",
    y = "Effect on FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  scale_x_continuous(limits = c(0, 2400), breaks = seq(0, 2400, by =600), 
                     labels = scales::number_format(accuracy = 0.01)) +
  # scale_y_continuous(name = "Effect on FVC", limits = c(-0.4, 0.6),
  #                    breaks = seq(-0.4, 0.6, by = 0.2),
  #                    labels = scales::number_format(accuracy = 0.1)) +
  theme(
    axis.text.y = element_text(color = "black", size = 15),  
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),  
    axis.title.y = element_text(size = 16),  
    legend.position = c(0.2, 0.2), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 14),  
    legend.background = element_blank()
  ) +
  scale_color_manual(values = urban_colors) +
  scale_fill_manual(values = urban_colors) 

print(effect_plot)

# 绘制直方图
hist_plot <- ggplot(merged_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill = urban), binwidth =  5,  alpha = 0.5) + 
  scale_x_continuous(limits = c(0, 2400), breaks = seq(0, 2400, by =600)) +
  scale_y_continuous(limits = c(0.0, 0.006),
                     breaks = seq(0.0, 0.006, by = 0.002),
                     labels = scales::number_format(accuracy = 0.001)) +
  labs(x = "Water deficit (mm)", 
       y = "",
       subtitle = paste0("")) +
  theme_classic() + 
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.title.x = element_text(size = 16), 
    plot.margin = margin(t=-1.32,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text.y = element_text(color = "black", size = 15),  
    axis.text.x = element_text(color = "black", size = 15),
    legend.position = "none") +
  scale_fill_manual(values = urban_colors_hist) 

print(hist_plot)

# 合并
effect_waterdeficit <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(effect_waterdeficit)



# 绘制边际效应图
# Figure S4a
marginal_plot <- ggplot(marginal_df, aes(x = WaterDeficit, y = Marginal_Effect, color = Region)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha = 0.2, color = NA
  ) +
  labs(
    x = "WaterDeficit",
    y = "Marginal effect on FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  theme_classic() +
  scale_x_continuous(limits = c(0, 2400), breaks = seq(0, 2400, by =600)) +
  scale_y_continuous(limits = c(-0.0004, 0.0001),
                     breaks = seq(-0.0004, 0.0001, by = 0.0002),
                     labels = scales::number_format(accuracy = 0.0001)) +
  
  theme(
    axis.text.y = element_text(color = "black", size = 15),  
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),  
    axis.title.y = element_text(size = 16), 
    legend.position = c(0.3,0.9), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 14),  
    legend.background = element_blank()
  ) +
  scale_color_manual(values = urban_colors) +
  scale_fill_manual(values = urban_colors)

print(marginal_plot)

marginal_waterdeficit <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))

print(marginal_waterdeficit)




################################################### Rainfall Seasonalityindex ################################################## 
cols_to_center <- c('WaterDeficit', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
merged_dt_Seasonalityindex <- copy(merged_dt)
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[merged_dt_Seasonalityindex$Precipitation < 1000,]
merged_dt_Seasonalityindex <- merged_dt_Seasonalityindex[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                         .SDcols = cols_to_center]
merged_dt_Seasonalityindex$urban <- factor(merged_dt_Seasonalityindex$urban, levels = c("out", "in"))

# linear regression
# Table 1 -- model 2
model_urban <- felm(FVC ~ 1 + WaterDeficit +
                      Seasonalityindex +
                      # Seasonalityindex_2 +
                      AverageTemperature +
                      TemperatureRange +
                      WindSpeed +
                      Elevation +
                      Latitude +
                      Precipitation +
                      SoilMoisture +
                      SoilPH + 
                      Longitude + 
                      
                      Seasonalityindex:urban +
                      # Seasonalityindex_2:urban +
                      
                      GDP_per_capita_PPP + HDI + Population +
                      ImperviousSurface + HumanSettlement + CityArea
                    |  year + CountryID, data = merged_dt_Seasonalityindex)  
summary(model_urban)

#####################保存模型参数#######################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_out.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

# nonlinear regression
# Table 1 -- model 4
model_urban <- felm(FVC ~ 1 + WaterDeficit +
                      Seasonalityindex +
                      Seasonalityindex_2 +
                      AverageTemperature +
                      TemperatureRange +
                      WindSpeed +
                      Elevation +
                      Latitude +
                      Precipitation +
                      SoilMoisture +
                      SoilPH + 
                      Longitude + 
                      
                      Seasonalityindex:urban +
                      Seasonalityindex_2:urban +
                      
                      GDP_per_capita_PPP + HDI + Population +
                      ImperviousSurface + HumanSettlement + CityArea
                    |  year + CountryID, data = merged_dt_Seasonalityindex)  
summary(model_urban)

#####################保存模型参数#######################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_out2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

######################################### 绘制效应图与边际效应图 #######################################
coeffs <- coef(model_urban)
beta1 <- coeffs["Seasonalityindex"]             
beta2 <- coeffs["Seasonalityindex_2"]             
beta1_city <- coeffs["Seasonalityindex:urbanin"]  
beta2_city <- coeffs["Seasonalityindex_2:urbanin"]
beta1_total <- beta1 + beta1_city
beta2_total <- beta2 + beta2_city
P_min <- quantile(merged_dt_Seasonalityindex$Seasonalityindex, 0.001, na.rm = TRUE)
P_max <- quantile(merged_dt_Seasonalityindex$Seasonalityindex, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)  

# 计算效应值
effect_non_urban <- beta1 * P_values + beta2 * (P_values^2)  
effect_urban <- beta1_total * P_values + beta2_total * (P_values^2)  
vcov_matrix <- vcovHC(model_urban, type = "HC1", cluster = "group", group = merged_dt_Seasonalityindex$CountryID)
cov_vars_non_urban <- c("Seasonalityindex", "Seasonalityindex_2")
vcov_non_urban <- vcov_matrix[cov_vars_non_urban, cov_vars_non_urban]
d_beta1 <- P_values  
d_beta2 <- P_values^2  
var_effect_non_urban <- (d_beta1^2) * vcov_non_urban["Seasonalityindex", "Seasonalityindex"] +
  (d_beta2^2) * vcov_non_urban["Seasonalityindex_2", "Seasonalityindex_2"] +
  2 * d_beta1 * d_beta2 * vcov_non_urban["Seasonalityindex", "Seasonalityindex_2"]
se_effect_non_urban <- sqrt(var_effect_non_urban)
ci_effect_non_urban_lower <- effect_non_urban - 1.96 * se_effect_non_urban
ci_effect_non_urban_upper <- effect_non_urban + 1.96 * se_effect_non_urban
cov_vars_urban <- c("Seasonalityindex", "Seasonalityindex_2", "Seasonalityindex:urbanin", "Seasonalityindex_2:urbanin")
vcov_urban <- vcov_matrix[cov_vars_urban, cov_vars_urban]
d_beta1_u <- P_values  
d_beta2_u <- P_values^2  
var_effect_urban <- (d_beta1_u^2) * vcov_urban["Seasonalityindex", "Seasonalityindex"] +
  (d_beta2_u^2) * vcov_urban["Seasonalityindex_2", "Seasonalityindex_2"] +
  (d_beta1_u^2) * vcov_urban["Seasonalityindex:urbanin", "Seasonalityindex:urbanin"] +
  (d_beta2_u^2) * vcov_urban["Seasonalityindex_2:urbanin", "Seasonalityindex_2:urbanin"] +
  2 * d_beta1_u^2 * vcov_urban["Seasonalityindex", "Seasonalityindex:urbanin"] +  
  2 * d_beta2_u^2 * vcov_urban["Seasonalityindex_2", "Seasonalityindex_2:urbanin"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["Seasonalityindex", "Seasonalityindex_2"] + 
  2 * d_beta1_u * d_beta2_u * vcov_urban["Seasonalityindex", "Seasonalityindex_2:urbanin"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["Seasonalityindex:urbanin", "Seasonalityindex_2"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["Seasonalityindex:urbanin", "Seasonalityindex_2:urbanin"]  
se_effect_urban <- sqrt(var_effect_urban)
ci_effect_urban_lower <- effect_urban - 1.96 * se_effect_urban
ci_effect_urban_upper <- effect_urban + 1.96 * se_effect_urban

# 计算边际效应
marginal_non_urban <- beta1 + 2 * beta2 * P_values  
marginal_urban <- beta1_total + 2 * beta2_total * P_values  
d_marg_beta1 <- 1 
d_marg_beta2 <- 2 * P_values  
var_marginal_non_urban <- (d_marg_beta1^2) * vcov_non_urban["Seasonalityindex", "Seasonalityindex"] +
  (d_marg_beta2^2) * vcov_non_urban["Seasonalityindex_2", "Seasonalityindex_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_non_urban["Seasonalityindex", "Seasonalityindex_2"]
se_marginal_non_urban <- sqrt(var_marginal_non_urban)
ci_marginal_non_urban_lower <- marginal_non_urban - 1.96 * se_marginal_non_urban
ci_marginal_non_urban_upper <- marginal_non_urban + 1.96 * se_marginal_non_urban
d_marg_beta1_u <- 1  
d_marg_beta2_u <- 2 * P_values  
var_marginal_urban <- (d_marg_beta1_u^2) * vcov_urban["Seasonalityindex", "Seasonalityindex"] +
  (d_marg_beta2_u^2) * vcov_urban["Seasonalityindex_2", "Seasonalityindex_2"] +
  (d_marg_beta1_u^2) * vcov_urban["Seasonalityindex:urbanin", "Seasonalityindex:urbanin"] +
  (d_marg_beta2_u^2) * vcov_urban["Seasonalityindex_2:urbanin", "Seasonalityindex_2:urbanin"] +
  2 * d_marg_beta1_u^2 * vcov_urban["Seasonalityindex", "Seasonalityindex:urbanin"] +  
  2 * d_marg_beta2_u^2 * vcov_urban["Seasonalityindex_2", "Seasonalityindex_2:urbanin"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Seasonalityindex", "Seasonalityindex_2"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Seasonalityindex", "Seasonalityindex_2:urbanin"] + 
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Seasonalityindex:urbanin", "Seasonalityindex_2"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Seasonalityindex:urbanin", "Seasonalityindex_2:urbanin"]  
se_marginal_urban <- sqrt(var_marginal_urban)
ci_marginal_urban_lower <- marginal_urban - 1.96 * se_marginal_urban
ci_marginal_urban_upper <- marginal_urban + 1.96 * se_marginal_urban

# 效应值结果
effect_df <- rbind(
  data.frame(
    Seasonalityindex = P_values,
    Region = "Non-urban",
    Effect = effect_non_urban,
    SE = se_effect_non_urban,
    CI_Lower = ci_effect_non_urban_lower,
    CI_Upper = ci_effect_non_urban_upper
  ),
  data.frame(
    Seasonalityindex = P_values,
    Region = "Urban",
    Effect = effect_urban,
    SE = se_effect_urban,
    CI_Lower = ci_effect_urban_lower,
    CI_Upper = ci_effect_urban_upper
  )
)

# 边际效应结果
marginal_df <- rbind(
  data.frame(
    Seasonalityindex = P_values,
    Region = "Non-urban",
    Marginal_Effect = marginal_non_urban,
    SE = se_marginal_non_urban,
    CI_Lower = ci_marginal_non_urban_lower,
    CI_Upper = ci_marginal_non_urban_upper
  ),
  data.frame(
    Seasonalityindex = P_values,
    Region = "Urban",
    Marginal_Effect = marginal_urban,
    SE = se_marginal_urban,
    CI_Lower = ci_marginal_urban_lower,
    CI_Upper = ci_marginal_urban_upper
  )
)

# 绘制效应值图 
# Figure 2b
effect_plot <- ggplot(effect_df, aes(x = Seasonalityindex, y = Effect, color = Region)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha = 0.2, color = NA
  ) +
  theme_classic() + 
  labs(
    x = "Seasonality index",
    y = "Effect on FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02), 
                     labels = scales::number_format(accuracy = 0.01)) +
  # scale_y_continuous(name = "Effect on FVC", limits = c(-0.4, 0.6),
  #                    breaks = seq(-0.4, 0.6, by = 0.2),
  #                    labels = scales::number_format(accuracy = 0.1)) +
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

#  绘制直方图
hist_plot <- ggplot(merged_dt_Seasonalityindex, aes(x = Seasonalityindex)) +
  geom_histogram(aes(y = ..density.., fill = urban), binwidth =  0.0001,  alpha = 0.5) + 
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02)) +
  # scale_y_continuous(name = "  ", limits = c(0.0, 0.006),
  #                    breaks = seq(0.0, 0.006, by = 0.002),
  #                    labels = scales::number_format(accuracy = 0.001)) +
  
  labs(x = "Seasonality index",  
       y = "",
       subtitle = paste0("")) +
  theme_classic() + 
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.title.x = element_text(size = 16),  
    plot.margin = margin(t=-1.32,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text.y = element_text(color = "black", size = 15),  
    axis.text.x = element_text(color = "black", size = 15),
    legend.position = "none") +
  scale_fill_manual(values = urban_colors_hist)

print(hist_plot)

# 合并
effect_seasonalityindex <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(effect_seasonalityindex)


# 绘制边际效应图
# Figure S4b
marginal_plot <- ggplot(marginal_df, aes(x = Seasonalityindex, y = Marginal_Effect, color = Region)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
    alpha = 0.2, color = NA
  ) +
  labs(
    x = "Seasonalityindex",
    y = "Marginal effect on FVC",
    color = "",
    fill = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  theme_classic() +
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02)) +
  # scale_y_continuous(name = "Marginal effects on FVC", limits = c(-0.0004, 0.0004),
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


marginal_seasonalityindex <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(marginal_seasonalityindex)



################################################### 城市内异质性分析 ################################################## 

##################################################### continent ####################################################
################################################### water deficit ################################################## 

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



################################################### water deficit ################################################## 
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
urban_in_dt_WaterDeficit <- copy(urban_in_dt)
urban_in_dt_WaterDeficit <- urban_in_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                     .SDcols = cols_to_center]
urban_in_dt_WaterDeficit$continent <- factor(urban_in_dt_WaterDeficit$continent, levels = c( "Asia",
                                                                                             "NorthAmerica",
                                                                                             "Africa",
                                                                                             "Europe",
                                                                                             "Oceania",
                                                                                             "SouthAmerica"))
# linear regression
# Table S3 -- model 1
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         # WaterDeficit_2 +
                         Seasonalityindex +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Latitude +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         WaterDeficit:continent +
                         # WaterDeficit_2:continent +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea 
                       |  year + CountryID, data = urban_in_dt_WaterDeficit) 
summary(model_urban_in)

#####################保存模型参数#######################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_continent.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

# nonlinear regression
# Table S3 -- model 5
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         WaterDeficit_2 +
                         Seasonalityindex +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Latitude +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         WaterDeficit:continent +
                         WaterDeficit_2:continent +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea 
                       |  year + CountryID, data = urban_in_dt_WaterDeficit) 
summary(model_urban_in)

#####################保存模型参数#######################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_continent2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

######################################### 绘制效应图与边际效应图 #######################################
beta1 <- coef(model_urban_in)["WaterDeficit"]               
beta2 <- coef(model_urban_in)["WaterDeficit_2"]             
continents <- c("NorthAmerica","Africa",  "Europe",  "Oceania", "SouthAmerica") 
beta1_continent <- sapply(continents, function(cont) {
  coef(model_urban_in)[paste0("WaterDeficit:continent", cont)]
})
beta2_continent <- sapply(continents, function(cont) {
  coef(model_urban_in)[paste0("WaterDeficit_2:continent", cont)]
})
beta1_total <- beta1 + beta1_continent
beta2_total <- beta2 + beta2_continent
names(beta1_total) <- continents
names(beta2_total) <- continents
P_min <- quantile(urban_in_dt_WaterDeficit$WaterDeficit, 0.001, na.rm = TRUE) 
P_max <- quantile(urban_in_dt_WaterDeficit$WaterDeficit, 0.999, na.rm = TRUE) 
P_values <- seq(P_min, P_max, length.out = 200)  

# 计算效应值
effect_base <- beta1 * P_values + beta2 * (P_values^2)
effect_continent <- lapply(continents, function(cont) {
  beta1_total[[cont]] * P_values + beta2_total[[cont]] * (P_values^2)
})
names(effect_continent) <- continents
vcov_matrix <- vcovHC(model_urban_in, type = "HC1", cluster = "group",
                      group = urban_in_dt_WaterDeficit$CountryID)
cov_vars_base <- c("WaterDeficit", "WaterDeficit_2")
vcov_base <- vcov_matrix[cov_vars_base, cov_vars_base]
d_beta1 <- P_values  
d_beta2 <- P_values^2  
var_effect_base <- (d_beta1^2) * vcov_base["WaterDeficit", "WaterDeficit"] +
  (d_beta2^2) * vcov_base["WaterDeficit_2", "WaterDeficit_2"] +
  2 * d_beta1 * d_beta2 * vcov_base["WaterDeficit", "WaterDeficit_2"]
se_effect_base <- sqrt(var_effect_base)
ci_effect_base_lower <- effect_base - 1.96 * se_effect_base
ci_effect_base_upper <- effect_base + 1.96 * se_effect_base
var_effect_continent <- list()
se_effect_continent <- list()
ci_effect_continent_lower <- list()
ci_effect_continent_upper <- list()
for (cont in continents) {
  cov_vars <- c("WaterDeficit", "WaterDeficit_2", 
                paste0("WaterDeficit:continent", cont), 
                paste0("WaterDeficit_2:continent", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  var_effect <- (d_beta1^2) * vcov_cont["WaterDeficit", "WaterDeficit"] +
    (d_beta2^2) * vcov_cont["WaterDeficit_2", "WaterDeficit_2"] +
    (d_beta1^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_beta2^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_beta1^2 * vcov_cont["WaterDeficit", cov_vars[3]] +  
    2 * d_beta2^2 * vcov_cont["WaterDeficit_2", cov_vars[4]] +  
    2 * d_beta1 * d_beta2 * vcov_cont["WaterDeficit", "WaterDeficit_2"] +  
    2 * d_beta1 * d_beta2 * vcov_cont["WaterDeficit", cov_vars[4]] +  
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], "WaterDeficit_2"] +  
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
var_marginal_base <- (d_marg_beta1^2) * vcov_base["WaterDeficit", "WaterDeficit"] +
  (d_marg_beta2^2) * vcov_base["WaterDeficit_2", "WaterDeficit_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_base["WaterDeficit", "WaterDeficit_2"]
se_marginal_base <- sqrt(var_marginal_base)
ci_marginal_base_lower <- marginal_base - 1.96 * se_marginal_base
ci_marginal_base_upper <- marginal_base + 1.96 * se_marginal_base
var_marginal_continent <- list()
se_marginal_continent <- list()
ci_marginal_continent_lower <- list()
ci_marginal_continent_upper <- list()
for (cont in continents) {
  cov_vars <- c("WaterDeficit", "WaterDeficit_2",
                paste0("WaterDeficit:continent", cont),
                paste0("WaterDeficit_2:continent", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  d_marg_beta1_u <- 1  
  d_marg_beta2_u <- 2 * P_values  
  var_marginal <- (d_marg_beta1_u^2) * vcov_cont["WaterDeficit", "WaterDeficit"] +
    (d_marg_beta2_u^2) * vcov_cont["WaterDeficit_2", "WaterDeficit_2"] +
    (d_marg_beta1_u^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_marg_beta2_u^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_marg_beta1_u^2 * vcov_cont["WaterDeficit", cov_vars[3]] +  
    2 * d_marg_beta2_u^2 * vcov_cont["WaterDeficit_2", cov_vars[4]] + 
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["WaterDeficit", "WaterDeficit_2"] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["WaterDeficit", cov_vars[4]] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], "WaterDeficit_2"] + 
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], cov_vars[4]]  
  var_marginal_continent[[cont]] <- var_marginal
  se_marginal_continent[[cont]] <- sqrt(var_marginal)
  ci_marginal_continent_lower[[cont]] <- marginal_continent[[cont]] - 1.96 * se_marginal_continent[[cont]]
  ci_marginal_continent_upper[[cont]] <- marginal_continent[[cont]] + 1.96 * se_marginal_continent[[cont]]
}


# 效应值结果
effect_df <- data.frame(
  WaterDeficit = P_values,
  Continent = "基准组",
  Effect = effect_base,
  SE = se_effect_base,
  CI_Lower = ci_effect_base_lower,
  CI_Upper = ci_effect_base_upper
)
for (cont in continents) {
  temp_df <- data.frame(
    WaterDeficit = P_values,
    Continent = cont,
    Effect = effect_continent[[cont]],
    SE = se_effect_continent[[cont]],
    CI_Lower = ci_effect_continent_lower[[cont]],
    CI_Upper = ci_effect_continent_upper[[cont]]
  )
  effect_df <- rbind(effect_df, temp_df)
}

# 边际效应结果
marginal_df <- data.frame(
  WaterDeficit = P_values,
  Continent = "基准组",
  Marginal_Effect = marginal_base,
  SE = se_marginal_base,
  CI_Lower = ci_marginal_base_lower,
  CI_Upper = ci_marginal_base_upper
)
for (cont in continents) {
  temp_df <- data.frame(
    WaterDeficit = P_values,
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


# 计算各大洲的实际最值
continent_limits <- tibble(
  Continent = c("Asia", "Africa", "North America", "South America", "Europe", "Oceania"),
  min_WD = c(
    min(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$continent == "Asia", ]$WaterDeficit, na.rm = TRUE),
    min(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$continent == "Africa", ]$WaterDeficit, na.rm = TRUE),
    min(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$continent == "NorthAmerica", ]$WaterDeficit, na.rm = TRUE),
    min(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$continent == "SouthAmerica", ]$WaterDeficit, na.rm = TRUE),
    min(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$continent == "Europe", ]$WaterDeficit, na.rm = TRUE),
    min(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$continent == "Oceania", ]$WaterDeficit, na.rm = TRUE)
  ),
  max_WD = c(
    max(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$continent == "Asia", ]$WaterDeficit, na.rm = TRUE),
    max(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$continent == "Africa", ]$WaterDeficit, na.rm = TRUE),
    max(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$continent == "NorthAmerica", ]$WaterDeficit, na.rm = TRUE),
    max(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$continent == "SouthAmerica", ]$WaterDeficit, na.rm = TRUE),
    max(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$continent == "Europe", ]$WaterDeficit, na.rm = TRUE),
    max(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$continent == "Oceania", ]$WaterDeficit, na.rm = TRUE)
  )
)

# 按各大洲最值截取
effect_df <- effect_df %>%
  left_join(continent_limits, by = "Continent") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)
marginal_df <- marginal_df %>%
  left_join(continent_limits, by = "Continent") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)

# 绘制效应值图 
continent_colors <- c(
  "Asia" = "#3498db",
  "Africa" = "#95a5a6",
  "North America" = "#2ecc71",
  "South America" = "#9b59b6",
  "Europe" = "#e74c3c",
  "Oceania" = "#f39c12"
)
effect_df$Continent <- factor(effect_df$Continent, levels = c("Asia", "Africa", "North America", "South America", "Europe", "Oceania"))

# 效应值图
# Figure 2c
effect_plot <- ggplot(effect_df, aes(x = WaterDeficit, y = Effect, color = Continent)) +
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
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.3, 0.1),
                     breaks = seq(-0.3, 0.1, by = 0.1),
                     labels = scales::number_format(accuracy = 0.1)) +
  theme(
    axis.text.y = element_text(color = "black", size = 13),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.35, 0.4),
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),
    legend.background = element_blank(),
  ) +
  scale_color_manual(values = continent_colors) +
  scale_fill_manual(values = continent_colors)

print(effect_plot)

# 直方图
hist_plot <- ggplot(urban_in_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill=continent), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600), name = "") +
  scale_y_continuous(limits = c(0, 0.023), breaks = seq(0, 0.023, by = 0.01)) +
  labs(y = "") +
  theme_classic() +
  theme(
    plot.margin = margin(t=-0.75,b=0.5,l=0.15,r=0.2, unit = "cm"),
    axis.text = element_text(color = "black", size = 13),
    axis.title = element_text(size = 14),
    legend.position = "none"
  ) + 
  scale_fill_manual(values = c("Asia" = "#3498db",
                               "Africa" = "#95a5a6",
                               "NorthAmerica" = "#2ecc71",
                               "SouthAmerica" = "#9b59b6",
                               "Europe" = "#e74c3c",
                               "Oceania" = "#f39c12"
  ))

print(hist_plot)

# 合并
effect_continent_waterdeficit <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_continent_waterdeficit)



# 绘制边际效应图 
# Figure S4c
marginal_df$Continent <- factor(marginal_df$Continent, levels = c("Asia", "Africa", "North America", "South America", "Europe", "Oceania"))

marginal_plot <- ggplot(marginal_df, aes(x = WaterDeficit, y = Marginal_Effect, color = Continent)) +
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
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.0004, 0.0007), 
                     breaks = seq(-0.0004, 0.0007, by = 0.0004),
                     labels = scales::number_format(accuracy = 0.0001)) +
  theme(
    axis.text.y = element_text(color = "black", size = 12),  
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),  
    axis.title.y = element_text(size = 12),  
    legend.position = c(0.45, 0.8), 
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),  
    legend.background = element_blank()
  ) +
  scale_color_manual(values = continent_colors) +
  scale_fill_manual(values = continent_colors)

print(marginal_plot)

marginal_continent_waterdeficit <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_continent_waterdeficit)




##################################################### continent ####################################################
################################################### Rainfall Seasonalityindex ######################################
cols_to_center <- c('WaterDeficit', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
urban_in_dt_Seasonalityindex <- copy(urban_in_dt)
urban_in_dt_Seasonalityindex <- urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$Precipitation < 1000,]
urban_in_dt_Seasonalityindex <- urban_in_dt_Seasonalityindex[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                             .SDcols = cols_to_center]
urban_in_dt_Seasonalityindex$continent <- factor(urban_in_dt_Seasonalityindex$continent, levels = c( "Asia",
                                                                                                     "NorthAmerica",
                                                                                                     "Africa",
                                                                                                     "Europe",
                                                                                                     "Oceania",
                                                                                                     "SouthAmerica"))
# linear regression
# Table S3 -- model 2
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         Seasonalityindex +
                         # Seasonalityindex_2 +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Latitude +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         Seasonalityindex:continent +
                         # Seasonalityindex_2:continent +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea 
                       |  year + CountryID, data = urban_in_dt_Seasonalityindex) 
summary(model_urban_in)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_continent.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

# nonlinear regression
# Table S3 -- model 6
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         Seasonalityindex +
                         Seasonalityindex_2 +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Latitude +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         Seasonalityindex:continent +
                         Seasonalityindex_2:continent +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea 
                       |  year + CountryID, data = urban_in_dt_Seasonalityindex) 
summary(model_urban_in)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_continent2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)


######################################### 绘制效应图与边际效应图 #######################################
beta1 <- coef(model_urban_in)["Seasonalityindex"] 
beta2 <- coef(model_urban_in)["Seasonalityindex_2"]            
continents <- c("NorthAmerica","Africa",  "Europe",  "Oceania", "SouthAmerica") 
beta1_continent <- sapply(continents, function(cont) {
  coef(model_urban_in)[paste0("Seasonalityindex:continent", cont)]
})
beta2_continent <- sapply(continents, function(cont) {
  coef(model_urban_in)[paste0("Seasonalityindex_2:continent", cont)]
})
beta1_total <- beta1 + beta1_continent
beta2_total <- beta2 + beta2_continent
names(beta1_total) <- continents
names(beta2_total) <- continents
P_min <- quantile(urban_in_dt_Seasonalityindex$Seasonalityindex, 0.001, na.rm = TRUE)
P_max <- quantile(urban_in_dt_Seasonalityindex$Seasonalityindex, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)  

# 计算效应值
effect_base <- beta1 * P_values + beta2 * (P_values^2)
effect_continent <- lapply(continents, function(cont) {
  beta1_total[[cont]] * P_values + beta2_total[[cont]] * (P_values^2)
})
names(effect_continent) <- continents
vcov_matrix <- vcovHC(model_urban_in, type = "HC1", cluster = "group",
                      group = urban_in_dt_Seasonalityindex$CountryID)
cov_vars_base <- c("Seasonalityindex", "Seasonalityindex_2")
vcov_base <- vcov_matrix[cov_vars_base, cov_vars_base]
d_beta1 <- P_values  
d_beta2 <- P_values^2  
var_effect_base <- (d_beta1^2) * vcov_base["Seasonalityindex", "Seasonalityindex"] +
  (d_beta2^2) * vcov_base["Seasonalityindex_2", "Seasonalityindex_2"] +
  2 * d_beta1 * d_beta2 * vcov_base["Seasonalityindex", "Seasonalityindex_2"]
se_effect_base <- sqrt(var_effect_base)
ci_effect_base_lower <- effect_base - 1.96 * se_effect_base
ci_effect_base_upper <- effect_base + 1.96 * se_effect_base
var_effect_continent <- list()
se_effect_continent <- list()
ci_effect_continent_lower <- list()
ci_effect_continent_upper <- list()
for (cont in continents) {
  cov_vars <- c("Seasonalityindex", "Seasonalityindex_2", 
                paste0("Seasonalityindex:continent", cont), 
                paste0("Seasonalityindex_2:continent", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  var_effect <- (d_beta1^2) * vcov_cont["Seasonalityindex", "Seasonalityindex"] +
    (d_beta2^2) * vcov_cont["Seasonalityindex_2", "Seasonalityindex_2"] +
    (d_beta1^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_beta2^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_beta1^2 * vcov_cont["Seasonalityindex", cov_vars[3]] +  
    2 * d_beta2^2 * vcov_cont["Seasonalityindex_2", cov_vars[4]] +  
    2 * d_beta1 * d_beta2 * vcov_cont["Seasonalityindex", "Seasonalityindex_2"] +  
    2 * d_beta1 * d_beta2 * vcov_cont["Seasonalityindex", cov_vars[4]] +  
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], "Seasonalityindex_2"] +  
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
var_marginal_base <- (d_marg_beta1^2) * vcov_base["Seasonalityindex", "Seasonalityindex"] +
  (d_marg_beta2^2) * vcov_base["Seasonalityindex_2", "Seasonalityindex_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_base["Seasonalityindex", "Seasonalityindex_2"]
se_marginal_base <- sqrt(var_marginal_base)
ci_marginal_base_lower <- marginal_base - 1.96 * se_marginal_base
ci_marginal_base_upper <- marginal_base + 1.96 * se_marginal_base
var_marginal_continent <- list()
se_marginal_continent <- list()
ci_marginal_continent_lower <- list()
ci_marginal_continent_upper <- list()
for (cont in continents) {
  cov_vars <- c("Seasonalityindex", "Seasonalityindex_2",
                paste0("Seasonalityindex:continent", cont),
                paste0("Seasonalityindex_2:continent", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  d_marg_beta1_u <- 1 
  d_marg_beta2_u <- 2 * P_values 
  var_marginal <- (d_marg_beta1_u^2) * vcov_cont["Seasonalityindex", "Seasonalityindex"] +
    (d_marg_beta2_u^2) * vcov_cont["Seasonalityindex_2", "Seasonalityindex_2"] +
    (d_marg_beta1_u^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_marg_beta2_u^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_marg_beta1_u^2 * vcov_cont["Seasonalityindex", cov_vars[3]] +  
    2 * d_marg_beta2_u^2 * vcov_cont["Seasonalityindex_2", cov_vars[4]] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["Seasonalityindex", "Seasonalityindex_2"] + 
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["Seasonalityindex", cov_vars[4]] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], "Seasonalityindex_2"] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], cov_vars[4]] 
  var_marginal_continent[[cont]] <- var_marginal
  se_marginal_continent[[cont]] <- sqrt(var_marginal)
  ci_marginal_continent_lower[[cont]] <- marginal_continent[[cont]] - 1.96 * se_marginal_continent[[cont]]
  ci_marginal_continent_upper[[cont]] <- marginal_continent[[cont]] + 1.96 * se_marginal_continent[[cont]]
}

# 效应值结果
effect_df <- data.frame(
  Seasonalityindex = P_values,
  Continent = "基准组",
  Effect = effect_base,
  SE = se_effect_base,
  CI_Lower = ci_effect_base_lower,
  CI_Upper = ci_effect_base_upper
)
for (cont in continents) {
  temp_df <- data.frame(
    Seasonalityindex = P_values,
    Continent = cont,
    Effect = effect_continent[[cont]],
    SE = se_effect_continent[[cont]],
    CI_Lower = ci_effect_continent_lower[[cont]],
    CI_Upper = ci_effect_continent_upper[[cont]]
  )
  effect_df <- rbind(effect_df, temp_df)
}

# 边际效应结果
marginal_df <- data.frame(
  Seasonalityindex = P_values,
  Continent = "基准组",
  Marginal_Effect = marginal_base,
  SE = se_marginal_base,
  CI_Lower = ci_marginal_base_lower,
  CI_Upper = ci_marginal_base_upper
)
for (cont in continents) {
  temp_df <- data.frame(
    Seasonalityindex = P_values,
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

# 计算各大洲的实际最值 
continent_limits <- tibble(
  Continent = c("Asia", "Africa", "North America", "South America", "Europe", "Oceania"),
  min_WD = c(
    min(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$continent == "Asia", ]$Seasonalityindex, na.rm = TRUE),
    min(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$continent == "Africa", ]$Seasonalityindex, na.rm = TRUE),
    min(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$continent == "NorthAmerica", ]$Seasonalityindex, na.rm = TRUE),
    min(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$continent == "SouthAmerica", ]$Seasonalityindex, na.rm = TRUE),
    min(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$continent == "Europe", ]$Seasonalityindex, na.rm = TRUE),
    min(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$continent == "Oceania", ]$Seasonalityindex, na.rm = TRUE)
  ),
  max_WD = c(
    max(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$continent == "Asia", ]$Seasonalityindex, na.rm = TRUE),
    max(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$continent == "Africa", ]$Seasonalityindex, na.rm = TRUE),
    max(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$continent == "NorthAmerica", ]$Seasonalityindex, na.rm = TRUE),
    max(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$continent == "SouthAmerica", ]$Seasonalityindex, na.rm = TRUE),
    max(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$continent == "Europe", ]$Seasonalityindex, na.rm = TRUE),
    max(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$continent == "Oceania", ]$Seasonalityindex, na.rm = TRUE)
  )
)
print(continent_limits)

# 按各大洲最值截取
effect_df <- effect_df %>%
  left_join(continent_limits, by = "Continent") %>%
  filter(Seasonalityindex >= min_WD & Seasonalityindex <= max_WD) %>%
  select(-min_WD, -max_WD)
marginal_df <- marginal_df %>%
  left_join(continent_limits, by = "Continent") %>%
  filter(Seasonalityindex >= min_WD & Seasonalityindex <= max_WD) %>%
  select(-min_WD, -max_WD)


# 效应值图
# Figure 2d
effect_plot <- ggplot(effect_df, aes(x = Seasonalityindex, y = Effect, color = Continent)) +
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
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.14, 0.1),
                     breaks = seq(-0.14, 0.1, by = 0.07),
                     labels = scales::number_format(accuracy = 0.01)) +
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

# 直方图
hist_plot <- ggplot(urban_in_dt_Seasonalityindex, aes(x = Seasonalityindex)) +
  geom_histogram(aes(y = ..density.., fill=continent), binwidth = 0.0001, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02), name = "") +
  labs(y = "") +
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
effect_continent_seasonalityindex <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_continent_seasonalityindex)


# 绘制边际效应图 
# Figure S4d
marginal_plot <- ggplot(marginal_df, aes(x = Seasonalityindex, y = Marginal_Effect, color = Continent)) +
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
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-12, 12),
                     breaks = seq(-12, 12, by = 6),
                     labels = scales::number_format(accuracy = 1)) +
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

marginal_continent_seasonalityindex <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_continent_seasonalityindex)




################################################### 城市内异质性分析 ################################################## 
##################################################### climate zone ################################################## 
################################################### water deficit ######################################
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
urban_in_dt_WaterDeficit <- copy(urban_in_dt)
urban_in_dt_WaterDeficit <- urban_in_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                     .SDcols = cols_to_center]
urban_in_dt_WaterDeficit$ClimateZone <- factor(urban_in_dt_WaterDeficit$ClimateZone, levels = c("Temperate",
                                                                                                "Tropical",
                                                                                                "Arid",
                                                                                                "Cold"))
# linear regression
# Table S4 -- model 1
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         # WaterDeficit_2 +
                         Seasonalityindex +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Latitude +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         WaterDeficit:ClimateZone +
                         # WaterDeficit_2:ClimateZone +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea 
                       |  year + CountryID, data = urban_in_dt_WaterDeficit)  
summary(model_urban_in)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_climate.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

# nonlinear regression
# Table S4 -- model 5
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         WaterDeficit_2 +
                         Seasonalityindex +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Latitude +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         WaterDeficit:ClimateZone +
                         WaterDeficit_2:ClimateZone +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea 
                       |  year + CountryID, data = urban_in_dt_WaterDeficit)  
summary(model_urban_in)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_climate2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

######################################### 绘制效应图与边际效应图 #######################################
beta1 <- coef(model_urban_in)["WaterDeficit"]           
beta2 <- coef(model_urban_in)["WaterDeficit_2"]             
continents <- c("Tropical","Arid",  "Cold") 
beta1_continent <- sapply(continents, function(cont) {
  coef(model_urban_in)[paste0("WaterDeficit:ClimateZone", cont)]
})
beta2_continent <- sapply(continents, function(cont) {
  coef(model_urban_in)[paste0("WaterDeficit_2:ClimateZone", cont)]
})
beta1_total <- beta1 + beta1_continent
beta2_total <- beta2 + beta2_continent
names(beta1_total) <- continents
names(beta2_total) <- continents
P_min <- quantile(urban_in_dt_WaterDeficit$WaterDeficit, 0.001, na.rm = TRUE)
P_max <- quantile(urban_in_dt_WaterDeficit$WaterDeficit, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200) 

# 计算效应值
effect_base <- beta1 * P_values + beta2 * (P_values^2)
effect_continent <- lapply(continents, function(cont) {
  beta1_total[[cont]] * P_values + beta2_total[[cont]] * (P_values^2)
})
names(effect_continent) <- continents
vcov_matrix <- vcovHC(model_urban_in, type = "HC1", cluster = "group",
                      group = urban_in_dt_WaterDeficit$CountryID)
cov_vars_base <- c("WaterDeficit", "WaterDeficit_2")
vcov_base <- vcov_matrix[cov_vars_base, cov_vars_base]
d_beta1 <- P_values  
d_beta2 <- P_values^2  
var_effect_base <- (d_beta1^2) * vcov_base["WaterDeficit", "WaterDeficit"] +
  (d_beta2^2) * vcov_base["WaterDeficit_2", "WaterDeficit_2"] +
  2 * d_beta1 * d_beta2 * vcov_base["WaterDeficit", "WaterDeficit_2"]
se_effect_base <- sqrt(var_effect_base)
ci_effect_base_lower <- effect_base - 1.96 * se_effect_base
ci_effect_base_upper <- effect_base + 1.96 * se_effect_base
var_effect_continent <- list()
se_effect_continent <- list()
ci_effect_continent_lower <- list()
ci_effect_continent_upper <- list()

for (cont in continents) {
  cov_vars <- c("WaterDeficit", "WaterDeficit_2", 
                paste0("WaterDeficit:ClimateZone", cont), 
                paste0("WaterDeficit_2:ClimateZone", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  var_effect <- (d_beta1^2) * vcov_cont["WaterDeficit", "WaterDeficit"] +
    (d_beta2^2) * vcov_cont["WaterDeficit_2", "WaterDeficit_2"] +
    (d_beta1^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_beta2^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_beta1^2 * vcov_cont["WaterDeficit", cov_vars[3]] +  
    2 * d_beta2^2 * vcov_cont["WaterDeficit_2", cov_vars[4]] +  
    2 * d_beta1 * d_beta2 * vcov_cont["WaterDeficit", "WaterDeficit_2"] + 
    2 * d_beta1 * d_beta2 * vcov_cont["WaterDeficit", cov_vars[4]] +  
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], "WaterDeficit_2"] +  
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
var_marginal_base <- (d_marg_beta1^2) * vcov_base["WaterDeficit", "WaterDeficit"] +
  (d_marg_beta2^2) * vcov_base["WaterDeficit_2", "WaterDeficit_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_base["WaterDeficit", "WaterDeficit_2"]
se_marginal_base <- sqrt(var_marginal_base)
ci_marginal_base_lower <- marginal_base - 1.96 * se_marginal_base
ci_marginal_base_upper <- marginal_base + 1.96 * se_marginal_base
var_marginal_continent <- list()
se_marginal_continent <- list()
ci_marginal_continent_lower <- list()
ci_marginal_continent_upper <- list()
for (cont in continents) {
  cov_vars <- c("WaterDeficit", "WaterDeficit_2",
                paste0("WaterDeficit:ClimateZone", cont),
                paste0("WaterDeficit_2:ClimateZone", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  d_marg_beta1_u <- 1  
  d_marg_beta2_u <- 2 * P_values  
  var_marginal <- (d_marg_beta1_u^2) * vcov_cont["WaterDeficit", "WaterDeficit"] +
    (d_marg_beta2_u^2) * vcov_cont["WaterDeficit_2", "WaterDeficit_2"] +
    (d_marg_beta1_u^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_marg_beta2_u^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_marg_beta1_u^2 * vcov_cont["WaterDeficit", cov_vars[3]] +  
    2 * d_marg_beta2_u^2 * vcov_cont["WaterDeficit_2", cov_vars[4]] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["WaterDeficit", "WaterDeficit_2"] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["WaterDeficit", cov_vars[4]] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], "WaterDeficit_2"] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], cov_vars[4]]  
  
  var_marginal_continent[[cont]] <- var_marginal
  se_marginal_continent[[cont]] <- sqrt(var_marginal)
  ci_marginal_continent_lower[[cont]] <- marginal_continent[[cont]] - 1.96 * se_marginal_continent[[cont]]
  ci_marginal_continent_upper[[cont]] <- marginal_continent[[cont]] + 1.96 * se_marginal_continent[[cont]]
}


# 效应值结果
effect_df <- data.frame(
  WaterDeficit = P_values,
  Climate = "基准组",
  Effect = effect_base,
  SE = se_effect_base,
  CI_Lower = ci_effect_base_lower,
  CI_Upper = ci_effect_base_upper
)
for (cont in continents) {
  temp_df <- data.frame(
    WaterDeficit = P_values,
    Climate = cont,
    Effect = effect_continent[[cont]],
    SE = se_effect_continent[[cont]],
    CI_Lower = ci_effect_continent_lower[[cont]],
    CI_Upper = ci_effect_continent_upper[[cont]]
  )
  effect_df <- rbind(effect_df, temp_df)
}

# 边际效应结果
marginal_df <- data.frame(
  WaterDeficit = P_values,
  Climate = "基准组",
  Marginal_Effect = marginal_base,
  SE = se_marginal_base,
  CI_Lower = ci_marginal_base_lower,
  CI_Upper = ci_marginal_base_upper
)
for (cont in continents) {
  temp_df <- data.frame(
    WaterDeficit = P_values,
    Climate = cont,
    Marginal_Effect = marginal_continent[[cont]],
    SE = se_marginal_continent[[cont]],
    CI_Lower = ci_marginal_continent_lower[[cont]],
    CI_Upper = ci_marginal_continent_upper[[cont]]
  )
  marginal_df <- rbind(marginal_df, temp_df)
}

effect_df$Climate <- factor(effect_df$Climate, 
                            levels = c("基准组", continents),
                            labels = c("Temperate", "Tropical", "Arid",  "Cold")) 
marginal_df$Climate <- factor(marginal_df$Climate,
                              levels = c("基准组", continents),
                              labels = c("Temperate", "Tropical", "Arid",  "Cold"))


# 计算各气候区的实际最值 
climate_limits <- tibble(
  Climate = c("Temperate", "Tropical", "Arid",  "Cold"),
  min_WD = c(
    min(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$ClimateZone == "Temperate", ]$WaterDeficit, na.rm = TRUE),
    min(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$ClimateZone == "Tropical", ]$WaterDeficit, na.rm = TRUE),
    min(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$ClimateZone == "Arid", ]$WaterDeficit, na.rm = TRUE),
    min(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$ClimateZone == "Cold", ]$WaterDeficit, na.rm = TRUE)
  ),
  max_WD = c(
    max(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$ClimateZone == "Temperate", ]$WaterDeficit, na.rm = TRUE),
    max(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$ClimateZone == "Tropical", ]$WaterDeficit, na.rm = TRUE),
    max(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$ClimateZone == "Arid", ]$WaterDeficit, na.rm = TRUE),
    max(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$ClimateZone == "Cold", ]$WaterDeficit, na.rm = TRUE)
  )
)
print(climate_limits)

# 按各气候区最值截取
effect_df <- effect_df %>%
  left_join(climate_limits, by = "Climate") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)
marginal_df <- marginal_df %>%
  left_join(climate_limits, by = "Climate") %>%
  filter(WaterDecifit >= min_WD & WaterDecifit <= max_WD) %>%
  select(-min_WD, -max_WD)


# 绘制效应值图 
# Figure 2e
climate_colors <- c(
  "Temperate" = "#2ecc71",
  "Tropical" = "#e74c3c",
  "Arid" = "#95a5a6",
  "Cold" = "#3498db"
)
effect_df$Climate <- factor(effect_df$Climate, levels = c( "Tropical", "Temperate",  "Cold","Arid"))
effect_plot <- ggplot(effect_df, aes(x = WaterDeficit, y = Effect, color = Climate)) +
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
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.1, 0.2),
                     breaks = seq(-0.1, 0.2, by = 0.1),
                     labels = scales::number_format(accuracy = 0.1)) +
  theme(
    axis.text.y = element_text(color = "black", size = 13),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.35, 0.85),
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),
    legend.background = element_blank()
  ) +
  scale_color_manual(values = climate_colors) +
  scale_fill_manual(values = climate_colors)

print(effect_plot)

# 直方图
hist_plot <- ggplot(urban_in_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill=ClimateZone), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600), name = "") +
  scale_y_continuous(limits = c(0, 0.01), breaks = seq(0, 0.01, by =0.004)) +
  labs(y = "") +
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
effect_climate_waterdeficit <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_climate_waterdeficit)



# 绘制边际效应图 
# Figure S4e
marginal_df$Climate <- factor(marginal_df$Climate, levels = c( "Tropical", "Temperate",  "Cold","Arid"))

marginal_plot <- ggplot(marginal_df, aes(x = WaterDeficit, y = Marginal_Effect, color = Climate)) +
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
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.0003, 0.0003),
                     breaks = seq(-0.0003, 0.0003, by = 0.0003),
                     labels = scales::number_format(accuracy = 0.0001)) +
  theme(
    axis.text.y = element_text(color = "black", size = 12),  
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),  
    axis.title.y = element_text(size = 12),  
    legend.position = c(0.35,0.85), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),
    legend.background = element_blank()
  ) +
  scale_color_manual(values = climate_colors) +
  scale_fill_manual(values = climate_colors)

print(marginal_plot)

marginal_climate_waterdeficit <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_climate_waterdeficit)




##################################################### climate zone ####################################################
################################################### Rainfall Seasonalityindex ######################################
cols_to_center <- c('WaterDeficit', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
urban_in_dt_Seasonalityindex <- copy(urban_in_dt)
urban_in_dt_Seasonalityindex <- urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$Precipitation < 1000,]
urban_in_dt_Seasonalityindex <- urban_in_dt_Seasonalityindex[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                             .SDcols = cols_to_center]
urban_in_dt_Seasonalityindex$ClimateZone <- factor(urban_in_dt_Seasonalityindex$ClimateZone, levels = c("Temperate",
                                                                                                        "Tropical",
                                                                                                        "Arid",
                                                                                                        "Cold"))
# linear regression
# Table S4 -- model 2
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         Seasonalityindex +
                         # Seasonalityindex_2 +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Latitude +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         Seasonalityindex:ClimateZone +
                         # Seasonalityindex_2:ClimateZone +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea 
                       |  year + CountryID, data = urban_in_dt_Seasonalityindex)  
summary(model_urban_in)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_climate.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

# nonlinear regression
# Table S4 -- model 6
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         Seasonalityindex +
                         Seasonalityindex_2 +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Latitude +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         Seasonalityindex:ClimateZone +
                         Seasonalityindex_2:ClimateZone +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea 
                       |  year + CountryID, data = urban_in_dt_Seasonalityindex)  
summary(model_urban_in)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_climate2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)

######################################### 绘制效应图与边际效应图 #######################################
beta1 <- coef(model_urban_in)["Seasonalityindex"]            
beta2 <- coef(model_urban_in)["Seasonalityindex_2"]   
continents <- c("Tropical","Arid",  "Cold") 
beta1_continent <- sapply(continents, function(cont) {
  coef(model_urban_in)[paste0("Seasonalityindex:ClimateZone", cont)]
})
beta2_continent <- sapply(continents, function(cont) {
  coef(model_urban_in)[paste0("Seasonalityindex_2:ClimateZone", cont)]
})
beta1_total <- beta1 + beta1_continent
beta2_total <- beta2 + beta2_continent
names(beta1_total) <- continents
names(beta2_total) <- continents
P_min <- quantile(urban_in_dt_Seasonalityindex$Seasonalityindex, 0.001, na.rm = TRUE)
P_max <- quantile(urban_in_dt_Seasonalityindex$Seasonalityindex, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)  

# 计算效应值
effect_base <- beta1 * P_values + beta2 * (P_values^2)
effect_continent <- lapply(continents, function(cont) {
  beta1_total[[cont]] * P_values + beta2_total[[cont]] * (P_values^2)
})
names(effect_continent) <- continents
vcov_matrix <- vcovHC(model_urban_in, type = "HC1", cluster = "group",
                      group = urban_in_dt_Seasonalityindex$CountryID)
cov_vars_base <- c("Seasonalityindex", "Seasonalityindex_2")
vcov_base <- vcov_matrix[cov_vars_base, cov_vars_base]
d_beta1 <- P_values 
d_beta2 <- P_values^2  
var_effect_base <- (d_beta1^2) * vcov_base["Seasonalityindex", "Seasonalityindex"] +
  (d_beta2^2) * vcov_base["Seasonalityindex_2", "Seasonalityindex_2"] +
  2 * d_beta1 * d_beta2 * vcov_base["Seasonalityindex", "Seasonalityindex_2"]
se_effect_base <- sqrt(var_effect_base)
ci_effect_base_lower <- effect_base - 1.96 * se_effect_base
ci_effect_base_upper <- effect_base + 1.96 * se_effect_base
var_effect_continent <- list()
se_effect_continent <- list()
ci_effect_continent_lower <- list()
ci_effect_continent_upper <- list()
for (cont in continents) {
  cov_vars <- c("Seasonalityindex", "Seasonalityindex_2", 
                paste0("Seasonalityindex:ClimateZone", cont), 
                paste0("Seasonalityindex_2:ClimateZone", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  var_effect <- (d_beta1^2) * vcov_cont["Seasonalityindex", "Seasonalityindex"] +
    (d_beta2^2) * vcov_cont["Seasonalityindex_2", "Seasonalityindex_2"] +
    (d_beta1^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_beta2^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_beta1^2 * vcov_cont["Seasonalityindex", cov_vars[3]] +  
    2 * d_beta2^2 * vcov_cont["Seasonalityindex_2", cov_vars[4]] +  
    2 * d_beta1 * d_beta2 * vcov_cont["Seasonalityindex", "Seasonalityindex_2"] +  
    2 * d_beta1 * d_beta2 * vcov_cont["Seasonalityindex", cov_vars[4]] + 
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], "Seasonalityindex_2"] +  
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
var_marginal_base <- (d_marg_beta1^2) * vcov_base["Seasonalityindex", "Seasonalityindex"] +
  (d_marg_beta2^2) * vcov_base["Seasonalityindex_2", "Seasonalityindex_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_base["Seasonalityindex", "Seasonalityindex_2"]
se_marginal_base <- sqrt(var_marginal_base)
ci_marginal_base_lower <- marginal_base - 1.96 * se_marginal_base
ci_marginal_base_upper <- marginal_base + 1.96 * se_marginal_base
var_marginal_continent <- list()
se_marginal_continent <- list()
ci_marginal_continent_lower <- list()
ci_marginal_continent_upper <- list()
for (cont in continents) {
  cov_vars <- c("Seasonalityindex", "Seasonalityindex_2",
                paste0("Seasonalityindex:ClimateZone", cont),
                paste0("Seasonalityindex_2:ClimateZone", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  d_marg_beta1_u <- 1  
  d_marg_beta2_u <- 2 * P_values  
  var_marginal <- (d_marg_beta1_u^2) * vcov_cont["Seasonalityindex", "Seasonalityindex"] +
    (d_marg_beta2_u^2) * vcov_cont["Seasonalityindex_2", "Seasonalityindex_2"] +
    (d_marg_beta1_u^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_marg_beta2_u^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_marg_beta1_u^2 * vcov_cont["Seasonalityindex", cov_vars[3]] +  
    2 * d_marg_beta2_u^2 * vcov_cont["Seasonalityindex_2", cov_vars[4]] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["Seasonalityindex", "Seasonalityindex_2"] + 
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["Seasonalityindex", cov_vars[4]] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], "Seasonalityindex_2"] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], cov_vars[4]] 
  var_marginal_continent[[cont]] <- var_marginal
  se_marginal_continent[[cont]] <- sqrt(var_marginal)
  ci_marginal_continent_lower[[cont]] <- marginal_continent[[cont]] - 1.96 * se_marginal_continent[[cont]]
  ci_marginal_continent_upper[[cont]] <- marginal_continent[[cont]] + 1.96 * se_marginal_continent[[cont]]
}


# 效应值结果
effect_df <- data.frame(
  Seasonalityindex = P_values,
  Climate = "基准组",
  Effect = effect_base,
  SE = se_effect_base,
  CI_Lower = ci_effect_base_lower,
  CI_Upper = ci_effect_base_upper
)

for (cont in continents) {
  temp_df <- data.frame(
    Seasonalityindex = P_values,
    Climate = cont,
    Effect = effect_continent[[cont]],
    SE = se_effect_continent[[cont]],
    CI_Lower = ci_effect_continent_lower[[cont]],
    CI_Upper = ci_effect_continent_upper[[cont]]
  )
  effect_df <- rbind(effect_df, temp_df)
}

# 边际效应结果
marginal_df <- data.frame(
  Seasonalityindex = P_values,
  Climate = "基准组",
  Marginal_Effect = marginal_base,
  SE = se_marginal_base,
  CI_Lower = ci_marginal_base_lower,
  CI_Upper = ci_marginal_base_upper
)

for (cont in continents) {
  temp_df <- data.frame(
    Seasonalityindex = P_values,
    Climate = cont,
    Marginal_Effect = marginal_continent[[cont]],
    SE = se_marginal_continent[[cont]],
    CI_Lower = ci_marginal_continent_lower[[cont]],
    CI_Upper = ci_marginal_continent_upper[[cont]]
  )
  marginal_df <- rbind(marginal_df, temp_df)
}

effect_df$Climate <- factor(effect_df$Climate, 
                            levels = c("基准组", continents),
                            labels = c("Temperate", "Tropical", "Arid",  "Cold"))
marginal_df$Climate <- factor(marginal_df$Climate,
                              levels = c("基准组", continents),
                              labels = c("Temperate", "Tropical", "Arid",  "Cold"))

# 计算各气候区的实际最值 
climate_limits <- tibble(
  Climate = c("Temperate", "Tropical", "Arid",  "Cold"),
  min_WD = c(
    min(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$ClimateZone == "Temperate", ]$Seasonalityindex, na.rm = TRUE),
    min(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$ClimateZone == "Tropical", ]$Seasonalityindex, na.rm = TRUE),
    min(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$ClimateZone == "Arid", ]$Seasonalityindex, na.rm = TRUE),
    min(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$ClimateZone == "Cold", ]$Seasonalityindex, na.rm = TRUE)
  ),
  max_WD = c(
    max(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$ClimateZone == "Temperate", ]$Seasonalityindex, na.rm = TRUE),
    max(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$ClimateZone == "Tropical", ]$Seasonalityindex, na.rm = TRUE),
    max(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$ClimateZone == "Arid", ]$Seasonalityindex, na.rm = TRUE),
    max(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$ClimateZone == "Cold", ]$Seasonalityindex, na.rm = TRUE)
  )
)

print(climate_limits)

# 按各气候区最值截取
effect_df <- effect_df %>%
  left_join(climate_limits, by = "Climate") %>%
  filter(Seasonalityindex >= min_WD & Seasonalityindex <= max_WD) %>%
  select(-min_WD, -max_WD)
marginal_df <- marginal_df %>%
  left_join(climate_limits, by = "Climate") %>%
  filter(Seasonalityindex >= min_WD & Seasonalityindex <= max_WD) %>%
  select(-min_WD, -max_WD)

# 效应值图
# Figure 2f
effect_plot <- ggplot(effect_df, aes(x = Seasonalityindex, y = Effect, color = Climate)) +
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
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02), 
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

# 直方图
hist_plot <- ggplot(urban_in_dt_Seasonalityindex, aes(x = Seasonalityindex)) +
  geom_histogram(aes(y = ..density.., fill=ClimateZone), binwidth = 0.0001, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02), name = "") +
  scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, by =200)) +
  labs(y = "") +
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
effect_climate_seasonalityindex <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_climate_seasonalityindex)


# 绘制边际效应图 
# Figure S4f
marginal_plot <- ggplot(marginal_df, aes(x = Seasonalityindex, y = Marginal_Effect, color = Climate)) +
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
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-1, 1),
                     breaks = seq(-1, 1, by = 0.5),
                     labels = scales::number_format(accuracy = 0.1)) +
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

marginal_climate_seasonalityindex <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_climate_seasonalityindex)







################################################### 城市内异质性分析 ################################################## 
##################################################### development level ################################################## 
################################################### water deficit ######################################
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
urban_in_dt_WaterDeficit <- copy(urban_in_dt)
urban_in_dt_WaterDeficit <- urban_in_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                     .SDcols = cols_to_center]
urban_in_dt_WaterDeficit$DevelopmentLevel <- factor(urban_in_dt_WaterDeficit$DevelopmentLevel, levels = c("developing", "developed"))

# linear regression
# Table S5 -- model 1
model_urban <- felm(FVC ~ 1 + WaterDeficit +
                      # WaterDeficit_2 +
                      Seasonalityindex +
                      AverageTemperature +
                      TemperatureRange +
                      WindSpeed +
                      Elevation +
                      Latitude +
                      Precipitation +
                      SoilMoisture +
                      SoilPH + 
                      Longitude + 
                      
                      WaterDeficit:DevelopmentLevel +
                      # WaterDeficit_2:DevelopmentLevel +
                      
                      GDP_per_capita_PPP + HDI + Population + 
                      ImperviousSurface + HumanSettlement + CityArea   
                    |  year + CountryID, data = urban_in_dt_WaterDeficit)      

# 输出结果
summary(model_urban)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_level.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  


# nonlinear regression
# Table S5 -- model 5
model_urban <- felm(FVC ~ 1 + WaterDeficit +
                      WaterDeficit_2 +
                      Seasonalityindex +
                      AverageTemperature +
                      TemperatureRange +
                      WindSpeed +
                      Elevation +
                      Latitude +
                      Precipitation +
                      SoilMoisture +
                      SoilPH + 
                      Longitude + 
                      
                      WaterDeficit:DevelopmentLevel +
                      WaterDeficit_2:DevelopmentLevel +
                      
                      GDP_per_capita_PPP + HDI + Population + 
                      ImperviousSurface + HumanSettlement + CityArea   
                    |  year + CountryID, data = urban_in_dt_WaterDeficit)      

# 输出结果
summary(model_urban)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_level2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

######################################### 绘制效应图与边际效应图 #######################################
coeffs <- coef(model_urban)
beta1 <- coeffs["WaterDeficit"]               
beta2 <- coeffs["WaterDeficit_2"]             
beta1_city <- coeffs["WaterDeficit:DevelopmentLeveldeveloped"]  
beta2_city <- coeffs["WaterDeficit_2:DevelopmentLeveldeveloped"]
beta1_total <- beta1 + beta1_city
beta2_total <- beta2 + beta2_city
P_min <- quantile(urban_in_dt_WaterDeficit$WaterDeficit, 0.001, na.rm = TRUE)
P_max <- quantile(urban_in_dt_WaterDeficit$WaterDeficit, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)  # 生成100个等间隔点

# 计算效应值
effect_non_urban <- beta1 * P_values + beta2 * (P_values^2)  
effect_urban <- beta1_total * P_values + beta2_total * (P_values^2)  
vcov_matrix <- vcovHC(model_urban, type = "HC1", cluster = "group", group = urban_in_dt_WaterDeficit$CountryID)
cov_vars_non_urban <- c("WaterDeficit", "WaterDeficit_2")
vcov_non_urban <- vcov_matrix[cov_vars_non_urban, cov_vars_non_urban]
d_beta1 <- P_values  
d_beta2 <- P_values^2  
var_effect_non_urban <- (d_beta1^2) * vcov_non_urban["WaterDeficit", "WaterDeficit"] +
  (d_beta2^2) * vcov_non_urban["WaterDeficit_2", "WaterDeficit_2"] +
  2 * d_beta1 * d_beta2 * vcov_non_urban["WaterDeficit", "WaterDeficit_2"]
se_effect_non_urban <- sqrt(var_effect_non_urban)
ci_effect_non_urban_lower <- effect_non_urban - 1.96 * se_effect_non_urban
ci_effect_non_urban_upper <- effect_non_urban + 1.96 * se_effect_non_urban
cov_vars_urban <- c("WaterDeficit", "WaterDeficit_2", "WaterDeficit:DevelopmentLeveldeveloped", "WaterDeficit_2:DevelopmentLeveldeveloped")
vcov_urban <- vcov_matrix[cov_vars_urban, cov_vars_urban]
d_beta1_u <- P_values  
d_beta2_u <- P_values^2  
var_effect_urban <- (d_beta1_u^2) * vcov_urban["WaterDeficit", "WaterDeficit"] +
  (d_beta2_u^2) * vcov_urban["WaterDeficit_2", "WaterDeficit_2"] +
  (d_beta1_u^2) * vcov_urban["WaterDeficit:DevelopmentLeveldeveloped", "WaterDeficit:DevelopmentLeveldeveloped"] +
  (d_beta2_u^2) * vcov_urban["WaterDeficit_2:DevelopmentLeveldeveloped", "WaterDeficit_2:DevelopmentLeveldeveloped"] +
  2 * d_beta1_u^2 * vcov_urban["WaterDeficit", "WaterDeficit:DevelopmentLeveldeveloped"] +  
  2 * d_beta2_u^2 * vcov_urban["WaterDeficit_2", "WaterDeficit_2:DevelopmentLeveldeveloped"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["WaterDeficit", "WaterDeficit_2"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["WaterDeficit", "WaterDeficit_2:DevelopmentLeveldeveloped"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["WaterDeficit:DevelopmentLeveldeveloped", "WaterDeficit_2"] +
  2 * d_beta1_u * d_beta2_u * vcov_urban["WaterDeficit:DevelopmentLeveldeveloped", "WaterDeficit_2:DevelopmentLeveldeveloped"]  
se_effect_urban <- sqrt(var_effect_urban)
ci_effect_urban_lower <- effect_urban - 1.96 * se_effect_urban
ci_effect_urban_upper <- effect_urban + 1.96 * se_effect_urban

# 计算边际效应 
marginal_non_urban <- beta1 + 2 * beta2 * P_values  
marginal_urban <- beta1_total + 2 * beta2_total * P_values  
d_marg_beta1 <- 1  
d_marg_beta2 <- 2 * P_values  
var_marginal_non_urban <- (d_marg_beta1^2) * vcov_non_urban["WaterDeficit", "WaterDeficit"] +
  (d_marg_beta2^2) * vcov_non_urban["WaterDeficit_2", "WaterDeficit_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_non_urban["WaterDeficit", "WaterDeficit_2"]
se_marginal_non_urban <- sqrt(var_marginal_non_urban)
ci_marginal_non_urban_lower <- marginal_non_urban - 1.96 * se_marginal_non_urban
ci_marginal_non_urban_upper <- marginal_non_urban + 1.96 * se_marginal_non_urban
d_marg_beta1_u <- 1  
d_marg_beta2_u <- 2 * P_values  
var_marginal_urban <- (d_marg_beta1_u^2) * vcov_urban["WaterDeficit", "WaterDeficit"] +
  (d_marg_beta2_u^2) * vcov_urban["WaterDeficit_2", "WaterDeficit_2"] +
  (d_marg_beta1_u^2) * vcov_urban["WaterDeficit:DevelopmentLeveldeveloped", "WaterDeficit:DevelopmentLeveldeveloped"] +
  (d_marg_beta2_u^2) * vcov_urban["WaterDeficit_2:DevelopmentLeveldeveloped", "WaterDeficit_2:DevelopmentLeveldeveloped"] +
  2 * d_marg_beta1_u^2 * vcov_urban["WaterDeficit", "WaterDeficit:DevelopmentLeveldeveloped"] +  
  2 * d_marg_beta2_u^2 * vcov_urban["WaterDeficit_2", "WaterDeficit_2:DevelopmentLeveldeveloped"] + 
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["WaterDeficit", "WaterDeficit_2"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["WaterDeficit", "WaterDeficit_2:DevelopmentLeveldeveloped"] + 
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["WaterDeficit:DevelopmentLeveldeveloped", "WaterDeficit_2"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["WaterDeficit:DevelopmentLeveldeveloped", "WaterDeficit_2:DevelopmentLeveldeveloped"]  
se_marginal_urban <- sqrt(var_marginal_urban)
ci_marginal_urban_lower <- marginal_urban - 1.96 * se_marginal_urban
ci_marginal_urban_upper <- marginal_urban + 1.96 * se_marginal_urban

# 效应值结果
effect_df <- rbind(
  data.frame(
    WaterDeficit = P_values,
    Region = "Developing",
    Effect = effect_non_urban,
    SE = se_effect_non_urban,
    CI_Lower = ci_effect_non_urban_lower,
    CI_Upper = ci_effect_non_urban_upper
  ),
  data.frame(
    WaterDeficit = P_values,
    Region = "Developed",
    Effect = effect_urban,
    SE = se_effect_urban,
    CI_Lower = ci_effect_urban_lower,
    CI_Upper = ci_effect_urban_upper
  )
)

# 边际效应结果
marginal_df <- rbind(
  data.frame(
    WaterDeficit = P_values,
    Region = "Developing",
    Marginal_Effect = marginal_non_urban,
    SE = se_marginal_non_urban,
    CI_Lower = ci_marginal_non_urban_lower,
    CI_Upper = ci_marginal_non_urban_upper
  ),
  data.frame(
    WaterDeficit = P_values,
    Region = "Developed",
    Marginal_Effect = marginal_urban,
    SE = se_marginal_urban,
    CI_Lower = ci_marginal_urban_lower,
    CI_Upper = ci_marginal_urban_upper
  )
)



# 计算各发展水平的实际最值
level_limits <- tibble(
  Region = c("Developing", "Developed"),
  min_WD = c(
    min(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$DevelopmentLevel == "developing", ]$WaterDeficit, na.rm = TRUE),
    min(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$DevelopmentLevel == "developed", ]$WaterDeficit, na.rm = TRUE)
  ),
  max_WD = c(
    max(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$DevelopmentLevel == "developing", ]$WaterDeficit, na.rm = TRUE),
    max(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$DevelopmentLevel == "developed", ]$WaterDeficit, na.rm = TRUE)
  )
)

print(level_limits)

# 按各发展水平最值截取
effect_df <- effect_df %>%
  left_join(level_limits, by = "Region") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)
marginal_df <- marginal_df %>%
  left_join(level_limits, by = "Region") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)



# 绘制效应值图 
# Figure 2g
level_colors <- c(
  "Developed" = "#0071bc", 
  "Developing" = "#d95218"
)
effect_plot <- ggplot(effect_df, aes(x = WaterDeficit, y = Effect, color = Region)) +
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
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.08, 0.07),
                     breaks = seq(-0.08, 0.07, by = 0.04),
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.text.y = element_text(color = "black", size = 13),  
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),  
    axis.title.y = element_text(size = 14),  
    legend.position = c(0.35, 0.9), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),  
    legend.background = element_blank()
  ) +
  scale_color_manual(values = level_colors) +
  scale_fill_manual(values = level_colors) 

print(effect_plot)



# 直方图
hist_plot <- ggplot(urban_in_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill = DevelopmentLevel), binwidth =  5,  alpha = 0.5) + 
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(0.0, 0.01),
                     breaks = seq(0.0, 0.01, by = 0.004),
                     labels = scales::number_format(accuracy = 0.001)) +
  
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
                                 "developing" = "#d95218")) #+

print(hist_plot)

# 合并
effect_level_waterdeficit <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(effect_level_waterdeficit)

# 绘制边际效应图 
# Figure S4g
marginal_plot <- ggplot(marginal_df, aes(x = WaterDeficit, y = Marginal_Effect, color = Region)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = CI_Lower, ymax = CI_Upper, fill = Region),
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
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600)) +
  scale_y_continuous(limits = c(-0.00008, 0.00007),
                     breaks = seq(-0.00008, 0.00007, by = 0.00004),
                     labels = scales::number_format(accuracy = 0.00001)) +
  
  theme(
    axis.text.y = element_text(color = "black", size = 12),  
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),  
    axis.title.y = element_text(size = 12),  
    legend.position = c(0.45,0.9), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),  
    legend.background = element_blank()
  ) +
  scale_color_manual(values = level_colors) +
  scale_fill_manual(values = level_colors)

print(marginal_plot)

marginal_level_waterdeficit <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(marginal_level_waterdeficit)


##################################################### development level ####################################################
################################################### Rainfall Seasonalityindex ######################################
cols_to_center <- c('WaterDeficit', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
urban_in_dt_Seasonalityindex <- copy(urban_in_dt)
urban_in_dt_Seasonalityindex <- urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$Precipitation < 1000,]
urban_in_dt_Seasonalityindex <- urban_in_dt_Seasonalityindex[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                             .SDcols = cols_to_center]
urban_in_dt_Seasonalityindex$DevelopmentLevel <- factor(urban_in_dt_Seasonalityindex$DevelopmentLevel, levels = c("developing", "developed"))

# linear regression
# Table S5 -- model 2
model_urban <- felm(FVC ~ 1 + WaterDeficit +
                      Seasonalityindex +
                      # Seasonalityindex_2 +
                      AverageTemperature +
                      TemperatureRange +
                      WindSpeed +
                      Elevation +
                      Latitude +
                      Precipitation +
                      SoilMoisture +
                      SoilPH + 
                      Longitude + 
                      
                      Seasonalityindex:DevelopmentLevel +
                      # Seasonalityindex_2:DevelopmentLevel +
                      
                      GDP_per_capita_PPP + HDI + Population + 
                      ImperviousSurface + HumanSettlement + CityArea  
                    |  year + CountryID, data = urban_in_dt_Seasonalityindex)      

# 输出结果
summary(model_urban)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_level.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  


# nonlinear regression
# Table S5 -- model 6
model_urban <- felm(FVC ~ 1 + WaterDeficit +
                      Seasonalityindex +
                      Seasonalityindex_2 +
                      AverageTemperature +
                      TemperatureRange +
                      WindSpeed +
                      Elevation +
                      Latitude +
                      Precipitation +
                      SoilMoisture +
                      SoilPH + 
                      Longitude + 
                      
                      Seasonalityindex:DevelopmentLevel +
                      Seasonalityindex_2:DevelopmentLevel +
                      
                      GDP_per_capita_PPP + HDI + Population + 
                      ImperviousSurface + HumanSettlement + CityArea  
                    |  year + CountryID, data = urban_in_dt_Seasonalityindex)      

# 输出结果
summary(model_urban)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_level2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

######################################### 绘制效应图与边际效应图 #######################################
coeffs <- coef(model_urban)
beta1 <- coeffs["Seasonalityindex"]               
beta2 <- coeffs["Seasonalityindex_2"]             
beta1_city <- coeffs["Seasonalityindex:DevelopmentLeveldeveloped"]  
beta2_city <- coeffs["Seasonalityindex_2:DevelopmentLeveldeveloped"]
beta1_total <- beta1 + beta1_city
beta2_total <- beta2 + beta2_city
P_min <- quantile(urban_in_dt_Seasonalityindex$Seasonalityindex, 0.001, na.rm = TRUE)
P_max <- quantile(urban_in_dt_Seasonalityindex$Seasonalityindex, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)  

# 计算效应值
effect_non_urban <- beta1 * P_values + beta2 * (P_values^2)  
effect_urban <- beta1_total * P_values + beta2_total * (P_values^2)  
vcov_matrix <- vcovHC(model_urban, type = "HC1", cluster = "group", group = urban_in_dt_Seasonalityindex$CountryID)
cov_vars_non_urban <- c("Seasonalityindex", "Seasonalityindex_2")
vcov_non_urban <- vcov_matrix[cov_vars_non_urban, cov_vars_non_urban]
d_beta1 <- P_values 
d_beta2 <- P_values^2  
var_effect_non_urban <- (d_beta1^2) * vcov_non_urban["Seasonalityindex", "Seasonalityindex"] +
  (d_beta2^2) * vcov_non_urban["Seasonalityindex_2", "Seasonalityindex_2"] +
  2 * d_beta1 * d_beta2 * vcov_non_urban["Seasonalityindex", "Seasonalityindex_2"]
se_effect_non_urban <- sqrt(var_effect_non_urban)
ci_effect_non_urban_lower <- effect_non_urban - 1.96 * se_effect_non_urban
ci_effect_non_urban_upper <- effect_non_urban + 1.96 * se_effect_non_urban
cov_vars_urban <- c("Seasonalityindex", "Seasonalityindex_2", "Seasonalityindex:DevelopmentLeveldeveloped", "Seasonalityindex_2:DevelopmentLeveldeveloped")
vcov_urban <- vcov_matrix[cov_vars_urban, cov_vars_urban]
d_beta1_u <- P_values  
d_beta2_u <- P_values^2  
var_effect_urban <- (d_beta1_u^2) * vcov_urban["Seasonalityindex", "Seasonalityindex"] +
  (d_beta2_u^2) * vcov_urban["Seasonalityindex_2", "Seasonalityindex_2"] +
  (d_beta1_u^2) * vcov_urban["Seasonalityindex:DevelopmentLeveldeveloped", "Seasonalityindex:DevelopmentLeveldeveloped"] +
  (d_beta2_u^2) * vcov_urban["Seasonalityindex_2:DevelopmentLeveldeveloped", "Seasonalityindex_2:DevelopmentLeveldeveloped"] +
  2 * d_beta1_u^2 * vcov_urban["Seasonalityindex", "Seasonalityindex:DevelopmentLeveldeveloped"] +  
  2 * d_beta2_u^2 * vcov_urban["Seasonalityindex_2", "Seasonalityindex_2:DevelopmentLeveldeveloped"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["Seasonalityindex", "Seasonalityindex_2"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["Seasonalityindex", "Seasonalityindex_2:DevelopmentLeveldeveloped"] + 
  2 * d_beta1_u * d_beta2_u * vcov_urban["Seasonalityindex:DevelopmentLeveldeveloped", "Seasonalityindex_2"] +  
  2 * d_beta1_u * d_beta2_u * vcov_urban["Seasonalityindex:DevelopmentLeveldeveloped", "Seasonalityindex_2:DevelopmentLeveldeveloped"]  
se_effect_urban <- sqrt(var_effect_urban)
ci_effect_urban_lower <- effect_urban - 1.96 * se_effect_urban
ci_effect_urban_upper <- effect_urban + 1.96 * se_effect_urban

# 计算边际效应 
marginal_non_urban <- beta1 + 2 * beta2 * P_values 
marginal_urban <- beta1_total + 2 * beta2_total * P_values  
d_marg_beta1 <- 1  
d_marg_beta2 <- 2 * P_values  
var_marginal_non_urban <- (d_marg_beta1^2) * vcov_non_urban["Seasonalityindex", "Seasonalityindex"] +
  (d_marg_beta2^2) * vcov_non_urban["Seasonalityindex_2", "Seasonalityindex_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_non_urban["Seasonalityindex", "Seasonalityindex_2"]
se_marginal_non_urban <- sqrt(var_marginal_non_urban)
ci_marginal_non_urban_lower <- marginal_non_urban - 1.96 * se_marginal_non_urban
ci_marginal_non_urban_upper <- marginal_non_urban + 1.96 * se_marginal_non_urban
d_marg_beta1_u <- 1  
d_marg_beta2_u <- 2 * P_values  
var_marginal_urban <- (d_marg_beta1_u^2) * vcov_urban["Seasonalityindex", "Seasonalityindex"] +
  (d_marg_beta2_u^2) * vcov_urban["Seasonalityindex_2", "Seasonalityindex_2"] +
  (d_marg_beta1_u^2) * vcov_urban["Seasonalityindex:DevelopmentLeveldeveloped", "Seasonalityindex:DevelopmentLeveldeveloped"] +
  (d_marg_beta2_u^2) * vcov_urban["Seasonalityindex_2:DevelopmentLeveldeveloped", "Seasonalityindex_2:DevelopmentLeveldeveloped"] +
  2 * d_marg_beta1_u^2 * vcov_urban["Seasonalityindex", "Seasonalityindex:DevelopmentLeveldeveloped"] +  
  2 * d_marg_beta2_u^2 * vcov_urban["Seasonalityindex_2", "Seasonalityindex_2:DevelopmentLeveldeveloped"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Seasonalityindex", "Seasonalityindex_2"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Seasonalityindex", "Seasonalityindex_2:DevelopmentLeveldeveloped"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Seasonalityindex:DevelopmentLeveldeveloped", "Seasonalityindex_2"] +  
  2 * d_marg_beta1_u * d_marg_beta2_u * vcov_urban["Seasonalityindex:DevelopmentLeveldeveloped", "Seasonalityindex_2:DevelopmentLeveldeveloped"]  
se_marginal_urban <- sqrt(var_marginal_urban)
ci_marginal_urban_lower <- marginal_urban - 1.96 * se_marginal_urban
ci_marginal_urban_upper <- marginal_urban + 1.96 * se_marginal_urban

# 效应值结果
effect_df <- rbind(
  data.frame(
    Seasonalityindex = P_values,
    Region = "Developing",
    Effect = effect_non_urban,
    SE = se_effect_non_urban,
    CI_Lower = ci_effect_non_urban_lower,
    CI_Upper = ci_effect_non_urban_upper
  ),
  data.frame(
    Seasonalityindex = P_values,
    Region = "Developed",
    Effect = effect_urban,
    SE = se_effect_urban,
    CI_Lower = ci_effect_urban_lower,
    CI_Upper = ci_effect_urban_upper
  )
)

# 边际效应结果
marginal_df <- rbind(
  data.frame(
    Seasonalityindex = P_values,
    Region = "Developing",
    Marginal_Effect = marginal_non_urban,
    SE = se_marginal_non_urban,
    CI_Lower = ci_marginal_non_urban_lower,
    CI_Upper = ci_marginal_non_urban_upper
  ),
  data.frame(
    Seasonalityindex = P_values,
    Region = "Developed",
    Marginal_Effect = marginal_urban,
    SE = se_marginal_urban,
    CI_Lower = ci_marginal_urban_lower,
    CI_Upper = ci_marginal_urban_upper
  )
)


# 计算各发展水平的实际最值 
level_limits <- tibble(
  Region = c("Developing", "Developed"),
  min_WD = c(
    min(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$DevelopmentLevel == "developing", ]$Seasonalityindex, na.rm = TRUE),
    min(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$DevelopmentLevel == "developed", ]$Seasonalityindex, na.rm = TRUE)
  ),
  max_WD = c(
    max(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$DevelopmentLevel == "developing", ]$Seasonalityindex, na.rm = TRUE),
    quantile(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$DevelopmentLevel == "developed", ]$Seasonalityindex, probs = 0.99, na.rm = TRUE)
  )
)

print(level_limits)

# 按各发展水平最值截取
effect_df <- effect_df %>%
  left_join(level_limits, by = "Region") %>%
  filter(Seasonalityindex >= min_WD & Seasonalityindex <= max_WD) %>%
  select(-min_WD, -max_WD)
marginal_df <- marginal_df %>%
  left_join(level_limits, by = "Region") %>%
  filter(Seasonalityindex >= min_WD & Seasonalityindex <= max_WD) %>%
  select(-min_WD, -max_WD)


# 绘制效应值图
# Figure 2h
effect_plot <- ggplot(effect_df, aes(x = Seasonalityindex, y = Effect, color = Region)) +
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
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.06, 0.04),
                     breaks = seq(-0.06, 0.04, by = 0.03),
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



# 直方图
hist_plot <- ggplot(urban_in_dt_Seasonalityindex, aes(x = Seasonalityindex)) +
  geom_histogram(aes(y = ..density.., fill = DevelopmentLevel), binwidth =  0.0001,  alpha = 0.5) + 
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02)) +
  scale_y_continuous(limits = c(0, 400), breaks = seq(0, 400, by =150)) +
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
effect_level_seasonalityindex <-ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(effect_level_seasonalityindex)

# 绘制边际效应图
# Figure S4h
marginal_plot <- ggplot(marginal_df, aes(x = Seasonalityindex, y = Marginal_Effect, color = Region)) +
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
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02)) +
  scale_y_continuous(limits = c(-1, 1.5),
                     breaks = seq(-1, 1.5, by = 1),
                     labels = scales::number_format(accuracy = 0.1)) +
  
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

marginal_level_seasonalityindex <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2,1))
print(marginal_level_seasonalityindex)


################################################### 城市内异质性分析 ################################################## 
##################################################### expansion rate ################################################## 
################################################### water deficit ######################################
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
urban_in_dt_WaterDeficit <- copy(urban_in_dt)
urban_in_dt_WaterDeficit <- urban_in_dt_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                     .SDcols = cols_to_center]
urban_in_dt_WaterDeficit$expansionrate <- factor(urban_in_dt_WaterDeficit$expansionrate, levels = c("Slow",   "Moderate",    "Rapid",    "Intense"))

# linear regression
# Table S6 -- model 1
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         # WaterDeficit_2 +
                         Seasonalityindex +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Latitude +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         WaterDeficit:expansionrate +
                         # WaterDeficit_2:expansionrate +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea
                       |  year + CountryID, data = urban_in_dt_WaterDeficit)  

summary(model_urban_in)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_rate.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  


# nonlinear regression
# Table S6 -- model 5
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         WaterDeficit_2 +
                         Seasonalityindex +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Latitude +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         WaterDeficit:expansionrate +
                         WaterDeficit_2:expansionrate +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         ImperviousSurface + HumanSettlement + CityArea
                       |  year + CountryID, data = urban_in_dt_WaterDeficit)  

summary(model_urban_in)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_rate2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

######################################### 绘制效应图与边际效应图 #######################################
beta1 <- coef(model_urban_in)["WaterDeficit"]               
beta2 <- coef(model_urban_in)["WaterDeficit_2"]             
rates <- c("Moderate",    "Rapid",    "Intense") 
beta1_rate <- sapply(rates, function(cont) {
  coef(model_urban_in)[paste0("WaterDeficit:expansionrate", cont)]
})
beta2_rate <- sapply(rates, function(cont) {
  coef(model_urban_in)[paste0("WaterDeficit_2:expansionrate", cont)]
})
beta1_total <- beta1 + beta1_rate
beta2_total <- beta2 + beta2_rate
names(beta1_total) <- rates
names(beta2_total) <- rates
P_min <- quantile(urban_in_dt_WaterDeficit$WaterDeficit, 0.001, na.rm = TRUE)
P_max <- quantile(urban_in_dt_WaterDeficit$WaterDeficit, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)  

# 计算效应值
effect_base <- beta1 * P_values + beta2 * (P_values^2)
effect_rate <- lapply(rates, function(cont) {
  beta1_total[[cont]] * P_values + beta2_total[[cont]] * (P_values^2)
})
names(effect_rate) <- rates
vcov_matrix <- vcovHC(model_urban_in, type = "HC1", cluster = "group",
                      group = urban_in_dt_WaterDeficit$CountryID)
cov_vars_base <- c("WaterDeficit", "WaterDeficit_2")
vcov_base <- vcov_matrix[cov_vars_base, cov_vars_base]
d_beta1 <- P_values  
d_beta2 <- P_values^2  
var_effect_base <- (d_beta1^2) * vcov_base["WaterDeficit", "WaterDeficit"] +
  (d_beta2^2) * vcov_base["WaterDeficit_2", "WaterDeficit_2"] +
  2 * d_beta1 * d_beta2 * vcov_base["WaterDeficit", "WaterDeficit_2"]
se_effect_base <- sqrt(var_effect_base)
ci_effect_base_lower <- effect_base - 1.96 * se_effect_base
ci_effect_base_upper <- effect_base + 1.96 * se_effect_base
var_effect_rate <- list()
se_effect_rate <- list()
ci_effect_rate_lower <- list()
ci_effect_rate_upper <- list()
for (cont in rates) {
  cov_vars <- c("WaterDeficit", "WaterDeficit_2", 
                paste0("WaterDeficit:expansionrate", cont), 
                paste0("WaterDeficit_2:expansionrate", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  var_effect <- (d_beta1^2) * vcov_cont["WaterDeficit", "WaterDeficit"] +
    (d_beta2^2) * vcov_cont["WaterDeficit_2", "WaterDeficit_2"] +
    (d_beta1^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_beta2^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_beta1^2 * vcov_cont["WaterDeficit", cov_vars[3]] +  
    2 * d_beta2^2 * vcov_cont["WaterDeficit_2", cov_vars[4]] +  
    2 * d_beta1 * d_beta2 * vcov_cont["WaterDeficit", "WaterDeficit_2"] +  
    2 * d_beta1 * d_beta2 * vcov_cont["WaterDeficit", cov_vars[4]] +  
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], "WaterDeficit_2"] + 
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
var_marginal_base <- (d_marg_beta1^2) * vcov_base["WaterDeficit", "WaterDeficit"] +
  (d_marg_beta2^2) * vcov_base["WaterDeficit_2", "WaterDeficit_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_base["WaterDeficit", "WaterDeficit_2"]
se_marginal_base <- sqrt(var_marginal_base)
ci_marginal_base_lower <- marginal_base - 1.96 * se_marginal_base
ci_marginal_base_upper <- marginal_base + 1.96 * se_marginal_base
var_marginal_rate <- list()
se_marginal_rate <- list()
ci_marginal_rate_lower <- list()
ci_marginal_rate_upper <- list()
for (cont in rates) {
  cov_vars <- c("WaterDeficit", "WaterDeficit_2",
                paste0("WaterDeficit:expansionrate", cont),
                paste0("WaterDeficit_2:expansionrate", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  d_marg_beta1_u <- 1  
  d_marg_beta2_u <- 2 * P_values  
  var_marginal <- (d_marg_beta1_u^2) * vcov_cont["WaterDeficit", "WaterDeficit"] +
    (d_marg_beta2_u^2) * vcov_cont["WaterDeficit_2", "WaterDeficit_2"] +
    (d_marg_beta1_u^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_marg_beta2_u^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_marg_beta1_u^2 * vcov_cont["WaterDeficit", cov_vars[3]] +  
    2 * d_marg_beta2_u^2 * vcov_cont["WaterDeficit_2", cov_vars[4]] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["WaterDeficit", "WaterDeficit_2"] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["WaterDeficit", cov_vars[4]] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], "WaterDeficit_2"] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], cov_vars[4]]  
  var_marginal_rate[[cont]] <- var_marginal
  se_marginal_rate[[cont]] <- sqrt(var_marginal)
  ci_marginal_rate_lower[[cont]] <- marginal_rate[[cont]] - 1.96 * se_marginal_rate[[cont]]
  ci_marginal_rate_upper[[cont]] <- marginal_rate[[cont]] + 1.96 * se_marginal_rate[[cont]]
}

# 效应值结果
effect_df <- data.frame(
  WaterDeficit = P_values,
  rate = "基准组",
  Effect = effect_base,
  SE = se_effect_base,
  CI_Lower = ci_effect_base_lower,
  CI_Upper = ci_effect_base_upper
)

for (cont in rates) {
  temp_df <- data.frame(
    WaterDeficit = P_values,
    rate = cont,
    Effect = effect_rate[[cont]],
    SE = se_effect_rate[[cont]],
    CI_Lower = ci_effect_rate_lower[[cont]],
    CI_Upper = ci_effect_rate_upper[[cont]]
  )
  effect_df <- rbind(effect_df, temp_df)
}

# 边际效应结果
marginal_df <- data.frame(
  WaterDeficit = P_values,
  rate = "基准组",
  Marginal_Effect = marginal_base,
  SE = se_marginal_base,
  CI_Lower = ci_marginal_base_lower,
  CI_Upper = ci_marginal_base_upper
)

for (cont in rates) {
  temp_df <- data.frame(
    WaterDeficit = P_values,
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




# 计算各发展速率的实际最值
rate_limits <- tibble(
  rate = c("Low",   "Middle",    "High",    "Very high"),
  min_WD = c(
    min(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$expansionrate == "Slow", ]$WaterDeficit, na.rm = TRUE),
    min(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$expansionrate == "Moderate", ]$WaterDeficit, na.rm = TRUE),
    min(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$expansionrate == "Rapid", ]$WaterDeficit, na.rm = TRUE),
    min(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$expansionrate == "Intense", ]$WaterDeficit, na.rm = TRUE)
  ),
  max_WD = c(
    max(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$expansionrate == "Slow", ]$WaterDeficit, na.rm = TRUE),
    max(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$expansionrate == "Moderate", ]$WaterDeficit, na.rm = TRUE),
    max(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$expansionrate == "Rapid", ]$WaterDeficit, na.rm = TRUE),
    max(urban_in_dt_WaterDeficit[urban_in_dt_WaterDeficit$expansionrate == "Intense", ]$WaterDeficit, na.rm = TRUE)
  )
)

print(rate_limits)

# 按各发展速率最值截取
effect_df <- effect_df %>%
  left_join(rate_limits, by = "rate") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)
marginal_df <- marginal_df %>%
  left_join(rate_limits, by = "rate") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)



# 绘制效应值图 
# Figure 2i
rate_colors <- c(
  "Low" = "#27ae60",      
  "Middle" = "#3498db",  
  "High" = "#e67e22",     
  "Very high" = "#c0392b"   
)
effect_df$rate <- factor(effect_df$rate, levels = c( "Low",   "Middle",    "High",    "Very high"))
effect_plot <- ggplot(effect_df, aes(x = WaterDeficit, y = Effect, color = rate)) +
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
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.08, 0.02),
                     breaks = seq(-0.08, 0.02, by = 0.04),
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.text.y = element_text(color = "black", size = 13),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.25, 0.4),
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),
    legend.background = element_blank()
  ) +
  scale_color_manual(values = rate_colors) +
  scale_fill_manual(values = rate_colors)

print(effect_plot)

# 直方图
rate_colors_hist <- c(
  "Slow" = "#27ae60",      
  "Moderate" = "#3498db",  
  "Rapid" = "#e67e22",     
  "Intense" = "#c0392b"    
)

hist_plot <- ggplot(urban_in_dt_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill=expansionrate), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600), name = "") +
  scale_y_continuous(limits = c(0, 0.025), breaks = seq(0, 0.025, by =0.01)) +
  labs(y = "") +
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
effect_rate_waterdeficit <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_rate_waterdeficit)


# 绘制边际效应图 
# Figure S4i
marginal_df$rate <- factor(marginal_df$rate, levels = c( "Low",   "Middle",    "High",    "Very high"))
marginal_plot <- ggplot(marginal_df, aes(x = WaterDeficit, y = Marginal_Effect, color = rate)) +
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
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.00009, 0.00003),
                     breaks = seq(-0.00009, 0.00003, by = 0.00003),
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

marginal_rate_waterdeficit <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_rate_waterdeficit)


##################################################### expansion rate ####################################################
################################################### Rainfall Seasonalityindex ######################################
cols_to_center <- c('WaterDeficit', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
urban_in_dt_Seasonalityindex <- copy(urban_in_dt)
urban_in_dt_Seasonalityindex <- urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$Precipitation < 1000, ]
urban_in_dt_Seasonalityindex <- urban_in_dt_Seasonalityindex[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                             .SDcols = cols_to_center]
urban_in_dt_Seasonalityindex$expansionrate <- factor(urban_in_dt_Seasonalityindex$expansionrate, levels = c( "Slow",   "Moderate",    "Rapid",    "Intense"))

# linear regression
# Table S6 --  model 2
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         Seasonalityindex +
                         # Seasonalityindex_2 +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Latitude +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         Seasonalityindex:expansionrate +
                         # Seasonalityindex_2:expansionrate +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         
                         ImperviousSurface + HumanSettlement + CityArea
                       | year + CountryID, data = urban_in_dt_Seasonalityindex)

summary(model_urban_in)


######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_rate.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  


# nonlinear regression
# Table S6 --  model 6
model_urban_in <- felm(FVC ~ 1 + WaterDeficit +
                         Seasonalityindex +
                         Seasonalityindex_2 +
                         AverageTemperature +
                         TemperatureRange +
                         WindSpeed +
                         Elevation +
                         Latitude +
                         Precipitation +
                         SoilMoisture +
                         SoilPH + 
                         Longitude + 
                         
                         Seasonalityindex:expansionrate +
                         Seasonalityindex_2:expansionrate +
                         
                         GDP_per_capita_PPP + HDI + Population + 
                         
                         ImperviousSurface + HumanSettlement + CityArea
                       | year + CountryID, data = urban_in_dt_Seasonalityindex)

summary(model_urban_in)


######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_rate2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  


######################################### 绘制效应图与边际效应图 #######################################
beta1 <- coef(model_urban_in)["Seasonalityindex"]              
beta2 <- coef(model_urban_in)["Seasonalityindex_2"]            
rates <- c(  "Moderate",    "Rapid",    "Intense")
beta1_rate <- sapply(rates, function(cont) {
  coef(model_urban_in)[paste0("Seasonalityindex:expansionrate", cont)]
})
beta2_rate <- sapply(rates, function(cont) {
  coef(model_urban_in)[paste0("Seasonalityindex_2:expansionrate", cont)]
})
beta1_total <- beta1 + beta1_rate
beta2_total <- beta2 + beta2_rate
names(beta1_total) <- rates
names(beta2_total) <- rates
P_min <- quantile(urban_in_dt_Seasonalityindex$Seasonalityindex, 0.001, na.rm = TRUE)
P_max <- quantile(urban_in_dt_Seasonalityindex$Seasonalityindex, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)  

# 计算效应值
effect_base <- beta1 * P_values + beta2 * (P_values^2)
effect_rate <- lapply(rates, function(cont) {
  beta1_total[[cont]] * P_values + beta2_total[[cont]] * (P_values^2)
})
names(effect_rate) <- rates
vcov_matrix <- vcovHC(model_urban_in, type = "HC1", cluster = "group", 
                      group = urban_in_dt_Seasonalityindex$CountryID)
cov_vars_base <- c("Seasonalityindex", "Seasonalityindex_2")
vcov_base <- vcov_matrix[cov_vars_base, cov_vars_base]
d_beta1 <- P_values 
d_beta2 <- P_values^2  
var_effect_base <- (d_beta1^2) * vcov_base["Seasonalityindex", "Seasonalityindex"] +
  (d_beta2^2) * vcov_base["Seasonalityindex_2", "Seasonalityindex_2"] +
  2 * d_beta1 * d_beta2 * vcov_base["Seasonalityindex", "Seasonalityindex_2"]
se_effect_base <- sqrt(var_effect_base)
ci_effect_base_lower <- effect_base - 1.96 * se_effect_base
ci_effect_base_upper <- effect_base + 1.96 * se_effect_base
var_effect_rate <- list()
se_effect_rate <- list()
ci_effect_rate_lower <- list()
ci_effect_rate_upper <- list()
for (cont in rates) {
  cov_vars <- c("Seasonalityindex", "Seasonalityindex_2", 
                paste0("Seasonalityindex:expansionrate", cont), 
                paste0("Seasonalityindex_2:expansionrate", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  var_effect <- (d_beta1^2) * vcov_cont["Seasonalityindex", "Seasonalityindex"] +
    (d_beta2^2) * vcov_cont["Seasonalityindex_2", "Seasonalityindex_2"] +
    (d_beta1^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_beta2^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_beta1^2 * vcov_cont["Seasonalityindex", cov_vars[3]] +  
    2 * d_beta2^2 * vcov_cont["Seasonalityindex_2", cov_vars[4]] +  
    2 * d_beta1 * d_beta2 * vcov_cont["Seasonalityindex", "Seasonalityindex_2"] + 
    2 * d_beta1 * d_beta2 * vcov_cont["Seasonalityindex", cov_vars[4]] +  
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], "Seasonalityindex_2"] +  
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
var_marginal_base <- (d_marg_beta1^2) * vcov_base["Seasonalityindex", "Seasonalityindex"] +
  (d_marg_beta2^2) * vcov_base["Seasonalityindex_2", "Seasonalityindex_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_base["Seasonalityindex", "Seasonalityindex_2"]
se_marginal_base <- sqrt(var_marginal_base)
ci_marginal_base_lower <- marginal_base - 1.96 * se_marginal_base
ci_marginal_base_upper <- marginal_base + 1.96 * se_marginal_base
var_marginal_rate <- list()
se_marginal_rate <- list()
ci_marginal_rate_lower <- list()
ci_marginal_rate_upper <- list()
for (cont in rates) {
  cov_vars <- c("Seasonalityindex", "Seasonalityindex_2",
                paste0("Seasonalityindex:expansionrate", cont),
                paste0("Seasonalityindex_2:expansionrate", cont))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  d_marg_beta1_u <- 1  
  d_marg_beta2_u <- 2 * P_values  
  var_marginal <- (d_marg_beta1_u^2) * vcov_cont["Seasonalityindex", "Seasonalityindex"] +
    (d_marg_beta2_u^2) * vcov_cont["Seasonalityindex_2", "Seasonalityindex_2"] +
    (d_marg_beta1_u^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_marg_beta2_u^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_marg_beta1_u^2 * vcov_cont["Seasonalityindex", cov_vars[3]] +  
    2 * d_marg_beta2_u^2 * vcov_cont["Seasonalityindex_2", cov_vars[4]] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["Seasonalityindex", "Seasonalityindex_2"] + 
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont["Seasonalityindex", cov_vars[4]] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], "Seasonalityindex_2"] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cont[cov_vars[3], cov_vars[4]]  
  var_marginal_rate[[cont]] <- var_marginal
  se_marginal_rate[[cont]] <- sqrt(var_marginal)
  ci_marginal_rate_lower[[cont]] <- marginal_rate[[cont]] - 1.96 * se_marginal_rate[[cont]]
  ci_marginal_rate_upper[[cont]] <- marginal_rate[[cont]] + 1.96 * se_marginal_rate[[cont]]
}

# 效应值结果
effect_df <- data.frame(
  Seasonalityindex = P_values,
  rate = "基准组",
  Effect = effect_base,
  SE = se_effect_base,
  CI_Lower = ci_effect_base_lower,
  CI_Upper = ci_effect_base_upper
)

for (cont in rates) {
  temp_df <- data.frame(
    Seasonalityindex = P_values,
    rate = cont,
    Effect = effect_rate[[cont]],
    SE = se_effect_rate[[cont]],
    CI_Lower = ci_effect_rate_lower[[cont]],
    CI_Upper = ci_effect_rate_upper[[cont]]
  )
  effect_df <- rbind(effect_df, temp_df)
}

# 边际效应结果
marginal_df <- data.frame(
  Seasonalityindex = P_values,
  rate = "基准组",
  Marginal_Effect = marginal_base,
  SE = se_marginal_base,
  CI_Lower = ci_marginal_base_lower,
  CI_Upper = ci_marginal_base_upper
)

for (cont in rates) {
  temp_df <- data.frame(
    Seasonalityindex = P_values,
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

# 计算各发展速率的实际最值 
rate_limits <- tibble(
  rate = c("Low",   "Middle",    "High",    "Very high"),
  min_WD = c(
    min(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$expansionrate == "Slow", ]$Seasonalityindex, na.rm = TRUE),
    min(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$expansionrate == "Moderate", ]$Seasonalityindex, na.rm = TRUE),
    min(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$expansionrate == "Rapid", ]$Seasonalityindex, na.rm = TRUE),
    min(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$expansionrate == "Intense", ]$Seasonalityindex, na.rm = TRUE)
  ),
  max_WD = c(
    max(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$expansionrate == "Slow", ]$Seasonalityindex, na.rm = TRUE),
    max(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$expansionrate == "Moderate", ]$Seasonalityindex, na.rm = TRUE),
    max(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$expansionrate == "Rapid", ]$Seasonalityindex, na.rm = TRUE),
    max(urban_in_dt_Seasonalityindex[urban_in_dt_Seasonalityindex$expansionrate == "Intense", ]$Seasonalityindex, na.rm = TRUE)
  )
)

print(rate_limits)

# 按各发展速率最值截取
effect_df <- effect_df %>%
  left_join(rate_limits, by = "rate") %>%
  filter(Seasonalityindex >= min_WD & Seasonalityindex <= max_WD) %>%
  select(-min_WD, -max_WD)
marginal_df <- marginal_df %>%
  left_join(rate_limits, by = "rate") %>%
  filter(Seasonalityindex >= min_WD & Seasonalityindex <= max_WD) %>%
  select(-min_WD, -max_WD)

# 效应值图
# Figure 2j
effect_plot <- ggplot(effect_df, aes(x = Seasonalityindex, y = Effect, color = rate)) +
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
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02), 
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.text.y = element_text(color = "black", size = 13),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",        
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),
  ) +
  scale_color_manual(values = rate_colors) +
  scale_fill_manual(values = rate_colors)

print(effect_plot)

# 直方图
hist_plot <- ggplot(urban_in_dt_Seasonalityindex, aes(x = Seasonalityindex)) +
  geom_histogram(aes(y = ..density.., fill=expansionrate), binwidth = 0.0005, alpha = 0.5) +
  scale_x_continuous(limits = c(0,0.1), breaks = seq(0,0.1, by =0.02), name = "") +
  scale_y_continuous(limits = c(0,550), breaks = seq(0,550, by =200)) +
  labs(y = "") +
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
effect_rate_seasonalityindex <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_rate_seasonalityindex)


# 绘制边际效应图
# Figure S4j
marginal_plot <- ggplot(marginal_df, aes(x = Seasonalityindex, y = Marginal_Effect, color = rate)) +
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
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-2, 0.6),
                     breaks = seq(-2, 0.6, by = 1),
                     labels = scales::number_format(accuracy = 1)) +
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

marginal_rate_seasonalityindex <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_rate_seasonalityindex)



################################################### 城市内异质性分析 ################################################## 
##################################################### country ################################################## 
################################################### water deficit ######################################
cols_to_center <- c('Seasonalityindex', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
urban_in_dt_urbanization_WaterDeficit <- copy(urban_in_dt)
urban_in_dt_urbanization_WaterDeficit <- urban_in_dt_urbanization_WaterDeficit[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                                               .SDcols = cols_to_center]
urban_in_dt_urbanization_WaterDeficit$country <- factor(urban_in_dt_urbanization_WaterDeficit$country, 
                                                        levels = c("United States of America",
                                                                   "Europe",
                                                                   "Canada",
                                                                   "China",
                                                                   "India",
                                                                   "Brazil", "other"))
# linear regression
# Table 7 -- model 1
model_urban_in_country <- felm(FVC ~ 1 + WaterDeficit +
                                 # WaterDeficit_2 +
                                 Seasonalityindex +
                                 AverageTemperature +
                                 TemperatureRange +
                                 WindSpeed +
                                 Elevation +
                                 Latitude +
                                 Precipitation +
                                 SoilMoisture +
                                 SoilPH + 
                                 Longitude + 
                                 
                                 WaterDeficit:country +
                                 # WaterDeficit_2:country +
                                 
                                 GDP_per_capita_PPP + HDI + Population + 
                                 
                                 ImperviousSurface + HumanSettlement + CityArea
                               | year, data = urban_in_dt_urbanization_WaterDeficit)

summary(model_urban_in_country)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_country.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in_country))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  


# nonlinear regression
# Table 7 -- model 5
model_urban_in_country <- felm(FVC ~ 1 + WaterDeficit +
                                 WaterDeficit_2 +
                                 Seasonalityindex +
                                 AverageTemperature +
                                 TemperatureRange +
                                 WindSpeed +
                                 Elevation +
                                 Latitude +
                                 Precipitation +
                                 SoilMoisture +
                                 SoilPH + 
                                 Longitude + 
                                 
                                 WaterDeficit:country +
                                 WaterDeficit_2:country +
                                 
                                 GDP_per_capita_PPP + HDI + Population + 
                                 
                                 ImperviousSurface + HumanSettlement + CityArea
                               | year, data = urban_in_dt_urbanization_WaterDeficit)

summary(model_urban_in_country)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_country2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in_country))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

######################################### 绘制效应图与边际效应图 #######################################
beta1 <- coef(model_urban_in_country)["WaterDeficit"]               
beta2 <- coef(model_urban_in_country)["WaterDeficit_2"]          
countries <- c("Europe", "Canada",  "China", "India", "Brazil", "other") 
beta1_country <- sapply(countries, function(cty) {
  coef(model_urban_in_country)[paste0("WaterDeficit:country", cty)]
})
beta2_country <- sapply(countries, function(cty) {
  coef(model_urban_in_country)[paste0("WaterDeficit_2:country", cty)]
})
beta1_total <- beta1 + beta1_country
beta2_total <- beta2 + beta2_country
names(beta1_total) <- countries
names(beta2_total) <- countries
P_min <- quantile(urban_in_dt_urbanization_WaterDeficit$WaterDeficit, 0.001, na.rm = TRUE)
P_max <- quantile(urban_in_dt_urbanization_WaterDeficit$WaterDeficit, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)  

# 计算效应值
effect_base <- beta1 * P_values + beta2 * (P_values^2)
effect_country <- lapply(countries, function(cty) {
  beta1_total[[cty]] * P_values + beta2_total[[cty]] * (P_values^2)
})
names(effect_country) <- countries
vcov_matrix <- vcovHC(model_urban_in_country, type = "HC1", cluster = "group",
                      group = urban_in_dt_urbanization_WaterDeficit$CountryID)
cov_vars_base <- c("WaterDeficit", "WaterDeficit_2")
vcov_base <- vcov_matrix[cov_vars_base, cov_vars_base]
d_beta1 <- P_values  
d_beta2 <- P_values^2  
var_effect_base <- (d_beta1^2) * vcov_base["WaterDeficit", "WaterDeficit"] +
  (d_beta2^2) * vcov_base["WaterDeficit_2", "WaterDeficit_2"] +
  2 * d_beta1 * d_beta2 * vcov_base["WaterDeficit", "WaterDeficit_2"]
se_effect_base <- sqrt(var_effect_base)
ci_effect_base_lower <- effect_base - 1.96 * se_effect_base
ci_effect_base_upper <- effect_base + 1.96 * se_effect_base
var_effect_country <- list()
se_effect_country <- list()
ci_effect_country_lower <- list()
ci_effect_country_upper <- list()

for (cty in countries) {
  cov_vars <- c("WaterDeficit", "WaterDeficit_2", 
                paste0("WaterDeficit:country", cty), 
                paste0("WaterDeficit_2:country", cty))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  var_effect <- (d_beta1^2) * vcov_cont["WaterDeficit", "WaterDeficit"] +
    (d_beta2^2) * vcov_cont["WaterDeficit_2", "WaterDeficit_2"] +
    (d_beta1^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_beta2^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_beta1^2 * vcov_cont["WaterDeficit", cov_vars[3]] +  
    2 * d_beta2^2 * vcov_cont["WaterDeficit_2", cov_vars[4]] +  
    2 * d_beta1 * d_beta2 * vcov_cont["WaterDeficit", "WaterDeficit_2"] +  
    2 * d_beta1 * d_beta2 * vcov_cont["WaterDeficit", cov_vars[4]] +  
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], "WaterDeficit_2"] +  
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
var_marginal_base <- (d_marg_beta1^2) * vcov_base["WaterDeficit", "WaterDeficit"] +
  (d_marg_beta2^2) * vcov_base["WaterDeficit_2", "WaterDeficit_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_base["WaterDeficit", "WaterDeficit_2"]
se_marginal_base <- sqrt(var_marginal_base)
ci_marginal_base_lower <- marginal_base - 1.96 * se_marginal_base
ci_marginal_base_upper <- marginal_base + 1.96 * se_marginal_base
var_marginal_country <- list()
se_marginal_country <- list()
ci_marginal_country_lower <- list()
ci_marginal_country_upper <- list()
for (cty in countries) {
  cov_vars <- c("WaterDeficit", "WaterDeficit_2",
                paste0("WaterDeficit:country", cty),
                paste0("WaterDeficit_2:country", cty))
  vcov_cty <- vcov_matrix[cov_vars, cov_vars]
  d_marg_beta1_u <- 1  
  d_marg_beta2_u <- 2 * P_values  
  var_marginal <- (d_marg_beta1_u^2) * vcov_cty["WaterDeficit", "WaterDeficit"] +
    (d_marg_beta2_u^2) * vcov_cty["WaterDeficit_2", "WaterDeficit_2"] +
    (d_marg_beta1_u^2) * vcov_cty[cov_vars[3], cov_vars[3]] +
    (d_marg_beta2_u^2) * vcov_cty[cov_vars[4], cov_vars[4]] +
    2 * d_marg_beta1_u^2 * vcov_cty["WaterDeficit", cov_vars[3]] +  
    2 * d_marg_beta2_u^2 * vcov_cty["WaterDeficit_2", cov_vars[4]] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cty["WaterDeficit", "WaterDeficit_2"] + 
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cty["WaterDeficit", cov_vars[4]] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cty[cov_vars[3], "WaterDeficit_2"] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cty[cov_vars[3], cov_vars[4]]  
  var_marginal_country[[cty]] <- var_marginal
  se_marginal_country[[cty]] <- sqrt(var_marginal)
  ci_marginal_country_lower[[cty]] <- marginal_country[[cty]] - 1.96 * se_marginal_country[[cty]]
  ci_marginal_country_upper[[cty]] <- marginal_country[[cty]] + 1.96 * se_marginal_country[[cty]]
}

# 效应值结果
effect_df <- data.frame(
  WaterDeficit = P_values,
  Country = "基准组",
  Effect = effect_base,
  SE = se_effect_base,
  CI_Lower = ci_effect_base_lower,
  CI_Upper = ci_effect_base_upper
)

for (cty in countries) {
  temp_df <- data.frame(
    WaterDeficit = P_values,
    Country = cty,
    Effect = effect_country[[cty]],
    SE = se_effect_country[[cty]],
    CI_Lower = ci_effect_country_lower[[cty]],
    CI_Upper = ci_effect_country_upper[[cty]]
  )
  effect_df <- rbind(effect_df, temp_df)
}

# 边际效应结果
marginal_df <- data.frame(
  WaterDeficit = P_values,
  Country = "基准组",
  Marginal_Effect = marginal_base,
  SE = se_marginal_base,
  CI_Lower = ci_marginal_base_lower,
  CI_Upper = ci_marginal_base_upper
)

for (cty in countries) {
  temp_df <- data.frame(
    WaterDeficit = P_values,
    Country = cty,
    Marginal_Effect = marginal_country[[cty]],
    SE = se_marginal_country[[cty]],
    CI_Lower = ci_marginal_country_lower[[cty]],
    CI_Upper = ci_marginal_country_upper[[cty]]
  )
  marginal_df <- rbind(marginal_df, temp_df)
}
effect_df$Country <- factor(effect_df$Country, 
                            levels = c("基准组", countries),
                            labels = c("United States", "Europe (excluding Russia)", "Canada",  "China", "India", "Brazil", "other")) 
marginal_df$Country <- factor(marginal_df$Country, 
                              levels = c("基准组", countries),
                              labels = c("United States", "Europe (excluding Russia)", "Canada",  "China", "India", "Brazil", "other")) 


# 计算各国家的实际最值
country_limits <- tibble(
  Country = c("United States", "Europe (excluding Russia)", "Canada",  "China", "India", "Brazil", "other"),
  min_WD = c(
    min(urban_in_dt_urbanization_WaterDeficit[country == "United States of America"]$WaterDeficit, na.rm = TRUE),
    min(urban_in_dt_urbanization_WaterDeficit[country == "Europe" ]$WaterDeficit, na.rm = TRUE),
    min(urban_in_dt_urbanization_WaterDeficit[country == "Canada" ]$WaterDeficit, na.rm = TRUE),
    min(urban_in_dt_urbanization_WaterDeficit[country == "China" ]$WaterDeficit, na.rm = TRUE),
    min(urban_in_dt_urbanization_WaterDeficit[country == "India" ]$WaterDeficit, na.rm = TRUE),
    min(urban_in_dt_urbanization_WaterDeficit[country == "Brazil" ]$WaterDeficit, na.rm = TRUE),
    min(urban_in_dt_urbanization_WaterDeficit[country == "other" ]$WaterDeficit, na.rm = TRUE)
  ),
  max_WD = c(
    max(urban_in_dt_urbanization_WaterDeficit[country == "United States of America" ]$WaterDeficit, na.rm = TRUE),
    max(urban_in_dt_urbanization_WaterDeficit[country == "Europe" ]$WaterDeficit, na.rm = TRUE),
    max(urban_in_dt_urbanization_WaterDeficit[country == "Canada"]$WaterDeficit, na.rm = TRUE),
    max(urban_in_dt_urbanization_WaterDeficit[country == "China"]$WaterDeficit, na.rm = TRUE),
    max(urban_in_dt_urbanization_WaterDeficit[country == "India"]$WaterDeficit, na.rm = TRUE),
    max(urban_in_dt_urbanization_WaterDeficit[country == "Brazil"]$WaterDeficit, na.rm = TRUE),
    max(urban_in_dt_urbanization_WaterDeficit[country == "other"]$WaterDeficit, na.rm = TRUE)
  )
)

print(country_limits)

# 按各国家最值截取
effect_df <- effect_df %>%
  left_join(country_limits, by = "Country") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)
marginal_df <- marginal_df %>%
  left_join(country_limits, by = "Country") %>%
  filter(WaterDeficit >= min_WD & WaterDeficit <= max_WD) %>%
  select(-min_WD, -max_WD)



# 绘制效应值图 
# Figure 2k
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
effect_plot <- ggplot(effect_df, aes(x = WaterDeficit, y = Effect, color = Country)) +
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
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.15, 0.55),
                     breaks = seq(-0.15, 0.55, by = 0.15),
                     labels = scales::number_format(accuracy = 0.01)) +
  theme(
    axis.text.y = element_text(color = "black", size = 13),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.5, 0.75),
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),
    legend.background = element_blank(),
  ) +
  scale_color_manual(values = country_colors) +
  scale_fill_manual(values = country_colors)

print(effect_plot)

# 直方图
hist_plot <- ggplot(urban_in_dt_urbanization_WaterDeficit, aes(x = WaterDeficit)) +
  geom_histogram(aes(y = ..density.., fill=country), binwidth = 5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600), name = "Water deficit (mm)") +
  scale_y_continuous(limits = c(0, 0.025), breaks = seq(0, 0.025, by =0.01)) +
  labs(y = "") +
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
effect_country_waterdeficit <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_country_waterdeficit)


# 绘制边际效应图
# Figure S4k
marginal_df$Country <- factor(marginal_df$Country, levels = c("United States", "Europe (excluding Russia)", "Canada", "China", "India", "Brazil" ,"other"))

marginal_plot <- ggplot(marginal_df, aes(x = WaterDeficit, y = Marginal_Effect, color = Country)) +
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
  scale_x_continuous(limits = c(0, 2430), breaks = seq(0, 2430, by =600), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-0.0005, 0.001),
                     breaks = seq(-0.0005, 0.001, by = 0.0005),
                     labels = scales::number_format(accuracy = 0.0001)) +
  theme(
    axis.text.y = element_text(color = "black", size = 12),  
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),  
    axis.title.y = element_text(size = 12),  
    legend.position = c(0.6,0.8), 
    legend.direction = "vertical", 
    legend.text = element_text(size = 12), legend.key.size = unit(0.3, 'cm'),  
    legend.background = element_blank()
  ) +
  scale_color_manual(values = country_colors) +
  scale_fill_manual(values = country_colors)

print(marginal_plot)

marginal_country_waterdeficit <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_country_waterdeficit)

##################################################### country ####################################################
################################################### Rainfall Seasonalityindex ######################################
cols_to_center <- c('WaterDeficit', 'AverageTemperature', 'TemperatureRange', 'WindSpeed', 'Elevation', 'Latitude', 
                    'Precipitation', 'SoilMoisture', 'SoilPH', 'Longitude', 'GDP_per_capita_PPP', 'HDI', 'Population', 
                    'ImperviousSurface', 'HumanSettlement', 'CityArea')  
urban_in_dt_urbanization_Seasonalityindex <- copy(urban_in_dt)
urban_in_dt_urbanization_Seasonalityindex <- urban_in_dt_urbanization_Seasonalityindex[urban_in_dt_urbanization_Seasonalityindex$Precipitation < 1000, ]
urban_in_dt_urbanization_Seasonalityindex <- urban_in_dt_urbanization_Seasonalityindex[, (cols_to_center) := lapply(.SD, function(x) x - mean(x, na.rm = TRUE)), 
                                                                                       .SDcols = cols_to_center]
urban_in_dt_urbanization_Seasonalityindex$country <- factor(urban_in_dt_urbanization_Seasonalityindex$country, 
                                                            levels = c( "United States of America",
                                                                        "Europe",
                                                                        "Canada",
                                                                        "China",
                                                                        "India",
                                                                        "Brazil", "other"))
# linear regression
# Table S7 -- model 2
model_urban_in_country <- felm(FVC ~ 1 + WaterDeficit +
                                 Seasonalityindex +
                                 # Seasonalityindex_2 +
                                 AverageTemperature +
                                 TemperatureRange +
                                 WindSpeed +
                                 Elevation +
                                 Latitude +
                                 Precipitation +
                                 SoilMoisture +
                                 SoilPH + 
                                 Longitude +
                                 
                                 Seasonalityindex:country +
                                 # Seasonalityindex_2:country +
                                 
                                 GDP_per_capita_PPP + HDI + Population + 
                                 
                                 ImperviousSurface + HumanSettlement + CityArea
                               | year, data = urban_in_dt_urbanization_Seasonalityindex)

summary(model_urban_in_country)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_country.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in_country))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  


# nonlinear regression
# Table S7 -- model 6
model_urban_in_country <- felm(FVC ~ 1 + WaterDeficit +
                                 Seasonalityindex +
                                 Seasonalityindex_2 +
                                 AverageTemperature +
                                 TemperatureRange +
                                 WindSpeed +
                                 Elevation +
                                 Latitude +
                                 Precipitation +
                                 SoilMoisture +
                                 SoilPH + 
                                 Longitude +
                                 
                                 Seasonalityindex:country +
                                 Seasonalityindex_2:country +
                                 
                                 GDP_per_capita_PPP + HDI + Population + 
                                 
                                 ImperviousSurface + HumanSettlement + CityArea
                               | year, 
                               data = urban_in_dt_urbanization_Seasonalityindex)

summary(model_urban_in_country)

######################################### 保存模型系数 #######################################
output_file <- "E:/2025/UrbanGreenSpaceWorks/Rcodes/urban_in_country2.txt"
con <- file(output_file, open = "a")
model_summary <- capture.output(summary(model_urban_in_country))
cat(model_summary, file = output_file, sep = "\n", append = TRUE)
close(con)  

######################################### 绘制效应图与边际效应图 #######################################
beta1 <- coef(model_urban_in_country)["Seasonalityindex"]              
beta2 <- coef(model_urban_in_country)["Seasonalityindex_2"]             
countries <- c("Europe", "Canada",  "China", "India", "Brazil", "other") 
beta1_country <- sapply(countries, function(cty) {
  coef(model_urban_in_country)[paste0("Seasonalityindex:country", cty)]
})
beta2_country <- sapply(countries, function(cty) {
  coef(model_urban_in_country)[paste0("Seasonalityindex_2:country", cty)]
})
beta1_total <- beta1 + beta1_country
beta2_total <- beta2 + beta2_country
names(beta1_total) <- countries
names(beta2_total) <- countries
P_min <- quantile(urban_in_dt_urbanization_Seasonalityindex$Seasonalityindex, 0.001, na.rm = TRUE)
P_max <- quantile(urban_in_dt_urbanization_Seasonalityindex$Seasonalityindex, 0.999, na.rm = TRUE)
P_values <- seq(P_min, P_max, length.out = 200)  

# 计算效应值
effect_base <- beta1 * P_values + beta2 * (P_values^2)
effect_country <- lapply(countries, function(cty) {
  beta1_total[[cty]] * P_values + beta2_total[[cty]] * (P_values^2)
})
names(effect_country) <- countries
vcov_matrix <- vcovHC(model_urban_in_country, type = "HC1", cluster = "group",
                      group = urban_in_dt_urbanization_Seasonalityindex$CountryID)
cov_vars_base <- c("Seasonalityindex", "Seasonalityindex_2")
vcov_base <- vcov_matrix[cov_vars_base, cov_vars_base]
d_beta1 <- P_values  
d_beta2 <- P_values^2  
var_effect_base <- (d_beta1^2) * vcov_base["Seasonalityindex", "Seasonalityindex"] +
  (d_beta2^2) * vcov_base["Seasonalityindex_2", "Seasonalityindex_2"] +
  2 * d_beta1 * d_beta2 * vcov_base["Seasonalityindex", "Seasonalityindex_2"]
se_effect_base <- sqrt(var_effect_base)
ci_effect_base_lower <- effect_base - 1.96 * se_effect_base
ci_effect_base_upper <- effect_base + 1.96 * se_effect_base
var_effect_country <- list()
se_effect_country <- list()
ci_effect_country_lower <- list()
ci_effect_country_upper <- list()

for (cty in countries) {
  cov_vars <- c("Seasonalityindex", "Seasonalityindex_2", 
                paste0("Seasonalityindex:country", cty), 
                paste0("Seasonalityindex_2:country", cty))
  vcov_cont <- vcov_matrix[cov_vars, cov_vars]
  var_effect <- (d_beta1^2) * vcov_cont["Seasonalityindex", "Seasonalityindex"] +
    (d_beta2^2) * vcov_cont["Seasonalityindex_2", "Seasonalityindex_2"] +
    (d_beta1^2) * vcov_cont[cov_vars[3], cov_vars[3]] +
    (d_beta2^2) * vcov_cont[cov_vars[4], cov_vars[4]] +
    2 * d_beta1^2 * vcov_cont["Seasonalityindex", cov_vars[3]] +  
    2 * d_beta2^2 * vcov_cont["Seasonalityindex_2", cov_vars[4]] +  
    2 * d_beta1 * d_beta2 * vcov_cont["Seasonalityindex", "Seasonalityindex_2"] +  
    2 * d_beta1 * d_beta2 * vcov_cont["Seasonalityindex", cov_vars[4]] +  
    2 * d_beta1 * d_beta2 * vcov_cont[cov_vars[3], "Seasonalityindex_2"] +  
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
var_marginal_base <- (d_marg_beta1^2) * vcov_base["Seasonalityindex", "Seasonalityindex"] +
  (d_marg_beta2^2) * vcov_base["Seasonalityindex_2", "Seasonalityindex_2"] +
  2 * d_marg_beta1 * d_marg_beta2 * vcov_base["Seasonalityindex", "Seasonalityindex_2"]
se_marginal_base <- sqrt(var_marginal_base)
ci_marginal_base_lower <- marginal_base - 1.96 * se_marginal_base
ci_marginal_base_upper <- marginal_base + 1.96 * se_marginal_base
var_marginal_country <- list()
se_marginal_country <- list()
ci_marginal_country_lower <- list()
ci_marginal_country_upper <- list()
for (cty in countries) {
  cov_vars <- c("Seasonalityindex", "Seasonalityindex_2",
                paste0("Seasonalityindex:country", cty),
                paste0("Seasonalityindex_2:country", cty))
  vcov_cty <- vcov_matrix[cov_vars, cov_vars]
  d_marg_beta1_u <- 1  
  d_marg_beta2_u <- 2 * P_values  
  var_marginal <- (d_marg_beta1_u^2) * vcov_cty["Seasonalityindex", "Seasonalityindex"] +
    (d_marg_beta2_u^2) * vcov_cty["Seasonalityindex_2", "Seasonalityindex_2"] +
    (d_marg_beta1_u^2) * vcov_cty[cov_vars[3], cov_vars[3]] +
    (d_marg_beta2_u^2) * vcov_cty[cov_vars[4], cov_vars[4]] +
    2 * d_marg_beta1_u^2 * vcov_cty["Seasonalityindex", cov_vars[3]] +  
    2 * d_marg_beta2_u^2 * vcov_cty["Seasonalityindex_2", cov_vars[4]] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cty["Seasonalityindex", "Seasonalityindex_2"] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cty["Seasonalityindex", cov_vars[4]] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cty[cov_vars[3], "Seasonalityindex_2"] +  
    2 * d_marg_beta1_u * d_marg_beta2_u * vcov_cty[cov_vars[3], cov_vars[4]]  
  var_marginal_country[[cty]] <- var_marginal
  se_marginal_country[[cty]] <- sqrt(var_marginal)
  ci_marginal_country_lower[[cty]] <- marginal_country[[cty]] - 1.96 * se_marginal_country[[cty]]
  ci_marginal_country_upper[[cty]] <- marginal_country[[cty]] + 1.96 * se_marginal_country[[cty]]
}


# 边际效应结果
marginal_df <- data.frame(
  Seasonalityindex = P_values,
  Country = "基准组",
  Marginal_Effect = marginal_base,
  SE = se_marginal_base,
  CI_Lower = ci_marginal_base_lower,
  CI_Upper = ci_marginal_base_upper
)

for (cty in countries) {
  temp_df <- data.frame(
    Seasonalityindex = P_values,
    Country = cty,
    Marginal_Effect = marginal_country[[cty]],
    SE = se_marginal_country[[cty]],
    CI_Lower = ci_marginal_country_lower[[cty]],
    CI_Upper = ci_marginal_country_upper[[cty]]
  )
  marginal_df <- rbind(marginal_df, temp_df)
}


# 效应值结果
effect_df <- data.frame(
  Seasonalityindex = P_values,
  Country = "基准组",
  Effect = effect_base,
  SE = se_effect_base,
  CI_Lower = ci_effect_base_lower,
  CI_Upper = ci_effect_base_upper
)

for (cty in countries) {
  temp_df <- data.frame(
    Seasonalityindex = P_values,
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



# 计算各国家的实际最值 
country_limits <- tibble(
  Country = c("United States", "Europe (excluding Russia)", "Canada",  "China", "India", "Brazil", "other"),
  min_WD = c(
    min(urban_in_dt_urbanization_Seasonalityindex[country == "United States of America"]$Seasonalityindex, na.rm = TRUE),
    min(urban_in_dt_urbanization_Seasonalityindex[country == "Europe" ]$Seasonalityindex, na.rm = TRUE),
    min(urban_in_dt_urbanization_Seasonalityindex[country == "Canada" ]$Seasonalityindex, na.rm = TRUE),
    min(urban_in_dt_urbanization_Seasonalityindex[country == "China" ]$Seasonalityindex, na.rm = TRUE),
    min(urban_in_dt_urbanization_Seasonalityindex[country == "India" ]$Seasonalityindex, na.rm = TRUE),
    min(urban_in_dt_urbanization_Seasonalityindex[country == "Brazil" ]$Seasonalityindex, na.rm = TRUE),
    min(urban_in_dt_urbanization_Seasonalityindex[country == "other" ]$Seasonalityindex, na.rm = TRUE)
  ),
  max_WD = c(
    max(urban_in_dt_urbanization_Seasonalityindex[country == "United States of America" ]$Seasonalityindex, na.rm = TRUE),
    max(urban_in_dt_urbanization_Seasonalityindex[country == "Europe" ]$Seasonalityindex, na.rm = TRUE),
    quantile(urban_in_dt_urbanization_Seasonalityindex[country == "Canada"]$Seasonalityindex, probs = 0.99, na.rm = TRUE),
    max(urban_in_dt_urbanization_Seasonalityindex[country == "China"]$Seasonalityindex, na.rm = TRUE),
    max(urban_in_dt_urbanization_Seasonalityindex[country == "India"]$Seasonalityindex, na.rm = TRUE),
    max(urban_in_dt_urbanization_Seasonalityindex[country == "Brazil"]$Seasonalityindex, na.rm = TRUE),
    max(urban_in_dt_urbanization_Seasonalityindex[country == "other"]$Seasonalityindex, na.rm = TRUE)
  )
)

print(country_limits)

# 按各国家最值截取
effect_df <- effect_df %>%
  left_join(country_limits, by = "Country") %>%
  filter(Seasonalityindex >= min_WD & Seasonalityindex <= max_WD) %>%
  select(-min_WD, -max_WD)
marginal_df <- marginal_df %>%
  left_join(country_limits, by = "Country") %>%
  filter(Seasonalityindex >= min_WD & Seasonalityindex <= max_WD) %>%
  select(-min_WD, -max_WD)

# 绘制效应值图
# Figure 2l
effect_plot <- ggplot(effect_df, aes(x = Seasonalityindex, y = Effect, color = Country)) +
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
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02), 
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

# 直方图
hist_plot <- ggplot(urban_in_dt_urbanization_Seasonalityindex, aes(x = Seasonalityindex)) +
  geom_histogram(aes(y = ..density.., fill=country), binwidth = 0.0002, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02), name = "Seasonality index") +
  scale_y_continuous(limits = c(0, 900), breaks = seq(0, 900, by =400)) +
  labs(y = "") +
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
effect_country_seasonalityindex <- ggarrange(effect_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(effect_country_seasonalityindex)


# 绘制边际效应图
# Figure S4l
marginal_plot <- ggplot(marginal_df, aes(x = Seasonalityindex, y = Marginal_Effect, color = Country)) +
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
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by =0.02), 
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(limits = c(-20, 15),
                     breaks = seq(-20, 15, by = 10),
                     labels = scales::number_format(accuracy = 1)) +
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

marginal_country_seasonalityindex <- ggarrange(marginal_plot, hist_plot, ncol = 1, align = "v", heights = c(2, 1))
print(marginal_country_seasonalityindex)


################################################合并图形###############################################################
# Figure 2: The effect of urbanization on the impacts of hydroclimatic stresses on UGS coverage
ii <- cowplot::plot_grid(effect_waterdeficit, nrow = 1, ncol = 1, labels = c("a"),               
                         label_size = 18, scale=1, label_x = 0.02, label_y = 1.00)

jj <- cowplot::plot_grid(effect_seasonalityindex, nrow = 1, ncol = 1, labels = c("b"),
                         label_size = 18, scale=1, label_x = 0.02, label_y = 1.00)

kk <- cowplot::plot_grid(effect_continent_waterdeficit, nrow = 1, ncol = 1, labels = c("c"),
                         label_size = 18, scale=1, label_x = -0.05, label_y = 1.15)

mm <- cowplot::plot_grid(effect_continent_seasonalityindex, nrow = 1, ncol = 1, labels = c("d"),
                         label_size = 18, scale=1, label_x = 0.02, label_y = 1.15)

oo <- cowplot::plot_grid(effect_climate_waterdeficit, nrow = 1, ncol = 1, labels = c("e"),
                         label_size = 18, scale=1, label_x = -0.05, label_y = 1.15)

pp <- cowplot::plot_grid(effect_climate_seasonalityindex, nrow = 1, ncol = 1, labels = c("f"),
                         label_size = 18, scale=1, label_x = 0.02, label_y = 1.15)

oo2 <- cowplot::plot_grid(effect_level_waterdeficit, nrow = 1, ncol = 1, labels = c("g"),
                          label_size = 18, scale=1, label_x = -0.05, label_y = 1.15)

pp2 <- cowplot::plot_grid(effect_level_seasonalityindex, nrow = 1, ncol = 1, labels = c("h"),
                          label_size = 18, scale=1, label_x = 0.02, label_y = 1.15)

oo3 <- cowplot::plot_grid(effect_rate_waterdeficit, nrow = 1, ncol = 1, labels = c("i"),
                          label_size = 18, scale=1, label_x = -0.05, label_y = 1.15)

pp3 <- cowplot::plot_grid(effect_rate_seasonalityindex, nrow = 1, ncol = 1, labels = c("j"),
                          label_size = 18, scale=1, label_x = 0.02, label_y = 1.15)

oo4 <- cowplot::plot_grid(effect_country_waterdeficit, nrow = 1, ncol = 1, labels = c("k"),
                          label_size = 18, scale=1, label_x = -0.05, label_y = 1.15)

pp4 <- cowplot::plot_grid(effect_country_seasonalityindex, nrow = 1, ncol = 1, labels = c("l"),
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
  theme(plot.margin = margin(25, 5, 5, 5))  

print(p_final)

ggsave("E:/2025/UrbanGreenSpaceWorks/Rcodes/Fig2.pdf", 
       plot = p_final, width = 300, height = 360, units = "mm",dpi=1000)



################################################合并图形###############################################################
# Figure S4: The marginal effect of urbanization on the impacts of hydroclimatic stresses on UGS coverage
ii <- cowplot::plot_grid(marginal_waterdeficit, nrow = 1, ncol = 1, labels = c("a"),               
                         label_size = 18, scale=1, label_x = 0.02, label_y = 1.00)

jj <- cowplot::plot_grid(marginal_seasonalityindex, nrow = 1, ncol = 1, labels = c("b"),
                         label_size = 18, scale=1, label_x = 0.02, label_y = 1.00)

kk <- cowplot::plot_grid(marginal_continent_waterdeficit, nrow = 1, ncol = 1, labels = c("c"),
                         label_size = 18, scale=1, label_x = -0.05, label_y = 1.15)

mm <- cowplot::plot_grid(marginal_continent_seasonalityindex, nrow = 1, ncol = 1, labels = c("d"),
                         label_size = 18, scale=1, label_x = 0.02, label_y = 1.15)

oo <- cowplot::plot_grid(marginal_climate_waterdeficit, nrow = 1, ncol = 1, labels = c("e"),
                         label_size = 18, scale=1, label_x = -0.05, label_y = 1.15)

pp <- cowplot::plot_grid(marginal_climate_seasonalityindex, nrow = 1, ncol = 1, labels = c("f"),
                         label_size = 18, scale=1, label_x = 0.02, label_y = 1.15)

oo2 <- cowplot::plot_grid(marginal_level_waterdeficit, nrow = 1, ncol = 1, labels = c("g"),
                          label_size = 18, scale=1, label_x = -0.05, label_y = 1.15)

pp2 <- cowplot::plot_grid(marginal_level_seasonalityindex, nrow = 1, ncol = 1, labels = c("h"),
                          label_size = 18, scale=1, label_x = 0.02, label_y = 1.15)

oo3 <- cowplot::plot_grid(marginal_rate_waterdeficit, nrow = 1, ncol = 1, labels = c("i"),
                          label_size = 18, scale=1, label_x = -0.05, label_y = 1.15)

pp3 <- cowplot::plot_grid(marginal_rate_seasonalityindex, nrow = 1, ncol = 1, labels = c("j"),
                          label_size = 18, scale=1, label_x = 0.02, label_y = 1.15)

oo4 <- cowplot::plot_grid(marginal_country_waterdeficit, nrow = 1, ncol = 1, labels = c("k"),
                          label_size = 18, scale=1, label_x = -0.05, label_y = 1.15)

pp4 <- cowplot::plot_grid(marginal_country_seasonalityindex, nrow = 1, ncol = 1, labels = c("l"),
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
  theme(plot.margin = margin(25, 5, 5, 5))  

print(p_final)

ggsave("E:/2025/UrbanGreenSpaceWorks/Rcodes/FigS4.pdf", 
       plot = p_final, width = 320, height = 380, units = "mm",dpi=1000)





