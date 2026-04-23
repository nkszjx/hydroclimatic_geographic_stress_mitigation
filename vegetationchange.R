

# Figure 1c-e

library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(scales)
library(showtext)

# Fig1_c
urban_in_dt <- fread("E:/2025/UrbanGreenSpaceWorks/Rcodes/dataset_in_urban.csv")
dim(urban_in_dt)

font_add("times", "times.ttf")
showtext_auto()

showtext_opts(dpi = 1000)

dt_clean <- urban_in_dt[FVC != 0, .(continent, year, FVC, CityArea)]

weighted_mean <- function(x, w) sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE)

continent_year_fvc <- dt_clean %>% 
  .[, .(weighted_fvc = weighted_mean(FVC, CityArea)), by = .(continent, year)] %>% 
  .[, year := as.integer(year)] %>% 
  .[, continent := case_when(
    continent == "NorthAmerica" ~ "North America",
    continent == "SouthAmerica" ~ "South America",
    TRUE ~ continent
  )] %>% 
  .[, continent := factor(continent, levels = c("Asia", "Africa", "North America", "Europe", "Oceania", "South America"))]

continent_colors <- c(
  "Asia" = "#3498db",
  "Africa" = "#95a5a6",
  "North America" = "#2ecc71",
  "South America" = "#9b59b6",
  "Europe" = "#e74c3c",
  "Oceania" = "#f39c12"
)

p1_3 <- ggplot(continent_year_fvc, aes(x = year, y = weighted_fvc, color = continent, group = continent)) +
  geom_line(linewidth = 0.5, alpha = 1) +
  geom_point(size = 1.5, stroke = 0, shape = 19) +
  scale_x_continuous(name = "Year", breaks = seq(2000, 2020, 5), expand = c(0.05,0.05)) +
  scale_y_continuous(name = "Average urban FVC", limits = c(0.2,0.5), breaks = seq(0.2,0.5,0.1), expand = c(0.02,0.02)) +
  scale_color_manual(values = continent_colors, name = "Continent",
                     labels = c("Asia", "Africa", "North\nAmerica", "Europe", "Oceania", "South\nAmerica")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 0.6, color = "black"),
    axis.text.x = element_text(size = 9.5, family = "times", color = "black"),
    axis.text.y = element_text(size = 9.5, family = "times", color = "black"),
    axis.title.x = element_text(size = 10, family = "times", color = "black", margin = margin(t=5)),
    axis.title.y = element_text(size = 10, family = "times", color = "black", margin = margin(r=5)),
    axis.ticks.length = unit(-0.1, "cm"),
    axis.ticks = element_line(color = "black", size = 0.4),
    legend.position = c(0.50, 0.85),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, family = "times", color = "black", lineheight = 0.7),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.x = unit(0.6, "cm"),
    legend.key = element_rect(fill=NA, colour=NA),
    legend.background = element_blank()
  ) +
  guides(color = guide_legend(nrow = 3, ncol =2))

p1_3

ggsave("E:/2025/UrbanGreenSpaceWorks/Rcodes/Fig1_c.png", 
       plot = p1_3, 
       width = 60, 
       height = 60, 
       units = "mm",
       dpi = 1000)



# Fig1_d
urban_2000 <- urban_in_dt[year == 2000] %>% select(continent, CityArea, FVC) %>% rename(area2000 = CityArea, FVC2000 = FVC)
urban_2005 <- urban_in_dt[year == 2005] %>% select(CityArea, FVC) %>% rename(area2005 = CityArea, FVC2005 = FVC)
urban_2010 <- urban_in_dt[year == 2010] %>% select(CityArea, FVC) %>% rename(area2010 = CityArea, FVC2010 = FVC)
urban_2015 <- urban_in_dt[year == 2015] %>% select(CityArea, FVC) %>% rename(area2015 = CityArea, FVC2015 = FVC)
urban_2020 <- urban_in_dt[year == 2020] %>% select(CityArea, FVC) %>% rename(area2020 = CityArea, FVC2020 = FVC)
urban_in_dt0 <- cbind(urban_2000,urban_2005,urban_2010,urban_2015,urban_2020)

urban_in_dt0$FVC2000_2020 <- (urban_in_dt0$FVC2020 - urban_in_dt0$FVC2000)/urban_in_dt0$FVC2000 *100  

urban_in_dt0 <- urban_in_dt0 %>%
  filter(
    !is.na(FVC2000_2020) & !is.infinite(FVC2000_2020)
  )

continent_fvc_20y <- urban_in_dt0 %>% 
  select(continent, area2000, FVC2000_2020) %>% 
  mutate(
    continent = case_when(
      continent == "NorthAmerica" ~ "North America",
      continent == "SouthAmerica" ~ "South America",
      TRUE ~ continent
    ),
    continent = factor(continent, levels = c("Asia", "Africa", "North America", "South America", "Europe", "Oceania"))
  )

continent_colors <- c(
  "Asia" = "#3498db",
  "Africa" = "#95a5a6",
  "North America" = "#2ecc71",
  "South America" = "#9b59b6",
  "Europe" = "#e74c3c",
  "Oceania" = "#f39c12"
)

font_add("times", "times.ttf")
showtext_auto()

showtext_opts(dpi = 1000)

p1_4 <- ggplot(continent_fvc_20y, aes(x = FVC2000_2020, fill = continent)) +
  geom_histogram(
    aes(y = ..density..),
    color = NA,
    linewidth = 0.2,
    binwidth = 2,
    alpha = 0.5
  ) +
  scale_fill_manual(values = continent_colors, name = "Continent") +
  scale_x_continuous(name = "Urban ΔFVC (%, 2000-2020)", expand = c(0.02, 0.02),limits = c(-100,200)) +
  scale_y_continuous(name = "Probability density", expand = c(0.02, 0.02)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9.5, family = "times", color = "black"),
    axis.text.y = element_text(size = 9.5, family = "times", color = "black"),
    axis.title = element_text(size = 10, family = "times", color = "black"),
    axis.ticks.length = unit(-0.1, "cm"),
    axis.ticks = element_line(color = "black", size = 0.4),
    axis.text.x.top = element_text(margin = margin(b = 12)),
    axis.text.y.right = element_text(margin = margin(l = 12)),
    legend.position = c(0.71,0.8),
    legend.text = element_text(size = 9, family = "times", hjust = 0, margin = margin(l = 0, unit = "cm")),
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(0.1, "cm"),
    legend.spacing.x = unit(0, "cm"),
    legend.key = element_rect(fill=NA, colour=NA),
    legend.background = element_blank(),
    legend.title = element_blank()
  )

p1_4

ggsave("E:/2025/UrbanGreenSpaceWorks/Rcodes/Fig1_d.png", 
       plot = p1_4, width = 60, height = 60, units = "mm",dpi=1000)



# Fig1_e
urban_2000 <- urban_in_dt[year == 2000] %>% select(continent, CityArea, FVC) %>% rename(area2000 = CityArea, FVC2000 = FVC)
urban_2005 <- urban_in_dt[year == 2005] %>% select(CityArea, FVC) %>% rename(area2005 = CityArea, FVC2005 = FVC)
urban_2010 <- urban_in_dt[year == 2010] %>% select(CityArea, FVC) %>% rename(area2010 = CityArea, FVC2010 = FVC)
urban_2015 <- urban_in_dt[year == 2015] %>% select(CityArea, FVC) %>% rename(area2015 = CityArea, FVC2015 = FVC)
urban_2020 <- urban_in_dt[year == 2020] %>% select(CityArea, FVC) %>% rename(area2020 = CityArea, FVC2020 = FVC)

urban_in_dt1 <- cbind(urban_2000,urban_2005,urban_2010,urban_2015,urban_2020)

unique(urban_in_dt1$continent)

urban_in_dt1$FVC2000_2005 <- (urban_in_dt1$FVC2005 - urban_in_dt1$FVC2000)/urban_in_dt1$FVC2000 *100
urban_in_dt1$FVC2005_2010 <- (urban_in_dt1$FVC2010 - urban_in_dt1$FVC2005)/urban_in_dt1$FVC2005 *100
urban_in_dt1$FVC2010_2015 <- (urban_in_dt1$FVC2015 - urban_in_dt1$FVC2010)/urban_in_dt1$FVC2010 *100
urban_in_dt1$FVC2015_2020 <- (urban_in_dt1$FVC2020 - urban_in_dt1$FVC2015)/urban_in_dt1$FVC2015 *100
urban_in_dt1$FVC2000_2020 <- (urban_in_dt1$FVC2020 - urban_in_dt1$FVC2000)/urban_in_dt1$FVC2000 *100
urban_in_dt1 <- urban_in_dt1 %>%
  filter(
    across(c(FVC2000_2005, FVC2005_2010, FVC2010_2015, FVC2015_2020, FVC2000_2020),
           ~ !is.na(.) & !is.infinite(.))
  )

weighted_mean <- function(x, w) sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
weighted_sd <- function(x, w) {
  w_mean <- weighted_mean(x, w)
  sqrt(sum(w * (x - w_mean)^2, na.rm = TRUE) / sum(w, na.rm = TRUE))
}

continent_period_stats <- urban_in_dt1 %>% 
  select(continent, area2000, FVC2000_2005, FVC2005_2010, FVC2010_2015, FVC2015_2020) %>% 
  pivot_longer(
    cols = c(FVC2000_2005, FVC2005_2010, FVC2010_2015, FVC2015_2020),
    names_to = "period",
    values_to = "FVC_5y_rate"
  ) %>% 
  mutate(
    period = case_when(
      period == "FVC2000_2005" ~ "2000-2005",
      period == "FVC2005_2010" ~ "2005-2010",
      period == "FVC2010_2015" ~ "2010-2015",
      period == "FVC2015_2020" ~ "2015-2020"
    ),
    period = factor(period, levels = c("2000-2005","2005-2010","2010-2015","2015-2020"))
  ) %>% 
  group_by(continent, period) %>% 
  summarise(
    mean_rate = weighted_mean(FVC_5y_rate, area2000),
    sd_rate = weighted_sd(FVC_5y_rate, area2000),
    .groups = "drop"
  )

continent_period_stats <- continent_period_stats %>% 
  mutate(
    continent = case_when(
      continent == "NorthAmerica" ~ "North America",
      continent == "SouthAmerica" ~ "South America",
      TRUE ~ continent
    )
  )

font_add("times", "times.ttf")
showtext_auto()

showtext_opts(dpi = 1000)

continent_period_stats$continent <- gsub("North America", "North\nAmerica", continent_period_stats$continent)
continent_period_stats$continent <- gsub("South America", "South\nAmerica", continent_period_stats$continent)
p1_5 <- ggplot(continent_period_stats, aes(x = continent, y = mean_rate, fill = period, color=period)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    linewidth = 0.3,
  ) +
  geom_errorbar(
    aes(ymin = mean_rate - sd_rate, ymax = mean_rate + sd_rate),
    position = position_dodge(width = 0.8),
    width = 0.3,
    linewidth = 0.3,
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  geom_vline(xintercept = seq(1.5, 5.5, 1), linetype = "dashed", color = "gray50", linewidth = 0.4) +
  scale_fill_brewer(palette = "Set2", name = "Period", direction = 1) +
  scale_color_brewer(palette = "Set2", name = "Period", direction = 1) +
  scale_y_continuous(name = "Urban ΔFVC (%, per 5-year)", limits = c(-50, 100),
                     breaks = seq(-50, 100, by = 50),
                     labels = scales::number_format(accuracy = 1)) +
  scale_x_discrete(limits = c("Asia", "Africa", "North\nAmerica", "South\nAmerica", "Europe", "Oceania")) +
  labs(title = NULL, x="") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9.5, family = "times", color = "black", lineheight = 0.7),
    axis.text.y = element_text(size = 9.5, family = "times", color = "black"),
    axis.title = element_text(size = 10, family = "times", color = "black"),
    axis.ticks.length = unit(-0.1, "cm"),
    axis.ticks = element_line(color = "black", size = 0.4),
    axis.text.x.top = element_text(margin = margin(b = 12)),
    axis.text.y.right = element_text(margin = margin(l = 12)),
    legend.position = c(0.2,0.85),
    legend.text = element_text(size = 9, family = "times", hjust = 0, margin = margin(l = 0.1, unit = "cm")),
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(0.1, "cm"),
    legend.spacing.x = unit(0, "cm"),
    legend.key = element_rect(fill=NA, colour=NA),
    legend.background = element_blank(),
    legend.title = element_blank()
  )

p1_5

ggsave("E:/2025/UrbanGreenSpaceWorks/Rcodes/Fig1_e.png", 
       plot = p1_5, width = 90, height = 62, units = "mm",dpi=1000)



