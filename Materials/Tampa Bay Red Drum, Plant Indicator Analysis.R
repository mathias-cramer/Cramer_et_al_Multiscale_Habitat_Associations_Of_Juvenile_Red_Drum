#Project Details_______________________________________________________________:

#Title: Multiscale Habitat Associations of Juvenile Red Drum (Sciaenops 
#ocellatus) in Tampa Bay, Florida Rivers

#Species indicator analysis, 0-50mm SL_________________________________________:

library(dplyr)
library(ggplot2)
library(indicspecies)
library(permute)
library(stringr)
library(ggtext)

#Grouping variable
q25 <- quantile(data$N_0_50_SL, 0.25, na.rm = TRUE)
q75 <- quantile(data$N_0_50_SL, 0.75, na.rm = TRUE)

group <- ifelse(data$N_0_50_SL <= q25, "Low",
                ifelse(data$N_0_50_SL >= q75, "High", NA))

keep <- !is.na(group)

#Species data
veg_data <- data %>% dplyr::select(48:69)

#Correct species name spelling
names(veg_data)[names(veg_data) == "Acrostichum_danaefolium"] <- "Acrostichum_danaeifolium"

veg_sub <- veg_data[keep, ]

#Remove species not present in either group
veg_sub <- veg_sub[, colSums(veg_sub > 0, na.rm = TRUE) > 0]

#Indicator analysis
set.seed(123)
ind_species <- multipatt(
  veg_sub,
  group[keep],
  control = how(nperm = 9999)
)

#Build results table
ind_vals <- data.frame(
  Species   = rownames(ind_species$sign),
  Indicator = ind_species$sign$stat,
  p.value   = ind_species$sign$p.value,
  High      = ind_species$sign$s.High,
  Low       = ind_species$sign$s.Low
)

#Benjamini-Hochberg false discovery rate adjustment
#Only taxa with a permutation-based p-value are included
ind_vals$p.adjusted <- NA_real_

testable <- !is.na(ind_vals$p.value)

ind_vals$p.adjusted[testable] <- p.adjust(
  ind_vals$p.value[testable],
  method = "BH"
)

#Keep NA p-values so untestable/both-group species remain visible
#ind_vals <- ind_vals %>% filter(!is.na(p.value))  # REMOVE / DO NOT USE

#Labels
# Italicize scientific names, but leave "sp." and "spp." in roman type
ind_vals <- ind_vals %>%
  mutate(
    Species_label = case_when(
      str_detect(Species, "_sp$") ~ paste0(
        "*",
        str_replace(str_remove(Species, "_sp$"), "_", " "),
        "* sp."
      ),
      str_detect(Species, "_spp$") ~ paste0(
        "*",
        str_replace(str_remove(Species, "_spp$"), "_", " "),
        "* spp."
      ),
      TRUE ~ paste0(
        "*",
        str_replace_all(Species, "_", " "),
        "*"
      )
    )
  )

#Direction / association
ind_vals$Association <- dplyr::case_when(
  ind_vals$High == 1 & ind_vals$Low == 1 ~ "Both",
  ind_vals$High == 1 & ind_vals$Low == 0 ~ "High",
  ind_vals$High == 0 & ind_vals$Low == 1 ~ "Low",
  TRUE ~ "None"
)

#Three-level indicator test status based on FDR-adjusted p-values
ind_vals$Indicator_Test <- dplyr::case_when(
  is.na(ind_vals$p.value) ~ "Not testable",
  ind_vals$p.adjusted < 0.05 ~ "Significant",
  TRUE ~ "Not significant"
)

#Lock order
ind_vals <- ind_vals %>%
  arrange(Indicator)

ind_vals$Species_label <- factor(
  ind_vals$Species_label,
  levels = ind_vals$Species_label
)

#Test-status symbols based on FDR-adjusted p-values
ind_vals$Test_Label <- dplyr::case_when(
  is.na(ind_vals$p.value) ~ "\u2020",       # dagger = Not testable
  ind_vals$p.adjusted < 0.05 ~ "*",         # asterisk = FDR significant
  TRUE ~ ""                                 # blank = Not significant
)

#Association
ind_vals$Association <- dplyr::case_when(
  ind_vals$High == 1 & ind_vals$Low == 1 ~ "Both",
  ind_vals$High == 1 & ind_vals$Low == 0 ~ "High",
  ind_vals$High == 0 & ind_vals$Low == 1 ~ "Low",
  TRUE ~ "None"
)

#Calculate mean percent cover separately for High and Low abundance groups
veg_group <- group[keep]

mean_cover_high <- colMeans(
  veg_sub[veg_group == "High", , drop = FALSE],
  na.rm = TRUE
)

mean_cover_low <- colMeans(
  veg_sub[veg_group == "Low", , drop = FALSE],
  na.rm = TRUE
)

cover_df <- data.frame(
  Species = names(mean_cover_high),
  MeanCover_High = mean_cover_high,
  MeanCover_Low = mean_cover_low
)

ind_vals <- ind_vals %>%
  left_join(cover_df, by = "Species")

#Plot
ggplot(ind_vals, aes(x = Species_label,
                     y = Indicator,
                     fill = Association)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.4) +
  geom_text(aes(label = Test_Label),
            hjust = -0.3,
            size = 5,
            fontface = "bold") +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c(
      "High" = "#f0a377",  # Okabe-Ito vermillion
      "Low"  = "#8baad1",  # Okabe-Ito blue
      "Both" = "#F7F7F7",  # Okabe-Ito white
      "None" = "grey70"
    ),
    name = "Association"
  ) +
  labs(
    x = NULL,
    y = "Indicator value",
    #title = "Indicator Species — Red Drum (0–50 mm SL)",
    caption = "\u2020 Associated with both High and Low groups; permutation-based p value not applicable"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.y = ggtext::element_markdown(
      size = 10,
      family = "sans"
    ),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0),
    plot.margin = margin(5.5, 30, 5.5, 5.5)
  ) +
  expand_limits(y = max(ind_vals$Indicator, na.rm = TRUE) * 1.12)

#750x750

u50_ind <- ind_vals
u50_ST <- data.frame(
  Species = rownames(ind_species$sign),
  Indicator = ind_species$sign$stat,
  p.value = ind_species$sign$p.value,
  p.adjusted = ind_vals$p.adjusted[
    match(rownames(ind_species$sign), ind_vals$Species)
  ],
  High = ind_species$sign$s.High,
  Low = ind_species$sign$s.Low,
  A_High = ind_species$A[, "High"],
  A_Low = ind_species$A[, "Low"],
  B_High = ind_species$B[, "High"],
  B_Low = ind_species$B[, "Low"],
  MeanCover_High = ind_vals$MeanCover_High[
    match(rownames(ind_species$sign), ind_vals$Species)
  ],
  MeanCover_Low = ind_vals$MeanCover_Low[
    match(rownames(ind_species$sign), ind_vals$Species)
  ]
)
u50_ST

#Species indicator analysis, 51-100mm SL_______________________________________:

#Grouping variable
q25 <- quantile(data$N_51_100_SL, 0.25, na.rm = TRUE)
q75 <- quantile(data$N_51_100_SL, 0.75, na.rm = TRUE)

group <- ifelse(data$N_51_100_SL <= q25, "Low",
                ifelse(data$N_51_100_SL >= q75, "High", NA))

keep <- !is.na(group)

#Species data
veg_data <- data %>% dplyr::select(48:69)

#Correct species name spelling
names(veg_data)[names(veg_data) == "Acrostichum_danaefolium"] <- "Acrostichum_danaeifolium"

veg_sub <- veg_data[keep, ]

#Remove species not present in either group
veg_sub <- veg_sub[, colSums(veg_sub > 0, na.rm = TRUE) > 0]

#Indicator analysis
set.seed(123)
ind_species <- multipatt(
  veg_sub,
  group[keep],
  control = how(nperm = 9999)
)

# Build results table
ind_vals <- data.frame(
  Species   = rownames(ind_species$sign),
  Indicator = ind_species$sign$stat,
  p.value   = ind_species$sign$p.value,
  High      = ind_species$sign$s.High,
  Low       = ind_species$sign$s.Low
)

#Benjamini-Hochberg false discovery rate adjustment
#Only taxa with a permutation-based p-value are included
ind_vals$p.adjusted <- NA_real_

testable <- !is.na(ind_vals$p.value)

ind_vals$p.adjusted[testable] <- p.adjust(
  ind_vals$p.value[testable],
  method = "BH"
)

#Keep NA p-values so untestable/both-group species remain visible
#ind_vals <- ind_vals %>% filter(!is.na(p.value))  # REMOVE / DO NOT USE

# Labels
# Italicize scientific names, but leave "sp." and "spp." in roman type
ind_vals <- ind_vals %>%
  mutate(
    Species_label = case_when(
      str_detect(Species, "_sp$") ~ paste0(
        "*",
        str_replace(str_remove(Species, "_sp$"), "_", " "),
        "* sp."
      ),
      str_detect(Species, "_spp$") ~ paste0(
        "*",
        str_replace(str_remove(Species, "_spp$"), "_", " "),
        "* spp."
      ),
      TRUE ~ paste0(
        "*",
        str_replace_all(Species, "_", " "),
        "*"
      )
    )
  )

#Direction / association
ind_vals$Association <- dplyr::case_when(
  ind_vals$High == 1 & ind_vals$Low == 1 ~ "Both",
  ind_vals$High == 1 & ind_vals$Low == 0 ~ "High",
  ind_vals$High == 0 & ind_vals$Low == 1 ~ "Low",
  TRUE ~ "None"
)

#Three-level indicator test status based on FDR-adjusted p-values
ind_vals$Indicator_Test <- dplyr::case_when(
  is.na(ind_vals$p.value) ~ "Not testable",
  ind_vals$p.adjusted < 0.05 ~ "Significant",
  TRUE ~ "Not significant"
)

#Lock order
ind_vals <- ind_vals %>%
  arrange(Indicator)

ind_vals$Species_label <- factor(
  ind_vals$Species_label,
  levels = ind_vals$Species_label
)

#Test-status symbols based on FDR-adjusted p-values
ind_vals$Test_Label <- dplyr::case_when(
  is.na(ind_vals$p.value) ~ "\u2020",       # dagger = Not testable
  ind_vals$p.adjusted < 0.05 ~ "*",         # asterisk = FDR significant
  TRUE ~ ""                                 # blank = Not significant
)

#Association
ind_vals$Association <- dplyr::case_when(
  ind_vals$High == 1 & ind_vals$Low == 1 ~ "Both",
  ind_vals$High == 1 & ind_vals$Low == 0 ~ "High",
  ind_vals$High == 0 & ind_vals$Low == 1 ~ "Low",
  TRUE ~ "None"
)

#Calculate mean percent cover separately for High and Low abundance groups
veg_group <- group[keep]

mean_cover_high <- colMeans(
  veg_sub[veg_group == "High", , drop = FALSE],
  na.rm = TRUE
)

mean_cover_low <- colMeans(
  veg_sub[veg_group == "Low", , drop = FALSE],
  na.rm = TRUE
)

cover_df <- data.frame(
  Species = names(mean_cover_high),
  MeanCover_High = mean_cover_high,
  MeanCover_Low = mean_cover_low
)

ind_vals <- ind_vals %>%
  left_join(cover_df, by = "Species")

#Plot
ggplot(ind_vals, aes(x = Species_label,
                     y = Indicator,
                     fill = Association)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.4) +
  geom_text(aes(label = Test_Label),
            hjust = -0.3,
            size = 5,
            fontface = "bold") +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c(
      "High" = "#f0a377",  # Okabe-Ito vermillion
      "Low"  = "#8baad1",  # Okabe-Ito blue
      "Both" = "#F7F7F7",  # Okabe-Ito white
      "None" = "grey70"
    ),
    name = "Association"
  ) +
  labs(
    x = NULL,
    y = "Indicator value",
    #title = "Indicator Species — Red Drum (51–100 mm SL)",
    caption = "\u2020 Associated with both High and Low groups; permutation-based p value not applicable"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.y = ggtext::element_markdown(
      size = 10,
      family = "sans"
    ),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0),
    plot.margin = margin(5.5, 30, 5.5, 5.5)
  ) +
  expand_limits(y = max(ind_vals$Indicator, na.rm = TRUE) * 1.12)

#750x750

u100_ind <- ind_vals
u100_ST <- data.frame(
  Species = rownames(ind_species$sign),
  Indicator = ind_species$sign$stat,
  p.value = ind_species$sign$p.value,
  p.adjusted = ind_vals$p.adjusted[
    match(rownames(ind_species$sign), ind_vals$Species)
  ],
  High = ind_species$sign$s.High,
  Low = ind_species$sign$s.Low,
  A_High = ind_species$A[, "High"],
  A_Low = ind_species$A[, "Low"],
  B_High = ind_species$B[, "High"],
  B_Low = ind_species$B[, "Low"],
  MeanCover_High = ind_vals$MeanCover_High[
    match(rownames(ind_species$sign), ind_vals$Species)
  ],
  MeanCover_Low = ind_vals$MeanCover_Low[
    match(rownames(ind_species$sign), ind_vals$Species)
  ]
)
u100_ST

#Summary
u50_ind
u100_ind

u50_ST
u100_ST

#END