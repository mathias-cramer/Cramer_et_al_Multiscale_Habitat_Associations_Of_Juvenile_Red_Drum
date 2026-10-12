#Project Details______________________________________________________________:

#Title: Multiscale Habitat Associations of Juvenile Red Drum (Sciaenops 
#ocellatus) in Tampa Bay, Florida Rivers

#Analysis 1 — What ecological characteristics differ among Shore Types?

library(dplyr)
library(FSA)

group_var <- "ShoreType"

continuous_vars <- c(
  "RvrWdt_",
  "ChnnlDp",
  "DstncU_",
  "TidlSlp",
  "BankSlp",
  "AvgBslV",
  "AvgBslL",
  "AvgBslW",
  "AvgCnpy",
  "AvgVgH_",
  "LvOystP",
  "OysterArea",
  "Tr_IC_P",
  "Tr_OM_P",
  "MC_IC_P",
  "MC_OM_P",
  "StartDepth",
  "temperature",
  "salinity",
  "dissolvedO2",
  "Mangrove",
  "SaltMarsh"
)

categorical_vars <- c(
  "River",
  "AvgInnd",
  "Bayou",
  "bottom",
  "VegeSummary"
)

#Kruskal-Wallis Tests

shore_kw <- lapply(continuous_vars, function(v){
  
  test <- kruskal.test(
    as.formula(paste(v, "~", group_var)),
    data = work_data
  )
  
  data.frame(
    Variable = v,
    ChiSq = test$statistic,
    df = test$parameter,
    p = test$p.value
  )
  
})

shore_kw <- bind_rows(shore_kw) %>%
  mutate(
    p_FDR = p.adjust(p, method="BH")
  ) %>%
  arrange(p)

shore_kw

#Significant Variables

sig_vars <- shore_kw %>%
  filter(p < 0.05) %>%
  pull(Variable)

sig_vars

#Dunn Tests

shore_dunn <- lapply(sig_vars, function(v){
  
  dunnTest(
    as.formula(paste(v, "~", group_var)),
    data = work_data,
    method = "bh"
  )$res
  
})

names(shore_dunn) <- sig_vars

shore_dunn

#Chi-square Tests

shore_cat <- lapply(categorical_vars, function(v){
  
  tab <- table(
    work_data[[group_var]],
    work_data[[v]]
  )
  
  test <- suppressWarnings(chisq.test(tab))
  
  data.frame(
    Variable = v,
    ChiSq = test$statistic,
    df = test$parameter,
    p = test$p.value
  )
  
})

shore_cat <- bind_rows(shore_cat) %>%
  mutate(
    p_FDR = p.adjust(p, method="BH")
  ) %>%
  arrange(p)

shore_cat

#Descriptive summaries
library(tidyr)

shore_summary <- work_data %>%
  dplyr::select(all_of(group_var), all_of(sig_vars)) %>%
  pivot_longer(
    cols = -all_of(group_var),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  group_by(
    Variable,
    !!sym(group_var)
  ) %>%
  summarise(
    n = sum(!is.na(Value)),
    Mean = mean(Value, na.rm=TRUE),
    SD = sd(Value, na.rm=TRUE),
    Median = median(Value, na.rm=TRUE),
    Q1 = quantile(Value,0.25,na.rm=TRUE),
    Q3 = quantile(Value,0.75,na.rm=TRUE),
    .groups="drop"
  )

shore_summary

#Cross-tabs for significant categorical variables

sig_cat <- shore_cat %>%
  filter(p < 0.05) %>%
  pull(Variable)

for(v in sig_cat){
  
  cat("\n")
  cat("=====================================\n")
  cat(v,"\n")
  cat("=====================================\n")
  
  print(table(
    work_data[[group_var]],
    work_data[[v]]
  ))
  
  cat("\nProportions\n")
  
  print(prop.table(
    table(
      work_data[[group_var]],
      work_data[[v]]
    ),
    margin=1
  ))
  
}

table(work_data$ShoreType, work_data$AvgInnd)

prop.table(
  table(work_data$ShoreType, work_data$AvgInnd),
  margin = 1
)

table(work_data$ShoreType, work_data$VegeSummary)

prop.table(
  table(work_data$ShoreType, work_data$VegeSummary),
  margin = 1
)

#Analysis 2 — What ecological characteristics differ among Bottom Types?

group_var <- "bottom"

continuous_vars <- c(
  "RvrWdt_",
  "ChnnlDp",
  "DstncU_",
  "TidlSlp",
  "BankSlp",
  "AvgBslV",
  "AvgBslL",
  "AvgBslW",
  "AvgCnpy",
  "AvgVgH_",
  "LvOystP",
  "OysterArea",
  "Tr_IC_P",
  "Tr_OM_P",
  "MC_IC_P",
  "MC_OM_P",
  "StartDepth",
  "temperature",
  "salinity",
  "dissolvedO2",
  "Mangrove",
  "SaltMarsh"
)

categorical_vars <- c(
  "River",
  "AvgInnd",
  "Bayou",
  "ShoreType",
  "VegeSummary"
)

#Kruskal-Wallis Tests

shore_kw <- lapply(continuous_vars, function(v){
  
  test <- kruskal.test(
    as.formula(paste(v, "~", group_var)),
    data = work_data
  )
  
  data.frame(
    Variable = v,
    ChiSq = test$statistic,
    df = test$parameter,
    p = test$p.value
  )
  
})

shore_kw <- bind_rows(shore_kw) %>%
  mutate(
    p_FDR = p.adjust(p, method="BH")
  ) %>%
  arrange(p)

shore_kw

#Significant Variables

sig_vars <- shore_kw %>%
  filter(p < 0.05) %>%
  pull(Variable)

sig_vars

#Dunn Tests

shore_dunn <- lapply(sig_vars, function(v){
  
  dunnTest(
    as.formula(paste(v, "~", group_var)),
    data = work_data,
    method = "bh"
  )$res
  
})

names(shore_dunn) <- sig_vars

shore_dunn

#Chi-square Tests

shore_cat <- lapply(categorical_vars, function(v){
  
  tab <- table(
    work_data[[group_var]],
    work_data[[v]]
  )
  
  test <- suppressWarnings(chisq.test(tab))
  
  data.frame(
    Variable = v,
    ChiSq = test$statistic,
    df = test$parameter,
    p = test$p.value
  )
  
})

shore_cat <- bind_rows(shore_cat) %>%
  mutate(
    p_FDR = p.adjust(p, method="BH")
  ) %>%
  arrange(p)

shore_cat

#Descriptive summaries
library(tidyr)

shore_summary <- work_data %>%
  dplyr::select(all_of(group_var), all_of(sig_vars)) %>%
  pivot_longer(
    cols = -all_of(group_var),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  group_by(
    Variable,
    !!sym(group_var)
  ) %>%
  summarise(
    n = sum(!is.na(Value)),
    Mean = mean(Value, na.rm=TRUE),
    SD = sd(Value, na.rm=TRUE),
    Median = median(Value, na.rm=TRUE),
    Q1 = quantile(Value,0.25,na.rm=TRUE),
    Q3 = quantile(Value,0.75,na.rm=TRUE),
    .groups="drop"
  )

shore_summary

#Cross-tabs for significant categorical variables

sig_cat <- shore_cat %>%
  filter(p < 0.05) %>%
  pull(Variable)

for(v in sig_cat){
  
  cat("\n")
  cat("=====================================\n")
  cat(v,"\n")
  cat("=====================================\n")
  
  print(table(
    work_data[[group_var]],
    work_data[[v]]
  ))
  
  cat("\nProportions\n")
  
  print(prop.table(
    table(
      work_data[[group_var]],
      work_data[[v]]
    ),
    margin=1
  ))
  
}

table(work_data$ShoreType, work_data$AvgInnd)

prop.table(
  table(work_data$ShoreType, work_data$AvgInnd),
  margin = 1
)

table(work_data$ShoreType, work_data$VegeSummary)

prop.table(
  table(work_data$ShoreType, work_data$VegeSummary),
  margin = 1
)

#Analysis 3 — What ecological characteristics differ by Bayou?

group_var <- "Bayou"

continuous_vars <- c(
  "RvrWdt_",
  "ChnnlDp",
  "DstncU_",
  "TidlSlp",
  "BankSlp",
  "AvgBslV",
  "AvgBslL",
  "AvgBslW",
  "AvgCnpy",
  "AvgVgH_",
  "LvOystP",
  "OysterArea",
  "Tr_IC_P",
  "Tr_OM_P",
  "MC_IC_P",
  "MC_OM_P",
  "StartDepth",
  "temperature",
  "salinity",
  "dissolvedO2",
  "Mangrove",
  "SaltMarsh"
)

categorical_vars <- c(
  "River",
  "AvgInnd",
  "bottom",
  "ShoreType",
  "VegeSummary"
)

#Kruskal-Wallis Tests

shore_kw <- lapply(continuous_vars, function(v){
  
  test <- kruskal.test(
    as.formula(paste(v, "~", group_var)),
    data = work_data
  )
  
  data.frame(
    Variable = v,
    ChiSq = test$statistic,
    df = test$parameter,
    p = test$p.value
  )
  
})

shore_kw <- bind_rows(shore_kw) %>%
  mutate(
    p_FDR = p.adjust(p, method="BH")
  ) %>%
  arrange(p)

shore_kw

#Significant Variables

sig_vars <- shore_kw %>%
  filter(p < 0.05) %>%
  pull(Variable)

sig_vars

#Variables significant at the unadjusted p < 0.05 level
sig_vars <- shore_kw %>%
  filter(p < 0.05) %>%
  pull(Variable)

#Wilcoxon rank-sum tests for Bayou: no vs yes
bayou_wilcox <- lapply(sig_vars, function(v) {
  
  test <- wilcox.test(
    as.formula(paste(v, "~", group_var)),
    data = work_data,
    exact = FALSE,
    conf.int = FALSE
  )
  
  data.frame(
    Variable = v,
    W = unname(test$statistic),
    p = test$p.value
  )
})

bayou_wilcox <- bind_rows(bayou_wilcox) %>%
  mutate(p_FDR = p.adjust(p, method = "BH")) %>%
  arrange(p)

bayou_wilcox

#Chi-square Tests

shore_cat <- lapply(categorical_vars, function(v){
  
  tab <- table(
    work_data[[group_var]],
    work_data[[v]]
  )
  
  test <- suppressWarnings(chisq.test(tab))
  
  data.frame(
    Variable = v,
    ChiSq = test$statistic,
    df = test$parameter,
    p = test$p.value
  )
  
})

shore_cat <- bind_rows(shore_cat) %>%
  mutate(
    p_FDR = p.adjust(p, method="BH")
  ) %>%
  arrange(p)

shore_cat

#Descriptive summaries
library(tidyr)

shore_summary <- work_data %>%
  dplyr::select(all_of(group_var), all_of(sig_vars)) %>%
  pivot_longer(
    cols = -all_of(group_var),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  group_by(
    Variable,
    !!sym(group_var)
  ) %>%
  summarise(
    n = sum(!is.na(Value)),
    Mean = mean(Value, na.rm=TRUE),
    SD = sd(Value, na.rm=TRUE),
    Median = median(Value, na.rm=TRUE),
    Q1 = quantile(Value,0.25,na.rm=TRUE),
    Q3 = quantile(Value,0.75,na.rm=TRUE),
    .groups="drop"
  )

shore_summary

#Cross-tabs for significant categorical variables

sig_cat <- shore_cat %>%
  filter(p < 0.05) %>%
  pull(Variable)

for(v in sig_cat){
  
  cat("\n")
  cat("=====================================\n")
  cat(v,"\n")
  cat("=====================================\n")
  
  print(table(
    work_data[[group_var]],
    work_data[[v]]
  ))
  
  cat("\nProportions\n")
  
  print(prop.table(
    table(
      work_data[[group_var]],
      work_data[[v]]
    ),
    margin=1
  ))
  
}

#Analysis 4 — How do 0-50mm abundances differ by River?
library(dplyr)
library(emmeans)
library(ggplot2)

names(work_data)

#1. Raw abundance summaries by river

river_raw_summary <- work_data %>%
  group_by(River) %>%
  summarise(
    n = sum(!is.na(N_0_50_SL)),
    Mean = mean(N_0_50_SL, na.rm = TRUE),
    SD = sd(N_0_50_SL, na.rm = TRUE),
    Median = median(N_0_50_SL, na.rm = TRUE),
    Q1 = quantile(N_0_50_SL, 0.25, na.rm = TRUE),
    Q3 = quantile(N_0_50_SL, 0.75, na.rm = TRUE),
    Min = min(N_0_50_SL, na.rm = TRUE),
    Max = max(N_0_50_SL, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean))

river_raw_summary

#2. Model-adjusted estimated marginal means

river_emm <- emmeans(U50_model, ~ River)

river_emm_table <- as.data.frame(river_emm) %>%
  arrange(desc(emmean)) %>%
  mutate(Rank = row_number()) %>%
  dplyr::select(Rank, everything())

river_emm_table

#Tukey-adjusted comparisons among rivers
pairs(river_emm, adjust = "tukey")

library(multcomp)

river_cld <- cld(
  river_emm,
  adjust = "tukey",
  Letters = letters
) %>%
  as.data.frame() %>%
  arrange(desc(emmean))

river_cld

#3. Plot adjusted river means

ggplot(river_emm_table,
       aes(x = reorder(River, emmean), y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.15
  ) +
  coord_flip() +
  labs(
    x = "River",
    y = "Estimated marginal mean",
    title = "Model-adjusted abundance of 0–50 mm Red Drum by river"
  ) +
  theme_classic()

#Analysis 5 — Why might Alafia have such lower abundances?

continuous_vars <- c(
  "RvrWdt_",
  "ChnnlDp",
  "DstncU_",
  "TidlSlp",
  "BankSlp",
  "AvgBslV",
  "AvgBslL",
  "AvgBslW",
  "AvgCnpy",
  "AvgVgH_",
  "LvOystP",
  "OysterArea",
  "RipRap",
  "AnySwll",
  "anySdmn",
  "Tr_IC_P",
  "Tr_OM_P",
  "MC_IC_P",
  "MC_OM_P",
  "StartDepth",
  "temperature",
  "salinity",
  "dissolvedO2",
  "Mangrove",
  "SaltMarsh"
)

cat_vars <- c(
  "AvgInnd",
  "ShoreType",
  "Bayou",
  "bottom",
  "VegeSummary"
)

#For continuous variables, use Kruskal-Wallis:
library(FSA)
river_kw <- lapply(continuous_vars, function(v){
  test <- kruskal.test(
    as.formula(paste(v,"~ River")),
    data = work_data
  )
  data.frame(
    Variable = v,
    ChiSq = test$statistic,
    df = test$parameter,
    p = test$p.value
  )
})

river_kw <- bind_rows(river_kw) %>%
  arrange(p)

river_kw

#Post hoc:
vars <- c(
  "salinity",
  "RvrWdt_",
  "ChnnlDp",
  "LvOystP",
  "temperature",
  "SaltMarsh",
  "MC_IC_P",
  "OysterArea"
)

library(rstatix)
for(v in vars){
  
  cat("\n====================================================\n")
  cat(v,"\n")
  cat("====================================================\n\n")
  
  # Summary statistics
  summary_table <-
    work_data %>%
    group_by(River) %>%
    summarise(
      n = sum(!is.na(.data[[v]])),
      Mean = mean(.data[[v]], na.rm = TRUE),
      Median = median(.data[[v]], na.rm = TRUE),
      SD = sd(.data[[v]], na.rm = TRUE),
      .groups="drop"
    ) %>%
    arrange(desc(Mean))
  
  print(summary_table)
  
  #Kruskal-Wallis
  kw <-
    kruskal_test(
      work_data,
      formula = as.formula(paste(v,"~ River"))
    )
  
  print(kw)
  
  #Dunn test
  dunn <-
    dunn_test(
      work_data,
      formula = as.formula(paste(v,"~ River")),
      p.adjust.method = "BH"
    )
  
  print(dunn)
}

#For categorical variables, use Chi-square:
river_cat <- lapply(cat_vars,function(v){
  tab <- table(work_data$River,
               work_data[[v]])
  test <- chisq.test(tab)
  data.frame(
    Variable=v,
    ChiSq=test$statistic,
    df=test$parameter,
    p=test$p.value
  )
})

river_cat <- bind_rows(river_cat) %>%
  arrange(p)

river_cat

#Analysis 6 — How do 0-50mm abundances differ by Shore Type?
library(dplyr)
library(emmeans)
library(ggplot2)

#1. Raw abundance summaries by shore type

shoretype_raw_summary <- work_data %>%
  group_by(ShoreType) %>%
  summarise(
    n = sum(!is.na(N_0_50_SL)),
    Mean = mean(N_0_50_SL, na.rm = TRUE),
    SD = sd(N_0_50_SL, na.rm = TRUE),
    Median = median(N_0_50_SL, na.rm = TRUE),
    Q1 = quantile(N_0_50_SL, 0.25, na.rm = TRUE),
    Q3 = quantile(N_0_50_SL, 0.75, na.rm = TRUE),
    Min = min(N_0_50_SL, na.rm = TRUE),
    Max = max(N_0_50_SL, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean))

shoretype_raw_summary

#2. Model-adjusted estimated marginal means

shoretype_emm <- emmeans(U50_model, ~ ShoreType)

shoretype_emm_table <- as.data.frame(shoretype_emm) %>%
  arrange(desc(emmean)) %>%
  mutate(Rank = row_number()) %>%
  dplyr::select(Rank, everything())

shoretype_emm_table

#Tukey-adjusted comparisons among shore types
shoretype_pairs <- pairs(
  shoretype_emm,
  adjust = "tukey"
)

shoretype_pairs

#Compact-letter display
shoretype_cld <- cld(
  shoretype_emm,
  adjust = "tukey",
  Letters = letters
) %>%
  as.data.frame() %>%
  mutate(.group = trimws(.group)) %>%
  arrange(desc(emmean))

shoretype_cld

#3. Plot adjusted shore-type means

ggplot(
  shoretype_emm_table,
  aes(x = reorder(ShoreType, emmean), y = emmean)
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.15
  ) +
  coord_flip() +
  labs(
    x = "Shore type",
    y = "Estimated marginal mean",
    title = "Model-adjusted abundance of 0–50 mm Red Drum by shore type"
  ) +
  theme_classic()

#Analysis 7 — How do 0-50mm abundances differ by Bottom Type?
library(dplyr)
library(emmeans)
library(ggplot2)
library(multcomp)

names(work_data)

#1. Raw abundance summaries by bottom type

bottom_raw_summary <- work_data %>%
  group_by(bottom) %>%
  summarise(
    n = sum(!is.na(N_0_50_SL)),
    Mean = mean(N_0_50_SL, na.rm = TRUE),
    SD = sd(N_0_50_SL, na.rm = TRUE),
    Median = median(N_0_50_SL, na.rm = TRUE),
    Q1 = quantile(N_0_50_SL, 0.25, na.rm = TRUE),
    Q3 = quantile(N_0_50_SL, 0.75, na.rm = TRUE),
    Min = min(N_0_50_SL, na.rm = TRUE),
    Max = max(N_0_50_SL, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean))

bottom_raw_summary

#2. Model-adjusted estimated marginal means

bottom_emm <- emmeans(U50_model, ~ bottom)

bottom_emm_table <- as.data.frame(bottom_emm) %>%
  arrange(desc(emmean)) %>%
  mutate(Rank = row_number()) %>%
  dplyr::select(Rank, everything())

bottom_emm_table

#Tukey-adjusted comparisons among bottom types
bottom_pairs <- pairs(
  bottom_emm,
  adjust = "tukey"
)

bottom_pairs

#Compact-letter display
bottom_cld <- cld(
  bottom_emm,
  adjust = "tukey",
  Letters = letters
) %>%
  as.data.frame() %>%
  mutate(.group = trimws(.group)) %>%
  arrange(desc(emmean))

bottom_cld

#3. Plot adjusted bottom-type means

ggplot(
  bottom_emm_table,
  aes(x = reorder(bottom, emmean), y = emmean)
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.15
  ) +
  coord_flip() +
  labs(
    x = "Bottom type",
    y = "Estimated marginal mean",
    title = "Model-adjusted abundance of 0–50 mm Red Drum by bottom type"
  ) +
  theme_classic()

#Analysis 8 — How do 51-100mm abundances differ by Bayou?
library(dplyr)
library(emmeans)
library(ggplot2)
library(multcomp)

names(work_data)

#1. Raw abundance summaries by bayou occurrence

bayou_raw_summary <- work_data %>%
  group_by(Bayou) %>%
  summarise(
    n = sum(!is.na(N_51_100_SL)),
    Mean = mean(N_51_100_SL, na.rm = TRUE),
    SD = sd(N_51_100_SL, na.rm = TRUE),
    Median = median(N_51_100_SL, na.rm = TRUE),
    Q1 = quantile(N_51_100_SL, 0.25, na.rm = TRUE),
    Q3 = quantile(N_51_100_SL, 0.75, na.rm = TRUE),
    Min = min(N_51_100_SL, na.rm = TRUE),
    Max = max(N_51_100_SL, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean))

bayou_raw_summary

#2. Model-adjusted estimated marginal means

bayou_emm <- emmeans(U100_model, ~ Bayou)

bayou_emm_table <- as.data.frame(bayou_emm) %>%
  arrange(desc(emmean)) %>%
  mutate(Rank = row_number()) %>%
  dplyr::select(Rank, everything())

bayou_emm_table

#Pairwise comparison between bayou categories
bayou_pairs <- pairs(
  bayou_emm,
  adjust = "tukey"
)

bayou_pairs

#Compact-letter display
bayou_cld <- cld(
  bayou_emm,
  adjust = "tukey",
  Letters = letters
) %>%
  as.data.frame() %>%
  mutate(.group = trimws(.group)) %>%
  arrange(desc(emmean))

bayou_cld

#3. Plot adjusted bayou means

ggplot(
  bayou_emm_table,
  aes(x = reorder(Bayou, emmean), y = emmean)
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.15
  ) +
  coord_flip() +
  labs(
    x = "Bayou occurrence",
    y = "Estimated marginal mean",
    title = "Model-adjusted abundance of 51–100 mm Red Drum by bayou occurrence"
  ) +
  theme_classic()

#Analysis 9 — Why might bayou sites have higher abundances?

continuous_vars <- c(
  "RvrWdt_",
  "ChnnlDp",
  "DstncU_",
  "TidlSlp",
  "BankSlp",
  "AvgBslV",
  "AvgBslL",
  "AvgBslW",
  "AvgCnpy",
  "AvgVgH_",
  "LvOystP",
  "OysterArea",
  "RipRap",
  "AnySwll",
  "anySdmn",
  "Tr_IC_P",
  "Tr_OM_P",
  "MC_IC_P",
  "MC_OM_P",
  "StartDepth",
  "temperature",
  "salinity",
  "dissolvedO2",
  "Mangrove",
  "SaltMarsh"
)

cat_vars <- c(
  "AvgInnd",
  "ShoreType",
  "River",
  "bottom",
  "VegeSummary"
)

#Continuous variables (Kruskal-Wallis)

library(FSA)
library(dplyr)

bayou_kw <- lapply(continuous_vars, function(v){
  
  test <- kruskal.test(
    as.formula(paste(v, "~ Bayou")),
    data = work_data
  )
  
  data.frame(
    Variable = v,
    ChiSq = test$statistic,
    df = test$parameter,
    p = test$p.value
  )
  
})

bayou_kw <- bind_rows(bayou_kw) %>%
  arrange(p)

bayou_kw

#Detailed summaries for significant variables

library(rstatix)

vars <- bayou_kw %>%
  filter(p < 0.05) %>%
  pull(Variable)

for(v in vars){
  
  cat("\n====================================================\n")
  cat(v, "\n")
  cat("====================================================\n\n")
  
  summary_table <-
    work_data %>%
    group_by(Bayou) %>%
    summarise(
      n = sum(!is.na(.data[[v]])),
      Mean = mean(.data[[v]], na.rm = TRUE),
      Median = median(.data[[v]], na.rm = TRUE),
      SD = sd(.data[[v]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Mean))
  
  print(summary_table)
  
  #Kruskal-Wallis (equivalent to Wilcoxon with two groups)
  
  kw <-
    kruskal_test(
      work_data,
      formula = as.formula(paste(v, "~ Bayou"))
    )
  
  print(kw)
  
  #Wilcoxon rank-sum test
  
  wilcox <-
    wilcox_test(
      work_data,
      formula = as.formula(paste(v, "~ Bayou"))
    )
  
  print(wilcox)
  
}

table(work_data$Bayou, work_data$River)

#Categorical variables: Chi-square or Fisher's exact test

bayou_cat <- lapply(cat_vars, function(v) {
  
  tab <- table(
    work_data$Bayou,
    work_data[[v]]
  )
  
  #Calculate expected counts only to determine which test is appropriate
  chi_check <- suppressWarnings(chisq.test(tab))
  
  if (any(chi_check$expected < 5)) {
    
    test <- fisher.test(tab)
    
    data.frame(
      Variable = v,
      Test = "Fisher's exact",
      Statistic = NA_real_,
      df = NA_real_,
      p = test$p.value
    )
    
  } else {
    
    test <- chisq.test(tab)
    
    data.frame(
      Variable = v,
      Test = "Chi-square",
      Statistic = unname(test$statistic),
      df = unname(test$parameter),
      p = test$p.value
    )
  }
})

bayou_cat <- bind_rows(bayou_cat) %>%
  arrange(p)

bayou_cat

#Which bottom types were concentrated in bayous?
bayou_bottom_table <- table(
  Bayou = work_data$Bayou,
  Bottom = work_data$bottom
)

bayou_bottom_table

#Row percentages: bottom-type composition within each Bayou category
round(prop.table(bayou_bottom_table, margin = 1) * 100, 1)

#Analysis 10 - Does the abundance-weighted center of the 51–100 mm size class
#occurs farther downstream than that of the 0–50 mm size class.

# Calculate as:
# 51–100 mm centroid minus 0–50 mm centroid
#
# Therefore:
# Negative difference = larger juveniles centered farther downstream
# Positive difference = larger juveniles centered farther upstream

library(dplyr)
library(ggplot2)

#1. Function to calculate an abundance-weighted mean

weighted_centroid <- function(distance, abundance) {
  
  keep <- is.finite(distance) &
    is.finite(abundance) &
    abundance >= 0
  
  distance <- distance[keep]
  abundance <- abundance[keep]
  
  # Return NA if no usable observations or no fish were recorded
  if (length(distance) == 0 || sum(abundance) == 0) {
    return(NA_real_)
  }
  
  weighted.mean(
    x = distance,
    w = abundance
  )
}

#2. Retain complete paired shoreline sites

bootstrap_data <- work_data %>%
  select(
    River,
    DstncU_,
    N_0_50_SL,
    N_51_100_SL
  ) %>%
  filter(
    !is.na(River),
    is.finite(DstncU_),
    is.finite(N_0_50_SL),
    is.finite(N_51_100_SL)
  )

# Confirm sample sizes by river
bootstrap_data %>%
  count(River)

#3. Calculate observed abundance-weighted centroids

centroid_0_50 <- weighted_centroid(
  distance = bootstrap_data$DstncU_,
  abundance = bootstrap_data$N_0_50_SL
)

centroid_51_100 <- weighted_centroid(
  distance = bootstrap_data$DstncU_,
  abundance = bootstrap_data$N_51_100_SL
)

observed_difference <- centroid_51_100 - centroid_0_50

observed_results <- data.frame(
  SizeClass = c("0–50 mm SL", "51–100 mm SL"),
  WeightedDistance_m = c(
    centroid_0_50,
    centroid_51_100
  ),
  WeightedDistance_km = c(
    centroid_0_50,
    centroid_51_100
  ) / 1000
)

print(observed_results)

cat(
  "\nObserved difference, 51–100 minus 0–50 mm:",
  round(observed_difference, 1),
  "m\n"
)

if (observed_difference < 0) {
  cat(
    "The abundance-weighted center of the larger size class was",
    round(abs(observed_difference), 1),
    "m farther downstream.\n"
  )
} else {
  cat(
    "The abundance-weighted center of the larger size class was",
    round(observed_difference, 1),
    "m farther upstream.\n"
  )
}

#4. Stratified paired-site bootstrap

#Sites are resampled with replacement within each river.
#Both size-class abundance values remain attached to the same shoreline site.

set.seed(123)

n_boot <- 10000

bootstrap_differences <- replicate(
  n = n_boot,
  expr = {
    
    boot_sample <- bootstrap_data %>%
      group_by(River) %>%
      group_modify(
        ~ .x[
          sample(
            seq_len(nrow(.x)),
            size = nrow(.x),
            replace = TRUE
          ),
        ]
      ) %>%
      ungroup()
    
    weighted_centroid(
      distance = boot_sample$DstncU_,
      abundance = boot_sample$N_51_100_SL
    ) -
      weighted_centroid(
        distance = boot_sample$DstncU_,
        abundance = boot_sample$N_0_50_SL
      )
  }
)

#Remove any unsuccessful bootstrap iterations
bootstrap_differences <- bootstrap_differences[
  is.finite(bootstrap_differences)
]

#5. Summarize bootstrap results

lower_95 <- quantile(
  bootstrap_differences,
  probs = 0.025
)

upper_95 <- quantile(
  bootstrap_differences,
  probs = 0.975
)

proportion_below_zero <- mean(
  bootstrap_differences < 0
)

proportion_above_zero <- mean(
  bootstrap_differences > 0
)

bootstrap_summary <- data.frame(
  Centroid_0_50_m = centroid_0_50,
  Centroid_51_100_m = centroid_51_100,
  ObservedDifference_m = observed_difference,
  BootstrapMeanDifference_m = mean(bootstrap_differences),
  Lower95_m = unname(lower_95),
  Upper95_m = unname(upper_95),
  ProportionBelowZero = proportion_below_zero,
  ProportionAboveZero = proportion_above_zero
)

print(bootstrap_summary)

cat(
  "\n95% bootstrap interval:",
  round(lower_95, 1),
  "to",
  round(upper_95, 1),
  "m\n"
)

cat(
  "Proportion of bootstrap differences below zero:",
  round(proportion_below_zero, 3),
  "\n"
)

#6. Annotated bootstrap plot

plot_annotation <- paste0(
  "Observed difference: ",
  round(observed_difference, 0),
  " m\n",
  "95% bootstrap interval: ",
  round(lower_95, 0),
  " to ",
  round(upper_95, 0),
  " m\n",
  round(proportion_below_zero * 100, 1),
  "% of bootstrap estimates < 0"
)

bootstrap_plot <- ggplot(
  data.frame(Difference = bootstrap_differences),
  aes(x = Difference)
) +
  geom_histogram(
    bins = 50
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = observed_difference,
    linewidth = 1
  ) +
  annotate(
    geom = "label",
    x = observed_difference,
    y = Inf,
    label = plot_annotation,
    vjust = 1.15,
    hjust = ifelse(observed_difference < 0, 0, 1),
    size = 3.5
  ) +
  labs(
    title = "Bootstrap comparison of abundance-weighted distance upstream",
    subtitle = paste0(
      format(length(bootstrap_differences), big.mark = ","),
      " stratified paired-site bootstrap samples"
    ),
    x = "Difference in abundance-weighted distance upstream (m)",
    y = "Bootstrap frequency",
    caption = paste0(
      "Difference calculated as 51–100 mm minus 0–50 mm. ",
      "Negative values indicate that larger juveniles were centered farther downstream. ",
      "Sites were resampled within rivers."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

bootstrap_plot

observed_difference
lower_95
upper_95
proportion_below_zero

#END
