#Project Details______________________________________________________________:

#Title: Multiscale Habitat Associations of Juvenile Red Drum (Sciaenops 
#ocellatus) in Tampa Bay, Florida Rivers

#Mangrove habitat follow-up analysis___________________________________________:

library(dplyr)
library(ggplot2)

#1. Create original High / Low abundance groups

#0–50 mm SL
q25_0_50 <- quantile(work_data$N_0_50_SL, 0.25, na.rm = TRUE)
q75_0_50 <- quantile(work_data$N_0_50_SL, 0.75, na.rm = TRUE)

#51–100 mm SL
q25_51_100 <- quantile(work_data$N_51_100_SL, 0.25, na.rm = TRUE)
q75_51_100 <- quantile(work_data$N_51_100_SL, 0.75, na.rm = TRUE)

work_data <- work_data %>%
  mutate(
    Group_0_50 = case_when(
      N_0_50_SL <= q25_0_50 ~ "Low",
      N_0_50_SL >= q75_0_50 ~ "High",
      TRUE ~ NA_character_
    ),
    Group_51_100 = case_when(
      N_51_100_SL <= q25_51_100 ~ "Low",
      N_51_100_SL >= q75_51_100 ~ "High",
      TRUE ~ NA_character_
    )
  )

#2. Subset to mangrove-dominated transects

mangrove_data <- work_data %>%
  filter(VegeSummary == "mangroves")

#Check total number of mangrove-dominated sites
nrow(mangrove_data)

#Check High / Low sample sizes for each size class
table(mangrove_data$Group_0_50, useNA = "ifany")
table(mangrove_data$Group_51_100, useNA = "ifany")

table(mangrove_data$River, mangrove_data$Group_0_50, useNA = "ifany")
table(mangrove_data$River, mangrove_data$Group_51_100, useNA = "ifany")

table(mangrove_data$AvgInnd, mangrove_data$Group_0_50, useNA = "ifany")
table(mangrove_data$AvgInnd, mangrove_data$Group_51_100, useNA = "ifany")

#3. Define habitat variables to examine

#Continuous variables
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

#Categorical variables
categorical_vars <- c(
  "River",
  "Bayou",
  "ShoreType",
  "bottom",
  "AvgInnd"
)

#PART A — HIGH VS LOW MANGROVE SITES___________________________________________:

#4. Function for continuous High vs Low comparisons

compare_continuous <- function(data, group_var, variables) {
  
  results <- lapply(variables, function(v) {
    
    temp <- data %>%
      select(all_of(c(group_var, v))) %>%
      filter(
        !is.na(.data[[group_var]]),
        !is.na(.data[[v]])
      )
    
    #Ensure both groups are represented
    if (length(unique(temp[[group_var]])) < 2) {
      return(NULL)
    }
    
    high <- temp[[v]][temp[[group_var]] == "High"]
    low  <- temp[[v]][temp[[group_var]] == "Low"]
    
    test <- wilcox.test(
      high,
      low,
      exact = FALSE
    )
    
    data.frame(
      Variable = v,
      N_High = length(high),
      N_Low = length(low),
      High_Median = median(high, na.rm = TRUE),
      Low_Median = median(low, na.rm = TRUE),
      High_Mean = mean(high, na.rm = TRUE),
      Low_Mean = mean(low, na.rm = TRUE),
      Difference_Median = median(high, na.rm = TRUE) -
        median(low, na.rm = TRUE),
      W = unname(test$statistic),
      p.value = test$p.value
    )
  })
  
  results <- bind_rows(results)
  
  #False-discovery-rate correction
  results$p.adjusted <- p.adjust(
    results$p.value,
    method = "BH"
  )
  
  results %>%
    arrange(p.value)
}

#5. Continuous-variable comparisons: 0–50 mm

results_0_50 <- compare_continuous(
  mangrove_data,
  "Group_0_50",
  continuous_vars
)

results_0_50

#6. Continuous-variable comparisons: 51–100 mm

results_51_100 <- compare_continuous(
  mangrove_data,
  "Group_51_100",
  continuous_vars
)

results_51_100

#7. Function for categorical High vs Low comparisons

set.seed(123)

compare_categorical <- function(data, group_var, variables) {
  
  results <- lapply(variables, function(v) {
    
    temp <- data %>%
      select(all_of(c(group_var, v))) %>%
      filter(
        !is.na(.data[[group_var]]),
        !is.na(.data[[v]])
      )
    
    tab <- table(
      temp[[group_var]],
      temp[[v]]
    )
    
    #Need at least 2 rows and 2 columns
    if (nrow(tab) < 2 || ncol(tab) < 2) {
      return(NULL)
    }
    
    test <- fisher.test(
      tab,
      simulate.p.value = TRUE,
      B = 9999
    )
    
    data.frame(
      Variable = v,
      N = sum(tab),
      p.value = test$p.value
    )
  })
  
  results <- bind_rows(results)
  
  results$p.adjusted <- p.adjust(
    results$p.value,
    method = "BH"
  )
  
  results %>%
    arrange(p.value)
}

#8. Categorical comparisons: 0–50 mm

categorical_0_50 <- compare_categorical(
  mangrove_data,
  "Group_0_50",
  categorical_vars
)

categorical_0_50

#9. Categorical comparisons: 51–100 mm

categorical_51_100 <- compare_categorical(
  mangrove_data,
  "Group_51_100",
  categorical_vars
)

categorical_51_100

#PART B — CONTINUOUS ABUNDANCE RELATIONSHIPS___________________________________:

#10. Spearman relationships across ALL mangrove-dominated sites

spearman_mangrove <- function(data, abundance_var, variables) {
  
  results <- lapply(variables, function(v) {
    
    temp <- data %>%
      select(all_of(c(abundance_var, v))) %>%
      filter(
        !is.na(.data[[abundance_var]]),
        !is.na(.data[[v]])
      )
    
    if (nrow(temp) < 3) {
      return(NULL)
    }
    
    test <- suppressWarnings(
      cor.test(
        temp[[abundance_var]],
        temp[[v]],
        method = "spearman",
        exact = FALSE
      )
    )
    
    data.frame(
      Variable = v,
      N = nrow(temp),
      rho = unname(test$estimate),
      p.value = test$p.value
    )
  })
  
  results <- bind_rows(results)
  
  results$p.adjusted <- p.adjust(
    results$p.value,
    method = "BH"
  )
  
  results %>%
    arrange(p.value)
}

#11. Spearman analysis: 0–50 mm

spearman_0_50 <- spearman_mangrove(
  mangrove_data,
  "N_0_50_SL",
  continuous_vars
)

spearman_0_50

#12. Spearman analysis: 51–100 mm

spearman_51_100 <- spearman_mangrove(
  mangrove_data,
  "N_51_100_SL",
  continuous_vars
)

spearman_51_100

#END
