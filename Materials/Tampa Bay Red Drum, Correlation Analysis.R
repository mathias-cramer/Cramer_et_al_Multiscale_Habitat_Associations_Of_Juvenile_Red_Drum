# Project Details______________________________________________________________

# Title: Multiscale Habitat Associations of Juvenile Red Drum (Sciaenops 
#ocellatus) in Tampa Bay, Florida Rivers

# Correlation Analysis_________________________________________________________

# Load packages
library(dplyr)
library(ggplot2)
library(Hmisc)
library(reshape2)

# Import final analysis dataset
data <- read.csv("Tampa_Bay_Red_Drum_Data.csv")

# Select continuous predictor and response variables included in the
# correlation analysis
corr_vars <- c(
  "N_0_50_SL",
  "N_51_100_SL",
  "BankSlp",
  "AvgCnpy",
  "ChnnlDp",
  "MC_IC_P",
  "MC_OM_P",
  "dissolvedO2",
  "DstncU_",
  "LvOystP",
  "Mangrove",
  "OysterArea",
  "RipRap",
  "RvrWdt_",
  "salinity",
  "SaltMarsh",
  "AnySwll",
  "anySdmn",
  "StartDepth",
  "AvgBslL",
  "AvgBslV",
  "AvgBslW",
  "temperature",
  "TidlSlp",
  "AvgVgH_"
)

# Confirm that all required variables are present
missing_vars <- setdiff(corr_vars, names(data))

if (length(missing_vars) > 0) {
  stop(
    "The following required variables are missing from the dataset: ",
    paste(missing_vars, collapse = ", ")
  )
}

# Create correlation dataset
corr_data <- data %>%
  dplyr::select(all_of(corr_vars)) %>%
  mutate(across(everything(), as.numeric))

# Calculate Spearman rank correlations and p-values
rcorr_res <- Hmisc::rcorr(
  as.matrix(corr_data),
  type = "spearman"
)

cor_mat <- rcorr_res$r
p_mat   <- rcorr_res$P

# View full correlation and p-value matrices
print(round(cor_mat, 3))
print(round(p_mat, 3))


# Format data for Figure 2______________________________________________________

# Convert correlation and p-value matrices to long format
cor_df <- reshape2::melt(cor_mat, na.rm = FALSE)
p_df   <- reshape2::melt(p_mat, na.rm = FALSE)

# Variable labels used in Figure 2
# Asterisks identify variables derived from fisheries-independent
# monitoring data collected concurrently with seining.
var_labels <- c(
  N_0_50_SL   = "0–50 mm Red Drum Abundance *",
  N_51_100_SL = "51–100 mm Red Drum Abundance *",
  BankSlp      = "Bank Slope",
  AvgCnpy      = "Canopy Cover",
  ChnnlDp      = "Channel Depth",
  MC_IC_P      = "Channel Inorganic Carbon",
  MC_OM_P      = "Channel Organic Matter",
  dissolvedO2  = "Dissolved Oxygen *",
  DstncU_      = "Distance Upstream",
  LvOystP      = "Live Oysters",
  Mangrove     = "Mangrove Cover",
  OysterArea   = "Oyster Area",
  RipRap       = "Riprap Cover",
  RvrWdt_      = "River Width",
  salinity     = "Salinity *",
  SaltMarsh    = "Salt Marsh Cover",
  AnySwll      = "Seawall Cover",
  anySdmn      = "Sediment Cover",
  StartDepth   = "Seine Depth *",
  AvgBslL      = "Basal Leaf Cover",
  AvgBslV      = "Basal Vegetation Cover",
  AvgBslW      = "Basal Wood Cover",
  temperature  = "Temperature *",
  TidlSlp      = "Tidal Slope",
  AvgVgH_      = "Vegetation Height"
)

# Apply display labels
cor_df$Var1 <- var_labels[as.character(cor_df$Var1)]
cor_df$Var2 <- var_labels[as.character(cor_df$Var2)]

p_df$Var1 <- var_labels[as.character(p_df$Var1)]
p_df$Var2 <- var_labels[as.character(p_df$Var2)]

# Preserve variable ordering
display_order <- unname(var_labels)

cor_df$Var1 <- factor(cor_df$Var1, levels = display_order)
cor_df$Var2 <- factor(cor_df$Var2, levels = display_order)

p_df$Var1 <- factor(p_df$Var1, levels = display_order)
p_df$Var2 <- factor(p_df$Var2, levels = display_order)

# Retain lower triangle, including diagonal
cor_df <- cor_df %>%
  filter(as.numeric(Var1) >= as.numeric(Var2))

p_df <- p_df %>%
  filter(as.numeric(Var1) >= as.numeric(Var2))

# Add p-values to correlation data
cor_df$p_value <- p_df$value

# Bold statistically significant correlations
cor_df$signif <- ifelse(
  !is.na(cor_df$p_value) & cor_df$p_value < 0.05,
  "bold",
  "plain"
)


# Report strong correlations___________________________________________________

# Correlations with |rho| >= 0.70, excluding the diagonal
strong_correlations <- cor_df %>%
  filter(
    Var1 != Var2,
    abs(value) >= 0.70
  ) %>%
  transmute(
    Variable_1 = as.character(Var1),
    Variable_2 = as.character(Var2),
    rho = round(value, 3),
    p = p_value
  ) %>%
  arrange(desc(abs(rho)))

print(strong_correlations)


# Figure 2_____________________________________________________________________

correlation_heatmap <- ggplot(
  cor_df,
  aes(x = Var1, y = Var2, fill = value)
) +
  geom_tile(color = "white") +
  geom_text(
    aes(
      label = sprintf("%.2f", value),
      fontface = signif
    ),
    color = "black",
    size = 4
  ) +
  geom_text(
    data = subset(cor_df, Var1 == Var2),
    aes(
      x = Var1,
      y = Var2,
      label = Var1
    ),
    fontface = "plain",
    nudge_x = -0.5,
    hjust = 1.05,
    vjust = 0.3,
    size = 5
  ) +
  scale_fill_gradient2(
    low = "#0072B2",
    mid = "#F7F7F7",
    high = "#D55E00",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish
  ) +
  scale_x_discrete(
    expand = expansion(add = c(9, 0))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 15,
      color = "black"
    ),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  coord_fixed()

correlation_heatmap

# END