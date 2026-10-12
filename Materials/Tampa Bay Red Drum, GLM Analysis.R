#Project Details______________________________________________________________:

#Title: Multiscale Habitat Associations of Juvenile Red Drum (Sciaenops 
#ocellatus) in Tampa Bay, Florida Rivers

#Data Exploration_____________________________________________________________:

#Load processing packages
library(FSA)
library(ggplot2)
library(tigris)
library(dplyr)
library(lubridate)
library(emmeans)

#Set working directory


#Import data
data <- read.csv("Tampa_Bay_Red_Drum_Data.csv")
headtail(data)

#Data Analysis________________________________________________________________:

#Feature selection
names(data)
work_data <- dplyr::select(data,   "River",
                            "RvrWdt_",
                            "ChnnlDp",
                            "DstncU_",
                            "TidlSlp",
                            "BankSlp",
                            "AvgBslV",
                            "AvgBslL",
                            "AvgBslW",
                            "AvgInnd",
                            "AvgCnpy",
                            "AvgVgH_",
                            "LvOystP",
                            "OysterArea",
                            "RipRap",
                            "AnySwll",
                            "anySdmn",
                            "ShoreType",
                            "Tr_IC_P",
                            "Tr_OM_P",
                            "MC_IC_P",
                            "MC_OM_P",
                            "Bayou",
                            "MidGran",
                            "TrnsctG",
                            "StartDepth",
                            "temperature",
                            "salinity",
                            "dissolvedO2",
                            "bottom",
                            "VegeSummary",
                            "Mangrove",
                            "SaltMarsh",
                            "N_0_50_SL",
                            "N_51_100_SL",
                            "N_101_150_SL"
                            )

#Confirm correct data types
work_data <- work_data %>%
  mutate(
    across(c(River, AvgInnd, ShoreType, Bayou, MidGran, TrnsctG, bottom, VegeSummary), ~ factor(as.character(.x))),
    across(-c(River, AvgInnd, ShoreType, Bayou, MidGran, TrnsctG, bottom, VegeSummary), as.numeric)
  )

#Check factor variables for dominance
table(work_data$River)
table(work_data$AvgInnd)
table(work_data$ShoreType)
table(work_data$Bayou)
table(work_data$MidGran)
table(work_data$TrnsctG)
table(work_data$bottom)
table(work_data$VegeSummary)
#MidGran and TrnsctG are highly dominant. Remove them.
work_data <- work_data %>%
  dplyr::select(-MidGran, -TrnsctG)

names(work_data)
headtail(work_data)

#Load analysis packages:
library(MASS)
library(car)

#-----------------------------------------------------------------------------#
#0-50 mm SL Model:
#-----------------------------------------------------------------------------#

#Subset to size class of interest
model_data <- work_data %>%
  dplyr::select(-N_51_100_SL, -N_101_150_SL)
names(model_data)

#Check the distribution
hist(model_data$N_0_50_SL, 
     main = "", 
     col = "lightblue", 
     xlab = "0-50 mm Redfish Count")

#Contains many zero values and is strongly right-skewed. Apply a Box-Cox transformation.
#Add 1 to accommodate zero values, then estimate the Box-Cox transformation
#parameter to reduce skew and stabilize variance.
x <- model_data$N_0_50_SL
bc <- powerTransform(x + 1)
model_data$N_0_50_SL <- bcPower(x + 1, bc$lambda)
bc$lambda

#Check the distribution again
hist(model_data$N_0_50_SL, 
     main = "", 
     col = "lightblue", 
     xlab = "0-50 mm Redfish Count")

#Check for multicollinearity
vif(glm(N_0_50_SL ~ ., data = model_data))
#Identify aliased coefficients and resolve
alias(glm(N_0_50_SL ~ ., data = model_data))
levels(model_data$ShoreType)
levels(model_data$VegeSummary)
model_data <- model_data |>
  dplyr::select(
    -RipRap, #RipRap, AnySwll, and anySdmn aliased. Remove and simply keep ShoreType.
    -AnySwll,
    -anySdmn,
    -Mangrove, #Mangrove and SaltMarsh aliased. Remove and simply keep VegeSummary.
    -SaltMarsh
  )
names(model_data)
#Try again:
vif(glm(N_0_50_SL ~ ., data = model_data))
#Remove if GVIF^(1/(2*Df)) > 5.
model_data <- model_data |> 
  dplyr::select(-DstncU_)
vif(glm(N_0_50_SL ~ ., data = model_data))
#Potential collinearity issues resolved

#Model via step-wise selection
library(MASS)
library(rsq)
full_model <- glm(N_0_50_SL ~ River + RvrWdt_ + ChnnlDp + TidlSlp + BankSlp + AvgBslV +
                    AvgBslL + AvgBslW + AvgInnd + AvgCnpy + AvgVgH_ + LvOystP + OysterArea +
                    ShoreType + Tr_IC_P + Tr_OM_P + MC_IC_P + MC_OM_P + Bayou + StartDepth + 
                    temperature + salinity + dissolvedO2 + bottom + VegeSummary,
                  data = model_data)

step_model <- stepAIC(full_model, direction = "both")
summary(step_model)
rsq(step_model, adj = TRUE)
Anova(step_model, type="II", test="LR")

drop1(step_model, test = "Chisq")
#High Impact: River, ChnnlDp, TidlSlp, ShoreType, bottom  
#Medium Impact: StartDepth, VegeSummary
#Low Impact: AvgBslW, LvOystP, OysterArea, temperature, dissolvedO2

#Compare the selected model with a reduced predictor set
GLM_model <- glm(N_0_50_SL ~  River + RvrWdt_ + ChnnlDp + TidlSlp + AvgBslW + AvgInnd + 
                   LvOystP + OysterArea + ShoreType + Tr_OM_P + MC_IC_P + MC_OM_P + 
                   StartDepth + temperature + dissolvedO2 + bottom + VegeSummary,
                 data = model_data)
GLM_model2 <- glm(N_0_50_SL ~  River + ChnnlDp + TidlSlp + AvgBslW + LvOystP + OysterArea + 
                    ShoreType + StartDepth + temperature + dissolvedO2 + bottom + VegeSummary,
                 data = model_data)
AIC(GLM_model, GLM_model2)
#Proceed with original model

#Return model diagnostics
U50_model <- GLM_model
summary(U50_model)
rsq(GLM_model, adj = TRUE)
Anova(U50_model, type="II", test="LR")
par(mfrow = c(2, 2))
plot(U50_model)
par(mfrow = c(1, 1))
AIC(U50_model)

# Order Table
coefs <- as.data.frame(summary(U50_model)$coefficients)
coefs$Variable <- rownames(coefs)
coefs <- coefs %>% arrange(`Pr(>|t|)`)
coefs

#Post hoc pairwise comparisons
pairs(emmeans(U50_model, "ShoreType"))
pairs(emmeans(U50_model, "bottom"))
pairs(emmeans(U50_model, "River"))
pairs(emmeans(U50_model, "AvgInnd"))

#Data Visualization:
library(MASS)
library(ggplot2)
library(gridExtra)
library(dplyr)

blue_orange_pal <- colorRampPalette(c(
  "#0072B2",  # blue
  "#F7F7F7",  # off-white
  "#D55E00"   # vermillion/orange
))

make_pred_plot <- function(model, varname, var_type = "continuous") {
  
  # Get model frame
  mf <- model.frame(model)
  
  # Set safe base values for factors
  riv_fixed <- levels(mf$River)[1]
  ind_fixed <- levels(mf$AvgInnd)[1]
  sho_fixed <- levels(mf$ShoreType)[1]
  bot_fixed <- levels(mf$bottom)[1]
  veg_fixed <- levels(mf$VegeSummary)[1]
  
  # Create base row of predictors
  base_row <- data.frame(
    River       = factor(riv_fixed, levels = levels(mf$River)),
    RvrWdt_     = mean(mf$RvrWdt_, na.rm = TRUE),
    ChnnlDp     = mean(mf$ChnnlDp, na.rm = TRUE),
    TidlSlp     = mean(mf$TidlSlp, na.rm = TRUE),
    AvgBslW     = mean(mf$AvgBslW, na.rm = TRUE),
    AvgInnd     = factor(ind_fixed, levels = levels(mf$AvgInnd)),
    LvOystP     = mean(mf$LvOystP, na.rm = TRUE),
    OysterArea  = mean(mf$OysterArea, na.rm = TRUE),
    ShoreType   = factor(sho_fixed, levels = levels(mf$ShoreType)),
    Tr_OM_P     = mean(mf$Tr_OM_P, na.rm = TRUE),
    MC_IC_P     = mean(mf$MC_IC_P, na.rm = TRUE),
    MC_OM_P     = mean(mf$MC_OM_P, na.rm = TRUE),
    StartDepth  = mean(mf$StartDepth, na.rm = TRUE),
    temperature = mean(mf$temperature, na.rm = TRUE),
    dissolvedO2 = mean(mf$dissolvedO2, na.rm = TRUE),
    bottom      = factor(bot_fixed, levels = levels(mf$bottom)),
    VegeSummary = factor(veg_fixed, levels = levels(mf$VegeSummary))
  )
  
  if (var_type == "continuous") {
    vals <- mf[[varname]]
    vals <- vals[!is.na(vals)]
    pred_df <- base_row[rep(1, 100), ]
    pred_df[[varname]] <- seq(from = min(vals), to = max(vals), length.out = 100)
  } else if (var_type == "factor") {
    # Use all levels of the factor variable
    levels_vec <- levels(mf[[varname]])
    pred_df <- base_row[rep(1, length(levels_vec)), ]
    pred_df[[varname]] <- factor(levels_vec, levels = levels_vec)
    
    #If factor, override base factor values
    if (varname %in% c("River", "AvgInnd", "ShoreType", "bottom", "VegeSummary")) {
      pred_df[[varname]] <- factor(levels_vec, levels = levels(mf[[varname]]))
    }
  }
  
  pred <- predict(model, newdata = pred_df, type = "response", se.fit = TRUE)
  pred_df$predicted <- pred$fit
  pred_df$conf.low  <- pred$fit - 1.96 * pred$se.fit
  pred_df$conf.high <- pred$fit + 1.96 * pred$se.fit
  
  # Publication-ready labels
  pretty_labels <- c(
    "River"       = "River",
    "RvrWdt_"     = "River width (m)",
    "ChnnlDp"     = "Channel depth (m)",
    "TidlSlp"     = "Tidal slope",
    "AvgBslW"     = "Basal wood cover",
    "AvgInnd"     = "Inundation (m)",
    "LvOystP"     = "Live oysters (%)",
    "OysterArea"  = "Oyster area",
    "ShoreType"   = "Shore type",
    "Tr_OM_P"     = "Shoreline organic matter",
    "MC_IC_P"     = "Channel inorganic carbon",
    "MC_OM_P"     = "Channel organic matter",
    "StartDepth"  = "Seine depth (m)",
    "temperature" = "Temperature",
    "dissolvedO2" = "Dissolved oxygen",
    "bottom"      = "Bottom type",
    "VegeSummary" = "Vegetation summary"
  )
  
  # Plot
  library(ggplot2)
  if (var_type == "continuous") {
    p <- ggplot(pred_df, aes(x = .data[[varname]], y = predicted)) +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.5,
        color = "grey50"
      ) +
      geom_ribbon(
        aes(ymin = conf.low, ymax = conf.high),
        fill = "#a4cae0",
        alpha = 0.3
      ) +
      geom_line(
        color = "#d55e00",
        linewidth = 1.2
      ) +
      labs(
        x = pretty_labels[varname],
        y = "Predicted Red Drum abundance (Box-Cox)"
      ) +
      theme_classic(base_size = 16)
  } else if (var_type == "factor") {
    p <- ggplot(
      pred_df,
      aes(
        x = .data[[varname]],
        y = predicted,
        fill = .data[[varname]]
      )
    ) +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.5,
        color = "grey50"
      ) +
      geom_errorbar(
        aes(
          ymin = conf.low,
          ymax = conf.high
        ),
        width = 0.15,
        linewidth = 0.8
      ) +
      geom_point(
        shape = 21,
        size = 3.5,
        color = "black",
        stroke = 0.7
      ) +
      labs(
        x = pretty_labels[varname],
        y = "Predicted Red Drum abundance (Box-Cox)"
      ) +
      theme_classic(base_size = 16) +
      theme(
        legend.position = "none"
      )
    
    if (varname == "River") {
      p <- p + scale_x_discrete(labels = c(
        AL = "Alafia",
        LMN = "Little Manatee",
        MN = "Manatee"
      )) +
        scale_fill_manual(values = setNames(
          blue_orange_pal(length(levels_vec)),
          levels_vec
        ))
    } else if (varname == "AvgInnd") {
      p <- p + scale_x_discrete(labels = c(
        "0" = "0",
        "1" = "1",
        "1.5" = "1–2",
        "2" = ">2"
      )) +
        scale_fill_manual(values = setNames(
          blue_orange_pal(length(levels_vec)),
          levels_vec
        ))
    } else if (varname == "ShoreType") {
      p <- p + scale_x_discrete(labels = c(
        "rip rap" = "Rip rap",
        "seawall" = "Seawall",
        "sediment" = "Sediment"
      )) +
        scale_fill_manual(values = setNames(
          blue_orange_pal(length(levels_vec)),
          levels_vec
        ))
    } else if (varname == "bottom") {
      p <- p + scale_x_discrete(labels = c(
        "Mud" = "Mud",
        "MudSan" = "Mud-sand",
        "MudSanStr" = "Mud-sand-structure",
        "MudStr" = "Mud-structure",
        "San" = "Sand",
        "SanStr" = "Sand-structure"
      )) 
    } else if (varname == "VegeSummary") {
      veg_levels <- c(
        "none",
        "freshwater marsh",
        "mangroves",
        "salt marsh",
        "upland vegetation"
      )
      p <- p + scale_x_discrete(
        limits = veg_levels,
        labels = c(
          "none" = "None",
          "freshwater marsh" = "Fresh marsh",
          "mangroves" = "Mangrove",
          "salt marsh" = "Salt marsh",
          "upland vegetation" = "Upland"
        )
      ) +
        scale_fill_manual(values = setNames(
          blue_orange_pal(length(veg_levels)),
          veg_levels
        ))
    }
  } 
  
  return(p)
}

#Plots
u50.1  <- make_pred_plot(U50_model, "River", var_type = "factor")
u50.2  <- make_pred_plot(U50_model, "RvrWdt_", var_type = "continuous")
u50.3  <- make_pred_plot(U50_model, "ChnnlDp", var_type = "continuous")
u50.4  <- make_pred_plot(U50_model, "TidlSlp", var_type = "continuous")
u50.5  <- make_pred_plot(U50_model, "AvgBslW", var_type = "continuous")
u50.6  <- make_pred_plot(U50_model, "AvgInnd", var_type = "factor")
u50.7  <- make_pred_plot(U50_model, "LvOystP", var_type = "continuous")
u50.8  <- make_pred_plot(U50_model, "OysterArea", var_type = "continuous")
u50.9  <- make_pred_plot(U50_model, "ShoreType", var_type = "factor")
u50.10  <- make_pred_plot(U50_model, "Tr_OM_P", var_type = "continuous")
u50.11  <- make_pred_plot(U50_model, "MC_IC_P", var_type = "continuous")
u50.12  <- make_pred_plot(U50_model, "MC_OM_P", var_type = "continuous")
u50.13  <- make_pred_plot(U50_model, "StartDepth", var_type = "continuous")
u50.14  <- make_pred_plot(U50_model, "temperature", var_type = "continuous")
u50.15  <- make_pred_plot(U50_model, "dissolvedO2", var_type = "continuous")
u50.16  <- make_pred_plot(U50_model, "bottom", var_type = "factor")
u50.17  <- make_pred_plot(U50_model, "VegeSummary", var_type = "factor")

#PvO Plot
# Extract the model frame used for fitting
model_data <- model.frame(U50_model)

# Predicted values on the response scale
model_data$predicted_CPUE <- predict(U50_model, type = "response")

# Observed vs predicted plot
pvo <- ggplot(model_data, aes(x = N_0_50_SL, y = predicted_CPUE)) +
  geom_point(alpha = 0.6, size = 3, color = "#a4cae0") +
  geom_smooth(method = "lm", color = "#d55e00", se = FALSE, linewidth = 1) +
  labs(x = "Observed Red Drum abundance (Box-Cox)",
       y = "Perdicted Red Drum abundance (Box-Cox)") +
  theme_classic(base_size = 16)

#Altogether:
grid.arrange(pvo, u50.1, u50.2, u50.3, u50.4, u50.5, u50.6, u50.7, u50.8, u50.9, u50.10,
             u50.11, u50.12, u50.13, u50.14, u50.15, u50.16, u50.17,
             ncol = 2, nrow = 9)

#Create a shared y-axis label
pvo <- pvo + ylab(NULL)
u50.1 <- u50.1 + ylab(NULL)
u50.2 <- u50.2 + ylab(NULL)
u50.3 <- u50.3 + ylab(NULL)
u50.4 <- u50.4 + ylab(NULL)
u50.5 <- u50.5 + ylab(NULL)
u50.6 <- u50.6 + ylab(NULL)
u50.7 <- u50.7 + ylab(NULL)
u50.8 <- u50.8 + ylab(NULL)
u50.9 <- u50.9 + ylab(NULL)
u50.10 <- u50.10 + ylab(NULL)
u50.11 <- u50.11 + ylab(NULL)
u50.12 <- u50.12 + ylab(NULL)
u50.13 <- u50.13 + ylab(NULL)
u50.14 <- u50.14 + ylab(NULL)
u50.15 <- u50.15 + ylab(NULL)
u50.16 <- u50.16 + ylab(NULL)
u50.17 <- u50.17 + ylab(NULL)

library(grid)
u50.visual <- grid.arrange(
  pvo, u50.1, u50.2, u50.3, u50.4, u50.5, u50.6, u50.7, u50.8, u50.9, u50.10,
  u50.11, u50.12, u50.13, u50.14, u50.15, u50.16, u50.17,
  ncol = 3,
  left = textGrob("Predicted Red Drum abundance (Box-Cox)", rot = 90, gp = gpar(fontsize = 16))
)
#3000x2500

pvo <- pvo +
  labs(tag = "(a)") +
  theme(
    plot.tag = element_text(face = "bold", size = 14)
  )

u50.6 <- u50.6 +
  labs(tag = "(b)") +
  theme(
    plot.tag = element_text(face = "bold", size = 14)
  )

u50.4 <- u50.4 +
  labs(tag = "(c)") +
  theme(
    plot.tag = element_text(face = "bold", size = 14)
  )

u50.1 <- u50.1 +
  labs(tag = "(d)") +
  theme(
    plot.tag = element_text(face = "bold", size = 14)
  )

u50.3 <- u50.3 +
  labs(tag = "(e)") +
  theme(
    plot.tag = element_text(face = "bold", size = 14)
  )

u50.13 <- u50.13 +
  labs(tag = "(f)") +
  theme(
    plot.tag = element_text(face = "bold", size = 14)
  )

library(cowplot)

aligned_plots <- align_plots(
  pvo,
  u50.6,
  u50.4,
  u50.1,
  u50.3,
  u50.13,
  align = "hv",
  axis = "tblr"
)

final_0_figure <- grid.arrange(
  grobs = aligned_plots,
  ncol = 2,
  left = textGrob(
    "Predicted Red Drum abundance (Box-Cox)",
    rot = 90,
    gp = gpar(fontsize = 16)
  )
)

#-----------------------------------------------------------------------------#
#51-100 mm SL Model:
#-----------------------------------------------------------------------------#

#Subset to size class of interest
model_data <- work_data %>%
  dplyr::select(-N_0_50_SL, -N_101_150_SL)
names(model_data)

#Check the distribution
hist(model_data$N_51_100_SL, 
     main = "", 
     col = "lightblue", 
     xlab = "51-100 mm Redfish Count")

#Contains many zero values and is strongly right-skewed. Apply a Box-Cox transformation.
#Add 1 to accommodate zero values, then estimate the Box-Cox transformation
#parameter to reduce skew and stabilize variance.
x <- model_data$N_51_100_SL
bc <- powerTransform(x + 1)
model_data$N_51_100_SL <- bcPower(x + 1, bc$lambda)

#Check the distribution again
hist(model_data$N_51_100_SL, 
     main = "", 
     col = "lightblue", 
     xlab = "51-100 mm Redfish Count")

#Model via step-wise selection
library(MASS)
full_model <- glm(N_51_100_SL ~ River + RvrWdt_ + ChnnlDp + TidlSlp + BankSlp + AvgBslV +
                    AvgBslL + AvgBslW + AvgInnd + AvgCnpy + AvgVgH_ + LvOystP + OysterArea +
                    ShoreType + Tr_IC_P + Tr_OM_P + MC_IC_P + MC_OM_P + Bayou + StartDepth + 
                    temperature + salinity + dissolvedO2 + bottom + VegeSummary,
                  data = model_data)

step_model <- stepAIC(full_model, direction = "both")
summary(step_model)
rsq(step_model, adj = TRUE)
Anova(step_model, type="II", test="LR")

drop1(step_model, test = "Chisq")
#High Impact: AvgInnd, Bayou, VegeSummary  
#Medium Impact: AvgVgH_
#Low Impact: ChnnlDp, BankSlp, LvOystP

#Compare the selected model with a reduced predictor set
GLM_model <- glm(N_51_100_SL ~  ChnnlDp + BankSlp + AvgInnd + AvgVgH_ + LvOystP + 
                   OysterArea + MC_IC_P + Bayou + temperature + dissolvedO2 + 
                   bottom + VegeSummary,
                 data = model_data)
GLM_model2 <- glm(N_51_100_SL ~  ChnnlDp + BankSlp + AvgInnd + AvgVgH_ + LvOystP + 
                    Bayou + VegeSummary,
                  data = model_data)
AIC(GLM_model, GLM_model2)
#Proceed with original model

#Return model diagnostics
U100_model <- GLM_model
summary(U100_model)
rsq(U100_model, adj = TRUE)
Anova(U100_model, type="II", test="LR")
par(mfrow = c(2, 2))
plot(U100_model)
par(mfrow = c(1, 1))
AIC(U100_model)

# Order Table
coefs <- as.data.frame(summary(U100_model)$coefficients)
coefs$Variable <- rownames(coefs)
coefs <- coefs %>% arrange(`Pr(>|t|)`)
coefs

# Post hoc pairwise comparisons
pairs(emmeans(U100_model, "AvgInnd"))
pairs(emmeans(U100_model, "bottom"))
pairs(emmeans(U100_model, "VegeSummary"))

#Data Visualization:
library(MASS)
library(ggplot2)
library(gridExtra)
library(dplyr)

make_pred_plot <- function(model, varname, var_type = "continuous") {
  
  # Get model frame
  mf <- model.frame(model)
  
  # Set safe base values for factors
  ind_fixed <- levels(mf$AvgInnd)[1]
  bay_fixed <- levels(mf$Bayou)[1]
  bot_fixed <- levels(mf$bottom)[1]
  veg_fixed <- levels(mf$VegeSummary)[1]
  
  # Create base row of predictors
  base_row <- data.frame(
    ChnnlDp     = mean(mf$ChnnlDp, na.rm = TRUE),
    BankSlp     = mean(mf$BankSlp, na.rm = TRUE),
    AvgInnd     = factor(ind_fixed, levels = levels(mf$AvgInnd)),
    AvgVgH_     = mean(mf$AvgVgH_, na.rm = TRUE),
    LvOystP     = mean(mf$LvOystP, na.rm = TRUE),
    OysterArea  = mean(mf$OysterArea, na.rm = TRUE),
    MC_IC_P     = mean(mf$MC_IC_P, na.rm = TRUE),
    Bayou       = factor(bay_fixed, levels = levels(mf$Bayou)),
    temperature = mean(mf$temperature, na.rm = TRUE),
    dissolvedO2 = mean(mf$dissolvedO2, na.rm = TRUE),
    bottom      = factor(bot_fixed, levels = levels(mf$bottom)),
    VegeSummary = factor(veg_fixed, levels = levels(mf$VegeSummary))
  )
  
  if (var_type == "continuous") {
    vals <- mf[[varname]]
    vals <- vals[!is.na(vals)]
    pred_df <- base_row[rep(1, 100), ]
    pred_df[[varname]] <- seq(from = min(vals), to = max(vals), length.out = 100)
  } else if (var_type == "factor") {
    # Use all levels of the factor variable
    levels_vec <- levels(mf[[varname]])
    pred_df <- base_row[rep(1, length(levels_vec)), ]
    pred_df[[varname]] <- factor(levels_vec, levels = levels_vec)
    
    #If factor, override base factor values
    if (varname %in% c("AvgInnd", "Bayou", "bottom", "VegeSummary")) {
      pred_df[[varname]] <- factor(levels_vec, levels = levels(mf[[varname]]))
    }
  }
  
  pred <- predict(model, newdata = pred_df, type = "response", se.fit = TRUE)
  pred_df$predicted <- pred$fit
  pred_df$conf.low  <- pred$fit - 1.96 * pred$se.fit
  pred_df$conf.high <- pred$fit + 1.96 * pred$se.fit
  
  # Publication-ready labels
  pretty_labels <- c(
    "ChnnlDp"     = "Channel depth (m)",
    "BankSlp"     = "Bank slope (rise:run)",
    "AvgInnd"     = "Inundation (m)",
    "AvgVgH_"     = "Vegetation height (m)",
    "LvOystP"     = "Live oysters (%)",
    "OysterArea"  = "Oyster area (m2)",
    "MC_IC_P"     = "Channel inorganic carbon (%)",
    "Bayou"       = "In embayment?",
    "temperature" = "Temperature (°C)",
    "dissolvedO2" = "Dissolved oxygen (mg/L)",
    "bottom"      = "Seine bottom type",
    "VegeSummary" = "Vegetation type"
  )
  
  # Plot
  library(ggplot2)
  if (var_type == "continuous") {
    
    p <- ggplot(pred_df, aes(x = .data[[varname]], y = predicted)) +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.5,
        color = "grey50"
      ) +
      geom_ribbon(
        aes(ymin = conf.low, ymax = conf.high),
        fill = "#a4cae0",
        alpha = 0.3
      ) +
      geom_line(
        color = "#d55e00",
        linewidth = 1.2
      ) +
      labs(
        x = pretty_labels[varname],
        y = "Predicted Red Drum abundance (Box-Cox)"
      ) +
      theme_classic(base_size = 16)
  } else if (var_type == "factor") {
    p <- ggplot(
      pred_df,
      aes(
        x = .data[[varname]],
        y = predicted,
        fill = .data[[varname]]
      )
    ) +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.5,
        color = "grey50"
      ) +
      geom_errorbar(
        aes(
          ymin = conf.low,
          ymax = conf.high
        ),
        width = 0.15,
        linewidth = 0.8
      ) +
      geom_point(
        shape = 21,
        size = 3.5,
        color = "black",
        stroke = 0.7
      ) +
      labs(
        x = pretty_labels[varname],
        y = "Predicted Red Drum abundance (Box-Cox)"
      ) +
      theme_classic(base_size = 16) +
      theme(
        legend.position = "none"
      )
    
    if (varname == "AvgInnd") {
      p <- p + scale_x_discrete(labels = c(
        "0" = "0",
        "1" = "1",
        "1.5" = "1–2",
        "2" = ">2"
      ))+
        scale_fill_manual(values = setNames(
          blue_orange_pal(length(levels_vec)),
          levels_vec
        ))
    } else if (varname == "Bayou") {
      p <- p + scale_x_discrete(labels = c(
        "no" = "No",
        "yes" = "Yes"
      )) +
        scale_fill_manual(values = setNames(
          blue_orange_pal(length(levels_vec)),
          levels_vec
        ))
    } else if (varname == "VegeSummary") {
      veg_levels <- c(
        "none",
        "freshwater marsh",
        "mangroves",
        "salt marsh",
        "upland vegetation"
      )
      p <- p + scale_x_discrete(
        limits = veg_levels,
        labels = c(
          "none" = "None",
          "freshwater marsh" = "Fresh marsh",
          "mangroves" = "Mangrove",
          "salt marsh" = "Salt marsh",
          "upland vegetation" = "Upland"
        )
      ) +
        scale_fill_manual(values = setNames(
          blue_orange_pal(length(veg_levels)),
          veg_levels
        ))
    }
  } 
  
  return(p)
}

#Plots
u100.1  <- make_pred_plot(U100_model, "ChnnlDp", var_type = "continuous")
u100.2  <- make_pred_plot(U100_model, "BankSlp", var_type = "continuous")
u100.3  <- make_pred_plot(U100_model, "AvgInnd", var_type = "factor")
u100.4  <- make_pred_plot(U100_model, "AvgVgH_", var_type = "continuous")
u100.5  <- make_pred_plot(U100_model, "LvOystP", var_type = "continuous")
u100.6  <- make_pred_plot(U100_model, "OysterArea", var_type = "continuous")
u100.7  <- make_pred_plot(U100_model, "MC_IC_P", var_type = "continuous")
u100.8  <- make_pred_plot(U100_model, "Bayou", var_type = "factor")
u100.9  <- make_pred_plot(U100_model, "temperature", var_type = "continuous")
u100.10  <- make_pred_plot(U100_model, "dissolvedO2", var_type = "continuous")
u100.11  <- make_pred_plot(U100_model, "bottom", var_type = "factor")
u100.12  <- make_pred_plot(U100_model, "VegeSummary", var_type = "factor")

#PvO Plot
# Extract the model frame used for fitting
model_data <- model.frame(U100_model)

# Predicted values on the response scale
model_data$predicted_CPUE <- predict(U100_model, type = "response")

# Observed vs predicted plot
pvo <- ggplot(model_data, aes(x = N_51_100_SL, y = predicted_CPUE)) +
  geom_point(alpha = 0.6, size = 3, color = "#a4cae0") +
  geom_smooth(method = "lm", color = "#d55e00", se = FALSE, linewidth = 1) +
  labs(x = "Observed Red Drum abundance (Box-Cox)",
       y = "Perdicted Red Drum abundance (Box-Cox)") +
  theme_classic(base_size = 16)

#Altogether:
grid.arrange(pvo, u100.1, u100.2, u100.3, u100.4, u100.5, u100.6, u100.7, u100.8, 
             u100.9, u100.10, u100.11, u100.12,
             ncol = 2, nrow = 7)

#Create a shared y-axis label
pvo <- pvo + ylab(NULL)
u100.1 <- u100.1 + ylab(NULL)
u100.2 <- u100.2 + ylab(NULL)
u100.3 <- u100.3 + ylab(NULL)
u100.4 <- u100.4 + ylab(NULL)
u100.5 <- u100.5 + ylab(NULL)
u100.6 <- u100.6 + ylab(NULL)
u100.7 <- u100.7 + ylab(NULL)
u100.8 <- u100.8 + ylab(NULL)
u100.9 <- u100.9 + ylab(NULL)
u100.10 <- u100.10 + ylab(NULL)
u100.11 <- u100.11 + ylab(NULL)
u100.12 <- u100.12 + ylab(NULL)

library(grid)
u100.visual <- grid.arrange(
  nullGrob(), pvo, nullGrob(), u100.1, u100.2, u100.3, u100.4, u100.5, u100.6, u100.7, u100.8, 
  u100.9, u100.10, u100.11, u100.12,
  ncol = 3,
  left = textGrob("Predicted Red Drum abundance (Box-Cox)", rot = 90, gp = gpar(fontsize = 16))
)

#3000x2083.33

library(ggplot2)

pvo <- pvo +
  labs(tag = "(a)") +
  theme(
    plot.tag = element_text(face = "bold", size = 14)
  )

u100.3 <- u100.3 +
  labs(tag = "(b)") +
  theme(
    plot.tag = element_text(face = "bold", size = 14)
  )

u100.12 <- u100.12 +
  labs(tag = "(c)") +
  theme(
    plot.tag = element_text(face = "bold", size = 14)
  )

u100.8 <- u100.8 +
  labs(tag = "(d)") +
  theme(
    plot.tag = element_text(face = "bold", size = 14)
  )

u100.9 <- u100.9 +
  labs(tag = "(e)") +
  theme(
    plot.tag = element_text(face = "bold", size = 14)
  )

u100.4 <- u100.4 +
  labs(tag = "(f)") +
  theme(
    plot.tag = element_text(face = "bold", size = 14)
  )

library(cowplot)

aligned_plots_51 <- align_plots(
  pvo,
  u100.3,
  u100.12,
  u100.8,
  u100.9,
  u100.4,
  align = "hv",
  axis = "tblr"
)

library(grid)
library(gridExtra)

final_51_figure <- grid.arrange(
  grobs = aligned_plots_51,
  ncol = 2,
  left = textGrob(
    "Predicted Red Drum abundance (Box-Cox)",
    rot = 90,
    gp = gpar(fontsize = 16)
  )
)

#Final figures:
grid.draw(final_0_figure)
grid.draw(final_51_figure)
#Poster: 1250x750

#END
