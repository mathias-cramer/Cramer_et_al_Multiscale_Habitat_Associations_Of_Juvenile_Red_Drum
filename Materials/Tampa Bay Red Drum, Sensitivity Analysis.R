#Project Details______________________________________________________________:

#Title: Multiscale Habitat Associations of Juvenile Red Drum (Sciaenops 
#ocellatus) in Tampa Bay, Florida Rivers

#Sensitivity analysis accounting for median distance between seine events and 
#their associated shoreline assessment site

#Data Exploration_____________________________________________________________:

#Load processing packages
library(FSA)
library(ggplot2)
library(tigris)
library(dplyr)
library(lubridate)
library(emmeans)
library(readr)

#Set working directory


#Import data
distance_data <- read_csv("Seine_Transect_Distance_Raw.csv")
headtail(distance_data)

distance_site <- distance_data |>
  group_by(River, SiteID) |>
  summarise(
    Median_Distance_m = median(NEAR_DIST, na.rm = TRUE),
    .groups = "drop"
  )
nrow(distance_site) #should be 51

#Set working directory


#Import data
data <- read.csv("Tampa_Bay_Red_Drum_Data.csv")
headtail(data)

#Calculate median seine-to-transect distance for each shoreline site
distance_site <- distance_data |>
  group_by(River, SiteID) |>
  summarise(
    Median_Distance_m = median(NEAR_DIST, na.rm = TRUE),
    n_seine_associations = n(),
    .groups = "drop"
  )

#Quality-control checks
nrow(distance_site)
summary(distance_site$Median_Distance_m)
summary(distance_site$n_seine_associations)
anyDuplicated(distance_site[c("River", "SiteID")])

#Inspect identifiers in both datasets
sort(unique(data$River))
sort(unique(distance_site$River))

head(data[c("River", "SiteID")])
head(distance_site[c("River", "SiteID")])

#Check for shoreline sites in the main dataset that do not match the distance table
unmatched_sites <- data |>
  distinct(River, SiteID) |>
  anti_join(
    distance_site,
    by = c("River", "SiteID")
  )
unmatched_sites

#Check for shoreline sites in the distance table that do not match the main dataset
unmatched_distance_sites <- distance_site |>
  anti_join(
    data |> distinct(River, SiteID),
    by = c("River", "SiteID")
  )
unmatched_distance_sites

#Merge site-level median distance into the original dataset
data_distance <- data |>
  left_join(
    distance_site,
    by = c("River", "SiteID")
  )

#Confirm successful merge
nrow(data)
nrow(data_distance)
sum(is.na(data_distance$Median_Distance_m))
summary(data_distance$Median_Distance_m)

#Inspect merged identifiers and distance
data_distance |>
  dplyr::select(
    River,
    SiteID,
    Median_Distance_m,
    n_seine_associations
  ) |>
  arrange(River, SiteID)

nrow(data)
nrow(data_distance)
unmatched_sites
unmatched_distance_sites
sum(is.na(data_distance$Median_Distance_m))
summary(data_distance$Median_Distance_m)

#-----------------------------------------------------------------------------#
#0-50 mm SL sensitivity analysis
#-----------------------------------------------------------------------------#

#Recreate modeling dataset
work_data_50 <- data_distance |>
  dplyr::select(
    River,
    RvrWdt_,
    ChnnlDp,
    DstncU_,
    TidlSlp,
    BankSlp,
    AvgBslV,
    AvgBslL,
    AvgBslW,
    AvgInnd,
    AvgCnpy,
    AvgVgH_,
    LvOystP,
    OysterArea,
    RipRap,
    AnySwll,
    anySdmn,
    ShoreType,
    Tr_IC_P,
    Tr_OM_P,
    MC_IC_P,
    MC_OM_P,
    Bayou,
    MidGran,
    TrnsctG,
    StartDepth,
    temperature,
    salinity,
    dissolvedO2,
    bottom,
    VegeSummary,
    Mangrove,
    SaltMarsh,
    N_0_50_SL,
    Median_Distance_m
  )

#Confirm correct data types
work_data_50 <- work_data_50 |>
  mutate(
    across(
      c(
        River,
        AvgInnd,
        ShoreType,
        Bayou,
        MidGran,
        TrnsctG,
        bottom,
        VegeSummary
      ),
      ~ factor(as.character(.x))
    ),
    across(
      -c(
        River,
        AvgInnd,
        ShoreType,
        Bayou,
        MidGran,
        TrnsctG,
        bottom,
        VegeSummary
      ),
      as.numeric
    )
  )

#Remove variables excluded in original analysis
work_data_50 <- work_data_50 |>
  dplyr::select(
    -MidGran,
    -TrnsctG,
    -RipRap,
    -AnySwll,
    -anySdmn,
    -Mangrove,
    -SaltMarsh,
    -DstncU_
  )

#Apply Box-Cox transformation exactly as in original analysis
x <- work_data_50$N_0_50_SL
bc_50 <- car::powerTransform(x + 1)
work_data_50$N_0_50_SL <- car::bcPower(
  x + 1,
  bc_50$lambda
)
bc_50$lambda

#Original final model
U50_original_check <- glm(
  N_0_50_SL ~
    River +
    RvrWdt_ +
    ChnnlDp +
    TidlSlp +
    AvgBslW +
    AvgInnd +
    LvOystP +
    OysterArea +
    ShoreType +
    Tr_OM_P +
    MC_IC_P +
    MC_OM_P +
    StartDepth +
    temperature +
    dissolvedO2 +
    bottom +
    VegeSummary,
  data = work_data_50
)
summary(U50_original_check)

car::Anova(
  U50_original_check,
  type = "II",
  test = "LR"
)
AIC(U50_original_check)
rsq::rsq(
  U50_original_check,
  adj = TRUE
)

#Sensitivity model including median seine-to-transect distance
U50_distance <- glm(
  N_0_50_SL ~
    River +
    RvrWdt_ +
    ChnnlDp +
    TidlSlp +
    AvgBslW +
    AvgInnd +
    LvOystP +
    OysterArea +
    ShoreType +
    Tr_OM_P +
    MC_IC_P +
    MC_OM_P +
    StartDepth +
    temperature +
    dissolvedO2 +
    bottom +
    VegeSummary +
    Median_Distance_m,
  data = work_data_50
)
summary(U50_distance)

car::Anova(
  U50_distance,
  type = "II",
  test = "LR"
)
AIC(U50_distance)
rsq::rsq(
  U50_distance,
  adj = TRUE
)

#Compare models directly
AIC(
  U50_original_check,
  U50_distance
)

anova(
  U50_original_check,
  U50_distance,
  test = "Chisq"
)

summary(U50_distance)$coefficients[
  "Median_Distance_m",
]

bc_50$lambda

summary(U50_distance)

car::Anova(
  U50_distance,
  type = "II",
  test = "LR"
)

rsq::rsq(
  U50_original_check,
  adj = TRUE
)

rsq::rsq(
  U50_distance,
  adj = TRUE
)

summary(U50_original_check)

#-----------------------------------------------------------------------------#
#51-100 mm SL sensitivity analysis
#-----------------------------------------------------------------------------#

#Recreate modeling dataset
work_data_100 <- data_distance |>
  dplyr::select(
    River,
    RvrWdt_,
    ChnnlDp,
    DstncU_,
    TidlSlp,
    BankSlp,
    AvgBslV,
    AvgBslL,
    AvgBslW,
    AvgInnd,
    AvgCnpy,
    AvgVgH_,
    LvOystP,
    OysterArea,
    RipRap,
    AnySwll,
    anySdmn,
    ShoreType,
    Tr_IC_P,
    Tr_OM_P,
    MC_IC_P,
    MC_OM_P,
    Bayou,
    MidGran,
    TrnsctG,
    StartDepth,
    temperature,
    salinity,
    dissolvedO2,
    bottom,
    VegeSummary,
    Mangrove,
    SaltMarsh,
    N_51_100_SL,
    Median_Distance_m
  )

#Confirm correct data types
work_data_100 <- work_data_100 |>
  mutate(
    across(
      c(
        River,
        AvgInnd,
        ShoreType,
        Bayou,
        MidGran,
        TrnsctG,
        bottom,
        VegeSummary
      ),
      ~ factor(as.character(.x))
    ),
    across(
      -c(
        River,
        AvgInnd,
        ShoreType,
        Bayou,
        MidGran,
        TrnsctG,
        bottom,
        VegeSummary
      ),
      as.numeric
    )
  )

#Remove variables excluded in original analysis
work_data_100 <- work_data_100 |>
  dplyr::select(
    -MidGran,
    -TrnsctG,
    -RipRap,
    -AnySwll,
    -anySdmn,
    -Mangrove,
    -SaltMarsh,
    -DstncU_
  )

#Apply Box-Cox transformation exactly as in original analysis
x <- work_data_100$N_51_100_SL
bc_100 <- car::powerTransform(x + 1)
work_data_100$N_51_100_SL <- car::bcPower(
  x + 1,
  bc_100$lambda
)
bc_100$lambda

#Original final model
U100_original_check <- glm(
  N_51_100_SL ~
    ChnnlDp +
    BankSlp +
    AvgInnd +
    AvgVgH_ +
    LvOystP +
    OysterArea +
    MC_IC_P +
    Bayou +
    temperature +
    dissolvedO2 +
    bottom +
    VegeSummary,
  data = work_data_100
)
summary(U100_original_check)

car::Anova(
  U100_original_check,
  type = "II",
  test = "LR"
)

AIC(U100_original_check)

rsq::rsq(
  U100_original_check,
  adj = TRUE
)

#Sensitivity model including median seine-to-transect distance
U100_distance <- glm(
  N_51_100_SL ~
    ChnnlDp +
    BankSlp +
    AvgInnd +
    AvgVgH_ +
    LvOystP +
    OysterArea +
    MC_IC_P +
    Bayou +
    temperature +
    dissolvedO2 +
    bottom +
    VegeSummary +
    Median_Distance_m,
  data = work_data_100
)
summary(U100_distance)

car::Anova(
  U100_distance,
  type = "II",
  test = "LR"
)

AIC(U100_distance)

rsq::rsq(
  U100_distance,
  adj = TRUE
)

#Compare models directly
AIC(
  U100_original_check,
  U100_distance
)

anova(
  U100_original_check,
  U100_distance,
  test = "Chisq"
)

summary(U100_distance)$coefficients[
  "Median_Distance_m",
]

bc_100$lambda

summary(U100_original_check)

summary(U100_distance)

car::Anova(
  U100_distance,
  type = "II",
  test = "LR"
)

AIC(
  U100_original_check,
  U100_distance
)

anova(
  U100_original_check,
  U100_distance,
  test = "Chisq"
)

#Final diagnostic check
car::vif(U50_distance)
car::vif(U100_distance)
rsq::rsq(U100_original_check, adj = TRUE)
rsq::rsq(U100_distance, adj = TRUE)

#END
