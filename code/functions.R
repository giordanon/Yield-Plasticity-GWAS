

library(tidyverse)
library(brms)
library(rstan)
library(mgcv)
library(tidybayes)
library(nlme)
library(lme4)
library(emmeans)
library(bayesQR)
library(FW)
library(extrafont)
library(nlraa)
rstan_options(auto_write = TRUE, silent = TRUE)
options(mc.cores = parallel::detectCores())
# Colors from MetBrewer
clrs <- MetBrewer::met.brewer("Java")

loadfonts(device = c("all", "pdf", "postscript", "win"))
`%nin%` <- Negate(`%in%`)
# Custom ggplot theme to make pretty plots
# Get the font at https://fonts.google.com/specimen/Jost
theme_nice <- function() {
  theme_minimal(base_family = "Jost") +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(family = "Jost", face = "bold"),
          axis.title = element_text(family = "Jost", size = 14),
          strip.text = element_text(family = "Jost", face = "bold",
                                    size = rel(1.2), hjust = 0),
          legend.position = "top",
          axis.line = element_line(),
          plot.background = element_rect(fill = "white", color = "white"),
          axis.text = element_text(size = 14),
          strip.background = element_rect(fill = "grey80", color = NA)
    )
}

azul = "#2c7c94"
rojo = "#a65852"

# Simulate Zaddocks growth stage

# Lollato Agronomy Journal 2020 paper

getGS <- function(acumGDD){
  set.seed(123)
  GDDZadocks <- function(acumGDD, a0, a1, a2, a3){
    a0 + a1 / (1+exp(-(acumGDD-a2)/a3))
  }
  
  n <- 1
  a0 <- rnorm(n, 24.74,2.8)
  a1 <- rnorm(n, 60, 3.9)
  a2 <- rnorm(n, 964.8, 20.4)
  a3 <- rnorm(n, 205.7, 20.2)
  GDDZadocks(acumGDD, a0, a1, a2, a3)
  
}


# This function returns the anthesis date given Lollato 2020 model and daymet data
anthesisDate <- function(data, tminCP){
  
  temp1 <- 
    data %>% 
    unnest(cols = c(daymet)) %>% 
    mutate(Tmean = (TEMP2MMIN + TEMP2MMAX)/2, 
           #Critical period growing degree units
           gmin = case_when(TEMP2MMIN >= tminCP ~ TEMP2MMIN, 
                            TRUE ~ tminCP),
           # Tmax threshold Growing Degrees.
           gmax = case_when(
             TEMP2MMAX >= tminCP ~ TEMP2MMAX,
             TEMP2MMAX <= tminCP ~ tminCP,
             TRUE ~ 99999),
           # Daily Growing Degree Units.
           gdu = case_when( ((gmin + gmax)/2) - tminCP <= tminCP ~ tminCP,
                            TRUE ~ ((gmin + gmax)/2) - tminCP) )
  
  temp2 <- 
    temp1[,c("TIMESTAMP","LOCATION",  "gdu")]%>% 
    group_by(LOCATION) %>% 
    mutate(GS = as.factor(round(getGS(cumsum(gdu)),0))) %>% 
    group_by(GS, LOCATION) %>% 
    summarise(ANTHESIS = mean(TIMESTAMP)) %>% 
    filter(GS == "55")
  
  out <- data %>% full_join(temp2[,2:3],by = join_by(LOCATION))
  
  return(out)
  
}

# Weather summaries
weather_summaries = function(data, tminCP = 4.5, tminGF = 0, 
                             minGDU_CP = 300, maxGCU_CP = 100, 
                             minGDU_GF = 100, maxGCU_GF = 600)
{
  
  
  VARS_SUM = c("PP", "Tmean", "Duration", "PQ", "VPD", "ETE_CP", "ETE_GF")
  
  out = 
    data %>% 
    unnest(cols = c(daymet)) %>% 
    mutate(Tmean = (TEMP2MMIN + TEMP2MMAX)/2, 
           #Critical period growing degree units
           gmin = case_when(TEMP2MMIN >= tminCP ~ TEMP2MMIN, 
                            TRUE ~ tminCP),
           # Tmax threshold Growing Degrees.
           gmax = case_when(
             TEMP2MMAX >= tminCP ~ TEMP2MMAX,
             TEMP2MMAX <= tminCP ~ tminCP,
             TRUE ~ 99999),
           # Daily Growing Degree Units.
           gdu_ant = case_when( ((gmin + gmax)/2) - tminCP <= tminCP ~ tminCP,
                                TRUE ~ ((gmin + gmax)/2) - tminCP),
           #Grain Filling period growing degree units
           gmin_gf = case_when(TEMP2MMIN >= tminGF ~ TEMP2MMIN, 
                               TRUE ~ tminGF), 
           gmax_gf = case_when(
             TEMP2MMAX >= tminGF ~ TEMP2MMAX,
             TEMP2MMAX <= tminGF ~ tminGF,
             TRUE ~ 999999),
           gdu_gf = case_when( ((gmin_gf + gmax_gf)/2) - tminGF <= tminGF ~ tminGF,
                               TRUE ~ ((gmin_gf + gmax_gf)/2) - tminGF) ) %>% 
    group_by(LOCATION, GENOTYPE) %>% 
    mutate(ANTHESIS = ANTHESIS + MATURITY) %>% 
    nest(weather = -group_cols()) %>% 
    mutate(weather = weather %>% map(~.x %>% mutate(condition = case_when(TIMESTAMP >= ANTHESIS ~ 1,T ~ 0)) %>%
                                       # FROSTS EVENTS CP - T = 0C
                                       group_by(condition) %>% 
                                       mutate(accum_gdu = ifelse(condition == 0,  rev(cumsum(rev(gdu_ant))), cumsum(gdu_ant)),
                                              PERIOD = case_when(condition == 0 & accum_gdu <= minGDU_CP | condition == 1 & accum_gdu <= maxGCU_CP ~ "CP",
                                                                 condition == 1 & accum_gdu >= minGDU_GF & accum_gdu <= maxGCU_GF ~ "GF", T~ NA_character_) ) %>% 
                                       drop_na(PERIOD) %>% 
                                       group_by(ANTHESIS, PERIOD) %>%
                                       mutate(.,
                                              ndays = n(), 
                                              # Extreme low temperature events can damage wheat stand
                                              #ETE_TIL = case_when(TEMP2MMIN <= ETE_TIL ~ 1, T ~ 0), 
                                              ETE_CP = case_when(TEMP2MMIN <= -2 ~ 1, T ~ 0), 
                                              # Heat Stress events in Grain Filling
                                              ETE_GF = case_when(TEMP2MMAX >= 30 ~ 1, T ~ 0), 
                                              ) %>% 
                                       summarise_at(vars(TEMP2MMIN, TEMP2MMAX, VPDEFAVG, Tmean, PRECIP, ndays, SR, gdu_ant, gdu_gf, ETE_CP, ETE_GF),
                                                    list(mean = ~mean(., na.rm=T), 
                                                         sum = ~sum(.) ))  %>% 
                                       mutate(PQ = SR_sum/ case_when(PERIOD == "CP"~ gdu_ant_sum,
                                                                     PERIOD == "GF"~ gdu_gf_sum), 
                                              #ETE_TIL = ETE_TIL_sum, 
                                              ETE_CP = ETE_CP_sum,
                                              ETE_GF = ETE_GF_sum,
                                              
                                              VPD = VPDEFAVG_mean,
                                              Duration = ndays_mean, 
                                              Tmean = Tmean_mean,
                                              PP = PRECIP_sum) %>% 
                                       dplyr::select(PERIOD, ANTHESIS, all_of(VARS_SUM)) %>%
                                       pivot_wider(names_from = PERIOD, values_from = c(PP:ncol(.)))  )) %>% 
    unnest(cols = weather) %>% 
    ungroup()
  return(out)
}
