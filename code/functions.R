
# Colors 
violeta <- '#75338a'
naranja <- '#de870d'
mix <- "#aa5d4c"

# Custom ggplot theme to make pretty plots
theme_clean <- function() {
  theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = rel(1), hjust = 0),
          strip.background = element_rect(fill = "grey80", color = NA),
          legend.title = element_text(face = "bold"))
}

## Mean square prediction error
MSPE <- function(model, testingData){
  
  yTest <- testingData[,2]
  yHat <- colMeans(rstan::extract(model, pars = "muTest")$muTest)
  
  MSPE <- sum((yTest - yHat)^2)/length(yHat)
  return(MSPE)
}

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
  
  
  VARS_SUM = c("PP", "Tmean", "Duration", "PQ", "VPD")
  
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
                                       mutate(ndays = n())
                                     %>% 
                                       summarise_at(vars(TEMP2MMIN, TEMP2MMAX, VPDEFAVG, Tmean, PRECIP, ndays, SR, gdu_ant, gdu_gf),
                                                    list(mean = ~mean(., na.rm=T), 
                                                         sum = ~sum(.) ))  %>% 
                                       mutate(PQ = SR_sum/ case_when(PERIOD == "CP"~ gdu_ant_sum,
                                                                     PERIOD == "GF"~ gdu_gf_sum), 
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



#Calculate Plasticity
plasticity3 <- function(data, genotype, variable){
  
  data1 = 
    data %>%
    group_by({{genotype}}) %>% 
    mutate_at(variable, list(VAR_G =  ~var(., na.rm=T))) %>% 
    ungroup() %>% 
    mutate_at(variable, list(VAR_E =  ~var(., na.rm=T))) %>% 
    mutate(name = VAR_G/VAR_E) %>% 
    dplyr::select({{genotype}}, name ) %>% 
    rename_at(vars(contains("name")), ~paste0("PP_", variable[[1]])) %>% 
    
    unique()
  
  out = full_join(data, data1, by = join_by({{genotype}}))
  return(out)
}

calculate_sd <- function(lb, ub, m) {
  # Z-Score for a 95% confidence interval
  z_score <- qnorm(0.95)
  ci <- ub - lb
  # Calculate standard deviation
  sd <- ci / (2 * z_score)
  return(sd)
}

# Formula for quantile regression equation
formula_quantile <- function(fit){
  b0 <- paste0(round(fixef(fit)[1,1],2))
  b1 <- paste0(round(fixef(fit)[2,1],2))
  pc <- paste0(round(fixef(fit)[3,1],2))
  
  pred <- linpred_draws(fit, 
                        newdata = data.frame(PP_GY = fixef(fit)[3,1]),
                        value = ".prediction")
  pc_mean <- round(mean(pred$.prediction),2)
  
  paste0("y = ", b0 ," + ", b1 ,"x if x < ", pc, "; y = ", pc_mean ," otherwise.")
}

