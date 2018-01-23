

###############################################
#Data obtantion
#get data from with MSKCC package 
require(cgdsr)
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

glioblastome_2013_id_sutdy = getCancerStudies(mycgds)[55,1]
glioblastome_2013_case_list = getCaseLists(mycgds, glioblastome_2013_id_sutdy)[2,1]
glioblastome_2013_clinical_data <-  getClinicalData(mycgds, glioblastome_2013_case_list)

#inspect dataframe
str(glioblastome_2013_clinical_data, no.list = T, vec.len = 2)

####################################################################
#Data Cleaning

#convert to lower case
names(glioblastome_2013_clinical_data) <- tolower(names(glioblastome_2013_clinical_data)) 

#convert missig values
convert_blank_to_na <- function(x) {
  if(!purrr::is_character(x)){
    warning('input variate not character - return original')
    return(x)
  } else {
    ifelse(x == '', NA, x)
  }
}
glio_clin_dat <- glioblastome_2013_clinical_data %>%
  dplyr::mutate_all(funs(convert_blank_to_na))

#inspect resulting dataframe
str(glio_clin_dat)

######################################################################
#Data Exploration
#Considering overall survival#
library(dplyr)

glio_clin_dat %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)

#filter unknown or negative survival times (os_monts < 0)

glio_clin_dat %>%
  filter(!is.na(os_status) & os_status != '') %>%
  filter(os_months < 0 | is.na(os_months)) %>%
  select(os_status, os_months) %>%
  head()

glio_clin_dat %>%
  filter(is.na(os_status) | os_status == '') %>%
  select(os_status, os_months) %>%
  str() 

#for now this observation will be remove from the analysis

glio_clin_dat <- glio_clin_dat %>%
  filter(!is.na(os_status) & os_status != '') %>%
  filter(os_months >= 0 & !is.na(os_months))
  
#Check 44 fewer obsrvations than original
assertthat::assert_that(nrow(glio_clin_dat) == nrow(glioblastome_2013_clinical_data) - 44)


########## Distribution of event times  ######################
library(ggplot2)
glio_clin_dat %>%
  ggplot(aes(x = os_months,
             group = os_status,
             colour = os_status,
             fill = os_status)) +
  geom_density(alpha = 0.5)

#KM curve
require(survival)
mle.surv <- survfit(Surv(os_months, os_deceased) ~ 1,
                    data = glio_clin_dat %>%
                      mutate(os_deceased = (os_status == "DECEASED")))
ggplot2::autoplot(mle.surv, conf.int = F) +
  ggtitle('KM survival for GGM Cohort')



############ Parametric Survival Model #####################


observed_data <- glio_clin_dat %>%
  filter(os_status == "DECEASED")

censored_data <- glio_clin_dat %>%
  filter(os_status != "DECEASED")

stan_data <- list(
  Nobs = nrow(observed_data),
  Ncen = nrow(censored_data),
  yobs = observed_data$os_months,
  yceb = censored_data$os_months
)
rm(censored_data)
rm(observed_data)

str(stan_data)

#Wraped in a function

gen_stan_data <- function(data) {
  observed_data <- data %>%
    dplyr::filter(os_status == 'DECEASED')
  
  censored_data <- data %>%
    dplyr::filter(os_status != 'DECEASED')
  
  stan_data <- list(
    Nobs = nrow(observed_data),
    Ncen = nrow(censored_data),
    yobs = observed_data$os_months,
    ycen = censored_data$os_months
  )
}

######### Setting intial values

gen_inits <- function() {
  list(
    alpha_raw = 0.01*rnorm(1),
    mu = rnorm(1)
  )
}

########################## Stan run ###########################
library(rstan)

stanfile <- "ClinicoGenomicBayesianModels/ClinicalGliobastomeParametricWithWeibull.stan"
#open stan file
if (interactive())
  file.edit(stanfile)

weibull_null_model <-  stan(stanfile,
              data = gen_stan_data(glio_clin_dat),
              chains = 4,
              iter = 1000,
              init = gen_inits
  )

####################### Checking convergence ###################

print(weibull_null_model) #(Check Rhat close to 1)

rstan::traceplot(weibull_null_model, 'lp__') #Review traceplot for log-posterior

rstan::traceplot(weibull_null_model, c('alpha','mu'), ncol = 1)    #Review traceplot for parameters of interest

if(interactive())
  shinystan::launch_shinystan(weibull_null_model)        #Launch shiny stan


#Review posterior distributions of parameters

pp_alpha <- rstan::extract(weibull_null_model, 'alpha')$alpha
pp_mu <- rstan::extract(weibull_null_model, 'mu')$mu


ggplot(data.frame(alpha = pp_alpha, mu = pp_mu)) + 
  geom_density(aes(x = alpha))+
  ggtitle('Posterior distribution of alpha')

ggplot(data.frame(alpha = pp_alpha, mu = pp_mu))+
  geom_density(aes(x = mu))+
  ggtitle('Posterior distribution of mu')
  
#Degree of correlation between alpha and mu
ggplot(data.frame(alpha = pp_alpha, mu = pp_mu))+
  geom_density2d(aes(x = alpha, y = mu))+
  ggtitle('Posterior distribution of alpha and mu')

######################### Posterior predicitive checks ###################################
#Simulate time to event data
#weibull_sim_data function takes two parameters (alpha and mu) as inputs and a desired sample size (n). 


weibull_sim_data <- function(alpha, mu, n) {
  
  data <- data.frame(surv_months = rweibull(n = n, alpha, exp(-(mu)/alpha)),
                     censor_months = rexp(n = n, rate = 1/100),
                     stringsAsFactors = F
  ) %>%
    dplyr::mutate(os_status = ifelse(surv_months < censor_months,
                                     'DECEASED', 'LIVING'
    ),
    os_months = ifelse(surv_months < censor_months,
                       surv_months, censor_months
    )
    )
  
  return(data)
}
#Censoring is "arbitrarily" rexp() , censoring is assumed to be noninformative.

######## Simulating data for each posterior draw #
test_n <- nrow(glio_clin_dat)
pp_newdata <- purrr::map2(.x = pp_alpha,
                          .y = pp_mu,
                          .f = ~weibull_sim_data(alpha = .x,
                                                 mu = .y,
                                                 n = test_n))

###### Plot time to event in the posterior draws compare to actual time to event
ggplot(pp_newdata %>%
         bind_rows() %>%
         mutate(type = 'posterior predicted values') %>%
         bind_rows(glio_clin_dat %>% mutate(type = 'actual data'))
       , aes(x = os_months, group = os_status, colour = os_status, fill = os_status))+ 
  geom_density(alpha = 0.5) +
  facet_wrap(~type, ncol = 1)

