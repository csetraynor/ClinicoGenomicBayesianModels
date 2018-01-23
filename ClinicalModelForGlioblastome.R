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


##########Distribution of event times######################
library(ggplot2)
library(ggfortify)
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
autoplot(mle.surv, conf.int = F) +
  ggtitle('KM survival for GGM Cohort')
