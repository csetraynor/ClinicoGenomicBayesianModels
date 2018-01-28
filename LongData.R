
#obtain unique times
times <- sort(unique(clinical_data$os_months))

#create long data with observed and censored times
longdata <- survival::survSplit(Surv(time = os_months, event = deceased) ~ . , 
                      cut = times, data = (clinical_data %>%
                        filter(!is.na(os_status) & os_status != '') %>%
                        filter(os_months > 0 & !is.na(os_months))  %>%
  mutate(deceased = os_status == "DECEASED")))

#Checkin the correctness of the new dataset
longdata %>% group_by(sample) %>%
  select(sample, tstart, os_months, deceased) %>%
  top_n(3, tstart) %>% arrange(desc(sample, tstart))

longdata %>% group_by(sample) %>%
  select(sample, tstart, os_months, deceased) %>%
  dplyr::slice(1:3) 

assertthat::assert_that((longdata %>% group_by(sample) %>% dplyr::n_groups()) == (clinical_data %>% filter(!is.na(os_status) & os_status != '') %>% filter(os_months > 0 & !is.na(os_months)) %>% group_by(sample) %>% dplyr::n_groups()))


#Create timepoint id
longdata <- longdata %>%
  group_by(sample) %>%
  mutate(t = seq(n())) %>%
  select(sample, t, tstart)

longdata %>% group_by(sample) %>%
  select(sample, tstart, t) %>%
  dplyr::slice(1:3) 

#generate data

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
