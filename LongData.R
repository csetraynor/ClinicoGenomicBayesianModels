
#create sample id
clinical_data <- rownames_to_column(clinical_data, var = "s") #sample id in numeric
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
  mutate(t = seq(n())) 

longdata %>% group_by(sample) %>%
  select(sample, tstart, t) %>%
  dplyr::slice(1:3) 

#Calculate the duration of each time point
longdata %>% 
  group_by(s) %>%
  mutate(t_dur = os_months - tstart) %>%
  select(s, t, tstart, os_months, t_dur) %>%
  View()

longdata <- longdata %>%
  group_by(sample) %>%
  mutate(t_dur = os_months - tstart) 

#Wrap in a function
#generate data

gen_stan_data <- function(data) {
  #create long data
  data <- tibble::rownames_to_column(data, var = "s")
  times <- sort(unique(data$os_months))
  longdata <- survival::survSplit(Surv(time = os_months, event = deceased) ~ . , 
                                  cut = times, data = (data %>%
                                                         filter(!is.na(os_status) & os_status != '') %>%
                                                         filter(os_months > 0 & !is.na(os_months))  %>%
                                                         mutate(deceased = os_status == "DECEASED")))
  #create time point id
  longdata <- longdata %>%
    group_by(sample) %>%
    mutate(t = seq(n())) 
  
  #calculate duration
  longdata <- longdata %>%
    group_by(sample) %>%
    mutate(t_dur = os_months - tstart) 
  
  stan_data <- list(
    N = nrow(longdata),
    S = dplyr::n_distinct(longdata$sample),
    "T" = length(times),
    s = as.numeric(longdata$s),
    t = longdata$t,
    event = longdata$deceased,
    t_obs = longdata$os_months,
    t_dur = longdata$t_dur
  )
}

stanfile <- "null_pem_survival_model.stan"
#open stan file
if (interactive())
  file.edit(stanfile)

pem_null_model <-  stan(stanfile,
                            data = gen_stan_data(clinical_data),
                            chains = 1,
                            iter = 5
)
