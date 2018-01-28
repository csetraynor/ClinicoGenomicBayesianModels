

times <- sort(unique(clinical_data$os_months)

longdata <- survival::survSplit(Surv(time = os_months, event = deceased) ~ . , 
                      cut = times, data = (clinical_data %>%
                        filter(!is.na(os_status) & os_status != '') %>%
                        filter(os_months > 0 & !is.na(os_months))  %>%
  mutate(deceased = os_status == "DECEASED")))

longdata %>% group_by(sample) %>%
  select(sample, tstart, os_months, deceased) %>%
  top_n(3, tstart) %>% arrange(desc(sample, tstart))

longdata %>% group_by(sample) %>%
  select(sample, tstart, os_months, deceased) %>%
  slice(1:3) 
