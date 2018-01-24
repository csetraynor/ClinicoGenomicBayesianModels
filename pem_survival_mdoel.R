stan_file <- system.file('stan', 'pem_survival_model.stan', package =  'biostan')

biostan::print_stan_file(stan_file)
if (interactive())
  file.edit(stan_file)


########----- generate dataset for Piecewise Exponential Model ------########
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