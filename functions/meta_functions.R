###############################
###  Functions to get I2  ####
###############################

# Instructions to hand calculate I2 are here https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
#Function to get I2
get_I2_overall <- function(.x, .y) { #.x is the dataset supplied to the rma.mv function, .y is the output of that function
  W <- diag(1/.x$var_dunb)
  X <- model.matrix(.y)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2 <- 100 * sum(.y$sigma2)/ (sum(.y$sigma2) + (.y$k-.y$p)/sum(diag(P)))
  return(I2)
}

# And getting it cleaned per a list of results

get_I2_overall_clean <- function(x, y){ #x is data that was mapped, y is meta results
  I2 <- data.frame(map2(x, y, ~get_I2_overall(.x, .y)))
  I2 %>% 
    pivot_longer(1:ncol(I2)) %>% 
    mutate(I2_overall = value) %>% 
    select(-value, -name)
}

# total variance attributed to between and within-level clusters separately
## From left to right in every model unless specified:
## Between sample heterogeneity
## Between condition heterogeneity
## Between outcome heterogeneity
## Hetereogeneity due to outcome type crossed between studies

get_I2_var_levels <- function(.x, .y) { #.x is the dataset supplied to the rma.mv function, .y is the output of that function
  W <- diag(1/.x$var_dunb)
  X <- model.matrix(.y)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2 <- round(100 * (.y$sigma2)/ (sum(.y$sigma2) + (.y$k-.y$p)/sum(diag(P))), digits = 2)
  return(I2)
}

# Less code for the same results as above and suggestion by the same person - the creator of metafor
## https://stats.stackexchange.com/questions/187197/computing-heterogeneity-assigned-to-random-factors-in-meta-analysis

get_I2_levels <- function(model){
  wi <- 1/model$vi
  k <- length(wi)
  s2m <- (k-1)*sum(wi)/(sum(wi)^2 - sum(wi^2))
  s2t <- s2m + sum(model$sigma2)
  round((model$sigma2 / s2t)*100, digits = 2)
}

# Get column for between study hetereogeneity to put in table
get_btw_study_I2 <- function(x){#x is output from function that gets I2 var levels per study
  between_study_I2 <- list(length = length(x))
  
  for (i in seq_along(x)){
    between_study_I2[[i]] <- x[[i]][1]
  }
  
  btw_study_I2 <- data.frame(between_study_I2) 
  btw_study_I2 %>% 
    pivot_longer(1:ncol(btw_study_I2)) %>% 
    mutate(I2_btw_study = value) %>% 
    select(-value, -name)
}

# Get column for within study heterogeneity
get_within_study_I2 <- function(x){#ne meta-analyses with three-level nesting structure, this woudl correspond to within-study variance (within-outcome variance)
  within_study_I2 <- list(length = length(x))
  
  for (i in seq_along(x)){
    within_study_I2[[i]] <- x[[i]][2]
  }
  
  within_I2 <- data.frame(within_study_I2) 
  within_I2 %>% 
    pivot_longer(1:ncol(within_I2)) %>% 
    mutate(I2_within_study = value) %>% 
    select(-value, -name)
}

### Getting values for crossed scale hetereogeneity into dataframe column to put in table later

get_crossed_study_I2 <- function(x){#x is output from function that gets I2 var levels per study
  crossed_study_I2 <- list(length = length(x))
  
  for (i in seq_along(x)){
    crossed_study_I2[[i]] <- x[[i]][3]
  }
  
  crossed_I2 <- data.frame(crossed_study_I2) 
  crossed_I2 %>% 
    pivot_longer(1:ncol(crossed_I2)) %>% 
    mutate(I2_crossed_scale = value) %>% 
    select(-value, -name)
}

# Function that subsets data into multiple datasets based on two filters & runs an rma.mv meta-analysis 
filtered_meta <- function(x, y, z){ #x is the overall dataset, y is the first set of filters, z is the second set of filters
  df <- list(length = length(y))
  
  for (i in seq_along(y)){
    df[[i]] <- x %>%
      dplyr::filter(pt_comparison == y[i] & 
                      outcome_type == z[i])
  }
  
  return(map(df, ~rma.mv(dunb,
                         var_dunb,
                         random = list( ~ 1 | sample_number_total/outcomes_within_sample_var,
                                        ~ 1 | outcome_scale_group),
                         data = .x)))
}


###############################
### Function to clean data ####
###############################

meta_clean_func <- function(converted_data){
  converted_data %>%  
    mutate(cohens_d = as.numeric(cohens_d),
           reverse = as.factor(reverse)) %>% 
    mutate(coded_cohens_d_var = (((n_pt + n_comparison)/(n_pt*n_comparison)) + ((cohens_d^2)/(2*(n_pt+n_comparison))))) %>% #calculated variance by hand for the d's we extracted from articles directly; the ones we converted also came with var calculations
    pivot_longer(c(cohens_d, d.x, d.y, yi, d_values)) %>% #  Wrangling all converted d's into one column
    mutate(d = value) %>% 
    filter(!is.na(d)) %>% 
    dplyr::select(-c(name, value)) %>% 
    pivot_longer(c(coded_cohens_d_var, var_d.x, var_d.y, vi, d_values_var)) %>% #wrangling all variances into one column
    mutate(var = value) %>% 
    filter(!is.na(var)) %>% 
    dplyr::select(-c(name, value)) %>% 
    mutate(df = (n_pt + n_comparison - 2), #getting degrees of freedom
           J = 1- (3/((4*df)-1)), #calculating hedges correction factor for d unbiased
           dunb = J*d, #d unbiased
           var_dunb = ((J^2)*var),  #variance for d unbiased
           lowCI_dunb = dunb-1.96*sqrt(var_dunb), #getting 95% CI's for d unbiased
           upCI_dunb  = dunb+1.96*sqrt(var_dunb)) %>% 
    mutate(dunb = ifelse(reverse == "yes", 
                         dunb*-1,
                         ifelse(reverse == "no",
                                dunb*1, NA))) %>% #reverse scored dunb that needed to be 
    dplyr::select(-c(DOI, Notes, outcome_description_from_coding, target_long_description,
                     weird_sample, F_score, t_score, effect_size_type,
                     effect_direction, p_value, mean_pt, sd_pt, se_pt, 
                     mean_comparison, sd_comparison,se_comparison, 
                     `Note about directionality of cohen's d`)) %>% 
    mutate(sample_number_total = as.factor(sample_number_total),
           pt_comparison = as.factor(pt_comparison),
           outcome_type = as.factor(outcome_type),
           target_ingroup_nonspecific = as.factor(target_ingroup_nonspecific),
           outcome_scale_group = as.factor(outcome_scale_group),
           target_out_minor = as.factor(target_out_minor),
           target_adversary = as.factor(target_adversary),
           target_emapthetic_need = as.factor(target_empathetic_need),
           between_within = as.factor(between_within),
           target_information = as.factor(target_information))
}

################################
### Extract output from RMA ####
################################

# Functions to extract values from rma.mv output; needed to help when cleaning my nested output
get_pval_rma <- function(x){
  x$pval
}

get_ci_ub_rma <- function(x){
  x$ci.ub
}

get_ci_lb_rma <- function(x){
  x$ci.lb
}

get_Q_rma <- function(x){
  x$QE
}

get_Q_p_rma <- function(x){
  x$QEp
}

# Functions to clean values extracted from rma.mv output

get_k <- function(x){
  x$k
}

get_k_clean <- function(x){
  k_list <- map(x, ~get_k(.x)) %>% 
    as.data.frame()
  k_list %>% 
    pivot_longer(1:ncol(k_list)) %>% 
    mutate(k = value) %>% 
    dplyr::select(-name, - value) %>% 
    round(digits = 3)
}

get_ests_clean <- function(x){
  est <- map(x, ~coef(.x)) %>% 
    as.data.frame()
  est %>% 
    pivot_longer(1:ncol(est)) %>% 
    mutate(estimate = value) %>% 
    dplyr::select(-name, - value) %>% 
    round(digits = 3)
}

get_pval_clean <- function(x){#x is data frame
  pvals <- map(x, ~get_pval_rma(.x)) %>% 
    as.data.frame()
  pvals %>% 
    pivot_longer(1:ncol(pvals)) %>% 
    mutate(est_pval = value) %>% 
    dplyr::select(-name, - value) %>% 
    round(digits = 3)  
}

get_ci_upper_clean <- function(x) {
  uppers <- map(x, ~get_ci_ub_rma(.x)) %>% 
    as.data.frame()
  uppers %>% 
    pivot_longer(1:ncol(uppers)) %>% 
    mutate(ci_upper = value) %>% 
    dplyr::select(-name, - value) %>% 
    round(digits = 3)
}

get_ci_lower_clean <- function(x) {
  lower <- map(x, ~get_ci_lb_rma(.x)) %>% 
    as.data.frame()
  lower %>% 
    pivot_longer(1:ncol(lower)) %>% 
    mutate(ci_upper = value) %>% 
    dplyr::select(-name, - value) %>% 
    round(digits = 3)
}

get_Q_clean <- function(x){
  Q <- map(x, ~get_Q_rma(.x)) %>% 
    as.data.frame()
  Q %>% 
    pivot_longer(1:ncol(Q)) %>% 
    mutate(Q = value) %>% 
    dplyr::select(-name, - value) %>% 
    round(digits = 3)
}

get_Q_p_clean <- function(x){
  Q_ps <- map(x, ~get_Q_p_rma(.x)) %>% 
    as.data.frame()
  Q_ps %>% 
    pivot_longer(1:ncol(Q_ps)) %>% 
    mutate(Q_pval = value) %>% 
    dplyr::select(-name, - value) %>% 
    round(digits = 4)
}

##### Function to print plots to pngs #####

ggpreview <- function(...) {
  fname <- tempfile(fileext = ".png")
  ggsave(filename = fname, ...)
  system2("open", fname)
  invisible(NULL)
}
