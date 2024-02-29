source(here::here("R/cv.cureitlasso.R"))
source(here::here("R/cureitlasso.R"))
source(here::here("R/simulasso.R"))
source(here::here("R/coxsplit.R"))

# Potentially change sample sizes and number of parameters
# We could try simulating scenarios with overlap between cure and cox

# Simulate Cure Data --------------------------------------------------------

# Test- Try simulating one train and test set


dat <- simulasso() # training data



dat_valid <- simulasso(id_cure = dat$id_cure,
          id_cox = dat$id_cox,
          coefs_cure = dat$coefs_cure,
          coefs_cox = dat$coefs_cox)


setdiff(dat$id_cox, dat_valid$id_cox)


# * Training- Simulate n Data Sets --------------------------------------------------------
sim_data_list <- list()
set.seed(999)

for (i in 1:100) {
  sim_data_list[[i]] = simulasso()
}

sim_data_list[[1]]$id_cure

# * Validation - Simulate n Validation Data ---------------------------------------------------

sim_valid_data_list <- list()

for (i in 1:100) {
  dat <- sim_data_list[[i]]
  sim_valid_data_list[[i]] = simulasso(id_cure = dat$id_cure,
                                       id_cox = dat$id_cox,
                                       coefs_cure = dat$coefs_cure,
                                       coefs_cox = dat$coefs_cox)
}

sim_data_list[[1]]$id_cure
sim_valid_data_list[[1]]$id_cure

save(sim_data_list, file = here::here("inst", 
                                      "cureit-simulation-results",
                                      "sim_data_list_v2.RData"))

save(sim_valid_data_list, file = here::here("inst", 
                                            "cureit-simulation-results",
                                            "sim_valid_data_list_v2.RData"))

# Fit Lasso  ------------------------------------------------------------

# test- A Single Lasso Fit 

dat <- sim_data_list[[1]]
test_single_run <- cureitlasso(
            t = dat$t,
            d = dat$d,
            x = dat$x,
            minmu.ratio=0.05, # minimum penalty value in the logistic model
            minlambda.ratio=0.05, # minimum penalty value in the cox model
            adaptive=FALSE,
            length.grid=8,
            tol=1e-2,
            maxiter=100,
            progress=FALSE)

# Run Cross Validation ----------------------------------------------------

# * Test a Single CV Fit ------
# To get a sense of what is in the output

dat_valid <- sim_data_list[[1]]

test_single_cv <- cv.cureitlasso(t = dat_valid$t,
                                 d = dat_valid$d,
                                 x = dat_valid$x,
                                 minmu.ratio=0.05, # minimum penalty value in the logistic model
                                 minlambda.ratio=0.05, # minimum penalty value in the cox model
                                 adaptive=FALSE,
                                 length.grid=10,

                                 tol=1e-2,
                                 maxiter=100,

                                 nfolds=5,
                                 ncore=5,
                                 seed=NULL,
                                 progress=FALSE)

test_list <- list()
test_list[[1]] <- test_single_cv

# * 100 Runs ----

cv_fits_list <- list()
start_time <- Sys.time()

for (i in 1:50) {

  dat_valid <- sim_data_list[[i]]
  
  cv_fits_list[[i]] <- cv.cureitlasso(t = dat_valid$t,
                                    d = dat_valid$d,
                                    x = dat_valid$x,
                                    minmu.ratio=0.05, # minimum penalty value in the logistic model
                                    minlambda.ratio=0.05, # minimum penalty value in the cox model
                                    adaptive=FALSE,
                                    length.grid=10,

                                    tol=1e-2,
                                    maxiter=100,

                                    nfolds=5,
                                    ncore=5,
                                    seed=NULL,
                                    progress=FALSE
  )


}

#  8.670569 mins for 3 data sets 
end_time <- Sys.time()
end_time - start_time


save(cv_fits_list, file = here::here("inst", "cureit-simulation-results", "cv_fits_50_v2.RData"))

# * 50 Runs ----

cv_fits_list2 <- list()

start_time <- Sys.time()

for (i in 1:50) {
  
  dat_valid <- sim_data_list[[i]]
  
  cv_fits_list2[[i]] <- cv.cureitlasso(t = dat_valid$t,
                                      d = dat_valid$d,
                                      x = dat_valid$x,
                                      minmu.ratio=0.05, # minimum penalty value in the logistic model
                                      minlambda.ratio=0.05, # minimum penalty value in the cox model
                                      adaptive=FALSE,
                                      length.grid=10,
                                      
                                      tol=1e-2,
                                      maxiter=100,
                                      
                                      nfolds=5,
                                      ncore=5,
                                      seed=NULL,
                                      progress=FALSE
  )
  
  
}

# 2.369856 hours for 50 fits
end_time <- Sys.time()
end_time - start_time


save(cv_fits_list2, file = here::here("inst", "cv_fits_50.RData"))

