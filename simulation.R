# data simulation

################################################

# simulate graded response data
num_items <- 10
num_persons <- 1000

item_discri <- runif(num_items, min = 1.5, max = 2.5)
item_diffi_avg <- runif(num_items, min = 0, max = 1.25)
item_diff <- matrix(c(item_diffi_avg-1, item_diffi_avg-0.3, item_diffi_avg+0.3, item_diffi_avg+1),num_items, 4)
itempar <- cbind(item_discri, item_diff)

theta <- sort(rnorm(num_persons, 0, 1))

response <- GRM_sim(theta, itempar)
