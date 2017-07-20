#packages to install

install.packages("microbenchmark")
install.packages("benchmarkme")
install.packages("profvis")

#set directory for download files
#with packrat, the file automatically download into the folder project 
#setwd("~/R/Download_Files")

url <- "http://s3.amazonaws.com/assets.datacamp.com/production/course_2987/datasets/movies.csv"
# <>  download.file(url, "movies.csv")
movies <- read.csv("movies.csv")


# save local copies
write.csv(movies, file = "movies.csv", row.names = FALSE)
saveRDS(movies, file = "movies.rds")


#------------------------------------------------------------------------The Art of Benchmarking

# ----------------------------------- R version

# Print the R version details using version
version

# Assign the variable `major` to the major component
major <-  version$major

# Assign the variable `minor` to the minor component
minor <- version$minor


# ----------------------------------- Comparing read times of CSV and RDS files

# How long does it take to read movies from CSV?
system.time(read.csv("movies.csv"))

# How long does it take to read movies from RDS?
system.time(readRDS("movies.rds"))


# ----------------------------------- Elapsed time

# Load the microbenchmark package
library("microbenchmark")

# Compare the two functions
compare <- microbenchmark(read.csv("movies.csv"), 
                          readRDS("movies.rds"), 
                          times = 10)

# Print compare
print(compare)


# ----------------------------------- DataCamp hardware

library("benchmarkme")

# Assign the variable `ram` to the amount of RAM on this machine
ram <- get_ram()
ram

# Assign the variable `cpu` to the cpu specs
cpu <- get_cpu()
cpu

# ----------------------------------- Benchmark DataCamp's machine

# Load the package
library(benchmarkme)

# Run the io benchmark
res <- benchmark_io(runs = 1, size = 5)

# Plot the results
plot(res)


#------------------------------------------------------------------------ Fine Tuning: Efficient Base R

# ----------------------------------- Timings - growing a vector

n <- 30000
# Slow code
growing <- function(n) {
  x <- NULL
  for (i in 1:n)
    x <- c(x, rnorm(1))
  x
}

# Use `<-` with system.time() to store the result as res_grow
system.time(res_grow <- growing(30000))

# ----------------------------------- Timings - pre-allocation 

n <- 30000
# Fast code
pre_allocate <- function(n) {
  x <- numeric(n) # Pre-allocate
  for (i in 1:n) 
    x[i] <- rnorm(1)
  x
}

# Use `<-` with system.time() to store the result as res_allocate
n <- 30000
system.time(res_allocate <- pre_allocate(n))


# ----------------------------------- Vectorized code: multiplication

x <- rnorm(10)
x2 <- numeric(length(x))
for (i in 1:10)
  x2[i] <- x[i] * x[i]

# Store your answer as x2_imp
x2_imp <- x*x

# ----------------------------------- Vectorized code: calculating a log-sum

# Initial code
n <- 100
total <- 0
x <- runif(n)
for (i in 1:n) 
  total <- total + log(x[i])

# Rewrite in a single line. Store the result in log_sum
log_sum <- sum(log(x))

# ----------------------------------- Data frames and matrices - column selection

mat = matrix(rnorm(1e5), ncol = 1000)
df = as.data.frame(mat)

# Which is faster, mat[, 1] or df[, 1]? 
microbenchmark(mat[,1], df[,1])

# ----------------------------------- Row timings

# Which is faster, mat[1, ] or df[1, ]? 
microbenchmark(mat[1,], df[1,])


#------------------------------------------------------------------------Diagnosing Problems: Code Profiling


# ----------------------------------- Profvis in action


# Load the profvis package
library(profvis)

# Profile the following code with the profvis function
profvis({
  # Load and select data
  movies <- movies[movies$Comedy == 1, ]
  
  # Plot data of interest
  plot(movies$year, movies$rating)
  
  # Loess regression line
  model <- loess(rating ~ year, data = movies)
  j <- order(movies$year)
  
  # Add a fitted line to the plot
  lines(movies$year[j], model$fitted[j], col = "red")
})     ## Remember the closing brackets!  


# ----------------------------------- Change the data frame to a matrix


# Load the microbenchmark package
library("microbenchmark")

# The previous data frame solution is defined
# d() Simulates 6 dices rolls
d <- function() {
  data.frame(
    d1 = sample(1:6, 3, replace = TRUE),
    d2 = sample(1:6, 3, replace = TRUE)
  )
}

# Complete the matrix solution
m <- function() {
  matrix(sample(1:6, 6, replace = TRUE), ncol = 2)
}

# Use microbenchmark to time m() and d()
microbenchmark(
  data.frame_solution = d(),
  matrix_solution     = m()
)

# -----------------------------------Calculating row sums

# Example data
rolls <- matrix(sample(1:6, 6, replace = TRUE), ncol = 2) 

# Define the previous solution 
app <- function(x) {
  apply(x, 1, sum)
}

# Define the new solution
r_sum <- function(x) {
  rowSums(x)
}

# Compare the methods
microbenchmark(
  app_sol = app(rolls),
  r_sum_sol = r_sum(rolls)
)


# ----------------------------------- Use && instead of &
# Example data

is_double = c(FALSE, TRUE, TRUE)


# Define the previous solution
move <- function(is_double) {
  if (is_double[1] & is_double[2] & is_double[3]) {
    current <- 11 # Go To Jail
  }
}

# Define the improved solution
improved_move <- function(is_double) {
  if (is_double[1] && is_double[2] && is_double[3]) {
    current <- 11 # Go To Jail
  }
}

## microbenchmark both solutions
microbenchmark(move, improved_move, times = 1e5)


#------------------------------------------------------------------------Turbo Charged Code: Parallel Programming Some problems can

# ----------------------------------- Can this loop run in parallel (1)?

total <- no_of_rolls <- 0
# Initialise
while (total < 10) {
  total <- total + sample(1:6, 1)
  if (total %% 2 == 0) {
    total <- 0
    # If even. Reset to 0
    no_of_rolls <- no_of_rolls + 1
  }
  no_of_rolls
}

#No: it's a sequential algorithm. The ith value depends on the previous value.

# ----------------------------------- Can this loop run in parallel (2)?

play <- function(){
  while (total < 10) {
    total <- total + sample(1:6, 1)
    if (total %% 2 == 0) {
      total <- 0
      # If even. Reset to 0
      no_of_rolls <- no_of_rolls + 1
    }
    no_of_rolls
  }
}

results <- numeric(100)
for (i in seq_along(results)) {
  results[i] <- play()
}

# Yes, this is embarrassingly parallel. We can simulate the games in any order.

# -----------------------------------Moving to parApply

library("parallel")

# Create a cluster using makeCluster().
# Do some work.
# Stop the cluster using stopCluster().
dd = matrix(rnorm(1000), ncol = 10)

# Determine the number of available cores
detectCores()

# Create a cluster via makeCluster
cl <- makeCluster(2)

# Parallelize this code
parApply(cl,dd, 2, median)

# Stop the cluster
stopCluster(cl)

# -----------------------------------The parallel package - parSapply

play <- function() {
  total <- no_of_rolls <- 0
  while (total < 10) {
    total <- total + sample(1:6, 1)
    
    # If even. Reset to 0
    if (total %% 2 == 0) total <- 0 
    no_of_rolls <- no_of_rolls + 1
  }
  no_of_rolls
}

# Create a cluster via makeCluster (2 cores)
cl <- makeCluster(2)

# Export the play() function to the cluster
clusterExport(cl, "play")

# Parallelize this code
res <- parSapply(cl, 1:100, function(i) play())

# Stop the cluster
stopCluster(cl)


# ----------------------------------- Timings parSapply()

# Set the number of games to play
no_of_games <- 1e5

## Time serial version
system.time(serial <- sapply(1:no_of_games, function(i) play()))

## Set up cluster
cl <- makeCluster(3)
clusterExport(cl, "play")

## Time parallel version
system.time(par <- parSapply(cl, 1:no_of_games, function(i) play()))

## Stop cluster
stopCluster(cl)


# use the microbechmark to get the time execution 

microbenchmark(serial, par, times = 1)
