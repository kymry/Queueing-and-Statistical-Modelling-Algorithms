# SMDE - Lab Assignment 2 - Kymry Burwell

install.packages("plyr")
install.packages("lmtest")
install.packages("FactoMineR")
install.packages("RCmdrPlugin.FactoMinerR")
install.packages("factoextra")
install.packages("car")
install.packages("gridExtra")
install.packages("randtests")
library(gridExtra)
library(car)
library(plyr)
library(lmtest)
library(ggplot2)
library(factoextra)
library(FactoMineR)
library(MASS)
library(randtests)

###### Random Number Generator #######

# Create a fuction for generting random numbers (multiplicative lagged fibonacci generator)
rm(list = ls())
fibonacci_RNG <- function(j = 10,k = 458 ,m = 2147483647,length) { 
    n <- 3000
    randomNumberList <- numeric(length)
    rnvector <- c(0:300000)
    
    for(i in 1:length){
        a <- rnvector[n-j] 
        b <- rnvector[n-k]
        rnvector[n] <- ((a*b) %% m)
        n <- n+1
        randomNumberList[i] <- rnvector[n-1]
    }
    return(randomNumberList)
}

# Generate a sequence of 100,000 random numbers
randomNumbers <- fibonacci_RNG(length = 100000)

# View a histogram of the random numbers to visually check for patterns
hist(randomNumbers)

# Test the random number generator for non-randomnes
difference.sign.test(randomNumbers)  # Tests conecutive pairs of random numbers and records the sign of the result
turning.point.test(randomNumbers) # Checks that the "turning points" in the data are normally distributed

###### Generate Data Set #######

# Create data frame of 10 factors and an answer 
rm(list = ls())
set.seed(23)
names <- paste(rep("Individual",2000),1:2000, sep="")
f1 <- rnorm(2000, mean = 50, sd = 10)
f2 <- rnorm(2000, mean = 20, sd = 5)
f3 <- rnorm(2000, mean = 35, sd = 7.5)
f4 <- rnorm(2000, mean = 40, sd = 10)
f5 <- rnorm(2000, mean = 5, sd = 1)
f6 <- (2*f1+4*f2+rnorm(2000,mean = 0, sd = 1))
f7 <- (3*f2+f5+rnorm(2000,mean = 2, sd = 1))
f8 <- (2*f3+f2+rnorm(2000,mean = 2, sd = .5))
f9 <- (5*f4+2*f5+rnorm(2000,mean = 0, sd = 2))
f10 <- (2*f1+2*f4+rnorm(200,mean = 1, sd = 1))
a1 <- (f5+f9+f6+f2+rnorm(2000,mean = 10, sd = 5))
df <- data.frame(names, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, a1)
df <- cbind(df,(sample(c(0,1), 2000, replace=T, prob = c(.8, .2)))) # Generate training/test data
colnames(df) <- c("names","factor_1","factor_2","factor_3","factor_4","factor_5","factor_6","factor_7","factor_8","factor_9","factor_10","answer","training_test")


####### Initial data exploration #######

# view the head of the data and some general information about the factors
head(df)
str(df)
summary(df)

# Plot each factor against the answer to visually see how they may be correlated
plot(df$factor_1,df$answer) 
plot(df$factor_2,df$answer)
plot(df$factor_3,df$answer)
plot(df$factor_4,df$answer)
plot(df$factor_5,df$answer)
plot(df$factor_6,df$answer)
plot(df$factor_7,df$answer)
plot(df$factor_8,df$answer)
plot(df$factor_9,df$answer)
plot(df$factor_10,df$answer)

# Test assumptions for PCA analysis
sapply(df, is.numeric) # Variables are all numeric
nrow(df)/ncol(df) # There are more observations than columns

# Perform PCA analysis
pca <- PCA(df[2:12], graph = TRUE)
fviz_eig(pca, addlabels = TRUE)

######## Define linear model #######

# Split data into training and test datasets
trainDF <- subset(df, training_test == 0) # training data
testDF <- subset(df, training_test ==1) # test data

# Compute mean of dependent variable (answer)
meanOfObservations <- mean(trainDF$answer) # Used to compare Residual Standard Error of linear model
meanOfObservations

# Generate multiple linear models
lm0 <- lm(answer ~ factor_1+factor_2+factor_3+factor_4+factor_5+factor_6+factor_7+factor_8+factor_9+factor_10, data = trainDF) # linear model 0
summary(lm0)
lm1 <- lm(answer ~ factor_9+factor_7+factor_1, data = trainDF) # linear model 1
summary(lm1)
lm2 <- lm(answer ~ factor_10+factor_6+factor_1, data = trainDF) # linear model 2
summary(lm2)
lm3 <- lm(answer ~ factor_4+factor_7+factor_1, data = trainDF) # linear model 3
summary(lm3)
lm4 <- lm(answer ~ factor_9+factor_6+factor_10, data = trainDF) # linear model 4
summary(lm4)
lm5 <- lm(answer ~ factor_4+factor_9+factor_10+factor_6, data = trainDF) # linear model 5
summary(lm5) 
lm6 <- lm(answer ~ factor_1+factor_2+factor_3+factor_4+factor_5, data = trainDF) # linear model 6
summary(lm6)

# Test linear model assumptions of "best" linear model (lm4)
dwtest(lm4, alternative = "two.sided") # independence of observations
shapiro.test(residuals(lm4)) # normality of the data
hist(residuals(lm4), freq=FALSE, main="Distribution of Residuals") # visual representation of normality test
lmtest::bptest(lm4) # homogeniety of variance
plot(fitted(lm4),residuals(lm4)) # visual representation of variance test

# Use test data set (created above) to run predictions using the "best" linear model (lm4)
predintDF <- predict(lm4, newdata = testDF, interval = "prediction")


######## Design Factorial Experiment #########

# Find min and max of each factor. This will be used as our high/low values in the factorial design
minF9 <- min(df$factor_9)
minF6 <- min(df$factor_6)
minF10 <- min(df$factor_10)
maxF9 <- max(df$factor_9)
maxF6 <- max(df$factor_6)
maxF10 <- max(df$factor_10)

# Create factorial design table
factorialDF <- expand.grid(factor_9 = c(minF9,maxF9), factor_6 = c(minF6,maxF6), factor_10 = c(minF10,maxF10))

# Generate answer using linear model and append to factorial design table (using lm4)
answerLM4 <- data.frame(predict(lm4, newdata = factorialDF, interval = "predict"))
factorialDF <- cbind(factorialDF,answerLM4$fit)
factorialDF <- cbind(factorialDF,c("mean","A","B","AB","C","AC","BC","ABC"))
colnames(factorialDF) <- c("factor_9","factor_6","factor_10","answer","ID")

# Calculate the size (range) of all confidence intervals of the experiment. 
ciExperiment <- data.frame(predict(lm4, newdata = factorialDF, interval="confidence"))
ciExperimentRanges <- ciExperiment$upr-ciExperiment$lwr
ciExperimentRanges

# Run Fisher-Yate algorithm - Please see excel file for computation. 


