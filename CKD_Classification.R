# Read original dataset
originalData <- read.csv("Data/chronic_kidney_disease_full.csv", header = TRUE)
colnames(originalData)

head(originalData, n = 5)

# Replace extraneous characters, whitespace and blankspace with NA values
replacedData <- read.csv("Data/chronic_kidney_disease_full.csv", 
	header = TRUE, na.strings = c("", " ", "?"))
head(replacedData, n = 1)

any(is.na(replacedData))

# Omit observations with NA entries
cleanedData <- na.omit(replacedData)
any(is.na(cleanedData))

str(cleanedData)

# Correct the variable types
formattedData <- cleanedData
formattedData$id <- as.integer(formattedData$id)
formattedData$X.sg. <- as.factor(formattedData$X.sg.)
formattedData$X.al. <- as.factor(formattedData$X.al.)
formattedData$X.su. <- as.factor(formattedData$X.su.)
formattedData$X.rbc. <- as.factor(formattedData$X.rbc.)
formattedData$X.pc. <- as.factor(formattedData$X.pc.)
formattedData$X.pcc. <- as.factor(formattedData$X.pcc.)
formattedData$X.ba. <- as.factor(formattedData$X.ba.)
formattedData$X.htn. <- as.factor(formattedData$X.htn.)
formattedData$X.dm. <- as.factor(formattedData$X.dm.)
formattedData$X.cad. <- as.factor(formattedData$X.cad.)
formattedData$X.appet. <- as.factor(formattedData$X.appet.)
formattedData$X.pe. <- as.factor(formattedData$X.pe.)
formattedData$X.ane. <- as.factor(formattedData$X.ane.)
formattedData$X.class. <- as.factor(formattedData$X.class.)
str(formattedData)

summary(formattedData)

# Extract continuous numerical variables
correlationMatrix <- data.frame(formattedData$X.age., formattedData$X.bp., 
	formattedData$X.bgr., formattedData$X.bu., formattedData$X.sc., 
		formattedData$X.sod., formattedData$X.pot., formattedData$X.hemo., 
			formattedData$X.pcv., formattedData$X.wbcc., formattedData$X.rbcc.)

# Abbreviate names
colnames(correlationMatrix) <- c("X.age.", "X.bp.", "X.bgr.", "X.bu.", "X.sc.", 
	"X.sod.", "X.pot.", "X.hemo.", "X.pcv.", "X.wbcc.", "X.rbcc")

# Correlation plot
plot(correlationMatrix, main = "Scatterplot of Correlation Matrix")

# Compute Pearson correlation coefficients and round r-values
round(cor(correlationMatrix, method = "pearson"), digits = 2)

abs(cor(correlationMatrix, method = "pearson")) > 0.7

table(formattedData$X.class.)

# Generate a row-randomized dataset
set.seed(42)
randomizedData <- formattedData[sample(nrow(formattedData)), ]

# Construct testing dataset
testingData <- randomizedData # Dataset is constructed by Complement Rule

# Construct training dataset
trainingData <- testingData[-c(1:157), ] # Copy column names & preserve type

for (i in 1:length(testingData$id)) # Extract first 22 CKD observations
{
	if ((testingData[i, ]$X.class. == "ckd") & 
		(sum(trainingData$X.class. == "ckd") < 22))
	{
		trainingData[nrow(trainingData) + 1, ] <- testingData[i, ]
		testingData <- testingData[-c(i), ]
	}
}

for (i in 1:length(testingData$id)) # Extract first 57 non-CKD observations 
{
	if ((testingData[i, ]$X.class. == "notckd") & 
		(sum(trainingData$X.class. == "notckd") < 57))
	{
		trainingData[nrow(trainingData) + 1, ] <- testingData[i, ]
		testingData <- testingData[-c(i), ]
	}
}

rm(i) # Clears the counter variable from the environment

# Reorder datasets
trainingData <- trainingData[order(trainingData$id), ]
testingData <- testingData[order(testingData$id), ]

# Generate the Null Model
nullModel <- glm(X.class. ~ 1, data = trainingData, family = "binomial")

# Generate the Full Model
fullModel <- glm(X.class. ~ X.age. + X.bp. + X.bgr. + X.bu. + X.sod. + X.pot. + 
	X.hemo. + X.pcv. + X.wbcc. + X.rbcc., data = trainingData, 
		family = "binomial")

# Run Forward Selection algorithm
forwardModel <- step(nullModel, direction = "forward", 
	scope = list(upper = fullModel, lower = ~1), trace = 0)
summary(forwardModel)

# Run the Backward Elimination algorithm
backwardModel <- step(fullModel, direction = "backward", trace = 0)
summary(backwardModel)

# Run Sequential Selection algorithm
sequentialModel <- step(nullModel, direction = "both", 
	scope = formula(fullModel), trace = 0)
summary(sequentialModel)

# Load the glmnet package
library(glmnet)

# k-fold cross-validation
lassoY <- trainingData$X.class.
lassoX <- data.matrix(trainingData[, colnames(trainingData)[2:25]])
lassoCVModel <- cv.glmnet(lassoX, lassoY, alpha = 1, family = "binomial")
lassoBestLambda <- lassoCVModel$lambda.min
lassoBestLambda

# Run LASSO Regression algorithm
lassoModel <- glmnet(lassoX, lassoY, alpha = 1, lambda = lassoBestLambda, 
	family = "binomial")
coef(lassoModel)

# k-fold cross-validation
ridgeY <- trainingData$X.class.
ridgeX <- data.matrix(trainingData[, colnames(trainingData)[2:25]])
ridgeCVModel <- cv.glmnet(ridgeX, ridgeY, alpha = 0, family = "binomial")
ridgeBestLambda <- ridgeCVModel$lambda.min
ridgeBestLambda

# Run Ridge Regression algorithm
ridgeModel <- glmnet(ridgeX, ridgeY, alpha = 0, lambda = ridgeBestLambda, 
	family = "binomial")
coef(ridgeModel)

table(formattedData$X.class.)

# Load car package
library(car)

# Multicollinearity Detection
vif(forwardModel)
vif(backwardModel)

# Cook's Distance Plot Backward Model
plot(cooks.distance(backwardModel), main = "Cook's Distance Values in the 
	 Backward Model (Cutoff = 1", xlab = "Observation #", ylab = "Di Value")

# Linearity Check for the Backward Model
backwardModelLogOdds <- backwardModel$linear.predictors
cor(trainingData$X.bgr., backwardModelLogOdds)
cor(trainingData$X.bu., backwardModelLogOdds)
cor(trainingData$X.wbcc., backwardModelLogOdds)


plot(trainingData$X.bgr., backwardModelLogOdds, 
	 main = "Backward Model Log Odds vs. Blood Glucose (X.bgr.) Scatterplot", 
		xlab = "Blood Glucose (mg/dL)", ylab = "Log Odds")
plot(trainingData$X.bu., backwardModelLogOdds,  
	 main = "Backward Model Log Odds vs. Blood Urea (X.bu.) Scatterplot", 
		xlab = "Blood Urea (mg/dL)", ylab = "Log Odds")
plot(trainingData$X.wbcc., backwardModelLogOdds,
	 main = "Backward Model Log Odds vs. White Blood Cell Count (X.wbcc.) 
		Scatterplot", xlab = "White Blood Cell Count (cells/mm^3)", 
			ylab = "Log Odds")

# Sample Size Check

table(trainingData$X.class.)

(10 * 4) / (22 / 79) # 4 regressors, 22 / 79 CKD-classified observations

# Load the dplyr package
library(dplyr)

# Run the Backward Model on the test dataset
backwardModelProbabilities <- backwardModel %>% predict(testingData, 
	type = "response")
backwardModelPredictedClasses <- ifelse(backwardModelProbabilities < 0.5, 
	"ckd", "notckd")

# Testing data matrix
testDataMatrix <- data.matrix(testingData[, 2:25])

# Run the LASSO Model on the test dataset
lassoModelPredictedProbabilities <- predict(lassoModel, s = lassoBestLambda, 
	newx = testDataMatrix, type = "response")
lassoModelPredictedClasses <- ifelse(lassoModelPredictedProbabilities < 0.5, 
	"ckd", "notckd")

ridgeModelPredictedProbabilities <- predict(ridgeModel, s = ridgeBestLambda, 
	newx = testDataMatrix, type = "response")
ridgeModelPredictedClasses <- ifelse(ridgeModelPredictedProbabilities < 0.5, 
	"ckd", "notckd")

# Load the caret package
library(caret)

confusionMatrix(table(backwardModelPredictedClasses, testingData$X.class.))

confusionMatrix(table(lassoModelPredictedClasses, testingData$X.class.))

confusionMatrix(table(ridgeModelPredictedClasses, testingData$X.class.))

knitr::purl("CKD_Classification_Paper.Rmd", documentation = 0)
