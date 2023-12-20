setwd(dir_ltla)

using("caTools")
using("gbm") # gradient boosting
using("xgboost")
using("pROC") # ROC

############################################################
#  Machine Learning analyses and ROC curves
############################################################

data_size = 2000000 # use only 2M data instead of 7M
split_ratio = 1/2 # half are training data, half validation

# reduce size of data
data <- read_csv("table_riskScores.csv.gz") %>%
  filter(!is.na(bg_rate_cases_app)) %>%
  slice_sample(n=data_size)# extract useful predictors
data <- data %>%
  summarise(
    individualId = 1:data_size,
    metadata = number_exposures,
    daily_rate_cases_app = bg_rate_cases_app,
    bin1_count = riskScore_1.1,
    bin2_count = riskScore_1.3,
    bin3_count = riskScore_1.6,
    bin4_count = riskScore_1.8,
    bin5_count = riskScore_2.2,
    bin6_count = riskScore_3.3,
    bin7_count = riskScore_4.4,
    bin8_count = riskScore_5.6,
    duration = duration,
    maxRiskScore = maxRiskScore,
    cumRiskScore = cumRiskScore,
    avgRiskScore = avgRiskScore,
    region = region,
    rural_urban_score = RuralUrbanScore,
    time_since_pandemic = diff_days(date,"2020-03-11"),
    peak_duration = peak_duration,
    peak_day = peak_day,
    result = positive
  )
rm(table_outcome_riskScores)

# returns numeric matrix
convertDataToInputMatrix = function(data){
  return(as.matrix(sapply(data, as.numeric)))
}

all_names = colnames(data)

risk_features = all_names[grepl( "riskScore*", all_names)]
basic_features = c(risk_features, c("duration", "cumRiskScore", "maxRiskScore", "avgRiskScore"
))
extra_features = c("daily_rate_cases_app", "region", "rural_urban_score", "time_since_pandemic", "peak_duration", "peak_day")

split = sample.split(data$result, SplitRatio = split_ratio)
train_data = subset(data, split == TRUE)
eval_data = subset(data, split == FALSE)
NROW(train_data)
NROW(eval_data)
rm(data)

# turn regions to numbers ordered by gross earnings
regions_to_numbers<-c(South_East=1,London=2,East_of_England=3,South_West=4,West_Midlands=5,North_West=6,East_Midlands=7,Wales=8,Yorkshire_and_The_Humber=9,North_East=10)
train_data <- train_data %>% mutate(region=regions_to_numbers[region])
eval_data <- eval_data %>% mutate(region=regions_to_numbers[region])

roclist = list()

# create ROCs using only basic individual features
vRisk<-c(100,120,140,160,200,300,400,500)

newRoc = roc(eval_data$result, eval_data$maxRiskScore/max(eval_data$maxRiskScore))
newRoc$description = "Max Risk Score"
roclist = append(roclist, list(maxRisk = newRoc ))

newRoc = roc(eval_data$result, eval_data$cumRiskScore/max(eval_data$cumRiskScore))
newRoc$description = "Cumulative Risk Score"
roclist = append(roclist, list(cumRisk = newRoc))

newRoc = roc(eval_data$result, eval_data$duration/max(eval_data$duration))
newRoc$description = "Cumulative Duration"
roclist = append(roclist, list(duration = newRoc))

newRoc = roc(eval_data$result, eval_data$avgRiskScore/max(eval_data$avgRiskScore))
newRoc$description = "Mean Risk Score"
roclist = append(roclist, list(meanRisk = newRoc))



# GBM with basic
fitGbm = gbm(result ~ .,  data = train_data[, c("result", basic_features)])
pred = predict(fitGbm, eval_data)
newRoc = roc(eval_data$result, pred)
newRoc$description = "GBM with basic features"
roclist = append(roclist, list(GBM = newRoc))
fitGbm = NULL #delete model to save space

# GBM with basic + extra
fitGbm = gbm(result ~ .,  data = train_data[, c("result", basic_features, extra_features)])
pred = predict(fitGbm, eval_data)
newRoc = roc(eval_data$result, pred)
newRoc$description = "GBM with all features"
roclist = append(roclist, list(GBM_Extra = newRoc))
fitGbm = NULL #delete model to save space


#xgboost with 10, 100, and 400 rounds (basic features)

basicTrainMatrix = convertDataToInputMatrix(train_data[, basic_features])
basicEvalMatrix = convertDataToInputMatrix(eval_data[, basic_features])

fitXg10 = xgboost(data = basicTrainMatrix , label = train_data$result, nrounds = 10)
pred = predict(fitXg10, newdata = basicEvalMatrix)
newRoc = roc(eval_data$result, pred)
newRoc$description = "XGB10 with basic features"
roclist = append(roclist, list(XGB10 = newRoc))
fitXg10 = NULL

fitXg100 = xgboost(data = basicTrainMatrix, label = train_data$result, nrounds = 100)
pred = predict(fitXg100, newdata = basicEvalMatrix)
newRoc = roc(eval_data$result, pred)
newRoc$description = "XGB100 with basic features"
roclist = append(roclist, list(XGB100 = newRoc))
fitXg100 = NULL

fitXg400 = xgboost(data = basicTrainMatrix, label = train_data$result, nrounds = 400)
pred = predict(fitXg400, newdata = basicEvalMatrix)
newRoc = roc(eval_data$result, pred)
newRoc$description = "XGB400 with basic features"
roclist = append(roclist, list(XGB400 = newRoc))
fitXg400 = NULL

#xgboost with 10, 100, and 400 rounds (basic and extra features)

extraTrainMatrix = convertDataToInputMatrix(train_data[, c(basic_features, extra_features)])
extraEvalMatrix = convertDataToInputMatrix(eval_data[, c(basic_features, extra_features)])

fitXg10 = xgboost(data = extraTrainMatrix , label = train_data$result, nrounds = 10)
pred = predict(fitXg10, newdata = extraEvalMatrix)
newRoc = roc(eval_data$result, pred)
newRoc$description = "XGB10 with all features"
roclist = append(roclist, list(XGB10_Extra = newRoc))
fitXg10 = NULL

fitXg100 = xgboost(data = extraTrainMatrix, label = train_data$result, nrounds = 100)
pred = predict(fitXg100, newdata = extraEvalMatrix)
newRoc = roc(eval_data$result, pred)
newRoc$description = "XGB100 with all features"
roclist = append(roclist, list(XGB100_Extra = newRoc))
fitXg100 = NULL

fitXg400 = xgboost(data = extraTrainMatrix, label = train_data$result, nrounds = 400)
pred = predict(fitXg400, newdata = extraEvalMatrix)
newRoc = roc(eval_data$result, pred)
newRoc$description = "XGB400 with all features"
roclist = append(roclist, list(XGB400_Extra = newRoc))
fitXg400 = NULL

#logistic regression (basic features)
fitGlm = glm(result ~ . , data = train_data[, c("result", basic_features)], family = "binomial")
pred = predict(fitGlm, eval_data[, basic_features], type = "response")
newRoc = roc(eval_data$result, pred)
newRoc$description = "Logistic Regression with basic features"
roclist = append(roclist, list(LogisticRegression = newRoc))
fitGlm = NULL

#logistic regression (basic and extra features)
fitGlm = glm(result ~ . , data = train_data[, c("result", basic_features, extra_features)], family = "binomial")
pred = predict(fitGlm, eval_data[, c(basic_features, extra_features)], type = "response")
newRoc = roc(eval_data$result, pred)
newRoc$description = "Logistic Regression with all features"
roclist = append(roclist, list(LogisticRegression_Extra = newRoc))
fitGlm = NULL

# XGB10 with basic features and some of the extra features

# "peak_duration", "peak_day"
features = c(basic_features, "peak_duration", "peak_day")
fitXg10 = xgboost(data = convertDataToInputMatrix(train_data[,features]) , label = train_data$result, nrounds = 10)
pred = predict(fitXg10, newdata = convertDataToInputMatrix(eval_data[,features]) )
newRoc = roc(eval_data$result, pred)
newRoc$description = "XGB10 with basic features, peak duration, and peak day"
roclist = append(roclist, list(XGB10_Timing = newRoc))
fitXg10 = NULL

#"time_since_pandemic"
features = c(basic_features, "time_since_pandemic")
fitXg10 = xgboost(data = convertDataToInputMatrix(train_data[,features]) , label = train_data$result, nrounds = 10)
pred = predict(fitXg10, newdata = convertDataToInputMatrix(eval_data[,features]) )
newRoc = roc(eval_data$result, pred)
newRoc$description = "XGB10 with basic features and date"
roclist = append(roclist, list(XGB10_Date = newRoc))
fitXg10 = NULL

# "region", "rural_urban_score"
features = c(basic_features, "region", "rural_urban_score")
fitXg10 = xgboost(data = convertDataToInputMatrix(train_data[,features]) , label = train_data$result, nrounds = 10)
pred = predict(fitXg10, newdata = convertDataToInputMatrix(eval_data[,features]) )
newRoc =  roc(eval_data$result, pred)
newRoc$description = "XGB10 with basic features and region"
roclist = append(roclist, list(XGB10_Location = newRoc))
fitXg10 = NULL

# "daily_rate_cases_app"
features = c(basic_features, "daily_rate_cases_app")
fitXg10 = xgboost(data = convertDataToInputMatrix(train_data[,features]) , label = train_data$result, nrounds = 10)
pred = predict(fitXg10, newdata = convertDataToInputMatrix(eval_data[,features]) )
newRoc = roc(eval_data$result, pred)
newRoc$description = "XGB10 with basic features and background cases"
roclist = append(roclist, list(XGB10_Background = newRoc))
fitXg10 = NULL

# duration+background risk
features = c("duration", "daily_rate_cases_app")
fitXg10 = xgboost(data = convertDataToInputMatrix(train_data[,features]) , label = train_data$result, nrounds = 10)
pred = predict(fitXg10, newdata = convertDataToInputMatrix(eval_data[,features]) )
newRoc = roc(eval_data$result, pred)
newRoc$description = "XGB10 with duration and background cases"
roclist = append(roclist, list(XGB10_duration_bg = newRoc))
fitXg10 = NULL

#runs with duration and date only

fitGbm = gbm(result ~ .,  data = train_data[, c("result", c("time_since_pandemic","duration"))])
pred = predict(fitGbm, eval_data)
newRoc = roc(eval_data$result, pred)
newRoc$description = "GBM with duration and date"
roclist = append(roclist, list(GBM_duration_date = newRoc))
fitGbm = NULL #delete model to save space

extraTrainMatrix = convertDataToInputMatrix(train_data[, c("time_since_pandemic","duration")])
extraEvalMatrix = convertDataToInputMatrix(eval_data[, c("time_since_pandemic","duration")])

fitXg10 = xgboost(data = extraTrainMatrix , label = train_data$result, nrounds = 10)
pred = predict(fitXg10, newdata = extraEvalMatrix)
newRoc = roc(eval_data$result, pred)
newRoc$description = "XGB10 with duration and date"
roclist = append(roclist, list(XGB10_duration_date = newRoc))
fitXg10 = NULL

fitXg100 = xgboost(data = extraTrainMatrix, label = train_data$result, nrounds = 100)
pred = predict(fitXg100, newdata = extraEvalMatrix)
newRoc = roc(eval_data$result, pred)
newRoc$description = "XGB100 with duration and date"
roclist = append(roclist, list(XGB100_duration_date = newRoc))
fitXg100 = NULL

fitXg400 = xgboost(data = extraTrainMatrix, label = train_data$result, nrounds = 400)
pred = predict(fitXg400, newdata = extraEvalMatrix)
newRoc = roc(eval_data$result, pred)
newRoc$description = "XGB400 with duration and date"
roclist = append(roclist, list(XGB400_duration_date = newRoc))
fitXg400 = NULL

#logistic regression (basic features)
fitGlm = glm(result ~ . , data = train_data[, c("result", c("time_since_pandemic","duration"))], family = "binomial")
pred = predict(fitGlm, eval_data[, c("time_since_pandemic","duration")], type = "response")
newRoc = roc(eval_data$result, pred)
newRoc$description = "Logistic regression with duration and date"
roclist = append(roclist, list(LogisticRegression_duration_date = newRoc))
fitGlm = NULL


#sort by auc
aucs = sapply(roclist, function(x){x$auc})
indices = sort(aucs, decreasing = TRUE, index.return = TRUE)$ix
roclist = roclist[indices]

# build AUC table
roundedAucs = sapply(roclist,function(x){round(x$auc, digits = 3)})
descriptions = sapply(roclist,function(x){x$description}) # assumes each ROC object has had a description added
table_AUCs = data.table(Name = names(roclist), Description = descriptions, AUC = roundedAucs)
SaveTable(table_AUCs)




# save ROC plot, and ROC list to RData file
logfile = paste0("OutputROCcurves.RData")
logfile = str_replace_all(logfile, " ", "-")
logfile = str_replace_all(logfile, ":", "-")
save(roclist, file = logfile)

##plot the ROCs
#rocplot = ggroc(roclist) + labs(colour = "Classifier") + theme_classic()
#print(rocplot)

# function to sort by AUC and add AUC to list names
sortAndConcatenate = function(roclist){
  aucs = sapply(roclist, function(x){x$auc})
  roundedAucs =round(aucs, digits = 3)
  names(roclist) = mapply(function(oldName, auc) {paste0(oldName, " (AUC=", auc, ")")}, 
                          names(roclist), roundedAucs )
  indices = sort(aucs, decreasing = TRUE, index.return = TRUE)$ix
  return(roclist[indices])
}


roclist_main <- list()
roclist_main[["cumulative risk score"]] <- roclist[["cumRisk"]]
roclist_main[["maximum risk score"]] <- roclist[["maxRisk"]]
roclist_main[["average risk score"]] <- roclist[["meanRisk"]]
#roclist_main[["maximum proximity score"]] <- roclist[["maxProximity"]]
#roclist_main[["average proximity score"]] <- roclist[["meanProximity"]]
roclist_main[["duration"]] <- roclist[["duration"]]
roclist_main[["best ML from app predictors \n (XGB10)"]] <- roclist[["XGB10"]]
roclist_main[["best ML from duration \n & background risk only"]] <- roclist[["XGB10_duration_bg"]]
roclist_main[["best ML from app predictors\n and extra predictors"]] <- roclist[["XGB10_Extra"]]

roclist_main = sortAndConcatenate(roclist_main)

#plot the ROCs for the main text
rocplot_main = ggroc(roclist_main, legacy.axes = TRUE) + labs(colour = "Classifier") + theme_classic() + 
  guides(color = guide_legend(byrow = TRUE)) + coord_cartesian(expand = FALSE) +
  ylab("Fraction of contacts testing positive\nclassified as at risk (sensitivity)") + xlab("Fraction of contacts not testing positive\nclassified as at risk (1-specificity)") + 
  theme(legend.position = c(0.75, 0.3),
        legend.text = element_text(size=9),
        legend.title = element_blank(),
        legend.spacing.y = unit(+1, "pt"),
        plot.margin = margin(t = 0.5,  # Top margin
                             r = 1.5,  # Right margin
                             b = 0.5,  # Bottom margin
                             l = 0.5,  # Left margin
                             unit = "cm"))  
SaveFigure(rocplot_main,height=12,width=14)

rm(roclist_main)
rm(rocplot_main)


roclist_basic <- list()
roclist_basic[["XGB10"]] <- roclist[["XGB10"]]
roclist_basic[["XGB100"]] <- roclist[["XGB100"]]
roclist_basic[["XGB400"]] <- roclist[["XGB400"]]
roclist_basic[["Logistic Regression"]] <- roclist[["LogisticRegression"]]
#roclist_basic[["Logistic Regression (Quadratic)"]] <- roclist[["LogisticRegressionQuadratic"]]
roclist_basic[["GBM"]] <- roclist[["GBM"]]

roclist_basic = sortAndConcatenate(roclist_basic)

#plot the ROCs for all algorithms with the basic features
rocplot_basic = ggroc(roclist_basic, legacy.axes = TRUE) + labs(colour = "Classifier") + theme_classic() + coord_cartesian(expand = FALSE) +
  ylab("Fraction of contacts testing positive\nclassified as at risk") + xlab("Fraction of contacts not testing positive\nclassified as at risk")
SaveFigure(rocplot_basic,height=10,width=15)

rm(roclist_basic)
rm(rocplot_basic)

roclist_extra <- list()
roclist_extra[["XGB10"]] <- roclist[["XGB10_Extra"]]
roclist_extra[["XGB100"]] <- roclist[["XGB100_Extra"]]
roclist_extra[["XGB400"]] <- roclist[["XGB400_Extra"]]
roclist_extra[["Logistic Regression"]] <- roclist[["LogisticRegression_Extra"]]
#roclist_extra[["Logistic Regression (Quadratic)"]] <- roclist[["LogisticRegressionQuadratic_Extra"]]
roclist_extra[["GBM"]] <- roclist[["GBM_Extra"]]

roclist_extra = sortAndConcatenate(roclist_extra)

#plot the ROCs for all algorithms with the basic and extra features
rocplot_extra = ggroc(roclist_extra, legacy.axes = TRUE) + labs(colour = "Classifier") + theme_classic() + coord_cartesian(expand = FALSE) +
  ylab("Fraction of contacts testing positive\nclassified as at risk") + xlab("Fraction of contacts not testing positive\nclassified as at risk")
SaveFigure(rocplot_extra,height=10,width=15)

rm(roclist_extra)
rm(rocplot_extra)

roclist_features <- list()
roclist_features[["Basic features only"]] <- roclist[["XGB10"]]
roclist_features[["Day of week and max duration"]] <- roclist[["XGB10_Timing"]]
roclist_features[["Date"]] <- roclist[["XGB10_Date"]]
roclist_features[["Region and rural/urban score"]] <- roclist[["XGB10_Location"]]
roclist_features[["Background infection rate"]] <- roclist[["XGB10_Background"]]
roclist_features[["All"]] <- roclist[["XGB10_Extra"]]

roclist_features = sortAndConcatenate(roclist_features)

#plot the ROCs for XGB10 with basic features and some of the extra features
rocplot_features = ggroc(roclist_features, legacy.axes = TRUE) + labs(colour = "Classifier") + theme_classic() + coord_cartesian(expand = FALSE) +
  ylab("Fraction of contacts testing positive\nclassified as at risk") + xlab("Fraction of contacts not testing positive\nclassified as at risk")
SaveFigure(rocplot_features,height=10,width=18)

rm(roclist_features)
rm(rocplot_features)

