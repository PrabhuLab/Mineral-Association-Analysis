# #Association Analysis
# 
library(readxl)
library(stringr)
library(dplyr)
library(arules)

#_______________________________________________

#predictions and query

# Subset: USA
mars_export <- read.csv("./data/mars_data/mars_export.csv")
USA_Rules<-readRDS("./data/association_data/Allmin_USA_localities_subset_cleaned.rds")
USA_Rules@info$support


library(data.table)

Location_List<-mars_export$sample_id #added 4/1

for(i in 1:length(Location_List)) {
  Loc_Min<-strsplit(mars_export$mars_minerals_all[mars_export$sample_id==Location_List[i]],',')
  
  filter_items <- intersect(Loc_Min[[1]], itemLabels(USA_Rules))
  ignored_items<- setdiff(Loc_Min[[1]], filter_items)
  print(Location_List[i])
  print(ignored_items)
}
# -------------- Start Looping Through MAA for each Location ------------------#

for(j in 1:length(Location_List)) { # added 4/1 to automate all, changed string-key calls to be j indexed
  
  Loc_Min<-strsplit(mars_export$mars_minerals_all[mars_export$sample_id==Location_List[j]],',')
  
  filter_items <- intersect(Loc_Min[[1]], itemLabels(USA_Rules))
  
  MinCombo<-combn(x = filter_items,m = 3)
  if (ncol(MinCombo)==0) { # Just adding this in case a location only has two minerals occurring
    next
  }
  # rules.sub.AP <- subset(USA_Rules, subset = ((lhs %pin% "Augite") | (lhs %pin% "Magnetite") | (lhs %pin% "Sanidine"))) # Checking if the union is working properly
  rules.sub.AP <- subset(USA_Rules, subset = lhs %ain% MinCombo[,1] & !rhs %in% filter_items) # This line needs to stay so there is something to union
  # rules.sub.AP2 <- subset(USA_Rules, subset = lhs %ain% MinCombo[,32] & !rhs %in% Loc_Min[[1]]) # Can be used to check if the union is working properly
  
  inspect(rules.sub.AP)
  
  
  for(i in 1:ncol(MinCombo))
  { 
    rules.sub.AP2 <- subset(USA_Rules, subset = lhs %ain% MinCombo[,i] & !rhs %in% filter_items)
    rules.sub.AP <- union(rules.sub.AP,rules.sub.AP2)
  }
  
  inspect(rules.sub.AP)
  rules_pred<-data.frame(lhs = labels(lhs(rules.sub.AP)),rhs = labels(rhs(rules.sub.AP)), 
                         rules.sub.AP@quality)
  
  # The old aggregation to be replaced one returning full rules for 1) highest confidence 2) highest lift
  # mars_predictions<-aggregate(x = rules_pred[,4:6], by = list(rules_pred$rhs), FUN = max) 
  
  rules_pred_copy <- data.table::copy(rules_pred) # Dataframe duplicate
  max_conf <- filter(rules_pred, confidence==max(confidence), .by=rhs) #all best confidences
  max_lift <- filter(rules_pred_copy, lift==max(lift), .by=rhs) # all best lifts
  mars_predictions <- rbind(max_conf, max_lift) # Put both results in one df
  mars_predictions <- sort_by(mars_predictions, mars_predictions[,2]) # Sorts by rhs
  
  mars_predictions$locality <- mars_export$sample_id[j] # Attaches the locality to the reading
  
  Loc_Min
  
  mars_export$mars_minerals_all[mars_export$sample_id==Location_List[j]]
  write.table(mars_predictions, './data/mars_data/mars_rules_pred.csv', append=TRUE)
  
}
