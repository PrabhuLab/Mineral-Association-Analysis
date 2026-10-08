# #Association Analysis
# 
library(readxl)
library(stringr)
library(dplyr)
# 
library(arules)


#predictions and query

# Subset: USA

mars_export <- read.csv("./data/mars_data/mars_export.csv")
USA_Rules<-readRDS("./data/association_data/Allmin_USA_localities_subset_cleaned.rds")
USA_Rules@info$support


library(data.table)

Location_List<-mars_export$sample_id #added 4/1

# Some quick code to see which minerals are ignored for each location, added 4/1

for(i in 1:length(Location_List)) {
  Loc_Min<-strsplit(mars_export$mars_minerals_all[mars_export$sample_id==Location_List[i]],',')
  
  filter_items <- intersect(Loc_Min[[1]], itemLabels(USA_Rules))
  ignored_items<- setdiff(Loc_Min[[1]], filter_items)
  print(Location_List[i])
  print(ignored_items)
}
# -------------- Start Looping Through MAA for each Location ------------------#
for(j in 1:length(Location_List)) { # added 4/1 to automate all, changed string-key calls to be j indexed
#for(j in 1:1) { # If want FOR ROCKNEST ONLY comment out previous line. Also remove # for Section A and change output file name:
# If want FOR ROCKNEST ONLY comment out previous line. Also remove # for Section A and change output file name:
    Loc_Min<-strsplit(mars_export$mars_minerals_all[mars_export$sample_id==Location_List[j]],',')
    
    filter_items <- intersect(Loc_Min[[1]], itemLabels(USA_Rules))
    for(k in 1:length(filter_items)) {
      removed_mineral<-filter_items[k]
      sublist<-filter_items[filter_items != removed_mineral]
      removed_mineral_as_list<-filter_items[filter_items == removed_mineral] # leaves one out
   
        
      MinCombo<-combn(x = sublist,m = 3)
      if (ncol(MinCombo)==0) { # Just adding this in case a location only has two minerals occurring
          next
      }
      
      rules.sub.AP <- subset(USA_Rules, subset = lhs %ain% MinCombo[,1] & rhs %in% removed_mineral_as_list)
      
      inspect(rules.sub.AP)
      
      
      for(i in 1:ncol(MinCombo))
      { 
          #rules.sub.AP2 <- subset(USA_Rules, subset = lhs %ain% MinCombo[,i] & !rhs %in% sublist)
          #TEMP LINE, remove and replace with above line for normal MAA
          rules.sub.AP2 <- subset(USA_Rules, subset = lhs %ain% MinCombo[,i] & rhs %in% removed_mineral_as_list)
          rules.sub.AP <- union(rules.sub.AP,rules.sub.AP2)
          #inspect(rules.sub.AP)
      }
      
      inspect(rules.sub.AP)
      # The following line breaks if there are no predictions because lhs and rhs return empty brackets, but @quality returns no values (support, conf, lift, etc)
      if(nrow(rules.sub.AP@quality)==0) {
          failed_attempt<-data.table(locality = mars_export$sample_id[j], mineral = removed_mineral)
          write.table(failed_attempt, './data/mars_data/failed_proof_of_concept.csv', append=TRUE)
          next
      }
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
      
      write.table(mars_predictions, './data/mars_data/proof_of_concept.csv', append=TRUE)
    }
}
