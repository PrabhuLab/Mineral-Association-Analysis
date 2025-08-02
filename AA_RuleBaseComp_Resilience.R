# Association Analysis Rules Comparison - Resilience
library(arules)

library('devtools')
install_github('brooksandrew/Rsenal',force=TRUE)
library('Rsenal')
?Rsenal::arulesApp

library('arules') 



data('AdultUCI') # data.frame
arulesApp(AdultUCI)

USA_Rules<-readRDS("~/Downloads/Allmin_USA_localities_subset_cleaned (1).rds")
VA_Rules<-readRDS("~/Downloads/VA_Rules0709.rds")


inspect(VA_Rules)

interesect_mins_without_uniqueLabels<-intersect(USA_Rules,VA_Rules)

inspect(interesect_mins_without_uniqueLabels)

all.equal(itemLabels(VA_Rules), itemLabels(USA_Rules))

itemLabels <- union(itemLabels(VA_Rules), itemLabels(USA_Rules))

VA_Rules_fixed <- new("rules",
                  lhs = recode(lhs(VA_Rules), itemLabels = itemLabels), 
                  rhs = recode(rhs(VA_Rules), itemLabels = itemLabels)
)

USA_Rules_fixed <- new("rules",
                  lhs = recode(lhs(USA_Rules), itemLabels = itemLabels), 
                  rhs = recode(rhs(USA_Rules), itemLabels = itemLabels)
)

itemLabels(VA_Rules_fixed)

union(USA_Rules_fixed, VA_Rules_fixed)

test_intersect<-intersect(USA_Rules_fixed, VA_Rules_fixed)

inspect(test_intersect)
inspect(tail(test_intersect))

set_of_minerals<-rowSums(as(rhs(VA_Rules), "ngCMatrix"))

# Generate full sets of LHS for both VA and USA

library('arules') 

VA_lhs<-as(lhs(VA_Rules), "list")
VA_rhs<-as(rhs(VA_Rules), "list")

VA_Rules_df<-data.frame(lhs = as.data.frame.character(VA_lhs),rhs= as.data.frame.character(VA_rhs))

USA_lhs<-as(lhs(USA_Rules), "list")
USA_rhs<-as(rhs(USA_Rules), "list")
USA_Rules_df<-data.frame(lhs = as.data.frame.character(USA_lhs),rhs= as.data.frame.character(USA_rhs))



VA_Rules_df_try<-VA_Rules_df[60:70,]

library(dplyr)

VA_Full_LHS_Set<-VA_Rules_df %>%
  group_by(VA_rhs) %>%
  summarise(combined_lhs = paste(unique(unlist(VA_lhs)), collapse = ", "))

USA_Full_LHS_Set<-USA_Rules_df %>%
  group_by(USA_rhs) %>%
  summarise(combined_lhs = paste(unique(unlist(USA_lhs)), collapse = ", "))


#____________________________________


#Processing the rules
write(USA_Rules, file = "U_Rules_Subset.csv", sep = ",")
write(rules_support_0_008_Full, file = "U_Rules_Full.csv", sep = ",")

#unlink("U_Rules_Subset.csv") # tidy up
#unlink("U_Rules_Full.csv") # tidy up



#rules_support_0_008 <- rules_support_0_008[!is.redundant(rules_support_0_008)]
rules.sorted = sort(VA_Rules, by="lift")
rules.sorted <- rules.sorted[!is.redundant(rules.sorted)]

rules.sorted.full = sort(USA_Rules, by="lift")
rules.sorted.full <- rules.sorted.full[!is.redundant(rules.sorted.full)]

write(rules.sorted, file = "U_Rules_Subset_SortedByLift.csv", sep = ",")
write(rules.sorted.full, file = "U_Rules_Full_SortedByLift.csv", sep = ",")

rules_df_subset<-DATAFRAME(rules.sorted)
rules_df_full<-DATAFRAME(rules.sorted.full)

#Remove {}
rules_df_subset$LHS<-gsub("{","",as.character(rules_df_subset$LHS),fixed = T)
rules_df_subset$LHS<-gsub("}","",as.character(rules_df_subset$LHS),fixed = T)
rules_df_subset$RHS<-gsub("{","",as.character(rules_df_subset$RHS),fixed = T)
rules_df_subset$RHS<-gsub("}","",as.character(rules_df_subset$RHS),fixed = T)

rules_df_full$LHS<-gsub("{","",as.character(rules_df_full$LHS),fixed = T)
rules_df_full$LHS<-gsub("}","",as.character(rules_df_full$LHS),fixed = T)
rules_df_full$RHS<-gsub("{","",as.character(rules_df_full$RHS),fixed = T)
rules_df_full$RHS<-gsub("}","",as.character(rules_df_full$RHS),fixed = T)

subset_RHS<-levels(as.factor(rules_df_subset$RHS))
levels(as.factor(rules_df_full$RHS))
levels(as.factor(rules_needed_RHS$RHS))

#rules_full_check<-rules_df_full
#Remove uncommon RHS rules since they cant be compared 
# These rules have no resilience
rules_needed_RHS<-rules_df_full[rules_df_full$RHS %in% subset_RHS,]

#rules_needed_RHS[duplicated(rules_needed_RHS[,1:2]),]

levels(as.factor(rules_needed_RHS$LHS))
levels(as.factor(rules_df_subset$LHS))

aaa<-rules_needed_RHS[rules_needed_RHS$LHS %in% levels(as.factor(rules_df_subset$LHS)),]

rules_needed_RHS$LHS_split<-strsplit(x = rules_needed_RHS$LHS,',')
#rules_df_full$LHS_split<-strsplit(x = rules_df_full$LHS,',')
rules_df_subset$LHS_split<-strsplit(x = rules_df_subset$LHS,',')

#aaa<-merge(x = rules_df_subset,y = rules_needed_RHS, by = "RHS")



rules_sub_vec<-list(rules_df_subset$LHS_split)
rules_sub_vec$RHS<-rules_df_subset$RHS

rules_full_vec<-list(rules_needed_RHS$LHS_split)
rules_full_vec$RHS<-rules_needed_RHS$RHS

compare_rules <- data.frame(LHS_Full = NA, LHS_Subset = NA, RHS = NA, Resilience = NA)
lhs_full <- NA
lhs_sub <- NA
rhs_all <- NA
resilience <- NA
j<-2

for(j in 1:length(rules_full_vec[[1]]))
{
  for(i in 1:length(rules_sub_vec[[1]]))  
  {
    if(rules_sub_vec[[2]][[i]]==rules_full_vec[[2]][[j]])
    {
      if(length(rules_full_vec[[1]][[j]]) >= length(rules_sub_vec[[1]][[i]]))
      {
        if(length(intersect(rules_sub_vec[[1]][[i]],rules_full_vec[[1]][[j]]))>1)
        {
          #print("Full set LHS")
          #print(rules_full_vec[[1]][j])
          lhs_full<-c(lhs_full,rules_full_vec[[1]][j])        
          #print("Subset LHS")
          lhs_sub<-c(lhs_sub,rules_sub_vec[[1]][i])        
          print(i)
          rhs_all<-c(rhs_all,rules_sub_vec[[2]][i])
          print(j)
          resilience<-c(resilience,(length(rules_full_vec[[1]][[j]])/length(rules_sub_vec[[1]][[i]]))*1/(length(intersect(rules_full_vec[[1]][[j]],rules_sub_vec[[1]][[i]])))*abs(rules_df_subset$lift[i]-rules_needed_RHS$lift[j]))
          
        }
      }
      
    }
  }
  #  if(j>20)
  #  {break}
}

compare_rules <- cbind(LHS_Full = lhs_full, LHS_Subset = lhs_sub, RHS = rhs_all, Resilience = resilience)

compare_rules <- as.data.frame(compare_rules)
compare_rules<-compare_rules[-1,]

compare_rules <- apply(compare_rules,2,as.character)

compare_rules[,4]<-as.numeric(compare_rules[,4])

write.csv(compare_rules,"~/USA_vs_VA_Rules_Stability.csv")

RHS_unique_iter<-unique(compare_rules[,3])


for(i in 1:length(RHS_unique_iter))
{
  print(i)
  print(RHS_unique_iter[i])
  print(max(compare_rules[which(compare_rules[,3]==RHS_unique_iter[i]),4]))
  print(min(compare_rules[which(compare_rules[,3]==RHS_unique_iter[i]),4]))  
}

max(compare_rules[which(compare_rules[,3]==RHS_unique_iter[i]),4])
min(compare_rules[which(compare_rules[,3]==RHS_unique_iter[i]),4])

length(intersect(rules_full_vec[[1]][[345]],rules_sub_vec[[1]][[3]]))

length(intersect(rules_sub_vec[[1]][[2]],rules_sub_vec[[1]][[3]]))


class(rules_sub_vec[[1]][[2]])

rules_sub_vec[[1]][[2]]


unique(compare_rules[,3])
