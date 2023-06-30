#Association Analysis

library(readxl)
MED_Ages_Locality <- read_excel("./data/association_data/MED_Ages_Locality.xlsx")
MED_export <- read.delim("./data/association_data/MED_export.csv")
#MED_export<-read.csv("~/Downloads/MED_export.csv")
library(stringr)
MED_export<-MED_export[which(str_count(MED_export$LLN,",")>2),]

#MED_export$MED_elements_all
MED_export_bottom_level<-MED_export[which(MED_export$bottom_level==1),]


# Matches Case

MED_U_subset<-MED_export[which(str_detect(MED_export$MED_elements_all,"U")),]

MED_USA_subset<-MED_export[which(str_detect(MED_export$LLN,"USA")),]

#Remove all localities for Mindat locality id = 0
MED_Ages_Locality_sub<-MED_Ages_Locality[which(MED_Ages_Locality$mindat_id!=0),]

MED_Age_merged<-merge(MED_export,MED_Ages_Locality_sub, by = "mindat_id")

MED_OlderThan2500Ma<-MED_Age_merged[which(MED_Age_merged$`Max Age`>2500),]

MED_YoungerThan540Ma<-MED_Age_merged[which(MED_Age_merged$`Max Age`<540),]

library(dplyr)
MED_from2500MaTo540Ma<-MED_Age_merged[which(between(MED_Age_merged$`Max Age`,540,2500)),]

MED_OlderThan2500Ma<-MED_OlderThan2500Ma[which(str_count(MED_OlderThan2500Ma$LLN,",")>2),]


#__________________________________________________________
# Run association analysis
U.associated.mins <- strsplit(as.character(MED_export_bottom_level$MED_minerals_all),',')

#length(U.associated.mins[[1]])
MED_export_bottom_level$mindat_url[which.max(lengths(U.associated.mins))]
which.max(lengths(U.associated.mins))
max(lengths(U.associated.mins))
#U_associated_localities_minerals <- data.frame(locality=rep(MED_OlderThan2500Ma$mindat_id,sapply(U.associated.mins, FUN=length)), mineral=unlist(U.associated.mins),stringsAsFactors = F)
U_associated_localities_minerals <- data.frame(locality=rep(MED_OlderThan2500Ma$mindat_id, sum(sapply(U.associated.mins, FUN=length))), mineral=unlist(U.associated.mins),stringsAsFactors = F)

# remove leading & trailing spaces from mineral names
U_associated_localities_minerals$mineral <- sapply(U_associated_localities_minerals$mineral,function(x) {sub("^\\s+|^\\s+$","",x)})

U_associated.loc.min <- U_associated_localities_minerals
U_associated.loc.min$locality <- as.character(U_associated.loc.min$locality)
U_associated.loc.min.mat <- as.data.frame.matrix(table(U_associated.loc.min[,1:2]))

#U_associated.loc.min.mat

U_associated_Min_Mat<-sapply(X = U_associated.loc.min.mat,FUN = function(x){replace(x, x > 0,1)})
#U_associated.loc.min.mat[which(U_associated.loc.min.mat != 0)] <- 1
dim(U_associated_Min_Mat)

U_associated_Min_Mat <- as.data.frame(U_associated_Min_Mat)
rownames(U_associated_Min_Mat) <- rownames(U_associated.loc.min.mat)

U_associated.Min.Mat<-as.data.frame(unclass(U_associated_Min_Mat))

#Logical Matrix Needed if only presence is to be considered. 
U_associated.Min.Mat<-sapply(U_associated.loc.min.mat, as.logical)

library(arules)
rules_support_0_008 <- apriori(data = U_associated.Min.Mat, parameter = list(supp = 0.002, conf = 0.7,target = "rules",maxtime = 0,maxlen = 4,minlen = 1))



saveRDS(rules_support_0_008, "Allmin_U_localities_subset.rds")
#rules_support_0_0001 <- apriori(data = U_associated.Min.Mat, parameter = list(supp = 0.0001, conf = 0.7,target = "rules",maxtime = 0,maxlen = 30))

#rules.sub.AP <- subset(rules_support_0_008, subset = ((lhs %pin% "Uraninite")))

inspect(head(rules_support_0_008))

#rules_support_0_008<-readRDS(file = "~/Allmin_OlderThan2500Ma_Rules.rds")
rules <- rules_support_0_008[!is.redundant(rules_support_0_008)]
saveRDS(rules, "Allmin_USA_localities_subset_cleaned.rds")

inspect(head(rules, n = 10, by = "lift"))
rules_by_lift <- sort(rules, by = "lift")


library(arulesViz)
plot(rules, method = "grouped")
plot(rules, method = "paracoord")
plot(rules, method = "graph")
plot(rules, method = "matrix3D")

summary(str_detect(MED_export$MED_minerals_all,"Cummingtonite")&str_detect(MED_export$MED_minerals_all,"Schreibersite")&str_detect(MED_export$MED_minerals_all,"Violarite")&str_detect(MED_export$MED_minerals_all,"Daubreelite",negate = TRUE))
summary(str_detect(MED_export$MED_minerals_all,"Schreibersite"))
summary(str_detect(MED_export$MED_minerals_all,"Violarite"))
summary(str_detect(MED_export$MED_minerals_all,"Daubreelite",negate = TRUE))

Pred<-MED_export[which(str_detect(MED_export$MED_minerals_all,"Siderophyllite")&str_detect(MED_export$MED_minerals_all,"Petalite")&str_detect(MED_export$MED_minerals_all,"Fluor-liddicoatite",negate = TRUE)),]

Pred<-MED_export[which(str_detect(MED_export$MED_minerals_all,"Opal")&str_detect(MED_export$MED_minerals_all,"Pyrolusite")&str_detect(MED_export$MED_minerals_all,"Halotrichite",negate = TRUE)),]



#_______________________________________________

#predictions and query

# Subset: USA
# Locality: Tecopa Basin 
# Mindat ID: 255486

USA_Rules<-readRDS("./data/association_data/Allmin_USA_localities_subset_cleaned.rds")
USA_Rules@info$support

Loc_Min<-strsplit(MED_export$MED_minerals_all[MED_export$mindat_id==255486],',')

MinCombo<-combn(x = Loc_Min[[1]],m = 3)


rules.sub.AP <- subset(USA_Rules, subset = ((lhs %pin% "Analcime") & (lhs %pin% "Calcite") & (lhs %pin% "Opal")))

rules.sub.AP <- subset(USA_Rules, subset = lhs %ain% MinCombo[,1] & !rhs %in% Loc_Min[[1]])

rules.sub.AP2 <- subset(USA_Rules, subset = lhs %ain% MinCombo[,32] & !rhs %in% Loc_Min[[1]])

inspect(rules.sub.AP)
inspect(rules.sub.AP2)

inspect(union(rules.sub.AP,rules.sub.AP2))


for(i in 1:ncol(MinCombo))
{
  rules.sub.AP2 <- subset(USA_Rules, subset = lhs %ain% MinCombo[,i] & !rhs %in% Loc_Min[[1]])
  rules.sub.AP <-union(rules.sub.AP,rules.sub.AP2)
  inspect(rules.sub.AP)
}


inspect(rules.sub.AP)

rules_pred<-data.frame(lhs = labels(lhs(rules.sub.AP)),rhs = labels(rhs(rules.sub.AP)), 
  rules.sub.AP@quality)

aggregate(x = rules_pred[,4:6], by = list(rules_pred$rhs), FUN = max)  

Loc_Min

MED_export$MED_minerals_all[MED_export$mindat_id==255486]

#_____________________________________________________________

# Mineral Assemblage: Rutherfordine, Andersonite, Schröechingerite
# Data subset: Geochemical (U)
# Predict locality for these rules. 

U_Rules<-readRDS("./data/association_data/Allmin_U_localities_subset.rds")
U_Rules@info

hist(U_Rules@quality$lift)

inspect(head(U_Rules, n = 5, by = "lift"))
plot(rules.sub.AP, method = "grouped")

inspect(U_Rules)[10]s

rules.sub.AP <- subset(U_Rules, subset = (rhs %pin% "Rutherfordine"))
inspect(rules.sub.AP)

rules.sub.AP <- subset(U_Rules, subset = (rhs %pin% "Andersonite"))
inspect(rules.sub.AP)
inspect(head(rules.sub.AP, n = 5, by = "lift"))
plot(rules.sub.AP, method = "grouped")


rules.sub.AP <- subset(U_Rules, subset = (rhs %in% "Bayleyite"))
rules.sub.AP <- subset(U_Rules, subset = (rhs %in% "Schrockingerite"))
rules.sub.AP <- subset(U_Rules, subset = (rhs %in% "Zippeite"))
inspect(rules.sub.AP)
inspect(head(rules.sub.AP, n = 5, by = "lift"))
plot(rules.sub.AP, method = "grouped")

mineral_assemblage <- c("Uraninite","Andersonite","Schrockingerite")
rules.sub.AP <- subset(U_Rules, subset = (lhs %oin% mineral_assemblage))
inspect(rules.sub.AP)

combn(mineral_assemblage,2)[,1]

rules.sub.AP <- subset(U_Rules, subset = (lhs %ain% combn(mineral_assemblage,2)[,3]))
inspect(rules.sub.AP)

rules.sub.AP <- subset(USA_Rules, subset = lhs %ain% MinCombo[,1] & !rhs %in% Loc_Min[[1]])

U_associated_Min_Mat[,which(U_associated_Min_Mat$Rutherfordine==1)]

rownames(U_associated_Min_Mat)[which(U_associated_Min_Mat$Saleeite==1 && U_associated_Min_Mat$Schoepite==1 && U_associated_Min_Mat$Torbernite==1 && U_associated_Min_Mat$Rutherfordine==0)]

which(U_associated_Min_Mat$Saleeite == 1 & U_associated_Min_Mat$Schoepite == 1 & U_associated_Min_Mat$Torbernite == 1 & U_associated_Min_Mat$Rutherfordine != 1)

library(data.table)
write.csv(MED_export[MED_export$MED_minerals_all %like% "Saleeite" & MED_export$MED_minerals_all %like% "Torbernite" & MED_export$MED_minerals_all %like% "Schoepite" & !MED_export$MED_minerals_all %like% "Rutherfordine", c(1,5,6)],file = "U_pred.csv")

write.csv(MED_export[MED_export$MED_minerals_all %like% "Bayleyite" & MED_export$MED_minerals_all %like% "Natrozippeite" & MED_export$MED_minerals_all %like% "Schrockingerite" & !MED_export$MED_minerals_all %like% "Andersonite", c(1,5,6)],file = "U_pred.csv")

write.csv(MED_export[!MED_export$MED_minerals_all %like% "Bayleyite" & MED_export$MED_minerals_all %like% "Natrozippeite" & MED_export$MED_minerals_all %like% "Schrockingerite" & MED_export$MED_minerals_all %like% "Carnotite", c(1,6)],file = "U_pred.csv")

MED_export[MED_export$MED_minerals_all %like% "Pyrite" & MED_export$MED_minerals_all %like% "Andersonite" & MED_export$MED_minerals_all %like% "Carnotite" & !MED_export$MED_minerals_all %like% "Schrockingerite", c(1,5,6)]

MED_export[MED_export$MED_minerals_all %like% "Uraninite" & MED_export$MED_minerals_all %like% "Andersonite" & !MED_export$MED_minerals_all %like% "Zippeite" & MED_export$MED_minerals_all %like% "Schrockingerite", c(1,5,6)]

write.csv(MED_export[MED_export$MED_minerals_all %like% "Uraninite" & MED_export$MED_minerals_all %like% "Andersonite" & !MED_export$MED_minerals_all %like% "Zippeite" & MED_export$MED_minerals_all %like% "Schrockingerite", c(1,6)],file = "U_pred.csv")

MED_export[!MED_export$MED_minerals_all %like% "Bayleyite" & MED_export$MED_minerals_all %like% "Natrozippeite" & MED_export$MED_minerals_all %like% "Schrockingerite" & MED_export$MED_minerals_all %like% "Carnotite", c(1,6)]
MED_export[!MED_export$MED_minerals_all %like% "Bayleyite" & MED_export$MED_minerals_all %like% "Chalcocite" & MED_export$MED_minerals_all %like% "Schrockingerite" & MED_export$MED_minerals_all %like% "Carnotite", c(1,6)]

#______________________________________________

# Query type: How have selected mineral associations changed through deep time? 
# Data subsets: Archean, Proterozoic, Phanerozoic Eons
# Minerals: Tephroite (Mn2+2SiO4), pyrite [Fe2+(S2)2-], akaganeite [(Fe3+,Ni2+)8(OH,O)16Cl1.25·nH2O]

Temporal_Rules<-readRDS("./data/association_data/Allmin_OlderThan2500Ma_Rules_cleaned.rds")
Temporal_Rules
plot(Temporal_Rules, method = "graph", control=list(main = "Top 100 mineral association rules by lift for Archean Eon"))

mineral_assemblage <- c("Altaite","Uraninite")
rules.sub.AP <- subset(Temporal_Rules, subset = (lhs %oin% mineral_assemblage & rhs %in% "Pyrite"))
inspect(rules.sub.AP)
rules.sub.AP <- subset(Temporal_Rules, subset = (rhs %pin% "Butlerite"))
inspect(rules.sub.AP)
plot(rules.sub.AP, method = "graph")


rules.sub.AP <- subset(Temporal_Rules, subset = (size(lhs) == 3 & rhs %in% "Pyrite"))

inspect(rules.sub.AP)
inspect(head(rules.sub.AP, n = 5, by = "lift"))
MED_export[!MED_export$MED_minerals_all %like% "Pyrite" & MED_export$MED_minerals_all %like% "Almandine" & MED_export$MED_minerals_all %like% "Pucherite" & MED_export$MED_minerals_all %like% "Molybdenite", c(1,6)]

plot(rules.sub.AP, method = "scatterplot")
plot(rules.sub.AP, method = "scatterplot", engine = "html")

rules.sub.AP <- subset(Temporal_Rules, subset = (rhs %pin% "Pyrite"))
inspect(rules.sub.AP)
plot(rules.sub.AP, method = "scatterplot")
plot(rules.sub.AP, method = "grouped")

rules.sub.AP <- subset(Temporal_Rules, subset = (size(lhs) == 3 & rhs %in% "Magnetite"))
inspect(rules.sub.AP)
inspect(head(rules.sub.AP, n = 5, by = "lift"))
MED_export[!MED_export$MED_minerals_all %like% "Magnetite" & MED_export$MED_minerals_all %like% "Antigorite" & MED_export$MED_minerals_all %like% "Chromite" & MED_export$MED_minerals_all %like% "Glaukosphaerite", c(1,6)]
plot(rules.sub.AP, method = "scatterplot", engine = "html")
plot(rules.sub.AP, method = "grouped")


Temporal_Rules<-readRDS("./data/association_data/Allmin_from2500MaTo540Ma_cleaned.rds")
Temporal_Rules@info
plot(Temporal_Rules, method = "graph", )
plot(Temporal_Rules, method = "grouped")

rules.sub.AP <- subset(Temporal_Rules, subset = (rhs %pin% "Butlerite"))
inspect(rules.sub.AP)

rules.sub.AP <- subset(Temporal_Rules, subset = (size(lhs) == 3 & rhs %in% "Pyrite"))
inspect(rules.sub.AP)
inspect(head(rules.sub.AP, n = 5, by = "lift"))
MED_export[!MED_export$MED_minerals_all %like% "Pyrite" & MED_export$MED_minerals_all %like% "Kaolinite" & MED_export$MED_minerals_all %like% "Crandallite" & MED_export$MED_minerals_all %like% "Wardite", c(1,6)]
plot(rules.sub.AP, method = "grouped")

rules.sub.AP <- subset(Temporal_Rules, subset = (size(lhs) == 3 & rhs %pin% "Magnetite"))
#inspect(rules.sub.AP)
inspect(head(rules.sub.AP, n = 5, by = "lift"))
MED_export[!MED_export$MED_minerals_all %like% "Magnetite" & MED_export$MED_minerals_all %like% "Aenigmatite" & MED_export$MED_minerals_all %like% "Ilmenite" & MED_export$MED_minerals_all %like% "Pyrrhotite", c(1,6)]

plot(rules.sub.AP, method = "grouped")

Temporal_Rules<-readRDS("./data/association_data/Allmin_YoungerThan540Ma_cleaned.rds")
Temporal_Rules@info

plot(USA_Rules, method = "graph",engine = "html")

rules.sub.AP <- subset(Temporal_Rules, subset = (rhs %pin% "Butlerite"))
inspect(rules.sub.AP)
write.csv(MED_export[!MED_export$MED_minerals_all %like% "Butlerite" & MED_export$MED_minerals_all %like% "Copiapite" & MED_export$MED_minerals_all %like% "Coquimbite" & MED_export$MED_minerals_all %like% "Fibroferrite", c(1,6)],file = "Butlerite_Archean_Localities.csv")


mineral_assemblage <- c("Altaite","Uraninite")
rules.sub.AP <- subset(Temporal_Rules, subset = (lhs %in% mineral_assemblage & rhs %in% "Pyrite"))
inspect(rules.sub.AP)

rules.sub.AP <- subset(Temporal_Rules, subset = (size(lhs) == 3 & rhs %in% "Pyrite"))
inspect(rules.sub.AP)
inspect(head(rules.sub.AP, n = 5, by = "lift"))
MED_export[!MED_export$MED_minerals_all %like% "Pyrite" & MED_export$MED_minerals_all %like% "Amesite" & MED_export$MED_minerals_all %like% "Andradite" & MED_export$MED_minerals_all %like% "Dolomite", c(1,6)]

plot(rules.sub.AP, method = "grouped", control = list(k=10,
                                                      lhs_label_items =5))

rules.sub.AP <- subset(Temporal_Rules, subset = (size(lhs) == 3 & rhs %in% "Magnetite"))
inspect(rules.sub.AP)
inspect(head(rules.sub.AP, n = 5, by = "lift"))
write.csv(MED_export[!MED_export$MED_minerals_all %like% "Magnetite" & MED_export$MED_minerals_all %like% "Chalcopyrite" & MED_export$MED_minerals_all %like% "Erlichmanite" & MED_export$MED_minerals_all %like% "Pentlandite", c(1,6)],file = "Phanerozoic_Magnetite.csv")

plot(rules.sub.AP, method = "grouped")

inspect(sort(Temporal_Rules, decreasing = T, by = "lift"))

inspect(head(Temporal_Rules, n = 20, by = "lift"))

Temporal


#____________________________________________________________________

library(scales)
rescale(U_Rules@quality$lift, from = c(0, max(U_Rules@quality$lift)), to = c(0, 100))

hist(rescale(U_Rules@quality$lift, from = c(0, max(U_Rules@quality$lift)), to = c(0, 100)))
hist(rescale(USA_Rules@quality$lift, from = c(0, max(USA_Rules@quality$lift)), to = c(0, 100)))

Temporal_Rules<-readRDS("./data/association_data/Allmin_OlderThan2500Ma_Rules_cleaned.rds")
hist(rescale(Temporal_Rules@quality$lift, from = c(0, max(Temporal_Rules@quality$lift)), to = c(0, 100)), xlab = "Normalized Lift",main = "Distribution of Lift for Archean")

Temporal_Rules<-readRDS("./data/association_data/Allmin_from2500MaTo540Ma_cleaned.rds")
hist(rescale(Temporal_Rules@quality$lift, from = c(0, max(Temporal_Rules@quality$lift)), to = c(0, 100)), xlab = "Normalized Lift",main = "Distribution of Lift for Proterozoic")
#hist(Temporal_Rules@quality$lift)

Temporal_Rules<-readRDS("./data/association_data/Allmin_YoungerThan540Ma_cleaned.rds")
hist(rescale(Temporal_Rules@quality$lift, from = c(0, max(Temporal_Rules@quality$lift)), to = c(0, 100)), xlab = "Normalized Lift",main = "Distribution of Lift for Phanerazoic")
#hist(Temporal_Rules@quality$lift)


#_________________________________________________________

# USA subset and search for critical minerals

USA_Rules<-readRDS("./data/association_data/Allmin_USA_localities_subset_cleaned.rds")

rules.sub.AP <- subset(USA_Rules, subset = (rhs %pin% "Monazite"))
inspect(rules.sub.AP)
inspect(head(rules.sub.AP, n = 50, by = "lift"))
# {Elbaite, Rutile, Sphalerite}                => {Monazite-(Ce)}

library(data.table)
write.csv(x = MED_USA_subset[MED_USA_subset$MED_minerals_all %like% "Sphalerite" & MED_USA_subset$MED_minerals_all %like% "Elbaite" & MED_USA_subset$MED_minerals_all %like% "Rutile" & !MED_USA_subset$MED_minerals_all %like% "Monazite-(Ce)", c(1,6)],file = "USA_Localities_For_Monazite-Ce.csv")



rules.sub.AP <- subset(USA_Rules, subset = (rhs %pin% "Spodumene"))
#inspect(rules.sub.AP)
inspect(head(rules.sub.AP, n = 60, by = "lift"))
#Beryl, Mitridatite, Pyrite
rules.sub.AP <- subset(USA_Rules, subset = (lhs %in% "Pyrite" & lhs %in% "Mitridatite" & rhs %pin% "Spodumene"))
write.csv(x = MED_USA_subset[MED_USA_subset$MED_minerals_all %like% "Beryl" & MED_USA_subset$MED_minerals_all %like% "Mitridatite" & MED_USA_subset$MED_minerals_all %like% "Pyrite" & !MED_USA_subset$MED_minerals_all %like% "Spodumene", c(1,6)],file = "USA_Localities_For_Spodumene.csv")


inspect(tail(rules.sub.AP))

rules.sub.AP <- subset(USA_Rules, subset = (rhs %pin% "Allanite"))
inspect(rules.sub.AP)
inspect(head(rules.sub.AP, n = 20, by = "lift"))
# {Azurite, Microcline, Thorite}          => {Allanite-(Ce)}

MED_USA_subset[MED_USA_subset$MED_minerals_all %like% "Azurite" & MED_USA_subset$MED_minerals_all %like% "Microcline" & MED_USA_subset$MED_minerals_all %like% "Thorite" & !MED_USA_subset$MED_minerals_all %like% "Allanite-(Ce)", c(1,6)]
MED_USA_subset[MED_USA_subset$MED_minerals_all %like% "Chalcopyrite" & MED_USA_subset$MED_minerals_all %like% "Thorite" & MED_USA_subset$MED_minerals_all %like% "Monazite-(Ce)" & !MED_USA_subset$MED_minerals_all %like% "Allanite-(Ce)", c(1,6)]


write.csv(x = MED_USA_subset[MED_USA_subset$MED_minerals_all %like% "Azurite" & MED_USA_subset$MED_minerals_all %like% "Microcline" & MED_USA_subset$MED_minerals_all %like% "Thorite" & !MED_USA_subset$MED_minerals_all %like% "Allanite-(Ce)", c(1,6)],file = "USA_Localities_For_Allanite-Ce.csv")
