# Forest & Kate Mars Mineral AA Generation

Mars_Min_AA_Data<-read.csv("./data/mars_data/mars_export_binary.csv") #FIXME File DNE

aa<-rownames(Mars_Min_AA_Data)

Mars_Min_AA_Data<-sapply(Mars_Min_AA_Data, as.logical)
rownames(Mars_Min_AA_Data)<-aa

Mars_Min_AA_Data<-as.data.frame(unclass(Mars_Min_AA_Data))
                                
library(arules)

rules_Mars_AA<-apriori(data = Mars_Min_AA_Data, parameter = list(supp = 0.003, conf = 0.7, target = "rules", maxtime = 0, minlen = 2, maxlen = 4))

inspect(tail(rules_Mars_AA))