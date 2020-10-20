library(tidyverse)
library(MSstats)
library(here)

load(here("data","DDA.proposed.rda"))
load(here("data","new_DDA.proposed.rda"))
load(here("data","DDA.inf.rda"))

mycomparison1<-matrix(c(-1,1,0,0),nrow=1)
mycomparison2<-matrix(c(-1,0,1,0),nrow=1)
mycomparison3<-matrix(c(-1,0,0,1),nrow=1)
mycomparison4<-matrix(c(0,-1,1,0),nrow=1)
mycomparison5<-matrix(c(0,-1,0,1),nrow=1)
mycomparison6<-matrix(c(0,0,-1,1),nrow=1)

mycomparison<-rbind(mycomparison1,mycomparison2, mycomparison3,mycomparison4,mycomparison5,mycomparison6)

row.names(mycomparison)<-c("C2-C1","C3-C1","C4-C1","C2-C3","C2-C4","C3-C4")

new_DDA.comparisons <- groupComparison(contrast.matrix = mycomparison, data = new_DDA.proposed)
DDA.inf.featureselected.comparisons <- groupComparison(contrast.matrix = mycomparison, data = DDA.inf)


Pc<-new_DDA.comparisons$ComparisonResult %>%
  filter(adj.pvalue<0.01) %>%
  filter(!Protein %in% spike_in) %>% 
  mutate(method="PeCorA",
         comparison= paste0(Protein,Label)) %>%
  
  select(Protein,Label,log2FC,comparison,method)


Fs<-DDA.inf.featureselected.comparisons$ComparisonResult %>%
  separate(Protein,into=c("Protein","specie"),sep="_") %>%
  filter(adj.pvalue<0.01) %>%
  filter(!Protein %in% spike_in) %>% 
  mutate(method="FeatureSelection",
         comparison= paste0(Protein,Label)) %>%
  select(Protein,Label,log2FC,comparison,method)

df<-inner_join(Fs,Pc,by="comparison") %>%
  filter_all(all_vars(!is.infinite(.))) %>% 
  select(log2FC.x,log2FC.y, method.x, method.y,comparison)

lm_fit <- lm(log2FC.x ~ log2FC.y, data=df)

df %>%
  ggplot(aes(x= log2FC.x,y=log2FC.y))+
  geom_point(color='blue', size=7,alpha=.5) +
  geom_abline(intercept=0, slope=lm_fit$coefficients[2], color='#2C3E50', size=1.1) +
  xlab("Log2 FC post Feature Selection")+
  ylab("Log2 FC post Pecora")+
  theme_minimal()+
  theme(axis.text.x = element_text(size =40),
        axis.text.y = element_text(size =40),
        axis.title.x = element_text(size =35),
        axis.title.y = element_text(size =35))
 #ggsave("scatter_Pecora.png")



