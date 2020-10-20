library(tidyverse)
library(MSstats)
library(here)


load(here("data","DDA.comparisons.rda"))
load(here("data","new_DDA.proposed.rda"))
load(here("data","DDA.inf.rda"))


SE_forPlots<-data.frame(name=c(rep("All",9045),rep("Proposed",8997),rep("Pecora",9018)),
value=c(DDA.comparisons$ComparisonResult$SE,
DDA.inf.featureselected.comparisons$ComparisonResult$SE,
new_DDA.comparisons$ComparisonResult$SE))

log2_forPlots<-data.frame(name=c(rep("All",9045),rep("Proposed",8997),rep("Pecora",9018)),
value=c(DDA.comparisons$ComparisonResult$log2FC,
        DDA.inf.featureselected.comparisons$ComparisonResult$log2FC,
        new_DDA.comparisons$ComparisonResult$log2FC))


SE_forPlots %>% ggplot(aes(x=name,y=value, fill=name))+
  geom_boxplot(size=1)+
coord_cartesian(ylim=c(0,0.3))+
  theme_minimal()+
scale_fill_manual(values=c("darkgrey", "darkgreen","darkorange")) +
theme(legend.position = "none") +
ylab("Standard error")+xlab("Method")+
  theme(text = element_text(size=30))
#ggsave(here::here("allplots", "Performance_eval_SE.png"))



log2_forPlots %>% ggplot(aes(x=name,y=value, fill=name))+
  geom_boxplot(size=1)+
  coord_cartesian(ylim=c(-0.5,0.5))+
  theme_minimal()+
  scale_fill_manual(values=c("darkgrey", "darkgreen","darkorange")) +
  theme(legend.position = "none") +
  ylab("log2 FC")+xlab("Method")+
  theme(text = element_text(size=30))
#ggsave(here::here("allplots", "Performance_eval_log2.png"))

