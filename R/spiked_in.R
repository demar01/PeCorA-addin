library(tidyverse)
library(here)
library(MSstats)
library(kableExtra)

## Table 1: spike-in proteins (https://pubmed.ncbi.nlm.nih.gov/27990823/)
spike_in<-c("sp|P44015|VAC2_YEAST", "sp|P55752|ISCB_YEAST",
  "sp|P44374|SFG2_YEAST", "sp|P44983|UTR6_YEAST",
  "sp|P44683|PGA4_YEAST", "sp|P55249|ZRT4_YEAST")
## tryptic digest spiked with the indicated amounts (in fmols) of tryptic digest of six individual proteins.
Sample_preparation <-
 matrix(c(65,55,15,2,
    55,15,2,65,
    15,2,65,55,
    2,65,55,15,
    11,0.6,10,500,
    10,500,11,0.6),
   nrow = 6, ncol = 4, byrow = TRUE) %>%
  as.data.frame()
rownames(Sample_preparation)<-spike_in

group_comparisons<-
  matrix(c(log2(55/65),log2(15/65),log2(2/65),log2(55/15),log2(55/2),log2(15/2),
           log2(15/55),log2(2/55),log2(65/55),log2(15/2),log2(15/65),log2(2/65),
           log2(2/15),log2(65/15),log2(55/15),log2(2/65),log2(2/55),log2(65/55),
           log2(65/2),log2(55/2),log2(15/2),log2(65/55),log2(65/15),log2(55/15),
           log2(0.6/11),log2(10/11),log2(500/11),log2(0.6/10),log2(0.6/500),log2(10/500),
          log2(500/10),log2(11/10),log2(0.6/10),log2(500/11),log2(500/0.6),log2(11/0.6)),
         nrow = 6, ncol = 6, byrow = TRUE) %>%
  as.data.frame()
rownames(group_comparisons)<-spike_in
colnames(group_comparisons)<-c("C2-C1","C3-C1","C4-C1","C2-C3","C2-C4","C3-C4")

real<-group_comparisons %>% tibble::rownames_to_column( "Protein") %>%
  pivot_longer(-Protein, names_to="Label",values_to="log2FC") %>%
  arrange(Protein,Label) %>%
  mutate(comparison=str_c(Protein,Label)) %>%
  mutate(log2FC=abs(log2FC))

spike_in %in% DDA.comparisons$ComparisonResult$Protein

load(here("data","DDA.proposed.rda"))
load(here("data","new_DDA.proposed.rda"))
load(here("data","DDA.inf.rda"))

#comparisons...
mycomparison1<-matrix(c(-1,1,0,0),nrow=1)
mycomparison2<-matrix(c(-1,0,1,0),nrow=1)
mycomparison3<-matrix(c(-1,0,0,1),nrow=1)
mycomparison4<-matrix(c(0,-1,1,0),nrow=1)
mycomparison5<-matrix(c(0,-1,0,1),nrow=1)
mycomparison6<-matrix(c(0,0,-1,1),nrow=1)

mycomparison<-rbind(mycomparison1,mycomparison2, mycomparison3,mycomparison4,mycomparison5,mycomparison6)

row.names(mycomparison)<-c("C2-C1","C3-C1","C4-C1","C2-C3","C2-C4","C3-C4")

DDA.comparisons <- groupComparison(contrast.matrix = mycomparison, data = DDA.proposed)
DDA.comparisons$ComparisonResult %>% dim()

DDA.comparisons$ComparisonResult %>% count(pvalue<0.05)


new_DDA.comparisons <- groupComparison(contrast.matrix = mycomparison, data = new_DDA.proposed)
new_DDA.comparisons$ComparisonResult %>% dim()

DDA.inf.featureselected.comparisons <- groupComparison(contrast.matrix = mycomparison, data = DDA.inf)
DDA.inf.featureselected.comparisons$ComparisonResult %>% dim()


estimated_all<-DDA.comparisons$ComparisonResult %>%
  filter(Protein %in% spike_in) %>% arrange(Protein) %>%
  select(Protein,Label,log2FC,adj.pvalue) %>%
  arrange(Protein,Label) %>%
  mutate(comparison=str_c(Protein,Label)) %>%
  mutate(log2FC=abs(log2FC))
spike_in %in% DDA.comparisons$ComparisonResult$Protein

proposed<-DDA.inf.featureselected.comparisons$ComparisonResult %>%
  filter(Protein %in% spike_in) %>% arrange(Protein) %>%
  select(Protein,Label,log2FC,adj.pvalue) %>%
  arrange(Protein,Label) %>%
  mutate(comparison=str_c(Protein,Label)) %>%
  mutate(log2FC=abs(log2FC))
spike_in %in% DDA.inf.featureselected.comparisons$ComparisonResult$Protein

pecora<-new_DDA.comparisons$ComparisonResult %>%
  filter(Protein %in% c("sp|P44015|VAC2","sp|P55752|ISCB", "sp|P44374|SFG2","sp|P44983|UTR6", "sp|P44683|PGA4","sp|P55249|ZRT4") ) %>% arrange(Protein) %>%
  select(Protein,Label,log2FC,adj.pvalue) %>%
  arrange(Protein,Label) %>%
  mutate(comparison=str_c(Protein,"_YEAST",Label)) %>%
  mutate(log2FC=abs(log2FC))
c("sp|P44015|VAC2","sp|P55752|ISCB", "sp|P44374|SFG2","sp|P44983|UTR6", "sp|P44683|PGA4","sp|P55249|ZRT4") %in%
  new_DDA.comparisons$ComparisonResult$Protein


a<-inner_join(real, estimated_all,by="comparison")  %>%
  dplyr::rename(real=log2FC.x) %>%
  dplyr::rename(estimated_all=log2FC.y)
fita <- lm(real ~ 0 + estimated_all, data = a) # Adding the 0 term tells the lm() to fit the line through the origin
format(summary(fita)$r.squared, digits = 2)
#0.884

b<-inner_join(real, proposed,by="comparison")  %>%
  dplyr:: rename(real=log2FC.x) %>%
  dplyr::rename(proposed=log2FC.y)
fitb <- lm(real ~ 0 + proposed, data = b) # Adding the 0 term tells the lm() to fit the line through the origin
format(summary(fitb)$r.squared, digits = 2)
#0.865

c<-inner_join(real, pecora,by="comparison")  %>%
  dplyr::rename(real=log2FC.x) %>%
  dplyr::rename(pecora=log2FC.y)
fitc <- lm(real ~ 0 + pecora, data = c) # Adding the 0 term tells the lm() to fit the line through the origin
format(summary(fitc)$r.squared, digits = 2)
#"0.874

data_real<-inner_join(real, estimated_all ,by="comparison") %>%
  dplyr::rename(real=log2FC.x) %>%
  dplyr::rename(estimated_all=log2FC.y)
  plot_all<- ggplot( data_real,aes(x = real , y = estimated_all))+
  geom_point(aes(color=adj.pvalue<0.05),alpha=0.5,size=7) +
  #geom_smooth(method="lm",level=0.95,se = TRUE) +
  geom_abline(intercept=0, slope=fita$coefficients[1], color='#2C3E50', size=1.1) +

    xlim(0,10)+    ylim(0,10)+
    theme_minimal()+
    scale_color_manual(values=c("black", "red")) +
    theme(legend.position = "none",
          axis.text.x = element_text(size =40),
          axis.text.y = element_text(size =40),
          axis.title.x = element_text(size =35),
          axis.title.y = element_text(size =35))+
    ylab("Estimated All log2FC") +
    xlab("True log2FC")

  # geom_text(label=data$comparison,nudge_x = 0.2, nudge_y = 0.4, size = 2,check_overlap = F )

data_proposed<-inner_join(real, proposed ,by="comparison") %>%
  dplyr::rename(real=log2FC.x) %>%
  dplyr::rename(proposed=log2FC.y)
  plot_proposed<-ggplot( data_proposed,aes(x = real , y = proposed))+
    geom_point(aes(color=adj.pvalue<0.05),alpha=0.5,size=7) +
   #geom_smooth(method="lm",level=0.95,se = TRUE) +
    geom_abline(intercept=0, slope=fitb$coefficients[1], color='#2C3E50', size=1.1) +

    xlim(0,10)+    ylim(0,10)+
    theme_minimal()+
    scale_color_manual(values=c("black", "red")) +
    theme(legend.position = "none",
          axis.text.x = element_text(size =40),
          axis.text.y = element_text(size =40),
          axis.title.x = element_text(size =35),
          axis.title.y = element_text(size =32))+
    ylab("Feature Selection log2FC") +
    xlab("True log2FC")

  data_pecora<-inner_join(real, pecora ,by="comparison") %>%
    dplyr::rename(real=log2FC.x) %>%
    dplyr::rename(pecora=log2FC.y)
  plot_pecora<-ggplot( data_pecora,aes(x = real , y = pecora))+
    geom_point(aes(color=adj.pvalue<0.05),alpha=0.5,size=7) +
   # geom_smooth(method="lm",level=0.95,se = TRUE) +
    geom_abline(intercept=0, slope=fitc$coefficients[1], color='#2C3E50', size=1.1) +

    xlim(0,10)+    ylim(0,10)+
    theme_minimal()+
    scale_color_manual(values=c("black", "red")) +
    theme(legend.position = "none",
          legend.title = element_text(size = 40),
          legend.text = element_text(size = 40),
          axis.text.x = element_text(size =40),
          axis.text.y = element_text(size =40),
          axis.title.x = element_text(size =35),
          axis.title.y = element_text(size =35))+
    ylab("PeCorA log2FC") +
    xlab("True log2FC")


  plot_all+ ggsave(here::here("plot_all.png"))
  plot_proposed+ ggsave(here::here("plot_proposed.png"))
  plot_pecora+ ggsave(here::here("plot_pecora.png"))

  ####counting TP####
  data_real %>%
    count(adj.pvalue <0.05)
  data_proposed %>%
    count(adj.pvalue <0.05)
  data_pecora %>%
    count(adj.pvalue <0.05)


table_Fig3c<-tribble(
  ~Comparisons, ~All,  ~Feature.Selection, ~PeCorA,
  "Diff. abundant",26,29,25 ,
  "Not diff. abundant", 10,7,11
)

table_Fig3c %>%
kable() %>%
kable_styling(font_size = 40)
#save_kable("Table_Fig3C.png")


