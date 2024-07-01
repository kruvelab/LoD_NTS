library(tidyverse)
library(readxl)
library(rcdk)
library(caret)
library(caTools)
library(Metrics)
library(RRF)
library(scales)
library(extrafont)
library(plotly)
library(cowplot)
library(ggplot2)
library(scales)
library(car)
library(caret)
library(caTools)
library(dplyr)
library(enviPat)
# install.packages("webchem")
library(gbm)
library(gsubfn)
library(janitor)
library(Metrics)
require(mgcv)
library(OrgMassSpecR)
library(plotly)
library(randomForest)
library(rcdk)
library(rcdklibs)
library(rgl)
library(rJava)
library(rjson)
library(RRF)
library(proto)
library(stringr)
library(webchem)
setwd("C:/transferred_data/amso1873/Desktop/Research/LoD new runs")
#------- Standard solutions-------
Alignment_solutions<- read.delim("Alignment_standard_solutions.txt", skip = 4)
Alignment_SN<- read.delim("SN_standard_solutions.txt", skip = 4)
Alignment_solutions <- Alignment_solutions %>% 
  filter(Adduct.type == "[M+H]+")  

suspect_list <- read_delim("C:/transferred_data/amso1873/Downloads/inclusion_list_Amina.csv", delim=",") %>% 
  filter(Comment!="Dazomet")  %>% filter(Comment!="Dimethyl phthalate H")  %>%
  filter(Comment!= "trans-ferulic acid") %>% filter(Comment!="Prometryn")  %>%
  filter(Comment!= "atenolol") %>% filter(Comment!="Imazalil")
suspect_list <- suspect_list %>% mutate(Average.Mz=round(suspect_list$`Mass [m/z]`, digits = 2))
Alignment_solutions <- Alignment_solutions %>% mutate(measured_mz=Alignment_solutions$Average.Mz)
Alignment_solutions$Average.Mz <-  round(Alignment_solutions$Average.Mz, digits = 2) 
colnames(suspect_list)[12] <- "Name"
Alignment_solutions <- Alignment_solutions %>% left_join(suspect_list)
Alignment_solutions_summary <- Alignment_solutions %>% filter(Name !="")
Alignment_solutions_summary <- Alignment_solutions_summary %>% select(Alignment.ID, Name, Average.Mz, Adduct.type, Average.Rt.min., MS.MS.spectrum, Blank_1_24012023, Blank_2_24012023, Blank_3_24012023, S1_1_24012023, S1_2_24012023, S1_3_24012023, S2_1_24012023, S2_2_24012023, S2_3_24012023, S3_1_24012023, S3_2_24012023, S3_3_24012023, S4_1_24012023, S4_2_24012023, S4_3_24012023, S5_1_24012023, S5_2_24012023, S5_3_24012035, S6_1_24012023, S6_2_24012023, S6_3_24012023, S7_1_24012023, S7_2_24012023, S7_3_24012023, 
                                                                      S8_1_24012023, S8_2_24012023, S8_3_24012023, S9_1_24012023, S9_2_24012023, S9_3_24012023, S10_1_24012023, S10_2_24012023, S10_3_24012023 , S11_1_24012023, S11_2_24012023, S11_3_24012023, S12_1_24012023, S12_2_24012023, S12_3_24012023, S13_1_24012023, S13_2_24012023, S13_3_24012023, S14_1_24012023, S14_2_24012023, S14_3_24012023, S15_1_24012023, S15_2_24012023, S15_3_24012023, S16_1_24012023, S16_2_24012023, S16_3_24012023 ) %>% 
  mutate(Diff_mz=abs(Alignment_solutions_summary$`Mass [m/z]`- Alignment_solutions_summary$measured_mz)*1000000/Alignment_solutions_summary$`Mass [m/z]`)
Alignment_solutions_summary <- Alignment_solutions_summary %>% filter(Diff_mz<6)
Alignment_solutions_summary_2 <- Alignment_solutions_summary %>%
  group_by(Alignment.ID, Name,  Average.Mz, Adduct.type, Average.Rt.min., MS.MS.spectrum) %>% 
  summarize(Blank=mean(Blank_1_24012023, Blank_2_24012023, Blank_3_24012023),
            S1=mean(S1_1_24012023, S1_2_24012023, S1_3_24012023), 
            S2=mean(S2_1_24012023, S2_2_24012023, S2_3_24012023),
            S3=mean(S3_1_24012023, S3_2_24012023, S3_3_24012023),
            S4=mean(S4_1_24012023, S4_2_24012023, S4_3_24012023),
            S5=mean(S5_1_24012023, S5_2_24012023, S5_3_24012035),
            S6=mean(S6_1_24012023, S6_2_24012023, S6_3_24012023),
            S7=mean(S7_1_24012023, S7_2_24012023, S7_3_24012023),
            S8=mean(S8_1_24012023, S8_2_24012023, S8_3_24012023),
            S9=mean(S9_1_24012023, S9_2_24012023, S9_3_24012023),
            S10=mean(S10_1_24012023, S10_2_24012023, S10_3_24012023),
            S11=mean(S11_1_24012023, S11_2_24012023, S11_3_24012023), 
            S12=mean(S12_1_24012023, S12_2_24012023, S12_3_24012023),
            S13=mean(S13_1_24012023, S13_2_24012023, S13_3_24012023),
            S14=mean(S14_1_24012023, S14_2_24012023, S14_3_24012023),
            S15=mean(S15_1_24012023, S15_2_24012023, S15_3_24012023),
            S16=mean(S16_1_24012023, S16_2_24012023, S16_3_24012023), 
            S_SD1=sd(c(S1_1_24012023, S1_2_24012023, S1_3_24012023)),
            S_SD2=sd(c(S2_1_24012023, S2_2_24012023, S2_3_24012023)),
            S_SD3=sd(c(S3_1_24012023, S3_2_24012023, S3_3_24012023)),
            S_SD4=sd(c(S4_1_24012023, S4_2_24012023, S4_3_24012023)),
            S_SD5=sd(c(S5_1_24012023, S5_2_24012023, S5_3_24012035)),
            S_SD6=sd(c(S6_1_24012023, S6_2_24012023, S6_3_24012023)),
            S_SD7=sd(c(S7_1_24012023, S7_2_24012023, S7_3_24012023)),
            S_SD8=sd(c(S8_1_24012023, S8_2_24012023, S8_3_24012023)),
            S_SD9=sd(c(S9_1_24012023, S9_2_24012023, S9_3_24012023)),
            S_SD10=sd(c(S10_1_24012023, S10_2_24012023, S10_3_24012023)),
            S_SD11=sd(c(S11_1_24012023, S11_2_24012023, S11_3_24012023)), 
            S_SD12=sd(c(S12_1_24012023, S12_2_24012023, S12_3_24012023)),
            S_SD13=sd(c(S13_1_24012023, S13_2_24012023, S13_3_24012023)),
            S_SD14=sd(c(S14_1_24012023, S14_2_24012023, S14_3_24012023)),
            S_SD15=sd(c(S15_1_24012023, S15_2_24012023, S15_3_24012023)),
            S_SD16=sd(c(S16_1_24012023, S16_2_24012023, S16_3_24012023)))%>% ungroup() 
Alignment_solutions_summary_3 <- Alignment_solutions_summary_2  %>% select(-c( Name, Average.Mz, Adduct.type, Average.Rt.min., MS.MS.spectrum))
Alignment_solutions_summary_3_SD <-Alignment_solutions_summary_3 %>% select(-c(Blank, S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16)) 
Alignment_solutions_summary_3_SD <-Alignment_solutions_summary_3_SD %>% gather(key = "Solution", value = "SD_Area", 2:17) 
Alignment_solutions_summary_3_SD$Solution[Alignment_solutions_summary_3_SD$Solution == 'S_SD1'] <- 'S1'
Alignment_solutions_summary_3_SD$Solution[Alignment_solutions_summary_3_SD$Solution == 'S_SD2'] <- 'S2'
Alignment_solutions_summary_3_SD$Solution[Alignment_solutions_summary_3_SD$Solution == 'S_SD3'] <- 'S3'
Alignment_solutions_summary_3_SD$Solution[Alignment_solutions_summary_3_SD$Solution == 'S_SD4'] <- 'S4'
Alignment_solutions_summary_3_SD$Solution[Alignment_solutions_summary_3_SD$Solution == 'S_SD5'] <- 'S5'
Alignment_solutions_summary_3_SD$Solution[Alignment_solutions_summary_3_SD$Solution == 'S_SD6'] <- 'S6'
Alignment_solutions_summary_3_SD$Solution[Alignment_solutions_summary_3_SD$Solution == 'S_SD7'] <- 'S7'
Alignment_solutions_summary_3_SD$Solution[Alignment_solutions_summary_3_SD$Solution == 'S_SD8'] <- 'S8'
Alignment_solutions_summary_3_SD$Solution[Alignment_solutions_summary_3_SD$Solution == 'S_SD9'] <- 'S9'
Alignment_solutions_summary_3_SD$Solution[Alignment_solutions_summary_3_SD$Solution == 'S_SD10'] <- 'S10'
Alignment_solutions_summary_3_SD$Solution[Alignment_solutions_summary_3_SD$Solution == 'S_SD11'] <- 'S11'
Alignment_solutions_summary_3_SD$Solution[Alignment_solutions_summary_3_SD$Solution == 'S_SD12'] <- 'S12'
Alignment_solutions_summary_3_SD$Solution[Alignment_solutions_summary_3_SD$Solution == 'S_SD13'] <- 'S13'
Alignment_solutions_summary_3_SD$Solution[Alignment_solutions_summary_3_SD$Solution == 'S_SD14'] <- 'S14'
Alignment_solutions_summary_3_SD$Solution[Alignment_solutions_summary_3_SD$Solution == 'S_SD15'] <- 'S15'
Alignment_solutions_summary_3_SD$Solution[Alignment_solutions_summary_3_SD$Solution == 'S_SD16'] <- 'S16'

Alignment_solutions_summary_3 <-Alignment_solutions_summary_3[1:18] 
Alignment_solutions_summary_3 <-Alignment_solutions_summary_3 %>% mutate(Ratio=S1/Blank)
Alignment_solutions_summary_3 <-Alignment_solutions_summary_3 %>% filter(Ratio>5)
Alignment_solutions_summary_3 <-Alignment_solutions_summary_3 %>% select(-Ratio)
Alignment_solutions_summary_3 <-Alignment_solutions_summary_3 %>%
  mutate(S1=S1-Blank) %>%
  mutate(S2=S2-Blank) %>%
  mutate(S3=S3-Blank) %>%
  mutate(S4=S4-Blank) %>%
  mutate(S5=S5-Blank) %>%
  mutate(S6=S6-Blank) %>%
  mutate(S7=S7-Blank) %>%
  mutate(S8=S8-Blank) %>%
  mutate(S9=S9-Blank) %>%
  mutate(S10=S10-Blank) %>%
  mutate(S11=S11-Blank) %>%
  mutate(S12=S12-Blank) %>%
  mutate(S13=S13-Blank) %>%
  mutate(S14=S14-Blank) %>%
  mutate(S15=S15-Blank) %>%
  mutate(S16=S16-Blank)
# Alignment_solutions_summary_3 <-Alignment_solutions_summary_3 %>% filter(Blank==0)
Alignment_solutions_summary_3 <-Alignment_solutions_summary_3 %>% select(-Blank)
Alignment_solutions_summary_3 <-Alignment_solutions_summary_3 %>%
  gather(key = "Solution", value = "Area", 2:17) 
Alignment_SN <- Alignment_SN %>%
  group_by(Alignment.ID) %>% 
  summarize(S1=mean(S1_1_24012023, S1_2_24012023, S1_3_24012023),
            S2=mean(S2_1_24012023, S2_2_24012023, S2_3_24012023),
            S3=mean(S3_1_24012023, S3_2_24012023, S3_3_24012023),
            S4=mean(S4_1_24012023, S4_2_24012023, S4_3_24012023),
            S5=mean(S5_1_24012023, S5_2_24012023, S5_3_24012035),
            S6=mean(S6_1_24012023, S6_2_24012023, S6_3_24012023),
            S7=mean(S7_1_24012023, S7_2_24012023, S7_3_24012023),
            S8=mean(S8_1_24012023, S8_2_24012023, S8_3_24012023),
            S9=mean(S9_1_24012023, S9_2_24012023, S9_3_24012023),
            S10=mean(S10_1_24012023, S10_2_24012023, S10_3_24012023),
            S11=mean(S11_1_24012023, S11_2_24012023, S11_3_24012023), 
            S12=mean(S12_1_24012023, S12_2_24012023, S12_3_24012023),
            S13=mean(S13_1_24012023, S13_2_24012023, S13_3_24012023),
            S14=mean(S14_1_24012023, S14_2_24012023, S14_3_24012023),
            S15=mean(S15_1_24012023, S15_2_24012023, S15_3_24012023),
            S16=mean(S16_1_24012023, S16_2_24012023, S16_3_24012023))%>% ungroup()
Alignment_SN<-Alignment_SN %>% 
  gather(key = "Solution", value = "SN", 2:17)
Alignment_solutions_summary_3 <-Alignment_solutions_summary_3 %>% left_join(Alignment_SN)
Alignment_solutions_summary_3 <-Alignment_solutions_summary_3 %>% left_join(Alignment_solutions_summary_3_SD)
Concentrations <- read_excel("C:/transferred_data/amso1873/Desktop/LOD mix 29122022 - Concentrations.xlsx", 
                             sheet = "Sheet2") %>% select(- `Cf (ug/L)`) %>% gather(key="Solution", value = "Conc", 2:17)
Concentrations <- Concentrations %>% left_join(suspect_list %>% select(c(Name, `Mass [m/z]`))) %>% mutate(Molar_mass = `Mass [m/z]`-1.007)
Names_list = Alignment_solutions_summary_2 %>% select(c(Name , Alignment.ID , Average.Mz, Average.Rt.min., MS.MS.spectrum))
Concentrations <- Concentrations %>% mutate(Conc=(Conc*10^(-6))/Molar_mass)
Alignment_solutions_summary_3 = Alignment_solutions_summary_3 %>% left_join(Names_list)
Concentrations <- Concentrations %>% select(c(Name, Solution, Conc)) %>% drop_na()
Alignment_solutions_summary_3 = Alignment_solutions_summary_3 %>% left_join(Concentrations)
Alignment_solutions_summary_3 = Alignment_solutions_summary_3 %>% drop_na()
suspect_list_al_id <- read_excel("C:/transferred_data/amso1873/Desktop/Research/LoD new runs/AligmentID_spiked_chemical.xlsx")
suspect_list_al_id <- suspect_list_al_id %>% select(Name, Alignment.ID, `Mass [m/z]`)
suspect_list_al_id <- suspect_list_al_id %>% drop_na()
suspect_list_al_id <- suspect_list_al_id %>% mutate(Correct=TRUE)
suspect_list_al_id$Alignment.ID <- as.integer(suspect_list_al_id$Alignment.ID)
Alignment_solutions_summary_3 <- Alignment_solutions_summary_3 %>% left_join(suspect_list_al_id)
Alignment_solutions_summary_3 <- Alignment_solutions_summary_3 %>% filter(Correct==TRUE)
Alignment_solutions_summary_3 <- Alignment_solutions_summary_3 %>% filter(Area!=0)
Alignment_solutions_summary_3 <- Alignment_solutions_summary_3 %>% filter(Name!="lufenuron")

linear_regression <- function(Signal, concentration, SN, SD_Area) {
  y = Signal
  x = concentration
  SN=SN
  slope = summary(lm(log(y) ~ log(x)))$coefficients[2]
  intercept = summary(lm(log(y) ~ log(x)))$coefficients[1]
  R2 = summary(lm(log(y) ~ log(x)))$r.squared
  residuals = (log(y) - (slope*log(x) +intercept))/log(y)*100
  lod1=min(x)
  if (SN[which(x==min(x))] < 11){
    lod2 = min(x)
  } else {
    lod2 = lod1*3/SN[which(x==min(x))]
  }
  lod4=SD_Area[which(x==min(x))]*3/slope
  max_residuals=max(abs(residuals))
  while (max_residuals>5 & length(y)>3){
    y=y[-length(y)]
    x=x[-length(x)]
    slope = summary(lm(log(y) ~ log(x)))$coefficients[2]
    intercept = summary(lm(log(y) ~ log(x)))$coefficients[1]
    R2 = summary(lm(log(y) ~ log(x)))$r.squared
    residuals = (log(y) - (slope*log(x) +intercept))/log(y)*100
    lod1=min(x)
    if (SN[which(x==min(x))] < 11){
      lod2 = min(x)
    } else {
      lod2 = lod1*3/SN[which(x==min(x))]
    }
    lod4=SD_Area[which(x==min(x))]*3/slope
    max_residuals=max(abs(residuals))
  }
  slope = summary(lm(y ~ x))$coefficients[2]
  intercept = summary(lm(y ~ x))$coefficients[1]
  R2 = summary(lm(y~ x))$r.squared
  residuals = (y - (slope*x +intercept))/y*100
  lod1=min(x)
  if (SN[which(x==min(x))] < 11){
    lod2 = min(x)
  } else {
    lod2 = lod1*3/SN[which(x==min(x))]
  }
  lod4=SD_Area[which(x==min(x))]*3/slope
  if  (length(x) >3) {
    y_lod3 = y[(length(y)-2): length(y)]
    x_lod3 = x[(length(x)-2): length(x)]
    slope_lod3 = summary(lm(log(y_lod3) ~ log(x_lod3)))$coefficients[2]
    intercept_lod3 = summary(lm(log(y_lod3) ~ log(x_lod3)))$coefficients[1]
    R2_lod3 = summary(lm(log(y_lod3) ~ log(x_lod3)))$r.squared
    residuals_lod3 = (log(y_lod3) - (slope_lod3*log(x_lod3) +intercept_lod3))/log(y_lod3)*100
    lod3=(3*sd(residuals_lod3))/slope_lod3}
  else {
    lod3=(3*sd(residuals))/slope}
  
  regression_parameters <- list("slope" = slope,
                                "intercept" = intercept,
                                "R2" = R2,
                                "lod1"=lod1,
                                "lod3"=lod3, 
                                "lod2"=lod2, 
                                "lod4"=lod4, 
                                "max_residuals"= max_residuals)
  return(regression_parameters)
}
Alignment_solutions_summary_3 <- Alignment_solutions_summary_3 %>% filter(Area>0)

LoD_summary<- Alignment_solutions_summary_3 %>%
  #combine the peak areas of different ions belonging to the same compound
  group_by(Name, Alignment.ID) %>%
  summarise(slope = linear_regression(Area, Conc, SN, SD_Area)$slope,
            residuals= linear_regression(Area, Conc, SN, SD_Area)$residuals,
            intercept = linear_regression(Area, Conc, SN, SD_Area)$intercept,
            RT = mean(Average.Rt.min.),  R2=linear_regression(Area, Conc, SN, SD_Area)$R2,
            LOD1=linear_regression(Area, Conc, SN, SD_Area)$lod1,
            LOD3=linear_regression(Area, Conc, SN, SD_Area)$lod3,
            LOD2=linear_regression(Area, Conc, SN, SD_Area)$lod2,
            LOD4=linear_regression(Area, Conc, SN, SD_Area)$lod4, 
            max_residuals=linear_regression(Area, Conc, SN, SD_Area)$max_residuals
  )
LoD_summary <- LoD_summary %>%
  mutate(LogLOD1=log10(LOD1)) %>%
  mutate(LogLOD2=log10(LOD2)) %>% 
  mutate(LogLOD4=log10(LOD4))  %>%
  mutate(LogLOD3=log10(LOD3)) %>%
  mutate(LogSlope=log10(slope))
LoD_summary <- LoD_summary %>% filter(R2>0.90) %>% filter(max_residuals<5)

setwd("C:/Users/amso1873/OneDrive - Kruvelab/Amina/manuscript 3/LoD new runs/plots for manuscript")
font <- choose_font("Raleway")
fontsize <- 12
basecolor <- "#7f7f7f"
my_theme <- theme(plot.background = element_blank(),
                  panel.background = element_blank(), axis.line=element_line( color=basecolor,  size=0.5),
                  legend.background = element_blank(),
                  legend.title = element_text(),
                  legend.position = c(0.9,0.8),
                  aspect.ratio = 1, text = element_text(family = font, size=fontsize, color=basecolor),
                  plot.title = element_text(color=basecolor, size=14, face="bold"),
                  legend.text=element_text(family = font, size=fontsize, color=basecolor),
                  axis.text=element_text(family = font, size=fontsize, color=basecolor),
                  panel.grid.major =  element_blank(),
                  panel.grid.minor=element_blank(),
                  panel.border = element_blank())
LoD_summary<- LoD_summary %>% filter(Name!="Acephate") %>% filter(Name!="Clotrimazole") %>%
  filter(Name!="Irbesartan") %>% filter(Name!="Irgarol")

p1 <- ggplot(data=LoD_summary) +
  geom_point(mapping = aes(x = slope, y =LOD1, text=Name), color= "#7f7f7f" , size = 4, alpha = 0.5) +
  # xlab(expression("Response factors (M"^ -1*")"))+
  # ylab("LoD"["Cut-Off"]~"(M)")+
  annotation_logticks(color = basecolor)+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p1+my_theme
# ggplotly(p1+my_theme)
# ggsave("slopes versus LOD-S from cut-off approach.svg",
#        width=10,
#        height=10,
#        units="cm")
cor.test(LoD_summary$LogLOD1, LoD_summary$LogSlope, method = "spearman")
# rho is -0.68
# p-value is below 2.2e-16

p2 <- ggplot(data=LoD_summary) +
  geom_point(mapping = aes(x = slope, y =LOD2, text=Name), color= "#7f7f7f" , size = 4, alpha = 0.5) +
  # xlab(expression("Response factors (M"^ -1*")"))+
  # ylab("LoD"["A2"]~"(M)")+
  # annotation_logticks(color = basecolor)+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p2+my_theme
ggplotly(p2+my_theme)
# ggsave("slopes versus LOD-S from SN approach.svg",
#        width=10,
#        height=10,
#        units="cm")
cor.test(LoD_summary$LogLOD2, LoD_summary$LogSlope, method = "spearman")
# rho is -0.59
# p-value is 8.813e-10

p3 <- ggplot(data=LoD_summary) +
  geom_point(mapping = aes(x = slope, y =LOD3), color= "#7f7f7f" , size = 4, alpha = 0.5) +
  xlab(expression("Response factors (M"^ -1*")"))+
  ylab("LoD"["A4"]~"(M)")+
  annotation_logticks(color = basecolor)+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p3+my_theme
# ggsave("slopes versus LOD-S from Residuals of CC.svg",
#        width=10,
#        height=10,
#        units="cm")
cor.test(LoD_summary$LogLOD4, LoD_summary$LogSlope, method = "spearman")
# rho is -0.51
# p-value is 1.813e-07
p4 <- ggplot(data=LoD_summary) +
  geom_point(mapping = aes(x = slope, y =LOD4, text=Name), color= "#7f7f7f" , size = 4, alpha = 0.5) +
  # xlab(expression("Response factors (M"^ -1*")"))+
  # ylab("LoD"["A3"]~"(M)")+
  annotation_logticks(color = basecolor)+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p4+my_theme
# ggplotly(p4+my_theme)
cor.test(LoD_summary$LogLOD3, LoD_summary$LogSlope, method = "spearman")
# ggsave("slopes versus LOD-S from SD of triplicates.svg",
#        width=10,
#        height=10,
#        units="cm")

#----- wastewater samples--------- 
setwd("C:/Users/amso1873/OneDrive - Kruvelab/Amina/manuscript 3/LoD new runs")
Alignment_solutions_WW<- read.delim("Alignment_WW_sample.txt", skip = 4)
Alignment_SN_WW<- read.delim("Alignment_WW_sample_SN.txt", skip = 4)
Alignment_solutions_WW <- Alignment_solutions_WW %>% 
  filter(Adduct.type == "[M+H]+") 

suspect_list <- read_delim("C:/transferred_data/amso1873/Downloads/inclusion_list_Amina.csv", delim=",") %>% 
  filter(Comment!="Dazomet")  %>% filter(Comment!="Dimethyl phthalate H")  %>%
  filter(Comment!= "trans-ferulic acid") %>% filter(Comment!="Prometryn")  %>%
  filter(Comment!= "atenolol") %>% filter(Comment!="Imazalil")
suspect_list <- suspect_list %>% mutate(Average.Mz=round(suspect_list$`Mass [m/z]`, digits = 2))
Alignment_solutions_WW <- Alignment_solutions_WW %>% mutate(measured_mz=Alignment_solutions_WW$Average.Mz)
Alignment_solutions_WW$Average.Mz <-  round(Alignment_solutions_WW$Average.Mz, digits = 2) 
colnames(suspect_list)[12] <- "Name"
Alignment_solutions_WW <- Alignment_solutions_WW %>% left_join(suspect_list)
Alignment_solutions_summary_WW <- Alignment_solutions_WW %>% filter(Name !="")
Alignment_solutions_summary_WW <- Alignment_solutions_summary_WW %>% 
  select(Alignment.ID, Name, Average.Mz, Adduct.type, Average.Rt.min., MS.MS.spectrum, WW_Blank_1_24012023, WW_Blank_2_24012023,
         WW_Blank_3_24012023, WW2_1_24012023, WW2_2_24012023, WW2_3_24012023, WW3_1_24012023, WW3_2_24012023, WW3_3_24012023, 
         WW4_1_24012023, WW4_2_24012023, WW4_3_24012023, WW5_1_24012023, WW5_2_24012023, WW5_3_24012035, WW6_1_24012023, 
         WW6_2_24012023, WW6_3_24012023, WW7_1_24012023, WW7_2_24012023, WW7_3_24012023, WW8_1_24012023, WW8_2_24012023,
         WW8_3_24012023, WW9_1_24012023, WW9_2_24012023, WW9_3_24012023, WW10_1_24012023, WW10_2_24012023, WW10_3_24012023 , 
         WW11_1_24012023, WW11_2_24012023, WW11_3_24012023, WW12_1_24012023, WW12_2_24012023, WW12_3_24012023, WW13_1_24012023,
         WW13_2_24012023, WW13_3_24012023, WW14_1_24012023, WW14_2_24012023, WW14_3_24012023, WW15_1_24012023, WW15_2_24012023,
         WW15_3_24012023, WW16_1_24012023, WW16_2_24012023, WW16_3_24012023 ) %>% 
  mutate(Diff_mz=abs(Alignment_solutions_summary_WW$`Mass [m/z]`- Alignment_solutions_summary_WW$measured_mz)*1000000/Alignment_solutions_summary_WW$`Mass [m/z]`)
Alignment_solutions_summary_WW <- Alignment_solutions_summary_WW %>% filter(Diff_mz<6)
Alignment_solutions_summary_WW_2 <- Alignment_solutions_summary_WW %>%
  group_by(Alignment.ID, Name,  Average.Mz, Adduct.type, Average.Rt.min., MS.MS.spectrum) %>% 
  summarize(Blank=mean(WW_Blank_1_24012023, WW_Blank_2_24012023,WW_Blank_3_24012023),
            WW2=mean(WW2_1_24012023, WW2_2_24012023, WW2_3_24012023),
            WW3=mean(WW3_1_24012023, WW3_2_24012023, WW3_3_24012023),
            WW4=mean(WW4_1_24012023, WW4_2_24012023, WW4_3_24012023),
            WW5=mean(WW5_1_24012023, WW5_2_24012023, WW5_3_24012035),
            WW6=mean(WW6_1_24012023, WW6_2_24012023, WW6_3_24012023),
            WW7=mean(WW7_1_24012023, WW7_2_24012023, WW7_3_24012023),
            WW8=mean(WW8_1_24012023, WW8_2_24012023, WW8_3_24012023),
            WW9=mean(WW9_1_24012023, WW9_2_24012023, WW9_3_24012023),
            WW10=mean(WW10_1_24012023, WW10_2_24012023, WW10_3_24012023),
            WW11=mean(WW11_1_24012023, WW11_2_24012023, WW11_3_24012023), 
            WW12=mean(WW12_1_24012023, WW12_2_24012023, WW12_3_24012023),
            WW13=mean(WW13_1_24012023, WW13_2_24012023, WW13_3_24012023),
            WW14=mean(WW14_1_24012023, WW14_2_24012023, WW14_3_24012023),
            WW15=mean(WW15_1_24012023, WW15_2_24012023, WW15_3_24012023),
            WW16=mean(WW16_1_24012023, WW16_2_24012023, WW16_3_24012023), 
            WW_SD2=sd(c(WW2_1_24012023, WW2_2_24012023, WW2_3_24012023)),
            WW_SD3=sd(c(WW3_1_24012023, WW3_2_24012023, WW3_3_24012023)),
            WW_SD4=sd(c(WW4_1_24012023, WW4_2_24012023, WW4_3_24012023)),
            WW_SD5=sd(c(WW5_1_24012023, WW5_2_24012023, WW5_3_24012035)),
            WW_SD6=sd(c(WW6_1_24012023, WW6_2_24012023, WW6_3_24012023)),
            WW_SD7=sd(c(WW7_1_24012023, WW7_2_24012023, WW7_3_24012023)),
            WW_SD8=sd(c(WW8_1_24012023, WW8_2_24012023, WW8_3_24012023)),
            WW_SD9=sd(c(WW9_1_24012023, WW9_2_24012023, WW9_3_24012023)),
            WW_SD10=sd(c(WW10_1_24012023, WW10_2_24012023, WW10_3_24012023)),
            WW_SD11=sd(c(WW11_1_24012023, WW11_2_24012023, WW11_3_24012023)), 
            WW_SD12=sd(c(WW12_1_24012023, WW12_2_24012023, WW12_3_24012023)),
            WW_SD13=sd(c(WW13_1_24012023, WW13_2_24012023, WW13_3_24012023)),
            WW_SD14=sd(c(WW14_1_24012023, WW14_2_24012023, WW14_3_24012023)),
            WW_SD15=sd(c(WW15_1_24012023, WW15_2_24012023, WW15_3_24012023)),
            WW_SD16=sd(c(WW16_1_24012023, WW16_2_24012023, WW16_3_24012023)))%>% ungroup() 
Alignment_solutions_summary_WW_3 <- Alignment_solutions_summary_WW_2  %>% select(-c( Name, Average.Mz, Adduct.type, Average.Rt.min., MS.MS.spectrum))
Alignment_solutions_summary_WW_3_SD <-Alignment_solutions_summary_WW_3 %>% select(-c(Blank, WW2, WW3, WW4, WW5, WW6, WW7, WW8, WW9, WW10, WW11, WW12, WW13, WW14, WW15, WW16)) 
Alignment_solutions_summary_WW_3_SD <-Alignment_solutions_summary_WW_3_SD %>% gather(key = "Solution", value = "SD_Area", 2:16) 
Alignment_solutions_summary_WW_3_SD$Solution[Alignment_solutions_summary_WW_3_SD$Solution == 'WW_SD2'] <- 'WW2'
Alignment_solutions_summary_WW_3_SD$Solution[Alignment_solutions_summary_WW_3_SD$Solution == 'WW_SD3'] <- 'WW3'
Alignment_solutions_summary_WW_3_SD$Solution[Alignment_solutions_summary_WW_3_SD$Solution == 'WW_SD4'] <- 'WW4'
Alignment_solutions_summary_WW_3_SD$Solution[Alignment_solutions_summary_WW_3_SD$Solution == 'WW_SD5'] <- 'WW5'
Alignment_solutions_summary_WW_3_SD$Solution[Alignment_solutions_summary_WW_3_SD$Solution == 'WW_SD6'] <- 'WW6'
Alignment_solutions_summary_WW_3_SD$Solution[Alignment_solutions_summary_WW_3_SD$Solution == 'WW_SD7'] <- 'WW7'
Alignment_solutions_summary_WW_3_SD$Solution[Alignment_solutions_summary_WW_3_SD$Solution == 'WW_SD8'] <- 'WW8'
Alignment_solutions_summary_WW_3_SD$Solution[Alignment_solutions_summary_WW_3_SD$Solution == 'WW_SD9'] <- 'WW9'
Alignment_solutions_summary_WW_3_SD$Solution[Alignment_solutions_summary_WW_3_SD$Solution == 'WW_SD10'] <- 'WW10'
Alignment_solutions_summary_WW_3_SD$Solution[Alignment_solutions_summary_WW_3_SD$Solution == 'WW_SD11'] <- 'WW11'
Alignment_solutions_summary_WW_3_SD$Solution[Alignment_solutions_summary_WW_3_SD$Solution == 'WW_SD12'] <- 'WW12'
Alignment_solutions_summary_WW_3_SD$Solution[Alignment_solutions_summary_WW_3_SD$Solution == 'WW_SD13'] <- 'WW13'
Alignment_solutions_summary_WW_3_SD$Solution[Alignment_solutions_summary_WW_3_SD$Solution == 'WW_SD14'] <- 'WW14'
Alignment_solutions_summary_WW_3_SD$Solution[Alignment_solutions_summary_WW_3_SD$Solution == 'WW_SD15'] <- 'WW15'
Alignment_solutions_summary_WW_3_SD$Solution[Alignment_solutions_summary_WW_3_SD$Solution == 'WW_SD16'] <- 'WW16'

Alignment_solutions_summary_WW_3 <-Alignment_solutions_summary_WW_3[1:17] 
# Alignment_solutions_summary_WW_3 <-Alignment_solutions_summary_WW_3 %>% filter(Blank==0)
Alignment_solutions_summary_WW_3 <-Alignment_solutions_summary_WW_3 %>% mutate(Ratio=WW2/Blank)
Alignment_solutions_summary_WW_3 <-Alignment_solutions_summary_WW_3 %>% filter(Ratio>5)
Alignment_solutions_summary_WW_3 <-Alignment_solutions_summary_WW_3 %>% select(-Ratio)
Alignment_solutions_summary_WW_3 <-Alignment_solutions_summary_WW_3 %>%
  mutate(WW2=WW2-Blank) %>%
  mutate(WW3=WW3-Blank) %>%
  mutate(WW4=WW4-Blank) %>%
  mutate(WW5=WW5-Blank) %>%
  mutate(WW6=WW6-Blank) %>%
  mutate(WW7=WW7-Blank) %>%
  mutate(WW8=WW8-Blank) %>%
  mutate(WW9=WW9-Blank) %>%
  mutate(WW10=WW10-Blank) %>%
  mutate(WW11=WW11-Blank) %>%
  mutate(WW12=WW12-Blank) %>%
  mutate(WW13=WW13-Blank) %>%
  mutate(WW14=WW14-Blank) %>%
  mutate(WW15=WW15-Blank) %>%
  mutate(WW16=WW16-Blank)
Alignment_solutions_summary_WW_3 <-Alignment_solutions_summary_WW_3 %>% select(-Blank)
Alignment_solutions_summary_WW_3 <-Alignment_solutions_summary_WW_3 %>%
  gather(key = "Solution", value = "Area", 2:16) 
Alignment_SN_WW <- Alignment_SN_WW %>%
  group_by(Alignment.ID) %>% 
  summarize(WW2=mean(WW2_1_24012023, WW2_2_24012023, WW2_3_24012023),
            WW3=mean(WW3_1_24012023, WW3_2_24012023, WW3_3_24012023),
            WW4=mean(WW4_1_24012023, WW4_2_24012023, WW4_3_24012023),
            WW5=mean(WW5_1_24012023, WW5_2_24012023, WW5_3_24012035),
            WW6=mean(WW6_1_24012023, WW6_2_24012023, WW6_3_24012023),
            WW7=mean(WW7_1_24012023, WW7_2_24012023, WW7_3_24012023),
            WW8=mean(WW8_1_24012023, WW8_2_24012023, WW8_3_24012023),
            WW9=mean(WW9_1_24012023, WW9_2_24012023, WW9_3_24012023),
            WW10=mean(WW10_1_24012023, WW10_2_24012023, WW10_3_24012023),
            WW11=mean(WW11_1_24012023, WW11_2_24012023, WW11_3_24012023), 
            WW12=mean(WW12_1_24012023, WW12_2_24012023, WW12_3_24012023),
            WW13=mean(WW13_1_24012023, WW13_2_24012023, WW13_3_24012023),
            WW14=mean(WW14_1_24012023, WW14_2_24012023, WW14_3_24012023),
            WW15=mean(WW15_1_24012023, WW15_2_24012023, WW15_3_24012023),
            WW16=mean(WW16_1_24012023, WW16_2_24012023, WW16_3_24012023))%>% ungroup()
Alignment_SN_WW<-Alignment_SN_WW %>% 
  gather(key = "Solution", value = "SN", 2:16)
Alignment_solutions_summary_WW_3 <-Alignment_solutions_summary_WW_3 %>% left_join(Alignment_SN_WW)
Alignment_solutions_summary_WW_3 <-Alignment_solutions_summary_WW_3 %>% left_join(Alignment_solutions_summary_WW_3_SD)
library(readxl)
Concentrations <- read_excel("C:/transferred_data/amso1873/Desktop/LOD mix 29122022 - Concentrations.xlsx", 
                             sheet = "Sheet1")
Concentrations <- Concentrations %>%select(- `Cf (ug/L)`)
Concentrations <- Concentrations %>% gather(key="Solution", value = "Conc", 2:17)
Concentrations <- Concentrations %>% left_join(suspect_list %>% select(c(Name, `Mass [m/z]`)))
Concentrations <- Concentrations %>% mutate(Molar_mass = `Mass [m/z]`-1.007)
Names_list =Alignment_solutions_summary_WW_2 %>% select(c(Name , Alignment.ID , Average.Mz, Average.Rt.min., MS.MS.spectrum))
Concentrations <- Concentrations %>% mutate(Conc=(Conc*10^(-6))/Molar_mass)
Alignment_solutions_summary_WW_3 = Alignment_solutions_summary_WW_3 %>% left_join(Names_list)
Concentrations <- Concentrations %>% select(c(Name, Solution, Conc)) %>% drop_na()
Alignment_solutions_summary_WW_3 = Alignment_solutions_summary_WW_3 %>% left_join(Concentrations)
Alignment_solutions_summary_WW_3 = Alignment_solutions_summary_WW_3 %>% drop_na()

Alignment_solutions_summary_3 <- Alignment_solutions_summary_3 %>% select(Name, Average.Rt.min.) %>% unique()
Alignment_solutions_summary_3 <- Alignment_solutions_summary_3 %>% mutate(Checked=TRUE)
Alignment_solutions_summary_3$Average.Rt.min.<- round(Alignment_solutions_summary_3$Average.Rt.min., digits = 1)
Alignment_solutions_summary_WW_3$Average.Rt.min.<- round(Alignment_solutions_summary_WW_3$Average.Rt.min., digits = 1)
Alignment_solutions_summary_WW_3 <- Alignment_solutions_summary_WW_3 %>% left_join(Alignment_solutions_summary_3)
Alignment_solutions_summary_WW_3 <- Alignment_solutions_summary_WW_3 %>% filter(Checked==TRUE)
Alignment_solutions_summary_ww_3 <- Alignment_solutions_summary_WW_3 %>% filter(Name!="lufenuron")
Alignment_solutions_summary_WW_3 <- Alignment_solutions_summary_WW_3 %>% ungroup()
Alignment_solutions_summary_WW_3 <- unique(Alignment_solutions_summary_WW_3) 
Alignment_solutions_summary_WW_3 <- Alignment_solutions_summary_WW_3 %>% filter(Area !=0)
Alignment_solutions_summary_WW_3 <- Alignment_solutions_summary_WW_3 %>% filter(Solution!="WW2") %>% filter(Solution!="WW3")
Alignment_solutions_summary_WW_3 <- Alignment_solutions_summary_WW_3 %>% filter(Area>0)
LoD_summary_WW <- Alignment_solutions_summary_WW_3 %>%
  #combine the peak areas of different ions belonging to the same compound
  group_by(Name, Alignment.ID) %>%
  summarise(slope = linear_regression(Area, Conc, SN, SD_Area)$slope,
            residuals= linear_regression(Area, Conc, SN, SD_Area)$residuals,
            intercept = linear_regression(Area, Conc, SN, SD_Area)$intercept,
            RT = mean(Average.Rt.min.),  R2=linear_regression(Area, Conc, SN, SD_Area)$R2,
            LOD1=linear_regression(Area, Conc, SN, SD_Area)$lod1, 
            LOD3=linear_regression(Area, Conc, SN, SD_Area)$lod3,
            LOD2=linear_regression(Area, Conc, SN, SD_Area)$lod2,
            LOD4=linear_regression(Area, Conc, SN, SD_Area)$lod4, 
            max_residuals=linear_regression(Area, Conc, SN, SD_Area)$max_residuals
  )

LoD_summary_WW <- LoD_summary_WW %>% 
  filter(R2>0.90) %>% 
  mutate(LogSlope=log10(slope)) %>%
  mutate(LogLOD1=log10(LOD1)) %>%
  mutate(LogLOD2=log10(LOD2)) %>%
  mutate(LogLOD3=log10(LOD3)) %>%
  mutate(LogLOD4=log10(LOD4)) %>% 
  filter(max_residuals<5)


# ------- Comparison of LoD accros different approaches-------
Comparison_LOD_SS <- LoD_summary %>% select(c(Name, LOD1, LOD2, LOD4))
Comparison_LOD_SS <- Comparison_LOD_SS %>% mutate(Ratio1_2 = LOD1/LOD2) %>%
  mutate(Ratio1_3 = LOD4/LOD1)

Comparison_LOD_WW <- LoD_summary_WW %>% select(c(Name, LOD1, LOD2, LOD4))
Comparison_LOD_WW <- Comparison_LOD_WW %>% mutate(Ratio1_2 = LOD1/LOD2) %>%
  mutate(Ratio1_3 = LOD4/LOD1)


p1_2 <- ggplot() +
  geom_point(mapping = aes(x = LoD_summary$LOD1, y = LoD_summary$LOD2), color= "#7f7f7f" , size = 4, alpha = 0.5) +
  geom_point(mapping = aes(x = LoD_summary_WW$LOD1, y = LoD_summary_WW$LOD2), color= "#0070C0" , size = 4, alpha = 0.5) +
  xlab(expression("LoD"["Cut-Off"]~"(M)"))+
  ylab(expression("LoD"["S/N Extrapoltion"]~"(M)"))+
  annotation_logticks(color = basecolor)+
  scale_x_log10(limits = c(1e-12, 1e-6), labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(limits = c(1e-12, 1e-6), labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p1_2+my_theme
ggsave("LOD1 vs LOD2.svg", height =4, width=4)
# ggplotly(p1_2+my_theme)
p1_3<- ggplot() +
  geom_point(mapping = aes(x = LoD_summary$LOD1, y =LoD_summary$LOD4), color= "#7f7f7f" , size = 4, alpha = 0.5) +
  geom_point(mapping = aes(x = LoD_summary_WW$LOD1, y =LoD_summary_WW$LOD4), color= "#0070C0" , size = 4, alpha = 0.5) +
  xlab(expression("LoD"["Cut-Off"]~"(M)"))+
  ylab(expression("LoD"["Standard deviation"]~"(M)"))+
  annotation_logticks(color = basecolor)+
  scale_x_log10(limits = c(1e-12, 1e-6), labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(limits = c(1e-12, 1e-6), labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p1_3+my_theme
ggsave("LOD1 vs LOD3.svg", height =4, width=4)
# ggplotly(p1_3+my_theme)
p1_4 <- ggplot() +
  geom_point(mapping = aes(x = LoD_summary$LOD1, y =LoD_summary$LOD3), color= "#7f7f7f" , size = 4, alpha = 0.5) +
  geom_point(mapping = aes(x = LoD_summary_WW$LOD1, y =LoD_summary_WW$LOD3), color= "#0070C0" , size = 4, alpha = 0.5) +
  xlab(expression("LoD"["Cut-Off"]~"(M)"))+
  ylab(expression("LoD"["Residuals"]~"(M)"))+
  annotation_logticks(color = basecolor)+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p1_4+my_theme
plot_grid(p1_2+my_theme , p1_3+my_theme, p1_4+my_theme, 
          labels =c("a) Cut-off vs S/N extrapolation ", "b) Cut-off vs Standard deviation ", "c) Cut-off vs Residuals"),
          align = "h",
          hjust =  0.1, vjust = 1.5 ,
          label_x = 0.3,
          label_y = 0.999,
          nrow = 1,
          label_size = 12,
          label_fontfamily = font,
          label_colour = basecolor,
          label_fontface = "plain")
ggsave("Figure 1 for manuscript 3.svg", height = 4, width=10)
LoD_summary_WW <- LoD_summary_WW %>% filter(Name!="Monocrotophos") %>%
  filter(Name!="Cefoperazone") %>% filter(Name!="Chlorpyrifos") %>% filter(Name!="Rimsulfuron") %>% 
  filter(Name!="Methamidophos") %>% filter(Name!="Nicosulfuron") %>% filter(Name!="Tetraethylammonium")  %>% 
  filter(Name!="4-dimethylaminopyridine") 

setwd("C:/Users/amso1873/OneDrive - Kruvelab/Amina/manuscript 3/LoD new runs/plots for manuscript")
p_W1 <- ggplot() +
  geom_point(mapping = aes(x = LoD_summary$slope, y =LoD_summary$LOD1, text=LoD_summary$Name),  color= "#7f7f7f" , size = 4, alpha = 0.5) +
  geom_point(mapping = aes(x = LoD_summary_WW$slope, y =LoD_summary_WW$LOD1, text=LoD_summary_WW$Name),  color= "#0070C0" , size = 4, alpha = 0.5) +
  xlab(expression("Response factors (M"^-1*")"))+
  ylab(expression("LoD"["Cut-off"]~"(M)"))+
  annotation_logticks(color = basecolor)+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p_W1+my_theme
cor.test(LoD_summary_WW$LogSlope, LoD_summary_WW$LogLOD1, method= "spearman")
# rho is -0.8113997 
# p-value is 2.2e-16
# ggplotly(p_W1+my_theme)
# ggsave("slopes versus LOD-S from cut-off approach - ww.svg",
#        width=10,
#        height=10,
#        units="cm")
p_W2 <- ggplot() +
  geom_point(mapping = aes(x = LoD_summary$slope, y =LoD_summary$LOD2, text=LoD_summary$Name),  color= "#7f7f7f" , size = 4, alpha = 0.5) +
  geom_point(mapping = aes(x = LoD_summary_WW$slope, y =LoD_summary_WW$LOD2, text=LoD_summary_WW$Name),  color= "#0070C0" , size = 4, alpha = 0.5) +
  xlab(expression("Response factors (M"^-1*")"))+
  ylab(expression("LoD"["S/N Extrapolation"]~"(M)"))+
  annotation_logticks(color = basecolor)+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p_W2+my_theme
# ggplotly(p_W2+my_theme)
# ggsave("slopes versus LOD-S from SN approach - WW.svg",
#        width=10,
#        height=10,
#        units="cm")
cor.test(LoD_summary_WW$LogSlope, LoD_summary_WW$LogLOD2, method= "spearman")
# rho is -0.8105339
# p-value is 2.2e-16

p_W3<- ggplot(data=) +
  geom_point(mapping = aes(x = LoD_summary$slope, y =LoD_summary$LOD3),  color= "#7f7f7f" , size = 4, alpha = 0.5) +
  geom_point(mapping = aes(x = LoD_summary_WW$slope, y =LoD_summary_WW$LOD3),  color= "#0070C0" , size = 4, alpha = 0.5) +
  xlab(expression("Response factors (M"^-1*")"))+
  ylab("LoD"["A4"]~"(M)")+
  annotation_logticks(color = basecolor)+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p_W3+my_theme
cor.test(LoD_summary_WW$LogSlope, LoD_summary_WW$LogLOD3, method= "spearman")
# ggsave("slopes versus LOD-S from approach 4 - WW.svg",
#        width=10,
#        height=10,
#        units="cm")
p_W4 <- ggplot() +
  geom_point(mapping = aes(x = LoD_summary$slope, y =LoD_summary$LOD4),  color= "#7f7f7f"  , size = 4, alpha = 0.5) +
  geom_point(mapping = aes(x = LoD_summary_WW$slope, y =LoD_summary_WW$LOD4),  color= "#0070C0" , size = 4, alpha = 0.5) +
  xlab(expression("Response factors (M"^ -1*")"))+
  ylab(expression("LoD"["Standard deviation"]~"(M)"))+
  annotation_logticks(color = basecolor)+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p_W4+my_theme
# ggplotly(p_W4+my_theme)
cor.test(LoD_summary_WW$LogSlope, LoD_summary_WW$LogLOD4, method= "spearman")
# rho is -0.7420635
# p-value is lower 2.2e-16
ggsave("slopes versus LOD-S from SD of triplicate approach - ww.svg",
       width=10,
       height=10,
       units="cm")

plot_grid(p_W1+my_theme , p_W2+my_theme, p_W4+my_theme, 
          labels =c("a) Cut-off ", "b) S/N extrapolation", "c) Standard deviation"),
          align = "h",
          hjust =  0.1, vjust = 1.5 ,
          label_x = 0.3,
          label_y = 0.999,
          nrow = 1,
          label_size = 12,
          label_fontfamily = font,
          label_colour = basecolor,
          label_fontface = "plain")
ggsave("Figure 2 for manuscript 3.svg", height = 4, width=10)
# ------- Predictions of IE-s-------- 
# -----with experimental retention times-----
SMILES <- read_excel("C:/transferred_data/amso1873/Desktop/SMILES_LIST_LoD_COMPOUNDS.xlsx")
DESC <-  read_delim('C:/transferred_data/amso1873/Desktop/descriptors_for_LoD_compounds.csv',
                    delim = ",",
                    col_names = TRUE,
                    trim_ws = TRUE)
DESC <- DESC %>% left_join(SMILES)
fn_viscosity <- function(organic,organic_modifier){
  viscosity <- case_when(
    organic_modifier == 1~ (-0.000103849885417527)*organic^2+0.00435719229180079*organic+0.884232851261593,
    organic_modifier == 2 ~ (-0.00035908)*organic^2+0.031972067*organic+0.90273943)
  return(viscosity)
}

#Pindpinevuse leidmine
fn_surface_tension <- function(organic,organic_modifier){
  surface_tension <- case_when(
    organic_modifier == 1 ~ 71.76-2.906*71.76*(organic/100)+(7.138*27.86+2.906*71.76-71.76)*(organic/100)^2+(27.86-7.138*27.86)*(organic/100)^3,
    organic_modifier == 2 ~ 71.76-2.245*71.76*(organic/100)+(5.625*22.12+2.245*71.76-71.76)*(organic/100)^2+(22.12-5.625*22.12)*(organic/100)^3)
  return(surface_tension)
}
#polaarsus indeksi leidmine
fn_polarity_index <- function(organic,organic_modifier){
  polarity_index <- case_when(
    organic_modifier == 1 ~ (organic/100)*5.1+((100-organic)/100)*10.2,
    organic_modifier == 2 ~ (organic/100)*5.1+((100-organic)/100)*10.2)
  return(polarity_index)
}
#calcualting org %
fn_organic_percentage <- function(eluent_parameters,ret_time){
  ApproxFun <- approxfun(x = eluent_parameters$time, y = eluent_parameters$B)
  organic <- ApproxFun(ret_time)
  return(organic)
}
setwd("C:/transferred_data/amso1873/OneDrive - Kruvelab/files")
descs_pos <-  read_rds("ESIpos_model_descs_191116.rds")
regressor_pos <- read_rds("ESIpos_model_191116.rds")
# #reading the gradient program for reversed phase LC
eluent_parameters <- read_delim('C:/transferred_data/amso1873/Desktop/eluent_LoD.csv',
                                delim = ",",
                                col_names = TRUE)
LoD_summary_WW <- LoD_summary_WW %>% left_join(DESC)%>% drop_na()

names(LoD_summary_WW)[names(LoD_summary_WW) == "RT"] <- "ret_time"
LoD_summary_WW<- data.frame(LoD_summary_WW ,"NH4"=0)
LoD_summary_WW$pH.aq. <- 2.7
LoD_summary_WW$Buffer <- as.factor(1)
LoD_summary_WW$organic_modifier <- as.factor(1)

names(LoD_summary_WW)[names(LoD_summary_WW) == "ASP.0"] <- "ASP-0"
names(LoD_summary_WW)[names(LoD_summary_WW) == "ASP.1"] <- "ASP-1"
names(LoD_summary_WW)[names(LoD_summary_WW) == "ASP.2"] <- "ASP-2"
names(LoD_summary_WW)[names(LoD_summary_WW) == "ASP.3"] <- "ASP-3"
names(LoD_summary_WW)[names(LoD_summary_WW) == "ASP.4"] <- "ASP-4"
names(LoD_summary_WW)[names(LoD_summary_WW) == "ASP.5"] <- "ASP-5"
names(LoD_summary_WW)[names(LoD_summary_WW) == "ASP.6"] <- "ASP-6"
names(LoD_summary_WW)[names(LoD_summary_WW) == "ASP.7"] <- "ASP-7"
names(LoD_summary_WW)[names(LoD_summary_WW) == "AVP.2"] <- "AVP-2"
names(LoD_summary_WW)[names(LoD_summary_WW) == "AVP.6"] <- "AVP-6"
names(LoD_summary_WW)[names(LoD_summary_WW) == "BCUTc.1h"] <- "BCUTc-1h"
names(LoD_summary_WW)[names(LoD_summary_WW) == "BCUTc.1l"] <- "BCUTc-1l"
names(LoD_summary_WW)[names(LoD_summary_WW) == "BCUTw.1h"] <- "BCUTw-1h"
names(LoD_summary_WW)[names(LoD_summary_WW) == "BCUTp.1l"] <- "BCUTp-1l"
names(LoD_summary_WW)[names(LoD_summary_WW) == "BCUTp.1h"] <- "BCUTp-1h"
names(LoD_summary_WW)[names(LoD_summary_WW) == "BCUTw.1l"] <- "BCUTw-1l"
names(LoD_summary_WW)[names(LoD_summary_WW) == "SCH.5"] <- "SCH-5"
names(LoD_summary_WW)[names(LoD_summary_WW) == "WTPT.4"] <- "WTPT-4"
names(LoD_summary_WW)[names(LoD_summary_WW) == "WTPT.5"] <- "WTPT-5"
LoD_summary_WW_2 <- LoD_summary_WW %>% 
  mutate(
    organic = fn_organic_percentage(eluent_parameters,ret_time),
    viscosity =  fn_viscosity(organic,organic_modifier),
    surface_tension = fn_surface_tension(organic,organic_modifier),
    polarity_index = fn_polarity_index(organic,organic_modifier)) %>%
  #though all columns are selected with "everything()" in the end, the order of the columns is changed
  dplyr::select( SMILES, ret_time , organic_modifier, organic, pH.aq., Buffer,  viscosity,surface_tension,polarity_index, descs_pos)
prediction_set_model_pos <- LoD_summary_WW_2 %>%
  mutate(logIE_pred = 0)

prediction <-  predict(regressor_pos, newdata = prediction_set_model_pos)
prediction <- prediction$aggregate
prediction_set_model_pos$logIE_pred <- prediction
prediction_set_model_pos <- prediction_set_model_pos %>%
  mutate(logIE_pred = prediction) %>%
  select(SMILES,logIE_pred, everything())
LoD_summary_WW <- LoD_summary_WW %>% left_join(prediction_set_model_pos)
LoD_summary_WW <- LoD_summary_WW %>% select(Name, SMILES, Alignment.ID, everything())
setwd("C:/Users/amso1873/OneDrive - Kruvelab/Amina/manuscript 3/LoD new runs/plots for manuscript")
isotopedistribution <- function(smiles){
  #convert SMILES to chemical formula
  molecule <- parse.smiles(smiles)[[1]]
  formula <- get.mol2formula(molecule,charge=0)
  formula <- formula@string
  
  # Chemical formula to isotope distribution
  data(isotopes)
  pattern<-isopattern(isotopes,
                      formula,
                      threshold=0.1,
                      charge = +1,
                      emass = 0.00054858,
                      plotit=FALSE,
                      algo=1)
  isotopes <- as.data.frame(pattern[[1]])
  isotope_dist <- as.numeric(sum(isotopes$abundance))
  return(isotope_dist)
}
LoD_summary_WW <- LoD_summary_WW %>% mutate(IC=isotopedistribution(LoD_summary_WW$SMILES))
LoD_summary_WW <- LoD_summary_WW %>% mutate(IE_Corrected_ex_RT= (10^(LoD_summary_WW$logIE_pred)/LoD_summary_WW$IC)*100)

# -----With predicted retention times----- 

# library(tidyverse)
# library(stringr)
# library(rjson)
# 
# 
# 
# SMILES_list<-as.data.frame(SMILES[2])
# PaDEL = function(SMILES_list) {
#   command = "java -jar descriptor-cli-0.1a-SNAPSHOT-all.jar"
#   SMILES = SMILES_list$SMILES[1]
#   command_final = paste(command, SMILES, sep =" ")
#   descs = calculation_by_SMILES(SMILES, command_final)
#   for (i in 2:189) {
#     SMILES = SMILES_list$SMILES[i]
#     command_final = paste(command, SMILES, sep =" ")
#     descs_this_smiles = calculation_by_SMILES(SMILES, command_final)
#     descs= descs %>%
#       bind_rows(descs_this_smiles)
#   }
#   
#   return(descs)
# }
# 
# 
# calculation_by_SMILES = function(smiles, command_final) {
#   tryCatch(
#     {
#       javaOutput = system(command_final, intern = TRUE)
#       output_string = str_split(javaOutput, pattern = " ")
#       descs_json = paste(output_string[1][[1]], output_string[2][[1]], output_string[3][[1]], output_string[4][[1]], sep = "")
#       descs_this_smiles = fromJSON(descs_json)
#       descs_this_smiles = data.frame(descs_this_smiles)
#       descs_this_smiles = descs_this_smiles %>%
#         mutate(SMILES = smiles)
#       return(descs_this_smiles)
#     },
#     error = function(e) {
#       return(tibble())
#     }
#   )
# }
# 
# 
# DESC <- PaDEL(SMILES_list)
# write_delim(DESC,
#             "descriptors_for_LoD_compounds.csv",
#             delim = ",")

LoD_summary_WW <- LoD_summary_WW %>% mutate(pH=2.7)
LoD_summary_WW <- LoD_summary_WW %>% mutate(Column=as.factor(2))
LoD_summary_WW<-LoD_summary_WW %>% mutate(Buffer=as.factor(1))
LoD_summary_WW <- LoD_summary_WW %>% mutate(Organic_modifier=as.factor(1))
LoD_summary_WW <- LoD_summary_WW %>% replace(is.na(.), 0)

library(rcdk)
fn_logP <- function(SMILES){
  mol <- parse.smiles(SMILES)[[1]]
  convert.implicit.to.explicit(mol)
  get.tpsa(mol)
  logP <- get.xlogp(mol)
  return(logP)
}
log_P_sum <-c()
for (i in 1:length(LoD_summary_WW$SMILES)){
  print(i)
  log_P_sum <- c(log_P_sum, fn_logP(LoD_summary_WW$SMILES[i]))
}
LoD_summary_WW <- LoD_summary_WW %>% ungroup() %>% mutate(log_P_sum) 
library("ranger")
regressor_MultiConditionRT <- read_rds("C:/Users/amso1873/OneDrive - Kruvelab/Amina/R code/regressorRF for RP 16082021.rds")
names(LoD_summary_WW)[names(LoD_summary_WW) == "ASP-6"] <- "ASP.6"
names(LoD_summary_WW)[names(LoD_summary_WW) == "ASP-7"] <- "ASP.7"
names(LoD_summary_WW)[names(LoD_summary_WW) == "SCH-5"] <- "SCH.5"
names(LoD_summary_WW)[names(LoD_summary_WW) == "WTPT-4"] <- "WTPT.4"
names(LoD_summary_WW)[names(LoD_summary_WW) == "WTPT-5"] <- "WTPT.5"
pred_RT=predict(regressor_MultiConditionRT, newdata = LoD_summary_WW) 
LoD_summary_WW <- LoD_summary_WW %>% mutate(pred_RT)
LoD_summary <- LoD_summary %>% left_join(SMILES)
dataset2 <-  read_delim('C:/Users/amso1873/OneDrive - Kruvelab/Amina/R code/Filter_RP_2010.csv',
                        delim = ",",
                        col_names = TRUE,
                        trim_ws = TRUE) 
dataset2 <- dataset2 %>%
  group_by(Compound_name, Column, Organic_modifier, pH, Buffer, SMILES ) %>%
  summarise(slope = mean(slope),
            RT_Miklos= mean(ret_time))  
dataset2 <- dataset2 %>% filter(pH==2.7)
dataset2 <- dataset2 %>% filter(Organic_modifier=="Acetonitrile")
dataset2 <- dataset2 %>% filter(Buffer=="Formic acid")
rt_for_projection_model <- LoD_summary %>% select(c(Name, SMILES,RT))
dataset2 <- dataset2 %>% left_join(rt_for_projection_model ) %>% drop_na()
p1_RT <- ggplot(data=dataset2) +
  geom_point(mapping = aes(x = RT_Miklos, y =RT), color= "#7f7f7f" , size = 4, alpha = 0.5) +
  xlab("RT on MCRT scale")+
  ylab("RT")+
  theme_bw()
p1_RT+my_theme
ggplotly(p1_RT)
library(mgcv)
regressor1 <- gam(RT~RT_Miklos, family = gaussian(), data=dataset2, method="GCV.Cp", optimizer = c("outer", "newton"))
names(LoD_summary_WW )[names(LoD_summary_WW ) == "pred_RT"] <- "RT_Miklos"
ret_time= predict(regressor1, newdata = LoD_summary_WW ) 
names(LoD_summary_WW)[names(LoD_summary_WW) == "ret_time"] <- "ret_time_exp"
LoD_summary_WW <- LoD_summary_WW %>% mutate(ret_time)
p_RT <- ggplot(data=LoD_summary_WW) +
  geom_point(mapping = aes(x =ret_time_exp, y =ret_time), color= "#7f7f7f" , size = 4, alpha = 0.5) +
  xlab("RT")+
  ylab("Predicted RT")+
  theme_bw()
p_RT+my_theme

LoD_summary_WW<- data.frame(LoD_summary_WW ,"NH4"=0)
LoD_summary_WW$pH.aq. <- 2.7
LoD_summary_WW$Buffer <- as.factor(1)
LoD_summary_WW$organic_modifier <- as.factor(1)
names(LoD_summary_WW)[names(LoD_summary_WW) == "ASP.0"] <- "ASP-0"
names(LoD_summary_WW)[names(LoD_summary_WW) == "ASP.1"] <- "ASP-1"
names(LoD_summary_WW)[names(LoD_summary_WW) == "ASP.2"] <- "ASP-2"
names(LoD_summary_WW)[names(LoD_summary_WW) == "ASP.3"] <- "ASP-3"
names(LoD_summary_WW)[names(LoD_summary_WW) == "ASP.4"] <- "ASP-4"
names(LoD_summary_WW)[names(LoD_summary_WW) == "ASP.5"] <- "ASP-5"
names(LoD_summary_WW)[names(LoD_summary_WW) == "ASP.6"] <- "ASP-6"
names(LoD_summary_WW)[names(LoD_summary_WW) == "ASP.7"] <- "ASP-7"
names(LoD_summary_WW)[names(LoD_summary_WW) == "AVP.2"] <- "AVP-2"
names(LoD_summary_WW)[names(LoD_summary_WW) == "AVP.6"] <- "AVP-6"
names(LoD_summary_WW)[names(LoD_summary_WW) == "BCUTc.1h"] <- "BCUTc-1h"
names(LoD_summary_WW)[names(LoD_summary_WW) == "BCUTc.1l"] <- "BCUTc-1l"
names(LoD_summary_WW)[names(LoD_summary_WW) == "BCUTw.1h"] <- "BCUTw-1h"
names(LoD_summary_WW)[names(LoD_summary_WW) == "BCUTp.1l"] <- "BCUTp-1l"
names(LoD_summary_WW)[names(LoD_summary_WW) == "BCUTp.1h"] <- "BCUTp-1h"
names(LoD_summary_WW)[names(LoD_summary_WW) == "BCUTw.1l"] <- "BCUTw-1l"
names(LoD_summary_WW)[names(LoD_summary_WW) == "SCH.5"] <- "SCH-5"
names(LoD_summary_WW)[names(LoD_summary_WW) == "WTPT.4"] <- "WTPT-4"
names(LoD_summary_WW)[names(LoD_summary_WW) == "WTPT.5"] <- "WTPT-5"
LoD_summary_WW <- LoD_summary_WW %>% select(-c(organic, viscosity, surface_tension, polarity_index, logIE_pred))
LoD_summary_WW_2 <- LoD_summary_WW %>% 
  mutate(
    organic = fn_organic_percentage(eluent_parameters,ret_time),
    viscosity =  fn_viscosity(organic,organic_modifier),
    surface_tension = fn_surface_tension(organic,organic_modifier),
    polarity_index = fn_polarity_index(organic,organic_modifier)) %>%
  #though all columns are selected with "everything()" in the end, the order of the columns is changed
  dplyr::select( SMILES, ret_time , organic_modifier, organic, pH.aq., Buffer,  viscosity,surface_tension,polarity_index, descs_pos)
prediction_set_model_pos_2 <- LoD_summary_WW_2 %>%
  mutate(logIE_pred = 0)

prediction_2 <-  predict(regressor_pos, newdata = prediction_set_model_pos_2)
prediction_2 <- prediction_2$aggregate
prediction_set_model_pos_2$logIE_pred <- prediction_2
prediction_set_model_pos_2 <- prediction_set_model_pos_2 %>%
  mutate(logIE_pred = prediction_2) %>%
  select(SMILES,logIE_pred, everything())
LoD_summary_WW <- LoD_summary_WW %>% left_join(prediction_set_model_pos_2)
LoD_summary_WW <- LoD_summary_WW %>% select(Name, SMILES, Alignment.ID, everything())
setwd("C:/Users/amso1873/OneDrive - Kruvelab/Amina/manuscript 3/LoD new runs/plots for manuscript")
isotopedistribution <- function(smiles){
  #convert SMILES to chemical formula
  molecule <- parse.smiles(smiles)[[1]]
  formula <- get.mol2formula(molecule,charge=0)
  formula <- formula@string
  
  # Chemical formula to isotope distribution
  data(isotopes)
  pattern<-isopattern(isotopes,
                      formula,
                      threshold=0.1,
                      charge = +1,
                      emass = 0.00054858,
                      plotit=FALSE,
                      algo=1)
  isotopes <- as.data.frame(pattern[[1]])
  isotope_dist <- as.numeric(sum(isotopes$abundance))
  return(isotope_dist)
}
LoD_summary_WW <- LoD_summary_WW %>% mutate(IC_pred_RT=isotopedistribution(LoD_summary_WW$SMILES))
LoD_summary_WW <- LoD_summary_WW %>% mutate(IE_Corrected_pred_RT= (10^(LoD_summary_WW$logIE_pred)/LoD_summary_WW$IC)*100)

p_ww_IE1 <- ggplot(data=LoD_summary_WW) +
  geom_point(mapping = aes(x = IE_Corrected_ex_RT , y =LOD1), color= "Red" , size = 4, alpha = 0.5) +
  geom_point(mapping = aes(x = IE_Corrected_pred_RT, y =LOD1), color= "Black" , size = 4, alpha = 0.5) +
  xlab("IE-s")+
  ylab("LoD"["A1"]~"(M)")+
  annotation_logticks(color = basecolor)+
  scale_x_log10( labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p_ww_IE1+my_theme
ggplotly(p_ww_IE1+my_theme)
cor.test(LoD_summary_WW$IE_Corrected_pred_RT, LoD_summary_WW$LogLOD1, method= "spearman")
# rho is -0.436273
# p-value is 0.003735
# ggplotly(p_ww_IE1+my_theme)
# ggsave("IE-s versus LOD-S from cut-off approach - ww.svg",
#        width=10,
#        height=10,
#        units="cm")

p_ww_IE2 <- ggplot(data=LoD_summary_WW) +
  geom_point(mapping = aes(x = IE_Corrected_pred_RT, y =LOD2, text=Name), color= "#0070C0" , size = 4, alpha = 0.5) +
  xlab("IE-s")+
  ylab("LoD"["A2"]~"(M)")+
  annotation_logticks(color = basecolor)+
  scale_x_log10( labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p_ww_IE2+my_theme
# ggplotly(p_ww_IE2+my_theme)
# ggsave("IE-s versus LOD-S from SN approach - ww.svg",
#        width=10,
#        height=10,
#        units="cm")
cor.test(LoD_summary_WW$IE_Corrected_pred_RT, LoD_summary_WW$LogLOD2, method= "spearman")
# rho is -0.3822108 
# p-value is 0.01189

p_ww_IE4 <- ggplot(data=LoD_summary_WW) +
  geom_point(mapping = aes(x = IE_Corrected_pred_RT, y =LOD4, text=Name ),color= "#0070C0" , size = 4, alpha = 0.5) +
  xlab("IE-s")+
  ylab("LoD"["A4"]~"(M)")+
  annotation_logticks(color = basecolor)+
  scale_x_log10( labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p_ww_IE4+my_theme
# ggplotly(p_ww_IE4+my_theme)
# ggsave("IE-s versus LOD-S from SD of triplicates approach - ww.svg",
#        width=10,
#        height=10,
#        units="cm")
cor.test(LoD_summary_WW$IE_Corrected_pred_RT, LoD_summary_WW$LogLOD4, method= "spearman")
# rho is -0.4341589 
# p-value is 0.003921
plot_grid(p_ww_IE1+my_theme , p_ww_IE2+my_theme,  p_ww_IE4+my_theme,
          labels =c("a) Cut-off approach", "b) Extrapolation of S/N",
                    "c) standard deviation-based approach"),
          align = "h",
          hjust =  0.1, vjust = 1.5 ,
          label_x = 0.3,
          label_y = 0.999,
          nrow = 1,
          label_size = 12,
          label_fontfamily = font,
          label_colour = basecolor,
          label_fontface = "plain")
ggsave("IE_vs_ww_Jaanus_model.svg", height = 4, width=12)
# LoD_summary <- LoD_summary %>% left_join(SMILES)

# ----- IE predictions for SS----------------
eluent_parameters <- read_delim('C:/transferred_data/amso1873/Desktop/eluent_LoD.csv',
                                delim = ",",
                                col_names = TRUE)
LoD_summary <- LoD_summary %>% left_join(DESC)%>% drop_na()
names(LoD_summary)[names(LoD_summary) == "RT"] <- "ret_time"
LoD_summary<- data.frame(LoD_summary ,"NH4"=0)
LoD_summary$pH.aq. <- 2.7
LoD_summary$Buffer <- as.factor(1)
LoD_summary$organic_modifier <- as.factor(1)

names(LoD_summary)[names(LoD_summary) == "ASP.0"] <- "ASP-0"
names(LoD_summary)[names(LoD_summary) == "ASP.1"] <- "ASP-1"
names(LoD_summary)[names(LoD_summary) == "ASP.2"] <- "ASP-2"
names(LoD_summary)[names(LoD_summary) == "ASP.3"] <- "ASP-3"
names(LoD_summary)[names(LoD_summary) == "ASP.4"] <- "ASP-4"
names(LoD_summary)[names(LoD_summary) == "ASP.5"] <- "ASP-5"
names(LoD_summary)[names(LoD_summary) == "ASP.6"] <- "ASP-6"
names(LoD_summary)[names(LoD_summary) == "ASP.7"] <- "ASP-7"
names(LoD_summary)[names(LoD_summary) == "AVP.2"] <- "AVP-2"
names(LoD_summary)[names(LoD_summary) == "AVP.6"] <- "AVP-6"
names(LoD_summary)[names(LoD_summary) == "BCUTc.1h"] <- "BCUTc-1h"
names(LoD_summary)[names(LoD_summary) == "BCUTc.1l"] <- "BCUTc-1l"
names(LoD_summary)[names(LoD_summary) == "BCUTw.1h"] <- "BCUTw-1h"
names(LoD_summary)[names(LoD_summary) == "BCUTp.1l"] <- "BCUTp-1l"
names(LoD_summary)[names(LoD_summary) == "BCUTp.1h"] <- "BCUTp-1h"
names(LoD_summary)[names(LoD_summary) == "BCUTw.1l"] <- "BCUTw-1l"
names(LoD_summary)[names(LoD_summary) == "SCH.5"] <- "SCH-5"
names(LoD_summary)[names(LoD_summary) == "WTPT.4"] <- "WTPT-4"
names(LoD_summary)[names(LoD_summary) == "WTPT.5"] <- "WTPT-5"
LoD_summary_2 <- LoD_summary %>%
  mutate(
    organic = fn_organic_percentage(eluent_parameters,ret_time),
    viscosity =  fn_viscosity(organic,organic_modifier),
    surface_tension = fn_surface_tension(organic,organic_modifier),
    polarity_index = fn_polarity_index(organic,organic_modifier)) %>%
  #though all columns are selected with "everything()" in the end, the order of the columns is changed
  dplyr::select( SMILES, ret_time , organic_modifier, organic, pH.aq., Buffer,  viscosity,surface_tension,polarity_index, descs_pos)
prediction_set_model_pos_SS <- LoD_summary_2 %>%
  mutate(logIE_pred = 0)

prediction_SS <-  predict(regressor_pos, newdata = prediction_set_model_pos_SS)
prediction_SS <- prediction_SS$aggregate
prediction_set_model_pos_SS$logIE_pred <- prediction_SS
prediction_set_model_pos_SS <- prediction_set_model_pos_SS %>%
  mutate(logIE_pred = prediction_SS) %>%
  select(SMILES,logIE_pred, everything())
LoD_summary <- LoD_summary %>% left_join(prediction_set_model_pos_SS)
LoD_summary <- LoD_summary %>% select(Name, SMILES, Alignment.ID, everything())
setwd("C:/Users/amso1873/OneDrive - Kruvelab/Amina/manuscript 3/LoD new runs/plots for manuscript")

LoD_summary <- LoD_summary %>% mutate(IC=isotopedistribution(LoD_summary$SMILES))
LoD_summary <- LoD_summary %>% mutate(IE_Corrected_ex_RT= (10^(LoD_summary$logIE_pred)/LoD_summary$IC)*100)
LoD_summary <- LoD_summary %>% left_join(DESC)%>% drop_na()
LoD_summary <- LoD_summary %>% mutate(pH=2.7)
LoD_summary <- LoD_summary %>% mutate(Column=as.factor(2))
LoD_summary<-LoD_summary %>% mutate(Buffer=as.factor(1))
LoD_summary <- LoD_summary %>% mutate(Organic_modifier=as.factor(1))
LoD_summary<- LoD_summary %>% replace(is.na(.), 0)


log_P_sum <-c()
for (i in 1:length(LoD_summary$SMILES)){
  print(i)
  log_P_sum <- c(log_P_sum, fn_logP(LoD_summary$SMILES[i]))
}
LoD_summary <- LoD_summary %>% ungroup() %>% mutate(log_P_sum) 

pred_RT=predict(regressor_MultiConditionRT, newdata = LoD_summary) 
LoD_summary<- LoD_summary %>% mutate(pred_RT)

names(LoD_summary )[names(LoD_summary) == "pred_RT"] <- "RT_Miklos"
ret_time= predict(regressor1, newdata = LoD_summary) 
names(LoD_summary)[names(LoD_summary) == "ret_time"] <- "ret_time_exp"
LoD_summary <- LoD_summary %>% mutate(ret_time)
p_RT <- ggplot(data=LoD_summary) +
  geom_point(mapping = aes(x = ret_time_exp, y =ret_time), color= "#7f7f7f" , size = 4, alpha = 0.5) +
  xlab("RT")+
  ylab("Predicted RT")+
  theme_bw()
p_RT+my_theme

LoD_summary<- data.frame(LoD_summary ,"NH4"=0)
LoD_summary$pH.aq. <- 2.7
LoD_summary$Buffer <- as.factor(1)
LoD_summary$organic_modifier <- as.factor(1)

names(LoD_summary)[names(LoD_summary) == "ASP.0"] <- "ASP-0"
names(LoD_summary)[names(LoD_summary) == "ASP.1"] <- "ASP-1"
names(LoD_summary)[names(LoD_summary) == "ASP.2"] <- "ASP-2"
names(LoD_summary)[names(LoD_summary) == "ASP.3"] <- "ASP-3"
names(LoD_summary)[names(LoD_summary) == "ASP.4"] <- "ASP-4"
names(LoD_summary)[names(LoD_summary) == "ASP.5"] <- "ASP-5"
names(LoD_summary)[names(LoD_summary) == "ASP.6"] <- "ASP-6"
names(LoD_summary)[names(LoD_summary) == "ASP.7"] <- "ASP-7"
names(LoD_summary)[names(LoD_summary) == "AVP.2"] <- "AVP-2"
names(LoD_summary)[names(LoD_summary) == "AVP.6"] <- "AVP-6"
names(LoD_summary)[names(LoD_summary) == "BCUTc.1h"] <- "BCUTc-1h"
names(LoD_summary)[names(LoD_summary) == "BCUTc.1l"] <- "BCUTc-1l"
names(LoD_summary)[names(LoD_summary) == "BCUTw.1h"] <- "BCUTw-1h"
names(LoD_summary)[names(LoD_summary) == "BCUTp.1l"] <- "BCUTp-1l"
names(LoD_summary)[names(LoD_summary) == "BCUTp.1h"] <- "BCUTp-1h"
names(LoD_summary)[names(LoD_summary) == "BCUTw.1l"] <- "BCUTw-1l"
names(LoD_summary)[names(LoD_summary) == "SCH.5"] <- "SCH-5"
names(LoD_summary)[names(LoD_summary) == "WTPT.4"] <- "WTPT-4"
names(LoD_summary)[names(LoD_summary) == "WTPT.5"] <- "WTPT-5"
LoD_summary <- LoD_summary %>% select(-c(organic, viscosity, surface_tension, polarity_index, logIE_pred))

LoD_summary_2 <- LoD_summary %>% 
  mutate(
    organic = fn_organic_percentage(eluent_parameters,ret_time),
    viscosity =  fn_viscosity(organic,organic_modifier),
    surface_tension = fn_surface_tension(organic,organic_modifier),
    polarity_index = fn_polarity_index(organic,organic_modifier)) %>%
  #though all columns are selected with "everything()" in the end, the order of the columns is changed
  dplyr::select( SMILES, ret_time , organic_modifier, organic, pH.aq., Buffer,  viscosity,surface_tension,polarity_index, descs_pos)
prediction_set_model_pos_SS_2 <- LoD_summary_2 %>%
  mutate(logIE_pred = 0)

prediction_SS_2 <-  predict(regressor_pos, newdata = prediction_set_model_pos_SS_2)
prediction_SS_2 <- prediction_SS_2$aggregate
prediction_set_model_pos_SS_2$logIE_pred <- prediction_SS_2
prediction_set_model_pos_SS_2 <- prediction_set_model_pos_SS_2 %>%
  mutate(logIE_pred = prediction_SS_2) %>%
  select(SMILES,logIE_pred, everything())
LoD_summary <- LoD_summary %>% left_join(prediction_set_model_pos_SS_2)
LoD_summary <- LoD_summary %>% select(Name, SMILES, Alignment.ID, everything())
setwd("C:/Users/amso1873/OneDrive - Kruvelab/Amina/manuscript 3/LoD new runs/plots for manuscript")
LoD_summary <- LoD_summary %>% mutate(IC=isotopedistribution(LoD_summary$SMILES))
LoD_summary <- LoD_summary %>% mutate(IE_Corrected_pred_RT= (10^(LoD_summary$logIE_pred)/LoD_summary$IC)*100)

p1_IE <- ggplot(data=) +
  geom_point(mapping = aes(x = LoD_summary$IE_Corrected_pred_RT , y = LoD_summary$LOD1),color= "#7f7f7f" , size = 4, alpha = 0.5) +
  geom_point(mapping = aes(x = LoD_summary_WW$IE_Corrected_pred_RT, y = LoD_summary_WW$LOD1),color= "#0070C0" , size = 4, alpha = 0.5) +
  xlab(expression(italic("IE")~"-s"))+
  ylab(expression("LoD"["Cut-Off"]~"(M)"))+
  annotation_logticks(color = basecolor)+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p1_IE+my_theme
ggplotly(p1_IE+my_theme)
ggsave("IE-s versus LOD-S from cut-off approach.svg",
       width=10,
       height=10,
       units="cm")
cor.test(LoD_summary$IE_Corrected_pred_RT, LoD_summary$LogLOD1, method= "spearman")
# rho is -0.4658591 
# p-value is 2.17e-05
cor.test(LoD_summary$IE_Corrected_ex_RT, LoD_summary$LogLOD1, method= "spearman")

# rho is  -0.46
# p-value is 2.123e-05

p2_IE <- ggplot() +
  geom_point(mapping = aes(x = LoD_summary$IE_Corrected_pred_RT , y = LoD_summary$LOD2),color= "#7f7f7f" , size = 4, alpha = 0.5) +
  geom_point(mapping = aes(x = LoD_summary_WW$IE_Corrected_pred_RT, y = LoD_summary_WW$LOD2),color= "#0070C0" , size = 4, alpha = 0.5) +  # xlab("IE-s")+
  xlab(expression(italic("IE")~"-s"))+
  ylab(expression("LoD"["S/N Extrapoltion"]~"(M)"))+
  annotation_logticks(color = basecolor)+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p2_IE+my_theme
ggplotly(p2_IE+my_theme)
ggsave("IE-s versus LOD-S from SN approach.svg",
       width=10,
       height=10,
       units="cm")
cor.test(LoD_summary$IE_Corrected_pred_RT, LoD_summary$LogLOD2, method= "spearman")
# rho is -0.3564916   
# p-value is 0.001447
cor.test(LoD_summary$IE_Corrected_ex_RT, LoD_summary$LogLOD2, method= "spearman")
# rho is -0.3663552  
# p-value is 0.001
p4_IE <- ggplot() +
  geom_point(mapping = aes(x = LoD_summary$IE_Corrected_pred_RT , y = LoD_summary$LOD4),color= "#7f7f7f" , size = 4, alpha = 0.5) +
  geom_point(mapping = aes(x = LoD_summary_WW$IE_Corrected_pred_RT, y = LoD_summary_WW$LOD4),color= "#0070C0" , size = 4, alpha = 0.5) +  # xlab("IE-s")+
  xlab(expression(italic("IE")~"-s"))+
  ylab(expression("LoD"["Standard deviation"]~"(M)"))+
  annotation_logticks(color = basecolor)+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p4_IE+my_theme
ggplotly(p4_IE+my_theme)
ggsave("IE-s versus LOD-S from SD of triplicates approach.svg",
       width=10,
       height=10,
       units="cm")
cor.test(LoD_summary$IE_Corrected_pred_RT, LoD_summary$LogLOD4, method= "spearman")
# rho is -0.3765222   
# p-value is 0.0007388
cor.test(LoD_summary$IE_Corrected_ex_RT, LoD_summary$LogLOD4, method= "spearman")
# rho is -0.3787984
# p-value is 0.0006826
library(cowplot)
plot_grid(p1_IE+my_theme , p2_IE+my_theme, p4_IE+my_theme, p1_IE+my_theme , p2_IE+my_theme, p4_IE+my_theme,
          labels =c("a)", " ", " ",  "b) Cut-off approach", "c) Extrapolation of S/N","d) Standard deviation"),
          
          align = "h",
          hjust =  0.1, vjust = 1.5 ,
          label_x = 0.3,
          label_y = 0.999,
          nrow = 2,
          label_size = 12,
          label_fontfamily = font,
          label_colour = basecolor,
          label_fontface = "plain")
ggsave("IE_vs_LoD_SS_Jaanus.svg", height = 8, width=12)
library(MASS)
LoD_summary <- LoD_summary %>% mutate(logIE_pred_RT= log10(IE_Corrected_pred_RT)) %>% mutate(logIE_exp_RT= log10(IE_Corrected_ex_RT))
LoD_summary_WW <- LoD_summary_WW %>% mutate(logIE_pred_RT= log10(IE_Corrected_pred_RT)) %>% mutate(logIE_exp_RT= log10(IE_Corrected_ex_RT))

LoD_regressor_SS_Pred_RT <- rlm(LogLOD1 ~ logIE_pred_RT , data= LoD_summary)

saveRDS(LoD_regressor_SS_Pred_RT ,  "LoD_regressor_SS_Pred_RT.rds")

LoD_regressor_WW_Pred_RT <- rlm(LogLOD1 ~ logIE_pred_RT, data= LoD_summary_WW)
saveRDS(LoD_regressor_WW_Pred_RT ,  "LoD_regressor_WW_Pred_RT.rds") 

LoD_regressor_SS_exp_RT <- rlm(LogLOD1 ~ logIE_exp_RT , data= LoD_summary)
saveRDS(LoD_regressor_SS_exp_RT,  "LoD_regressor_SS_exp_RT.rds")

LoD_regressor_WW_exp_RT <- rlm(LogLOD1 ~ logIE_exp_RT, data= LoD_summary_WW)
saveRDS(LoD_regressor_WW_exp_RT ,  "LoD_regressor_WW_exp_RT.rds") 

LoD_summary <- LoD_summary %>% mutate(pred_LoD=10^predict(LoD_regressor_SS_Pred_RT, newdata = LoD_summary) )
p1_SS_application <- ggplot(data=LoD_summary) +
  geom_point(mapping = aes(x = IE_Corrected_pred_RT , y = LOD1),color= "#7f7f7f" , size = 4, alpha = 0.5) +
  geom_line(mapping = aes(x = IE_Corrected_pred_RT , y = pred_LoD),color= "#7f7f7f" , size = 1, alpha = 0.5) +
  geom_line(mapping = aes(x = IE_Corrected_pred_RT , y = 10*pred_LoD), linetype = "dashed", color= "#7f7f7f" , size = 1, alpha = 0.5) +
  geom_line(mapping = aes(x = IE_Corrected_pred_RT , y = pred_LoD/10), linetype = "dashed", color= "#7f7f7f" , size = 1, alpha = 0.5) +
  xlab(expression(italic("IE")~"-s"))+
  ylab(expression("LoD"["Cut-Off"]~"(M)"))+
  annotation_logticks(color = basecolor)+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p1_SS_application+my_theme
ggsave("IE-s versus LOD-S from cut-off approach with robust regression.svg",
       width=10,
       height=10,
       units="cm")

# ------IE predictions using MS2Quant -------
# devtools::install_github("kruvelab/MS2Quant",
#                          ref="main",
#                          INSTALL_opts="--no-multiarch")
library(MS2Quant)
chemicals_SMILES = tibble(SMILES = SMILES$SMILES)

IE_pred = MS2Quant_predict_IE(chemicals_for_IE_prediction = chemicals_SMILES,
                              organic_modifier = "MeCN",
                              eluent = eluent_parameters,
                              pH_aq = 2.7)
data = IE_pred$chemicals_predicted_IEs
names(data)[names(data) == "pred_logIE"] <- "pred_logIE_MS2Quant"
LoD_summary <- LoD_summary %>% left_join(SMILES)
LoD_summary <- LoD_summary %>% left_join(data)

p1_IE_MS2Quant <- ggplot(data=LoD_summary) +
  geom_point(mapping = aes(x = pred_logIE_MS2Quant, y =LogLOD1, text=Name),color= "#0070C0" , size = 4, alpha = 0.5) +
  xlab("IE-s")+
  ylab("LoD-s")+
  # annotation_logticks(color = basecolor)+
  # scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  # scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p1_IE_MS2Quant+my_theme
ggplotly(p1_IE+my_theme)
LoD_regressor_SS_MS2QUANT <- rlm(LogLOD1 ~ pred_logIE_MS2Quant , data= LoD_summary)
saveRDS(LoD_regressor_SS_MS2QUANT,  "LoD_regressor_SS_MS2QUANT.rds")
ggsave("IE-s versus LOD-S from cut-off approach.svg",
       width=12,
       height=12,
       units="cm")
cor.test(LoD_summary$pred_logIE_MS2Quant, LoD_summary$LogLOD1, method= "spearman")
# rho is -0.2537585 
# p-value is 0.02062
p2_IE_MS2Quant <- ggplot(data=LoD_summary) +
  geom_point(mapping = aes(x = pred_logIE_MS2Quant, y =LogLOD2, text=Name ), color= "#0070C0" , size = 4, alpha = 0.5) +
  xlab("IE-s")+
  ylab("LoD-s")+
  # annotation_logticks(color = basecolor)+
  # scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  # scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p2_IE_MS2Quant+my_theme
ggplotly(p2_IE_MS2Quant+my_theme)
ggsave("IE-s versus LOD-S from SN approach.svg",
       width=12,
       height=12,
       units="cm")
cor.test(LoD_summary$pred_logIE_MS2Quant, LoD_summary$LogLOD2, method= "spearman")
# rho is -0.1821617 
# p-value is 0.09931

LoD_summary <- LoD_summary %>% mutate(LogLOD3=log10(LOD3))
p3_IE_MS2Quant <- ggplot(data=LoD_summary) +
  geom_point(mapping = aes(x = pred_logIE_MS2Quant, y =LogLOD3, text=Name), color= "#0070C0" , size = 4, alpha = 0.5) +
  xlab("IE-s")+
  ylab("LoD-s")+
  # annotation_logticks(color = basecolor)+
  # scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  # scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p3_IE_MS2Quant+my_theme
p4_IE_MS2Quant <- ggplot(data=LoD_summary) +
  geom_point(mapping = aes(x = pred_logIE_MS2Quant, y =LogLOD4, text=Name), color= "#0070C0" , size = 4, alpha = 0.5) +
  xlab("IE-s")+
  ylab("LoD-s")+
  # annotation_logticks(color = basecolor)+
  # scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  # scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p4_IE_MS2Quant+my_theme
ggplotly(p4_IE_MS2Quant+my_theme)
ggsave("IE-s versus LOD-S from SD of triplicates approach.svg",
       width=12,
       height=12,
       units="cm")
cor.test(LoD_summary$pred_logIE_MS2Quant, LoD_summary$LogLOD4, method= "spearman")
# rho is -0.1610037 
# p-value is 0.1459
setwd("C:/Users/amso1873/OneDrive - Kruvelab/Amina/manuscript 3/LoD new runs/plots for manuscript")
plot_grid(p1_IE_MS2Quant+my_theme , p2_IE_MS2Quant+my_theme, p3_IE_MS2Quant+my_theme, p4_IE_MS2Quant+my_theme,
          labels =c("a) Cut-off approach", "b) Extrapolation of S/N","c) Residuals approach",
                    "d) Triplicate approach"),
          align = "h",
          hjust =  0.1, vjust = 1.5 ,
          label_x = 0.3,
          label_y = 0.999,
          nrow = 2,
          label_size = 12,
          label_fontfamily = font,
          label_colour = basecolor,
          label_fontface = "plain")
ggsave("IE_vs_LoD_SS_MS2quant.svg", height = 7, width=7)

LoD_summary_WW <- LoD_summary_WW %>% left_join(data)
p1_IE_MS2Quant_WW <- ggplot(data=LoD_summary_WW) +
  geom_point(mapping = aes(x = pred_logIE_MS2Quant, y =LogLOD1, text=Name),color= "#0070C0" , size = 4, alpha = 0.5) +
  xlab("IE-s")+
  ylab("LoD-s")+
  # annotation_logticks(color = basecolor)+
  # scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  # scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p1_IE_MS2Quant_WW+my_theme
ggplotly(p1_IE_MS2Quant_WW+my_theme)
ggsave("IE-s versus LOD-S from cut-off approach.svg",
       width=12,
       height=12,
       units="cm")
cor.test(LoD_summary_WW$pred_logIE_MS2Quant, LoD_summary_WW$LogLOD1, method= "spearman")
# rho is -0.2254668
# p-value is 0.1232

p2_IE_MS2Quant_WW <- ggplot(data=LoD_summary_WW) +
  geom_point(mapping = aes(x = pred_logIE_MS2Quant, y =LogLOD2, text=Name ), color= "#0070C0" , size = 4, alpha = 0.5) +
  xlab("IE-s")+
  ylab("LoD-s")+
  # annotation_logticks(color = basecolor)+
  # scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  # scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p2_IE_MS2Quant_WW+my_theme
ggplotly(p2_IE_MS2Quant_WW+my_theme)

ggsave("IE-s versus LOD-S from SN approach.svg",
       width=12,
       height=12,
       units="cm")
cor.test(LoD_summary_WW$pred_logIE_MS2Quant, LoD_summary_WW$LogLOD2, method= "spearman")
# rho is -0.192792
# p-value is 0.1887

p3_IE_MS2Quant_WW <- ggplot(data=LoD_summary_WW) +
  geom_point(mapping = aes(x = pred_logIE_MS2Quant, y =LogLOD3, text=Name), color= "#0070C0" , size = 4, alpha = 0.5) +
  xlab("IE-s")+
  ylab("LoD-s")+
  # annotation_logticks(color = basecolor)+
  # scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  # scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p3_IE_MS2Quant_WW+my_theme
ggplotly(p3_IE_MS2Quant_WW+my_theme)
p4_IE_MS2Quant_WW <- ggplot(data=LoD_summary_WW) +
  geom_point(mapping = aes(x = pred_logIE_MS2Quant, y =LogLOD4, text=Name), color= "#0070C0" , size = 4, alpha = 0.5) +
  xlab("IE-s")+
  ylab("LoD-s")+
  # annotation_logticks(color = basecolor)+
  # scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  # scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p4_IE_MS2Quant_WW+my_theme
ggplotly(p4_IE_MS2Quant_WW+my_theme)
ggsave("IE-s versus LOD-S from SD of triplicates approach.svg",
       width=12,
       height=12,
       units="cm")
cor.test(LoD_summary_WW$pred_logIE_MS2Quant, LoD_summary_WW$LogLOD4, method= "spearman")
# rho is -0.1980026 
# p-value is 0.1768
plot_grid(p1_IE_MS2Quant_WW+my_theme , p2_IE_MS2Quant_WW+my_theme, p3_IE_MS2Quant_WW+my_theme, p4_IE_MS2Quant_WW+my_theme,
          labels =c("a) Cut-off approach", "b) Extrapolation of S/N","c) Residuals approach",
                    "d) Triplicate approach"),
          align = "h",
          hjust =  0.1, vjust = 1.5 ,
          label_x = 0.3,
          label_y = 0.999,
          nrow = 2,
          label_size = 12,
          label_fontfamily = font,
          label_colour = basecolor,
          label_fontface = "plain")
ggsave("IE_vs_LoD_ww_MS2quant.svg", height = 7, width=7)
# ------- Matrix effect  -------
Matrix_effect  <- LoD_summary %>% select(Name, slope)
names(Matrix_effect)[names(Matrix_effect) == "slope"] <- "slope_ss"
Matrix_effect <- Matrix_effect %>% left_join(LoD_summary_WW %>% select(Name, slope))
p_Matrix_effect<- ggplot(data=Matrix_effect) +
  geom_point(mapping = aes(x = slope_ss, y =slope, text=Name), color= "#800000" , size = 4, alpha = 0.5) +
  geom_line(mapping = aes(x = slope_ss, y = slope_ss),color= "#800000" , size = 1, alpha = 0.5) +
  xlab("Response factor "["standard solutions"]~"(M"^-1*")")+
  ylab("Response factor "["spiked wastewater"]~"(M"^-1*")")+
  # annotation_logticks(color = basecolor)+
  scale_x_log10(limits=c(10^14,10^17 ), labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(limits=c(10^14,10^17 ), labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()
p_Matrix_effect+my_theme
ggsave("Matrix effect.svg", height = 4, width=4)
Matrix_effect <- Matrix_effect %>% mutate(ME=Matrix_effect$slope/Matrix_effect$slope_ss*100)

# ------------NORMAN List--------

dataset <-  read_delim('C:/Users/amso1873/OneDrive - Kruvelab/Amina/New folder/Rcode/logP_final_list3.csv',
                       delim = ",",
                       trim_ws = TRUE)
dataset <- dataset %>% mutate(pH=2.7)
dataset <- dataset %>% mutate(Column=as.factor(2))
dataset <- dataset %>% mutate(Buffer=as.factor(1))
dataset <- dataset %>% mutate(Organic_modifier=as.factor(1))
dataset <- dataset %>% replace(is.na(.), 0)
library("ranger")
regressor_MultiConditionRT <- read_rds("C:/Users/amso1873/OneDrive - Kruvelab/Amina/R code/regressorRF for RP 16082021.rds")
names(dataset)[names(dataset) == "logP"] <- "log_P_sum"
pred_RT=predict(regressor_MultiConditionRT, newdata = dataset)
dataset <- dataset %>% mutate(pred_RT)

names(dataset)[names(dataset) == "pred_RT"] <- "RT_Miklos"
Trans_pred_RT= predict(regressor1, newdata = dataset)
dataset <- dataset %>% mutate(Trans_pred_RT)

names(dataset)[names(dataset) == "pH"] <- "pH.aq."
dataset <- data.frame(dataset ,"NH4"=0)
names(dataset)[names(dataset) == "Organic_modifier"] <- "organic_modifier"
names(dataset)[names(dataset) == "Trans_pred_RT"] <- "ret_time"
names(dataset)[names(dataset) == "ASP.3"] <- "ASP-3"
names(dataset)[names(dataset) == "ASP.4"] <- "ASP-4"
names(dataset)[names(dataset) == "ASP.5"] <- "ASP-5"
names(dataset)[names(dataset) == "ASP.6"] <- "ASP-6"
names(dataset)[names(dataset) == "ASP.7"] <- "ASP-7"
names(dataset)[names(dataset) == "AVP.2"] <- "AVP-2"
names(dataset)[names(dataset) == "AVP.6"] <- "AVP-6"
names(dataset)[names(dataset) == "BCUTc.1h"] <- "BCUTc-1h"
names(dataset)[names(dataset) == "BCUTc.1l"] <- "BCUTc-1l"
names(dataset)[names(dataset) == "BCUTp.1h"] <- "BCUTp-1h"
names(dataset)[names(dataset) == "BCUTp.1l"] <- "BCUTp-1l"
names(dataset)[names(dataset) == "BCUTw.1l"] <- "BCUTw-1l"
names(dataset)[names(dataset) == "SCH.5"] <- "SCH-5"
names(dataset)[names(dataset) == "WTPT.4"] <- "WTPT-4"
names(dataset)[names(dataset) == "WTPT.5"] <- "WTPT-5"


dataset_3 <- dataset %>%
  mutate(
    organic = fn_organic_percentage(eluent_parameters,ret_time),
    viscosity =  fn_viscosity(organic,organic_modifier),
    surface_tension = fn_surface_tension(organic,organic_modifier),
    polarity_index = fn_polarity_index(organic,organic_modifier)) %>%
  #though all columns are selected with "everything()" in the end, the order of the columns is changed
  dplyr::select( SMILES, ret_time , organic_modifier, organic, pH.aq., Buffer,  viscosity,surface_tension,polarity_index, descs_pos)
prediction_set_model_pos <- dataset_3 %>%
  mutate(logIE_pred = 0)

prediction <-  predict(regressor_pos, newdata = prediction_set_model_pos)
prediction <- prediction$aggregate
prediction_set_model_pos$logIE_pred <- prediction
prediction_set_model_pos <- prediction_set_model_pos %>%
  mutate(logIE_pred = prediction) %>%
  dplyr::select(SMILES,logIE_pred, everything())
dataset <- dataset %>% left_join(prediction_set_model_pos)
dataset_for_pca <- dataset %>% dplyr::select(-c(SMILES, ret_time, logIE_pred, polarity_index, viscosity, surface_tension, NH4, RT_Miklos, organic_modifier, Buffer, Column, MW, SPLIT, pH.aq., organic))
# Run lines 1606 until 1665
dataset_for_pca <- dataset_for_pca %>% dplyr::select(intersect(colnames(Saer_dataset), colnames(dataset_for_pca)))
dataset_for_pca <- dataset_for_pca %>%
  dplyr::select(-nearZeroVar(dataset_for_pca))

dataset_for_pca <- dataset_for_pca %>% dplyr::select(-c(bruto_formula))


#Calculating  the correlation between features
correlationMatrix <- cor(dataset_for_pca, use = "complete.obs")


# Finding the  highly corrected features (ideally >0.8)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.8)

dataset_for_pca <- dataset_for_pca %>%
  dplyr::select(-highlyCorrelated)

pca = prcomp(dataset_for_pca, center=TRUE, scale=TRUE)
dataset_used_for_pca <- dataset_for_pca
# write_delim(dataset_used_for_pca,
#             "dataset_used_for_pca.csv", delim=",")
predictions = predict(pca, dataset_for_pca)
dataset_pca_output <- as_tibble(pca$x)
loadings = as_tibble(pca$rotation)
dataset_for_pca <- dataset_for_pca %>% mutate(PC1=dataset_pca_output$PC1)
dataset_for_pca <- dataset_for_pca %>% mutate(PC2=dataset_pca_output$PC2)
names(dataset)[names(dataset) == "logIE_pred"] <- "logIE_pred_RT"
pred_LoD <- predict(LoD_regressor_SS_Pred_RT, newdata = dataset)
dataset <- dataset %>% mutate(pred_LoD)
dataset_for_pca <- dataset_for_pca %>% left_join(dataset)
dataset_for_pca <- dataset_for_pca %>% mutate(LoD=10^pred_LoD)
pca1 <-ggplot(data=) +
  geom_point(mapping = aes(x = dataset_for_pca$PC1, y = dataset_for_pca$PC2, color=dataset_for_pca$pred_LoD), size = 2, alpha = 1) +
  # geom_segment(data = loadings,
  #              mapping = aes(x = 0,
  #                            xend=PC1*150, # must multiply, since otherwise it just sits at 0
  #                            y = 0,
  #                            yend = PC2*150),     # must multiply, since otherwise it just sits at 0
  #              arrow = arrow(),
  #              color = "black") +
  # geom_text(mapping = aes( x = PC1 *175,    #since previously changed length of arrows, we need to do the same for the labels
  #                          y = PC2*175),    #since previously changed length of arrows, we need to do the same for the labels
  #           label = dataset_used_for_pca %>%
  #              colnames(),
  #           data = loadings,
  #            size = 3)+
  xlim(-38,20)+
  ylim(-25,30)+
  xlab("PC1 (8.7%)")+
  ylab("PC2 (6.7%)")+
  labs(color= "Predicted LoD")+
  theme_bw()
pca1+ my_theme
pca2 <-ggplot(data=) +
  geom_point(mapping = aes(x = dataset_for_pca$PC1, y = dataset_for_pca$PC2, color=dataset_for_pca$pred_LoD), size = 2, alpha = 1) +
  # geom_segment(data = loadings,
  #            mapping = aes(x = 0,
  #                          xend=PC1*150, # must multiply, since otherwise it just sits at 0
  #                          y = 0,
  #                          yend = PC2*150),     # must multiply, since otherwise it just sits at 0
  #            arrow = arrow(),
  #            color = "black") +
  # geom_text(mapping = aes( x = PC1 *175,    #since previously changed length of arrows, we need to do the same for the labels
  #                        y = PC2*175),    #since previously changed length of arrows, we need to do the same for the labels
  #         label = dataset_used_for_pca %>%
  #            colnames(),
  #         data = loadings,
  #          size = 3)+
  xlim(-38,20)+
  ylim(-25,30)+
  xlab("PC1 (8.7%)")+
  ylab("PC2 (6.7%)")+
  labs(color= "Predicted LoD")+
  theme_bw()
pca2+ my_theme
# plot_grid( pca1+ my_theme, pca2+my_theme,
#            labels =c("a) ", "b) "),
#            align = "h",
#            hjust =  0.1, vjust = 1.5 ,
#            label_x = 0.3,
#            label_y = 0.999,
#            nrow = 1,
#            label_size = 12,
#            label_fontfamily = font,
#            label_colour = basecolor,
#            label_fontface = "plain")
# ggsave("pca plot with predicted LoD with loadings and  scores separately.svg", height = 10, width= 20)
# setwd("C:/Users/amso1873/OneDrive - Kruvelab/Amina/manuscript 3/LoD new runs")
# dataset_for_pca <- dataset_for_pca %>% dplyr::select(-ret_time)
# write_delim(dataset_for_pca,
#             "pca_results_with_pred_LOD1.csv", delim=",")
dataset_for_pca_2 <-  read_delim('C:/Users/amso1873/OneDrive - Kruvelab/Amina/manuscript 3/LoD new runs/pca_results_with_pred_LOD1.csv',
                                 delim = ",",
                                 trim_ws = TRUE)


amenability <-  read_delim('C:/Users/amso1873/OneDrive - Kruvelab/Amina/manuscript 3/LoD new runs/susdat_2022-01-18-104316.csv',
                           delim = ",",
                           trim_ws = TRUE)
amenability <- amenability %>% dplyr::select(SMILES, `Prob_+ESI`)
# amenability <- amenability %>% filter(`Prob_+ESI` == 0 | `Prob_+ESI` ==1)
amenability <- amenability %>% mutate(Prob = case_when(`Prob_+ESI`<0.2 ~ "<0.2", 
                                                       `Prob_+ESI`<0.4 ~ "<0.4",
                                                       `Prob_+ESI`<0.6 ~ "<0.6", 
                                                       `Prob_+ESI`<0.8 ~ "<0.8", 
                                                       `Prob_+ESI`=0.8 & `Prob_+ESI`> 0.8 ~ ">0.8" ))
dataset_for_pca_2 <- dataset_for_pca_2 %>% left_join(amenability)
dataset_for_pca_2 <- dataset_for_pca_2 %>% drop_na()
dataset_for_pca_2$Prob <- as.factor(dataset_for_pca_2$Prob)
dataset_for_pca_2 <- dataset_for_pca_2 %>% mutate(LoD= 10^pred_LoD)
# dataset_for_pca_3 <- dataset_for_pca_2 %>% filter(Prob=="<0.2" |Prob==">0.8")
library("PupillometryR")
p_pca_2 <- ggplot(dataset_for_pca_2) +
  aes(x = Prob,
      y = LoD,
      fill = Prob) +
  geom_flat_violin(alpha = 0.5) +
  geom_boxplot(width = .25)+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_fill_manual(values = c("#7f7f7f","#CCCCCC", "#FFFFFF", "#99CCFF","#0070C0"))+
  ylab("Predicted LoD (M)")+
  xlab("Probability of amenability")+
  labs(fill = "Probability")+
  guides(fill = FALSE)  
p_pca_2+my_theme
wilcox.test(LoD~Prob, alternative="greater", data= dataset_for_pca_3)
p_pca_3 <- ggplot(dataset_for_pca_2) +
  aes(x = Prob,
      y = LoD,
      fill = Prob) +
  geom_flat_violin(alpha = 0.5) +
  geom_boxplot(width = .25)+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  ylab("Predicted LoD (M)")+
  xlab("Probability of amenability")+
  coord_flip()+ 
  labs(fill = "Probability")
p_pca_3+my_theme
ggsave("SI figure_LoD predictions of Norman list and compared to Probability of amenability in ESI.svg", height = 10, width=14)

# ------- unknown features -------

Unknwon_features <-  read_delim('C:/Users/amso1873/OneDrive - Kruvelab/Amina/manuscript 3/LoD new runs/IE-s for unknown features.csv',
                                       delim = ",",
                                       col_names = TRUE,
                                       trim_ws = TRUE)
Unknwon_features_FP <-  read_delim('C:/Users/amso1873/OneDrive - Kruvelab/Amina/manuscript 3/LoD new runs/Fingerprints for unknown features.csv',
                                delim = ",",
                                col_names = TRUE,
                                trim_ws = TRUE)
names(Unknwon_features)[names(Unknwon_features) == "pred_logIE"] <- "pred_logIE_MS2Quant"
LoD_regressor_SS_MS2QUANT <- read_rds("C:/Users/amso1873/OneDrive - Kruvelab/Amina/manuscript 3/LoD new runs/plots for manuscript/LoD_regressor_SS_MS2QUANT.rds")

pred_LoD <- predict(LoD_regressor_SS_MS2QUANT, newdata = Unknwon_features)
Unknwon_features <- Unknwon_features %>% mutate(pred_LoD)
Unknwon_features_filtered <- Unknwon_features %>% filter(pred_LoD > -7.54 )
Unknwon_features_filtered <- Unknwon_features_filtered %>% bind_rows(Unknwon_features %>% filter(pred_LoD < -10.007))
write_delim(Unknwon_features_filtered,
            "Unknwon_features_filtered_with_20_lowest_and_highest_LoD_LAtestversion.csv", delim=",")
# ----------Saer's dataset------------
Saer_dataset <-  read_delim('C:/Users/amso1873/OneDrive - Kruvelab/Amina/manuscript 3/LoD new runs/Descriptors fro saer dataset.csv',
                            delim = ",",
                            col_names = TRUE,
                            trim_ws = TRUE)
log_P_sum <-c()
for (i in 1:length(Saer_dataset$SMILES)){
  print(i)
  log_P_sum <- c(log_P_sum, fn_logP(Saer_dataset$SMILES[i]))
}
Saer_dataset <- Saer_dataset %>% ungroup() %>% mutate(log_P_sum) 
Saer_dataset <- Saer_dataset %>% mutate(pH=2.7) %>% 
  mutate(Column=as.factor(2)) %>% 
  mutate(Buffer=as.factor(1)) %>%
  mutate(Organic_modifier=as.factor(1))

Saer_dataset <- Saer_dataset %>% mutate(RT_Miklos=predict(regressor_MultiConditionRT, newdata = Saer_dataset))
Saer_dataset <- Saer_dataset %>% mutate(ret_time=predict(regressor1, newdata = Saer_dataset))
names(Saer_dataset)[names(Saer_dataset) == "pH"] <- "pH.aq."
Saer_dataset <- data.frame(Saer_dataset ,"NH4"=0)
names(Saer_dataset)[names(Saer_dataset) == "Organic_modifier"] <- "organic_modifier"
names(Saer_dataset)[names(Saer_dataset) == "ASP.3"] <- "ASP-3"
names(Saer_dataset)[names(Saer_dataset) == "ASP.4"] <- "ASP-4"
names(Saer_dataset)[names(Saer_dataset) == "ASP.5"] <- "ASP-5"
names(Saer_dataset)[names(Saer_dataset) == "ASP.6"] <- "ASP-6"
names(Saer_dataset)[names(Saer_dataset) == "ASP.7"] <- "ASP-7"
names(Saer_dataset)[names(Saer_dataset) == "AVP.2"] <- "AVP-2"
names(Saer_dataset)[names(Saer_dataset) == "AVP.6"] <- "AVP-6"
names(Saer_dataset)[names(Saer_dataset) == "BCUTc.1h"] <- "BCUTc-1h"
names(Saer_dataset)[names(Saer_dataset) == "BCUTc.1l"] <- "BCUTc-1l"
names(Saer_dataset)[names(Saer_dataset) == "BCUTp.1h"] <- "BCUTp-1h"
names(Saer_dataset)[names(Saer_dataset) == "BCUTp.1l"] <- "BCUTp-1l"
names(Saer_dataset)[names(Saer_dataset) == "BCUTw.1l"] <- "BCUTw-1l"
names(Saer_dataset)[names(Saer_dataset) == "SCH.5"] <- "SCH-5"
names(Saer_dataset)[names(Saer_dataset) == "WTPT.4"] <- "WTPT-4"
names(Saer_dataset)[names(Saer_dataset) == "WTPT.5"] <- "WTPT-5"
Saer_dataset_2 <- Saer_dataset %>% 
  mutate(
    organic = fn_organic_percentage(eluent_parameters,ret_time),
    viscosity =  fn_viscosity(organic,organic_modifier),
    surface_tension = fn_surface_tension(organic,organic_modifier),
    polarity_index = fn_polarity_index(organic,organic_modifier)) %>%
  #though all columns are selected with "everything()" in the end, the order of the columns is changed
  dplyr::select( SMILES, ret_time , organic_modifier, organic, pH.aq., Buffer,  viscosity,surface_tension,polarity_index, descs_pos)


prediction_set_model_pos <- Saer_dataset_2 %>%
  mutate(logIE_pred = 0)

prediction <-  predict(regressor_pos, newdata = prediction_set_model_pos)
prediction <- prediction$aggregate
prediction_set_model_pos$logIE_pred <- prediction
prediction_set_model_pos <- prediction_set_model_pos %>%
  mutate(logIE_pred = prediction) %>%
  dplyr::select(SMILES,logIE_pred, everything())
Saer_dataset_2 <- Saer_dataset_2 %>% left_join(prediction_set_model_pos)
names(Saer_dataset_2)[names(Saer_dataset_2) == "logIE_pred"] <- "logIE_pred_RT"
Saer_dataset_2 <- Saer_dataset_2 %>% mutate(Pred_LoD=predict(LoD_regressor_SS_Pred_RT, newdata = Saer_dataset_2))
Saer_dataset <- Saer_dataset%>% left_join(prediction_set_model_pos)
names(Saer_dataset)[names(Saer_dataset) == "logIE_pred_RT"] <- "logIE_pred"
project_Saer_dataset <- predict(pca, Saer_dataset)
project_Saer_dataset <- as_tibble(project_Saer_dataset)
Saer_dataset  <- Saer_dataset %>% 
  mutate(PC1=project_Saer_dataset$PC1) %>%
  mutate(PC2=project_Saer_dataset$PC2) %>%
  mutate(Pred_LoD =Saer_dataset_2$Pred_LoD)
project_Saer_dataset  <- project_Saer_dataset %>% unique()

p_saer <-ggplot() +
  geom_point(mapping = aes(x = dataset_for_pca$PC1, y = dataset_for_pca$PC2 ), color="#7f7f7f",  size = 2, alpha = 0.5) +
  geom_point(mapping = aes(x = Saer_dataset$PC1, y = Saer_dataset$PC2 ), color="#0070C0", size = 2, alpha = 0.5) +
  xlim(-38,20)+
  ylim(-25,30)+
  xlab("PC1 (8.7%)")+
  ylab("PC2 (6.7%)")+
  theme_bw()
p_saer+ my_theme
p_saer_LoD <-ggplot(data=Saer_dataset) +
  geom_point(mapping = aes(x = PC1, y = PC2, color=Pred_LoD ),  size = 2, alpha = 1) +
  xlim(-20,20)+
  ylim(-20,20)+
  xlab("PC1 (8.7%)")+
  ylab("PC2 (6.7%)")+
  labs(color= "Predicted LoD")+
  theme_bw()
p_saer_LoD+ my_theme


CID_SMILES <-  read_delim('C:/Users/amso1873/OneDrive - Kruvelab/Amina/manuscript 3/LoD new runs/All 2657 compounds from the papers.csv',
                            delim = ",",
                            col_names = TRUE,
                            trim_ws = TRUE)
CID_SMILES <- CID_SMILES  %>%  dplyr::select(c(SMILES, CID))
CID_Polarity <-  read_delim('C:/Users/amso1873/OneDrive - Kruvelab/Amina/manuscript 3/LoD new runs/All papers.csv',
                          delim = ",",
                          col_names = TRUE,
                          trim_ws = TRUE) %>%  dplyr::select(c(Polarity, CID))
CID_SMILES <- CID_SMILES %>% left_join(CID_Polarity) %>% unique()
Saer_dataset <- Saer_dataset %>% left_join(CID_SMILES) %>% unique()
polarity_p <- ggplot(Saer_dataset %>% filter(Polarity=="Only positive" |Polarity=="Only negative" )) +
  aes(x = Polarity,
      y =  10^Pred_LoD,
      fill = Polarity) +
  geom_flat_violin(alpha = 0.5) +
  geom_boxplot( width = .25)+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_fill_manual(values = c("#7f7f7f", "#0070C0"))+
  ylab("Predicted LoD (M)")+
  xlab("Polarity")+
  labs(fill = "Probability")+
  guides(fill = FALSE)  
polarity_p+my_theme
ggsave("Saer_dataset_with_polarity_LoD.svg", height = 10, width=10)
Norman_list_with_LoD <- dataset_for_pca %>% dplyr::select(c(SMILES, LoD))
Saer_list_with_LoD <- Saer_dataset %>% dplyr::select(c(SMILES, Pred_LoD)) %>% mutate(LoD=10^Pred_LoD) %>% unique()
plot_grid( p1_SS_application+my_theme, pca2+ my_theme, p_pca_2+my_theme, p_saer+ my_theme, p_saer_LoD+ my_theme , polarity_p+my_theme, 
           labels =c("a) ", "b) ",
                     "c) ", "d)", "e)", "f)"),
           align = "h",
           hjust =  0.1, vjust = 1.5 ,
           label_x = 0.3,
           label_y = 0.999,
           nrow = 2,
           label_size = 12,
           label_fontfamily = font,
           label_colour = basecolor,
           label_fontface = "plain")
ggsave("Final figure for manuscript.svg", height = 10, width=15)
LOD_SS <- LoD_summary %>% dplyr::select(c(Name, SMILES, Alignment.ID, slope, intercept, ret_time_exp, R2, LOD1, LOD2, LOD3, LOD4))
# write_delim(LOD_SS,
#             "LOD_SS.csv", delim=",")
LOD_ww <- LoD_summary_WW %>% dplyr::select(c(Name, SMILES, Alignment.ID, slope, intercept, ret_time_exp, R2, LOD1, LOD2, LOD3, LOD4))
# write_delim( LOD_ww,
#             " LOD_ww.csv", delim=",")
NORMAN_WITH_pred_LOD_export <- dataset_for_pca %>% dplyr::select(c(SMILES, pred_LoD))
write_delim( NORMAN_WITH_pred_LOD_export,
            "NORMAN_WITH_pred_LOD_export.csv", delim=",")
