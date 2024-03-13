########Reference to Qiaozhi's codes######

require(ggplot2)
require(grid)
require(gridExtra)
library(cowplot)
library(scales)
require(reshape)
require(egg)
library(data.table)
library(plyr)


# Load in data ------------------------------------------------------------
is <- read.csv('/mpp/2M_CSV/MPP_IS_20180529.csv',header = TRUE, fileEncoding='GB18030')
med_raw <- read.csv('/mpp/2M_CSV/MPP_ISMed_20180529.csv',header = TRUE, fileEncoding='latin1')
med <- med_raw[med_raw$mark=="InitialS",c("PID", "drugid")]
med$drugid <- as.character(med$drugid)
med$drugid[med$drugid==""] <- NA
med <- na.omit(med)
medsplit <- do.call("rbind",strsplit(med$drugid,"_"))
med3 <- data.frame(PID=med$PID,chem=as.numeric(medsplit[,1]))

drugbp <- read.csv("/mpp/codes/codes_xinyue/HTN_drug.csv",header=T,as.is=T, fileEncoding='GB18030')
drugbp <- drugbp[,c(1,2,4,5)]
names(drugbp) <- c("chem_id","chem_name","class_id","class_name")
drugbp$pid <- paste(drugbp$class_id,drugbp$chem_id,sep="_")
med4 <- med3[med3$chem %in% drugbp$chem_id,] ##from 399665 to 248949, select BP drugs
medWM <- med4[!duplicated(med4),] 
is$tx <- rep(FALSE,nrow(is))
is$tx[is$PID %in% medWM$PID] <- TRUE

baseline_data<- is

# Preprocess --------------------------------------------------------------
data <- subset(baseline_data, select = c(PID, IS_gender, is_AGE, is_income, IS_education, IS_occupation, IS_hukou, IS_Prov, IS_marriage, IS_SBP, IS_DBP,IS_smk,  IS_BMI,tx, IS_Ethnicity,IS_TC,IS_HDL,IS_TG,IS_LDL,IS_alc_freq
))
names(data) <- c("Patient_Id", "Sex.x", "Age", "Household_Income", "Education", "Occupation", "Hukou","Province", "Marital_Status","Systolic_Average", "Diastolic_Average",  
"Smoking_Currently",  "BMI","tx", "Ethnicity","TC","HDL","TG","LDL","Alcohol"
)

# Preprocess --------------------------------------------------------------
eth_tbl <- table(data$Ethnicity)
large_ethnicities <- names(eth_tbl[order(-eth_tbl)])[1:12]
data_clean <- na.omit(subset(data,Ethnicity %in% large_ethnicities  ))



data_clean$Hukou <- factor(data_clean$Hukou)
data_clean$Province <- factor(data_clean$Province)
data_clean$Sex.x <- factor(data_clean$Sex.x)
data_clean$Education <- factor(data_clean$Education)
data_clean$Household_Income <- factor(data_clean$Household_Income)
data_clean$Occupation <- factor(data_clean$Occupation)
data_clean$Marital_Status <- factor(data_clean$Marital_Status == 1)
data_clean$Smoking_Currently <- factor(data_clean$Smoking_Currently)

data_clean$tx <- factor(data_clean$tx)
data_clean$Alcohol <- factor(data_clean$Alcohol)
levels(data_clean$Sex.x) <- c("Male", "Female")
levels(data_clean$Occupation) <- c("Farmer", "Workers", "Administrators", "Admin. Clerk", "Technician", "Business", "Bus. Owner", "Military", "Others", "Retire", "Unemployed", "Housework", "Unknown", "Refuse")
levels(data_clean$Marital_Status) <- c("Unmarried", "Married")
levels(data_clean$Hukou) <- c("Rural", "Urban", "Unified","No Hukou")

levels(data_clean$Smoking_Currently) <- c("Smoking", "Not Smoking")
data_clean$Education[data_clean$Education==13] <- 12
data_clean$Education <- droplevels(data_clean$Education)
levels(data_clean$Education) <- c("Illiterate","<Primary", "Homeschool", "Elementary", 
                                  "Middle",  "High", "Vocational", "Associate",  "Bachelor",
                                  "Masters", "Doctoral", "No answer")
levels(data_clean$tx) <- c("No meds", "On meds")
levels(data_clean$Alcohol)<-c("Never","Less than once a month","2-4 times a month","2-3 times a week","More than 4 times a month","Unclear","Refuse")



eth_dict <- list('1'='汉族','2'='蒙古族','3'='回族','4'='藏族','5'='维吾尔族','6'='苗族','7'='彝族','8'='壮族','9'='布依族','10'='朝鲜族','11'='满族','12'='侗族','13'='瑶族','14'='白族','15'='土家族','16'='哈尼族','17'='哈萨克族','18'='傣族','19'='黎族','20'='傈僳族','21'='佤族','22'='畲族','23'='高山族','24'='拉祜族','25'='水族','26'='东乡族','27'='纳西族','28'='景颇族','29'='柯尔克孜族','30'='土族','31'='达斡尔族','32'='仫佬族','33'='羌族','34'='布朗族','35'='撒拉族','36'='毛南族','37'='仡佬族','38'='锡伯族','39'='阿昌族','40'='普米族','41'='塔吉克族','42'='怒族','43'='乌孜别克族','44'='俄罗斯族','45'='鄂温克族','46'='德昂族','47'='保安族','48'='裕固族','49'='京族','50'='塔塔尔族','51'='独龙族','52'='鄂伦春族','53'='赫哲族','54'='门巴族','55'='珞巴族','56'='基诺族')
ethnames <- unlist(eth_dict[levels(data_clean$Ethnicity)])

data_clean$Ethnicity  <- factor(data_clean$Ethnicity)
levels(data_clean$Ethnicity) <- c("Han", "Mongol", "Hui","Tibet", "Uygher", "Miao", "Yi","Zhuang", "Korean", "Man","Dong", "Tujia")

provdict <- list('11'='Beijing','12'='Tianjin','13'='Hebei','14'='Shanxi','15'='InnerMongol','21'='Liaoning','22'='Jilin','23'='Heilongjiang','31'='Shanghai','32'='Jiangsu','33'='Zhejiang','34'='Anhui','35'='Fujian','36'='Jiangxi','37'='Shandong','41'='Henan','42'='Hubei','43'='Hunan','44'='Guangdon','45'='Guangxi','46'='Hainan','50'='Chongqin','51'='Sichuan','52'='Guizhou','53'='Yunnan','54'='Tibet','61'='Shaanxi','62'='Gansu','63'='Qinghai','64'='Ningxia','65'='Xinjiang','71'='Taiwan','91'='HongKong','92'='Macao')
levels(data_clean$Province) =  unlist(provdict[levels(data_clean$Province)])

data_clean$Ethnicity  <- factor(data_clean$Ethnicity)
data_clean$Household_Income <- as.factor(data_clean$Household_Income)
incomedict <- list('1'='<10,000','2'='10,000- 25,000','3'='25,001- 50,000','4'='50,001- 100,000','5'='100,001- 200,000','6'='200,001- 300,000','7'='300,001- 600,000','8'='600,001 -1,500,000','9'='>1,500,000','10'='Unclear','11'='Refuse to answer','12'='少于5,000','13'='5,000-9,999','14'='10,000-19,999','15'='20,000-50,000','16'='>50,000')
levels(data_clean$Household_Income)<- unlist(incomedict[levels(data_clean$Household_Income)])
levels(data_clean$Household_Income) <- c(levels(data_clean$Household_Income), "<10k", "10k-50k", ">50k" )
data_clean$Household_Income[data_clean$Household_Income %in% c("<10,000", "少于5,000", '5,000-9,999')] <- "<10k"
data_clean$Household_Income[data_clean$Household_Income  %in% c("10,000- 25,000", "25,001- 50,000", "20,000-50,000", "10,000-19,999")] <- "10k-50k"
data_clean$Household_Income[data_clean$Household_Income  %in% c("50,001- 100,000","100,001- 200,000", "200,001- 300,000", "600,001 -1,500,000", ">1,500,000",">50,000" ,"300,001- 600,000" )] <- ">50k"
data_clean$Household_Income[data_clean$Household_Income == "Refuse to answer" ] <- "Unclear"
data_clean$Household_Income <- droplevels(data_clean$Household_Income)
#data_clean$TC <- data_clean$TC*18
#data_clean$TG <- data_clean$TG*18
data_clean$LDL <- data_clean$LDL*38.65
#data_clean$HDL <- data_clean$HDL*18
data_clean <- data_clean[which(data_clean$Age <=75 & data_clean$Age >=35),]
data_clean<- data_clean[which(data_clean$LDL<200 & data_clean$TC<200 & data_clean$HDL<200 ),]
data_clean$pp <- data_clean$Systolic_Average-data_clean$Diastolic_Average

####Figure 1
cor(as.numeric(data_clean$BMI),data_clean$Systolic_Average,method = "spearman")
cor(as.numeric(data_clean$Hx_Stroke),data_clean$Systolic_Average,method = "spearman")
cor(as.numeric(data_clean$Marital_Status),data_clean$Systolic_Average,method = "spearman")
cor(as.numeric(data_clean$Age),data_clean$Systolic_Average,method = "spearman")



cor(as.numeric(data_clean$BMI),data_clean$Diastolic_Average,method = "spearman")
cor(as.numeric(data_clean$Age),data_clean$Diastolic_Average,method = "spearman")
library(ggplot2)
title_theme =  theme(plot.title=element_text(face='plain'))

ggplot(data_clean, aes(Age, Diastolic_Average) ) +  stat_binhex() +geom_smooth()+theme(legend.position="none") + scale_fill_gradient(low="#ffe5e5", high = "#ff0000") + scale_y_continuous(name="Diastolic BP (mmHg)")+ xlab("Age") + theme_set(theme_gray(base_size = 22))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background =element_blank(), axis.line = element_line(colour = "black"),legend.position="none")


ggplot(data_clean, aes(Age, Systolic_Average) ) +  stat_binhex() +geom_smooth() +theme(legend.position="none") + scale_fill_gradient(low="#ffe5e5", high = "#ff0000") + scale_y_continuous(name="Systolic BP (mmHg)")+ xlab("Age") +theme_set(theme_gray(base_size = 22))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background =element_blank(), axis.line = element_line(colour = "black"),legend.position="none")

ggplot(subset(data_clean,pp>0), aes(Age, pp) ) +  stat_binhex() +geom_smooth() +theme(legend.position="none") + scale_fill_gradient(low="#ffe5e5", high = "#ff0000") + scale_y_continuous(name="Pulse Pressure (mmHg)")+ xlab("Age") + theme_set(theme_gray(base_size = 22))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background =element_blank(), axis.line = element_line(colour = "black"),legend.position="none")



ggplot(subset(data_clean, LDL>=40), aes(Age, LDL) ) +  stat_binhex() +geom_smooth() + ggtitle('c')+theme(legend.position="none") + scale_fill_gradient(low="#ffe5e5", high = "#ff0000") + scale_y_continuous(name="LDL (mg/dL)")+ xlab("Age") + title_theme

ggplot(subset(data_clean,BMI<50 &BMI>10), aes(Age, BMI) ) +  stat_binhex() +geom_smooth() + ggtitle('d')+theme(legend.position="none") + scale_fill_gradient(low="#ffe5e5", high = "#ff0000") + scale_y_continuous(name="BMI")+ xlab("Age") + title_theme
#####
require(egg)

data_clean_sample_idx <- 1:nrow(data_clean)
data_sub <- subset(data_clean[data_clean_sample_idx,],pp>0)

plot_titles <- c(`Sex.x`="Sex",  `Hukou`="Hukou", `Occupation`="Occupation", 
                 `Education`="Education", `Household_Income`="Household Income",  `Marital_Status`="Marital Status", `Province`="Province" ,`tx`="BP Meds" , `Ethnicity`="Ethnicity", `Smoking_Currently`="Currently Smoking") 


covariates <- c("Sex.x", "Household_Income", 
"Occupation", "Ethnicity", "Marital_Status",  "Hukou","Province","Education","Smoking_Currently", "tx")

df <- melt(data_sub[,c("Systolic_Average", "Age", covariates )], id =c("Systolic_Average", "Age"))


#950 x 1220
#SBP VS. AGE

theme_set(theme_cowplot(font_size=9))
g <-plot_grid_smooths(df,data_clean,plot_titles, 2, "AGE", "SBP", c(120, 160))
ggsave(g,filename = "figure2.pdf",width = 0.8*12,height=0.8*17)



#DBP VS. AGE
theme_set(theme_cowplot(font_size=10.5))
df <- melt(data_sub[,c("Diastolic_Average", "Age", covariates )], id =c("Diastolic_Average", "Age"))
g<- plot_grid_smooths_dbp(df,plot_titles, 2, "AGE", "DBP", c(60, 120))
ggsave(g,filename = "eFigure2.pdf",width = 0.8*12,height=0.8*17)


#BMI VS.AGE
plot_titles <- c(`Sex.x`="Sex",  `Hukou`="Hukou", `Occupation`="Occupation", 
                 `Education`="Education", `Household_Income`="Household Income",  `Marital_Status`="Marital Status", `Province`="Province"  , `Ethnicity`="Ethnicity", `Smoking_Currently`="Currently Smoking") 


covariates <- c("Sex.x", "Household_Income", 
                "Occupation", "Ethnicity", "Marital_Status",  "Hukou","Province","Education","Smoking_Currently")

theme_set(theme_cowplot(font_size=10.5))
df <- melt(data_sub[,c("BMI", "Age", covariates )], id =c("BMI", "Age"))
g<- plot_grid_smooths_bmi(df,plot_titles, 2, "AGE", "BMI", c(10, 50))

#PP VS. AGE
plot_titles <- c(`Sex.x`="Sex",  `Hukou`="Hukou", `Occupation`="Occupation", 
                 `Education`="Education", `Household_Income`="Household Income",  `Marital_Status`="Marital Status", `Province`="Province"  , `Ethnicity`="Ethnicity", `Smoking_Currently`="Currently Smoking") 


covariates <- c("Sex.x", "Household_Income", 
                "Occupation", "Ethnicity", "Marital_Status",  "Hukou","Province","Education","Smoking_Currently")

theme_set(theme_cowplot(font_size=10.5))
df <- melt(data_sub[,c("pp", "Age", covariates )], id =c("pp", "Age"))
g<- plot_grid_smooths_pp(df,plot_titles, 2, "AGE", "PP", c(35, 75))




#########LDL VS. AGE####

age_cat_1 <-  data_clean[which(data_clean$Age >=35 & data_clean$Age <45),]
age_cat_2 <-  data_clean[which(data_clean$Age >=45 & data_clean$Age <55),]
age_cat_3 <-  data_clean[which(data_clean$Age >=55 & data_clean$Age <65),]
age_cat_4 <-  data_clean[which(data_clean$Age >=65 & data_clean$Age <=75),]

p1<-round(100*sum(age_cat_1$LDL>=160)/dim(age_cat_1)[1],2)
p2<-round(100*sum(age_cat_2$LDL>=160)/dim(age_cat_2)[1],2)
p3<-round(100*sum(age_cat_3$LDL>=160)/dim(age_cat_3)[1],2)
p4<-round(100*sum(age_cat_4$LDL>=160)/dim(age_cat_4)[1],2)


#p5<-round(100*sum(age_cat_1$HDL<30)/dim(age_cat_1)[1],1)
#p6<-round(100*sum(age_cat_2$HDL<30)/dim(age_cat_2)[1],1)
#p7<-round(100*sum(age_cat_3$HDL<30)/dim(age_cat_3)[1],1)
#p8<-round(100*sum(age_cat_4$HDL<30)/dim(age_cat_4)[1],1)

Percentage <- c(p1,p2,p3,p4)
#group<- c(rep("LDL",4),rep("HDL",4))
age_cat<- c("35-44","45-54","55-64","65-75")


dat <- cbind(age_cat,Percentage)
dat<-as.data.frame(dat)
#dat1 = dat %>% group_by(Age) %>% mutate(position = rank(group))

ggplot(data=dat, aes(x=age_cat, y=Percentage))+geom_bar(stat="identity",fill="steelblue")+labs(x="Age category")+ggtitle("Proportion of Dyslipidaemia\n in each age groups")



##############SBP/DBP VS. AGE

p1<-round(100*sum(age_cat_1$Systolic_Average>=160)/dim(age_cat_1)[1],1)
p2<-round(100*sum(age_cat_2$Systolic_Average>=160)/dim(age_cat_2)[1],1)
p3<-round(100*sum(age_cat_3$Systolic_Average>=160)/dim(age_cat_3)[1],1)
p4<-round(100*sum(age_cat_4$Systolic_Average>=160)/dim(age_cat_4)[1],1)


p5<-round(100*sum(age_cat_1$Diastolic_Average>=100)/dim(age_cat_1)[1],1)
p6<-round(100*sum(age_cat_2$Diastolic_Average>=100)/dim(age_cat_2)[1],1)
p7<-round(100*sum(age_cat_3$Diastolic_Average>=100)/dim(age_cat_3)[1],1)
p8<-round(100*sum(age_cat_4$Diastolic_Average>=100)/dim(age_cat_4)[1],1)

p1<-round(100*sum(age_cat_1$Systolic_Average>=140 | age_cat_1$Diastolic_Average>=90)/dim(age_cat_1)[1],1)
p2<-round(100*sum(age_cat_2$Systolic_Average>=140 | age_cat_2$Diastolic_Average>=90)/dim(age_cat_2)[1],1)
p3<-round(100*sum(age_cat_3$Systolic_Average>=140 | age_cat_3$Diastolic_Average>=90)/dim(age_cat_3)[1],1)
p4<-round(100*sum(age_cat_4$Systolic_Average>=140 | age_cat_4$Diastolic_Average>=90)/dim(age_cat_4)[1],1)


p5<-round(100*sum(age_cat_1$Systolic_Average>=160 | age_cat_1$Diastolic_Average>=100)/dim(age_cat_1)[1],1)
p6<-round(100*sum(age_cat_2$Systolic_Average>=160 | age_cat_2$Diastolic_Average>=100)/dim(age_cat_2)[1],1)
p7<-round(100*sum(age_cat_3$Systolic_Average>=160 | age_cat_3$Diastolic_Average>=100)/dim(age_cat_3)[1],1)
p8<-round(100*sum(age_cat_4$Systolic_Average>=160 | age_cat_4$Diastolic_Average>=100)/dim(age_cat_4)[1],1)

Percentage <- c(p1,p2,p3,p4,p5,p6,p7,p8)
group<- c(rep("Stage1",4),rep("Stage2",4))
age_cat<- c("35-44","45-54","55-64","65-75")
Age <-rep(age_cat,2)

dat <- cbind(group,Age,Percentage)
dat<-as.data.frame(dat)
dat$Percentage <- lapply(dat$Percentage,function(x) as.numeric(as.character(x)))

dat$Percentage<-as.numeric(dat$Percentage)
dat$group<-as.character(dat$group)
dat$Age<-as.character(dat$Age)

ggplot(data=dat, aes(x=Age, y=Percentage, fill=group))+geom_bar(stat="identity", position=position_dodge())+labs(x="Age category",fill="Group")+ggtitle("Proportion of high blood pressure\n in each age groups")

####sd BMI vs. Age###
age_cat_1 <-  data_clean[which(data_clean$Age >=35 & data_clean$Age <45),]
age_cat_2 <-  data_clean[which(data_clean$Age >=45 & data_clean$Age <55),]
age_cat_3 <-  data_clean[which(data_clean$Age >=55 & data_clean$Age <65),]
age_cat_4 <-  data_clean[which(data_clean$Age >=65 & data_clean$Age <=75),]
b1<-round(mean(age_cat_1$BMI),1)
b2<-round(mean(age_cat_2$BMI),1)
b3<-round(mean(age_cat_3$BMI),1)
b4<-round(mean(age_cat_4$BMI),1)

bmi <- c(b1,b2,b3,b4)
age_cat<- c("35-44","45-54","55-64","65-75")

dat <- cbind(age_cat,bmi)
dat<-as.data.frame(dat)

ggplot(data=dat, aes(x=age_cat, y=bmi))+geom_bar(stat="identity",fill="steelblue")+labs(x="Age category")+ggtitle("Mean BMI in each age groups")


age <- seq(35,75,by=1)
for(i in 1:41){
  p[i] <- mean(data_clean$BMI[which(data_clean$Age==age[i])])
}

plot(age,p,main='Mean BMI for each age',ylab="BMI",xlab="Age")

m <- vector()
for(i in 1:41){
  m[i] <- mean(data_clean$LDL[which(data_clean$Age==age[i])])
}

plot(age,m,main='Mean LDL for each age',ylab="LDL",xlab="Age")
m <- vector()
for(i in 1:41){
  m[i] <- mean(data_clean$Systolic_Average[which(data_clean$Age==age[i])])
}

m <- vector()
for(i in 1:41){
  m[i] <- mean(data_clean$pp[which(data_clean$Age==age[i])])
}

plot(age,m,main='Mean SBP for each age',ylab="SBP",xlab="Age")

m <- vector()
for(i in 1:41){
  m[i] <- mean(data_clean$Diastolic_Average[which(data_clean$Age==age[i])])
}

plot(age,m,main='Mean DBP for each age',ylab="DBP",xlab="Age")

#######eFigure1#####

title_theme =  theme(plot.title=element_text(face='plain'))
#title_theme = theme(plot.title=element_text(family = "Helvetica", face = "bold", size = (15)))
require(cowplot)
theme_set(theme_cowplot(font_size=10.5))

require(scales)

f1_a <- ggplot(data_clean, aes(x=Systolic_Average)) + geom_histogram() + xlab("Systolic BP (mmHg)") + 
  background_grid(major="xy", minor="none")+ scale_y_continuous(labels=comma)+ ggtitle("(a)")  + title_theme
f1_b <- ggplot(subset(data_clean, Diastolic_Average< 200), aes(x=Diastolic_Average)) + geom_histogram() + xlab("Diastolic BP (mmHg)") + 
  background_grid(major="xy", minor="none")+ scale_y_continuous(labels=comma)+ ggtitle("(b)") + title_theme
f1_c <- ggplot(data_clean, aes(x=Age)) + geom_histogram() + xlab("Age (years)")+ 
  background_grid(major="xy", minor="none") + scale_y_continuous(labels=comma)+ ggtitle("(c)") + title_theme

g<-grid.arrange(f1_a,  f1_b, f1_c, ncol=3)
#ggsave(g,filename = "eFigure1.pdf",width = 8,height=3)


f1_d <- ggplot(data_clean, aes(x=TC)) + geom_histogram() + xlab("TC (mg/dl)") + 
  background_grid(major="xy", minor="none")+ scale_y_continuous(labels=comma)+ ggtitle("(d)")  + title_theme+xlim(0,200)
f1_e <- ggplot(data_clean, aes(x=HDL)) + geom_histogram() + xlab("HDL-C (mg/dl)") + 
  background_grid(major="xy", minor="none")+ scale_y_continuous(labels=comma)+ ggtitle("(e)") + title_theme+xlim(0,100)
f1_f <- ggplot(data_clean, aes(x=TG)) + geom_histogram() + xlab("TG (mg/dl)")+ 
  background_grid(major="xy", minor="none") + scale_y_continuous(labels=comma)+ ggtitle("(f)") + title_theme+xlim(0,200)
par(mar=c(2.1, 1.1, 2.1, 4.1))
f1_g <- ggplot(data_clean, aes(x=LDL)) + geom_histogram() + xlab("LDL-C (mg/dl)")+ 
  background_grid(major="xy", minor="none") + scale_y_continuous(labels=comma)+ ggtitle("(g)") + title_theme+xlim(0,200)
g2<-grid.arrange(f1_d,  f1_e, f1_f,f1_g, ncol=4)

########figure 3
lm_subgroup_by <- function(SA, Age,by_){
  byvec <- rep(FALSE,length(SA))
  byvec[by_] <- TRUE;
  outlm <- lm (SA[by_] ~ Age[by_])
  lower <- confint.default(outlm)[2,1]
  upper <- confint.default(outlm)[2,2]
  
  spear <- cor(SA[by_],Age[by_],method="spearman")
  c("coeflm"= outlm$coefficients[2], "lowerlm"= lower, "upperlm" =upper,spear=spear, "N"=length(by_));
}
subgroup_feats2 <- c("Sex.x", "Ethnicity", "Hukou", "Occupation", "Education", "Household_Income",  "Marital_Status", "Province", "tx")
ks <- 1:6;
#thresh <- 1E4
thresh <- 1E4
subgroup_results <- data.frame();
for (k in ks) {
  subsub_feats <- combn(subgroup_feats2, k)
  library(data.table)
  data_subgroups <- data_clean[,c("Systolic_Average", "Age", subgroup_feats2)]
  setDT(data_subgroups)
  for (i in 1:ncol(subsub_feats)){
    print(sprintf("%d/%d completed", i, ncol(subsub_feats)))
    subgroup_results_i <- data.frame(data_subgroups[,  as.list(lm_subgroup_by(data_clean$Systolic_Average, data_clean$Age,.I)),by = eval(subsub_feats[,i])])
    if (k > 1){
      subgroup_name <- apply( subgroup_results_i[,1:k], 1, paste, collapse=",")
    }else{
      subgroup_name <- subgroup_results_i[,1]
    }
    combined_df <- rbind.fill(data.frame(subgroup_results_i, group=subgroup_name))
    d =subset(combined_df, N>thresh)
    subgroup_results<- rbind.fill(subgroup_results,d)
  }
  
}



subgroup_results_dbp <- data.frame();
for (k in ks) {
  subsub_feats <- combn(subgroup_feats2, k)
  library(data.table)
  data_subgroups <- data_clean[,c("Diastolic_Average", "Age", subgroup_feats2)]
  setDT(data_subgroups)
  for (i in 1:ncol(subsub_feats)){
    print(sprintf("%d/%d completed", i, ncol(subsub_feats)))
    subgroup_results_i <- data.frame(data_subgroups[,  as.list(lm_subgroup_by(data_clean$Diastolic_Average, data_clean$Age,.I)),by = eval(subsub_feats[,i])])
    if (k > 1){
      subgroup_name <- apply( subgroup_results_i[,1:k], 1, paste, collapse=",")
    }else{
      subgroup_name <- subgroup_results_i[,1]
    }
    d =subset(rbind.fill(data.frame(subgroup_results_i, group=subgroup_name)), N>thresh)
    subgroup_results_dbp<- rbind.fill(subgroup_results_dbp,d)
  }
  
}
subgroup_results_pp <- data.frame();
for (k in ks) {
  subsub_feats <- combn(subgroup_feats2, k)
  library(data.table)
  data_subgroups <- data_clean[,c("pp", "Age", subgroup_feats2)]
  setDT(data_subgroups)
  for (i in 1:ncol(subsub_feats)){
    print(sprintf("%d/%d completed", i, ncol(subsub_feats)))
    subgroup_results_i <- data.frame(data_subgroups[,  as.list(lm_subgroup_by(data_clean$pp, data_clean$Age,.I)),by = eval(subsub_feats[,i])])
    if (k > 1){
      subgroup_name <- apply( subgroup_results_i[,1:k], 1, paste, collapse=",")
    }else{
      subgroup_name <- subgroup_results_i[,1]
    }
    d =subset(rbind.fill(data.frame(subgroup_results_i, group=subgroup_name)), N>thresh)
    subgroup_results_pp<- rbind.fill(subgroup_results_pp,d)
  }
  
}
write.csv(subgroup_results,"/mpp/Yahui Zhao/age project/sbp_subgroup.csv")
write.csv(subgroup_results_dbp,"/mpp/Yahui Zhao/age project/dbp_subgroup.csv")
write.csv(subgroup_results_pp,"/mpp/Yahui Zhao/age project/pp_subgroup.csv")
toplot <- cbind(x=1:nrow(subgroup_results), subgroup_results[order(-subgroup_results$N),])
f6_a <- ggplot(toplot, aes(x=coeflm.Age.by_.)) + geom_histogram(bins=40) + 
  xlab("Regression Coefficient") +ylab('Number of Subgroups') + ggtitle('(a)') 
f6_b <- ggplot(subset(toplot,!is.na(toplot$tx)), aes(x =coeflm.Age.by_. , fill=factor(tx))) + xlab("Regression Coefficients")  + 
  geom_density( position="identity", alpha=0.5 ) +labs (fill="Treatment")+ ggtitle('(b)')


toplot <- cbind(x=1:nrow(subgroup_results_dbp), subgroup_results_dbp[order(-subgroup_results_dbp$N),])
f7_a <- ggplot(toplot, aes(x=coeflm.Age.by_.)) + geom_histogram(bins=40) +ylab('Number of Subgroups')+ xlab("Regression Coefficients") + ggtitle('(c)')
f7_b <- ggplot(subset(toplot,!is.na(toplot$tx)), aes(x =coeflm.Age.by_. , fill=factor(tx))) + xlab("Subgroup Coefficients")  + geom_density( position="identity", alpha=0.5 ) +labs (fill="Treatment")+ ggtitle('(d)') 
g<-grid.arrange(f6_a, f6_b,f7_a,f7_b,ncol=2)


toplot <- cbind(x=1:nrow(subgroup_results_pp), subgroup_results_pp[order(-subgroup_results_pp$N),])
f8_a <- ggplot(toplot, aes(x=coeflm.Age.by_.)) + geom_histogram(bins=40) +ylab('Number of Subgroups')+ xlab("Regression Coefficients") + ggtitle('(e)')
f8_b <- ggplot(subset(toplot,!is.na(toplot$tx)), aes(x =coeflm.Age.by_. , fill=factor(tx))) + xlab("Subgroup Coefficients")  + geom_density( position="identity", alpha=0.5 ) +labs (fill="Treatment")+ ggtitle('(f)') 
