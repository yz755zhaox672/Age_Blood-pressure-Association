require(ggplot2)
require(grid)
require(gridExtra)
library(cowplot)
library(scales)
require(reshape)
require(egg)
library(data.table)
library(plyr)

#load('~/BP_Paper/baseline_data_gcl_1.7.Rdata')
load('/mpp/Qiaozhi/baseline_data_gcl_1.7.Rdata')
# Load in data ------------------------------------------------------------
is <- read.csv('/mpp/1.7M_CSV/MPP_IS_20170620.csv',header = TRUE, fileEncoding='GB18030')
med_raw <- read.csv('/mpp/1.7M_CSV/MPP_MED_20170620.csv',header = TRUE, fileEncoding='GB18030')
med <- med_raw[med_raw$mark=="InitialS",c("PID", "drugid")]
med$drugid <- as.character(med$drugid)
med$drugid[med$drugid==""] <- NA
med <- na.omit(med)
medsplit <- do.call("rbind",strsplit(med$drugid,"_"))
med3 <- data.frame(PID=med$PID,chem=as.numeric(medsplit[,1]))

drugbp <- read.csv("/mpp/codes/codes_xinyue/HTN_drug.csv",header=T,as.is=T, fileEncoding='GB18030')
drugbp <- drugbp[,c(1,2,4,5)]
names(drugbp) <- c("chem_id","chem_name","class_id","class_name")
drugbp$pid <- paste(drugbp$class_id,
                    drugbp$chem_id,sep="_")
med4 <- med3[med3$chem %in% drugbp$chem_id,] ##from 399665 to 248949, select BP drugs
medWM <- med4[!duplicated(med4),] 
is$tx <- rep(FALSE,nrow(is))
is$tx[is$PID %in% medWM$PID] <- TRUE

raw_baseline_data <- is
#save(raw_baseline_data, file="/mpp/Qiaozhi/baseline_data_gcl_1.7.Rdata")
# Load from RData file below ----------------------------------------------
#load("/mpp/Qiaozhi/baseline_data_gcl_1.7.Rdata")
baseline_data <- raw_baseline_data[!duplicated(raw_baseline_data$PID),]

#save(raw_baseline_data,baseline_data_merge, file = 'baseline_data_gcl.Rdata')
# Preprocess --------------------------------------------------------------
data <- subset(baseline_data, select = c(PID, IS_gender, is_AGE, is_income, IS_education, IS_occupation, IS_hukou, 
                                         IS_Prov, IS_marriage, 
                                          IS_Stroke_yes, IS_SBP, IS_DBP,
                                          IS_smk,  IS_BMI,
                                         tx, IS_Ethnicity
))
names(data) <- c("Patient_Id", "Sex.x", "Age", "Household_Income", "Education", "Occupation", "Hukou",
                 "Province", "Marital_Status", 
                 "Hx_Stroke","Systolic_Average", "Diastolic_Average",  
                 "Smoking_Currently",  "BMI",
                  "tx", "Ethnicity"
)
eth_tbl <- table(data$Ethnicity)
large_ethnicities <- names(eth_tbl[order(-eth_tbl)])[1:12]
data_clean <- na.omit(subset(data,  Age <=80 &Age >= 35& Ethnicity %in% large_ethnicities  ))
#testdata <- omit.na(data_clean[,c("Sex.x", "Age_cat2", "Ethnicity",  "Hukou", "Occupation", "Education", "Household_Income",  "Marital_Status", "Province", "tx")))


sprintf("Removed %f individuals not in major ethnicities ", nrow(data) - sum(data$Ethnicity %in% large_ethnicities))
sprintf("Went from %d to %d after filtering, %f ", nrow(data), nrow(data_clean), nrow(data_clean)/nrow(data))

data_clean$Age_cat <- factor(cut(data_clean$Age, breaks = quantile(data_clean$Age,probs = seq(0,1,0.2)) ))
data_clean$Age_cat2 <- factor(cut(data_clean$Age, breaks = c(35, 50, 60, 70, 80, 88),include.lowest = TRUE))

data_clean$Hukou <- factor(data_clean$Hukou)
data_clean$Province <- factor(data_clean$Province)
data_clean$Sex.x <- factor(data_clean$Sex.x)
data_clean$Education <- factor(data_clean$Education)
data_clean$Household_Income <- factor(data_clean$Household_Income)
data_clean$Occupation <- factor(data_clean$Occupation)
data_clean$Marital_Status <- factor(data_clean$Marital_Status == 1)
data_clean$Smoking_Currently <- factor(data_clean$Smoking_Currently)
data_clean$Hx_Stroke <- factor(data_clean$Hx_Stroke)
data_clean$tx <- factor(data_clean$tx)

levels(data_clean$Sex.x) <- c("Male", "Female")
#levels(data_clean$Occupation) <- c("Farmer", "Retire", "Unemployment", "Housework", "Unknown", "Refuse", "Worker", "Manager", "Admin", "Technician", "Serv. Ind.", "Bus. Owner", "Military", "Other")
levels(data_clean$Occupation) <- c("Farmer", "Workers", "Administrators", "Admin. Clerk", "Technician", "Business", "Bus. Owner", "Military", "Others", "Retire", "Unemployed", "Housework", "Unknown", "Refuse")
levels(data_clean$Marital_Status) <- c("Unmarried", "Married")
levels(data_clean$Hukou) <- c("Rural", "Urban", "Unified","No Hukou")
levels(data_clean$Hx_Stroke) <- c("No history", "Has history")
levels(data_clean$Smoking_Currently) <- c("Smoking", "Not Smoking")
data_clean$Education[data_clean$Education==13] <- 12
data_clean$Education <- droplevels(data_clean$Education)
#levels(data_clean$Education) <- c("Illiterate", "Masters", "Doctoral", "No answer",  "<Primary", "Homeschool", "Elementary",  "Middle",  "High", "Vocational", "Associate",  "Bachelor")
levels(data_clean$Education) <- c("Illiterate","<Primary", "Homeschool", "Elementary", 
                                   "Middle",  "High", "Vocational", "Associate",  "Bachelor",
                                  "Masters", "Doctoral", "No answer")
levels(data_clean$tx) <- c("No meds", "On meds")




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

# overview --------------------------------------------------------------

nrow(data_clean)
mean(data_clean$Age)
sd(data_clean$Age)
aggregate(data_clean$BMI, by=list(data_clean$Sex.x), FUN = mean)
mean(data_clean$BMI)
sd(data_clean$BMI)
mean(data_clean$Systolic_Average)
sd(data_clean$Systolic_Average)
mean(data_clean$Diastolic_Average)
sd(data_clean$Diastolic_Average)

mean(data_clean$BMI < 18.5)
mean(data_clean$BMI >= 18.5 & data_clean$BMI < 25)
mean(data_clean$BMI >= 25 & data_clean$BMI < 30)
mean(data_clean$BMI >= 30)

mean(data_clean$Systolic_Average >= 140 | data_clean$Diastolic_Average >= 90)
mean(data_clean$Systolic_Average >= 160 | data_clean$Diastolic_Average >= 100)

# Spearman's --------------------------------------------------------------

cor(as.numeric(data_clean$BMI),data_clean$Systolic_Average,method = "spearman")
cor(as.numeric(data_clean$Hx_Stroke),data_clean$Systolic_Average,method = "spearman")
cor(as.numeric(data_clean$Marital_Status),data_clean$Systolic_Average,method = "spearman")
cor(as.numeric(data_clean$Age),data_clean$Systolic_Average,method = "spearman")



cor(as.numeric(data_clean$BMI),data_clean$Diastolic_Average,method = "spearman")
cor(as.numeric(data_clean$Age),data_clean$Diastolic_Average,method = "spearman")

ggplot( subset(data_clean, BMI < 50 & BMI > 10), aes(BMI, Diastolic_Average) ) +  stat_binhex() +geom_smooth() + ggtitle('2a')+
 theme(legend.position="none") + scale_fill_gradient(low="#ffe5e5", high = "#ff0000") + scale_y_continuous(name="Systolic BP (mmHg)")+ xlab("Body Mass Index (kg/m^2)") + title_theme


ggplot( subset(data_clean, BMI < 50 & BMI > 10), aes(Age, Diastolic_Average) ) +  stat_binhex() +geom_smooth() + ggtitle('2a')+
 theme(legend.position="none") + scale_fill_gradient(low="#ffe5e5", high = "#ff0000") + scale_y_continuous(name="Systolic BP (mmHg)")+ xlab("Body Mass Index (kg/m^2)") + title_theme


temp <- subset(data_clean,Household_Income!="Unclear")
cor(as.numeric(temp$Household_Income),temp$Systolic_Average,method = "spearman")


# Figure 1 & 2 --------------------------------------------------------------
title_theme =  theme(plot.title=element_text(face='plain'))
#title_theme = theme(plot.title=element_text(family = "Helvetica", face = "bold", size = (15)))
theme_set(theme_cowplot(font_size=10.5))
require(cowplot)
require(scales)

f1_a <- ggplot(data_clean, aes(x=Systolic_Average)) + geom_histogram() + xlab("Systolic BP (mmHg)") + 
  background_grid(major="xy", minor="none")+ scale_y_continuous(labels=comma)+ ggtitle("(a)")  + title_theme
f1_b <- ggplot(subset(data_clean, Diastolic_Average< 200), aes(x=Diastolic_Average)) + geom_histogram() + xlab("Diastolic BP (mmHg)") + 
  background_grid(major="xy", minor="none")+ scale_y_continuous(labels=comma)+ ggtitle("(b)") + title_theme
f1_c <- ggplot(subset(data_clean,BMI < 50 ), aes(x=BMI)) + geom_histogram() + xlab("BMI (kg/m^2)")+ 
  background_grid(major="xy", minor="none") + scale_y_continuous(labels=comma)+ ggtitle("(c)") + title_theme
g<-grid.arrange(f1_a,  f1_b, f1_c, ncol=3)
ggsave(g,filename = "eFigure1.pdf",width = 8,height=3)

f1_d <- ggplot( subset(data_clean, BMI < 50 & BMI > 10), aes(BMI, Systolic_Average) ) +  stat_binhex() +geom_smooth() + ggtitle('(a)')+
 theme(legend.position="none") + scale_fill_gradient(low="#ffe5e5", high = "#ff0000") + scale_y_continuous(name="Systolic BP (mmHg)")+ xlab("Body Mass Index (kg/m^2)") + title_theme
f1_e <- ggplot( subset(data_clean, BMI < 50 & BMI > 10 & Diastolic_Average < 200), aes(BMI, Diastolic_Average) ) + ggtitle('(b)') +
  stat_binhex() +geom_smooth() + theme(legend.position="none") + scale_fill_gradient(low="#ffe5e5", high = "#ff0000") + scale_y_continuous(name="Diastolic BP (mmHg)")+ xlab("Body Mass Index (kg/m^2)") + title_theme

g <- grid.arrange(f1_d,  f1_e, ncol=2)
ggsave(g,filename = "figure1.pdf",width = 6.5,height=3)

summary(lm(Systolic_Average ~ BMI + I(BMI^2) + I(BMI^3), data = subset(data_clean, BMI > 18.5 & BMI< 30)))

# Figure 3 GAMS --------------------------------------------------------------------
require(egg)

data_clean_sample_idx <- 1:nrow(data_clean)
  data_sub <- subset(data_clean[data_clean_sample_idx,], BMI < 35 & BMI > 18.5 )

plot_titles <- c(`Sex.x`="Sex",  `Hukou`="Hukou", `Occupation`="Occupation", 
                 `Education`="Education", `Household_Income`="Household Income",  `Marital_Status`="Marital Status", 
                 `Province`="Province" ,`tx`="BP Meds" , `Ethnicity`="Ethnicity", `Age_cat`="Age", `Age_cat2`="Age", `Hx_Stroke`="History of Stroke",`Smoking_Currently`="Currently Smoking") 


#covariates <- c("Sex.x", "Age_cat2", "Household_Income", 
                #"Occupation", "Ethnicity", "Marital_Status",  "Hukou","Province","Education",  "tx")

covariates <- c("Sex.x", "Age_cat2", "Household_Income", 
                "Occupation", "Ethnicity", "Marital_Status",  "Hukou","Province","Education","Smoking_Currently","Hx_Stroke",  "tx")

df <- melt(data_sub[,c("Systolic_Average", "BMI", covariates )], id =c("Systolic_Average", "BMI"))

#950 x 1220
theme_set(theme_cowplot(font_size=10.5))
g <-plot_grid_smooths(df, data_clean,plot_titles, 2, "BMI", "SBP", c(120, 160))
ggsave(g,filename = "figure2.pdf",width = 0.8*12,height=0.8*17)







# Figure 4 GAMS --------------------------------------------------------------------
theme_set(theme_cowplot(font_size=10.5))
df <- melt(data_sub[,c("Diastolic_Average", "BMI", covariates )], id =c("Diastolic_Average", "BMI"))
g<- plot_grid_smooths_dbp(df, plot_titles, 2, "BMI", "DBP", c(60, 120))
ggsave(g,filename = "eFigure2.pdf",width = 0.8*12,height=0.8*17)

# Table 1 --------------------------------------------------------------------

# subgroup_feats <- c("All", "Sex.x", "Age_cat2", "Ethnicity",  "Hukou", "Occupation", "Education", "Household_Income",  "Marital_Status", "Province", "tx")
# for (var in subgroup_feats) {
#   cat(sprintf("\n%s\n", plot_titles[var]))
#   for (lvl in levels(data_clean_dum[,var])) {
#     data_sub_sub <- data_clean_dum[ data_clean_dum[,var] == lvl,]
#     n_sub_sub <- nrow(data_sub_sub);
#     
#       sbp_ave <- mean(data_sub_sub$Systolic_Average)
#       sbp_sd <- sd(data_sub_sub$Systolic_Average)
#       sbp_str <- sprintf('%.1f (%.1f)', sbp_ave, sbp_sd)
#       
#       dbp_ave <- mean(data_sub_sub$Diastolic_Average)
#       dbp_sd <- sd(data_sub_sub$Diastolic_Average)
#       dbp_str <- sprintf('%.1f (%.1f)', dbp_ave, dbp_sd)
#       n_str <- format(n_sub_sub,big.mark=",", scientific = FALSE)
#       
#       cat(sprintf("%s=%s=%s=%s\n", lvl,n_str, sbp_str, dbp_str ))
#   }
# }

# Table 2 Subgroup slopes  ---------------------------------------------------------------

predictors_htn <- c("BMI", "Age",  'Sex.x', 'Education', 'Ethnicity','Household_Income', 'Hukou', 'Occupation', 'Marital_Status', 'Smoking_Currently', "Hx_Stroke", "Province","tx");


data_clean_dum <- cbind(data_clean, "All" =factor(rep(1,nrow(data_clean))))
subgroup_feats <- c("All", covariates)
for (var in subgroup_feats) {
  cat(sprintf("\n%s\n", plot_titles[var]))
  for (lvl in levels(data_clean_dum[,var])) {
    data_sub_sub <- data_clean_dum[ data_clean_dum[,var] == lvl,]
    n_sub_sub <- nrow(data_sub_sub);
    if (n_sub_sub > 5E3){
      varlvl <- paste("BMI:",paste(var,lvl, sep = ""), sep="")
      
      outlm <- lm (Systolic_Average ~ BMI, data = data_sub_sub)
      if (var=="Age_cat2"){
        predictors_htn2 <- predictors_htn[!("Age" ==predictors_htn) ]
      }else{
        predictors_htn2 <- predictors_htn[!(var ==predictors_htn) ]
      }
      
      outlm_adj <- lm (Systolic_Average ~ ., data = data_sub_sub[,c("Systolic_Average", predictors_htn2)])
      coef <- coef_confint ( outlm, "BMI")
      coef_pm <- coef_confint ( outlm_adj, "BMI")
      
      outlm_dbp <- lm (Diastolic_Average ~ BMI, data = data_sub_sub)
      outlm_adj_dbp <- lm (Diastolic_Average ~ ., data = data_sub_sub[,c("Diastolic_Average", predictors_htn2)])
      coef_dbp <- coef_confint ( outlm_dbp, "BMI")
      coef_pm_dbp <- coef_confint ( outlm_adj_dbp, "BMI")
      
      spear_sbp <- cor(data_sub_sub$Systolic_Average, data_sub_sub$BMI,method = "spearman")
      spear_dbp <- cor(data_sub_sub$Diastolic_Average, data_sub_sub$BMI,method = "spearman")
      
      spear_sbp_confint <- spearman_confint(spear_sbp,nrow(data_sub_sub))
      spear_dbp_confint <- spearman_confint(spear_dbp,nrow(data_sub_sub))
      
      
      #cat(sprintf("%s=%s=%s=%s=%s=%s\n", lvl, format(n_sub_sub,big.mark=",", scientific = FALSE), coef, coef_pm,
                   #coef_dbp, coef_pm_dbp))
      
      
      sbp_ave <- mean(data_sub_sub$Systolic_Average)
      sbp_sd <- sd(data_sub_sub$Systolic_Average)
      sbp_str <- sprintf('%.1f (%.1f)', sbp_ave, sbp_sd)
      
      dbp_ave <- mean(data_sub_sub$Diastolic_Average)
      dbp_sd <- sd(data_sub_sub$Diastolic_Average)
      dbp_str <- sprintf('%.1f (%.1f)', dbp_ave, dbp_sd)
      bp_str <- sprintf("%s/%s\n", sbp_str, dbp_str )
      bp2_str <- sprintf("%.0f/%.0f (%.0f/%.0f)", sbp_ave, dbp_ave, sbp_sd,dbp_sd )
      n_str <- format(n_sub_sub,big.mark=",", scientific = FALSE)
      
      cat(sprintf("%s=%s=%s=%s=%s=%s=%s=%s=%s\n", lvl,n_str , bp2_str,coef, coef_pm,spear_sbp_confint,
                   coef_dbp, coef_pm_dbp,spear_dbp_confint))
    }
  }
}


# All combos --------------------------------------------------------------






#Using by in a cooler way
lm_subgroup_by <- function(SA, BMI,by_){
  byvec <- rep(FALSE,length(SA))
  byvec[by_] <- TRUE;
  outlm <- lm (SA[by_] ~ BMI[by_])
  lower <- confint.default(outlm)[2,1]
  upper <- confint.default(outlm)[2,2]
  
  spear <- cor(SA[by_],BMI[by_],method="spearman")
  c("coeflm"= outlm$coefficients[2], "lowerlm"= lower, "upperlm" =upper,spear=spear, "N"=length(by_));
}

subgroup_feats2 <- c("Sex.x", "Age_cat2",  "Ethnicity", "Hukou", "Occupation", "Education", "Household_Income",  "Marital_Status", "Province", "tx")
#subgroup_feats2 <- covariates

# How many subgroups are there?

for (i in 1:length(subgroup_feats2)){
 print(sprintf("%s: %d", subgroup_feats2[i], length(levels(data_clean[,subgroup_feats2[i]]))) )
}

num_groups <- 0;
for (k in 1:6) {
  subsub_feats <- combn(subgroup_feats2, k)
  for (j in 1:ncol(subsub_feats)){
    print(sprintf("%d/%d completed", j, ncol(subsub_feats)))
    vars <- paste(plot_titles[subsub_feats[,j]],collapse = ", ")
    mults <- sapply( subsub_feats[,j], function(x) length(levels(data_clean[,x])))
    #cat(sprintf("%s: %s = %d\n",vars, paste(mults, collapse = "*"),prod(mults) )  )
    num_groups <- prod(mults) + num_groups;
  }
}
print(num_groups)


ks <- 1:6;
#thresh <- 1E4
thresh <- 1E4
subgroup_results <- data.frame();
for (k in ks) {
  subsub_feats <- combn(subgroup_feats2, k)
  library(data.table)
  data_subgroups <- data_clean[,c("Systolic_Average", "BMI", subgroup_feats2)]
  setDT(data_subgroups)
  for (i in 1:ncol(subsub_feats)){
    print(sprintf("%d/%d completed", i, ncol(subsub_feats)))
    subgroup_results_i <- data.frame(data_subgroups[,  as.list(lm_subgroup_by(data_clean$Systolic_Average, data_clean$BMI,.I)),by = eval(subsub_feats[,i])])
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
#print(subgroup_results)
nrow(subgroup_results)
print(num_groups)

# subgroup_results <- subset(subgroup_results, N>1E4)

#For checking
#summary(lm(Systolic_Average ~ BMI, data = subset(data_clean,Sex.x =="Female" & Age_cat == "(52,59]")))
#save(subgroup_results,file="~/BP_Paper/subgroup_results.RData")
summary(lm(Systolic_Average ~ BMI, data = subset(data_clean,Household_Income != "<5k RMB"& Province == "Tibet" & tx == "No meds" & Ethnicity == "Tibet")))

# All combos for Diastolic_Average--------------------------------------------------------------
#Using by in a cooler way
subgroup_results_dbp <- data.frame();
for (k in ks) {
  subsub_feats <- combn(subgroup_feats2, k)
  library(data.table)
  data_subgroups <- data_clean[,c("Diastolic_Average", "BMI", subgroup_feats2)]
  setDT(data_subgroups)
  for (i in 1:ncol(subsub_feats)){
    print(sprintf("%d/%d completed", i, ncol(subsub_feats)))
    subgroup_results_i <- data.frame(data_subgroups[,  as.list(lm_subgroup_by(data_clean$Diastolic_Average, data_clean$BMI,.I)),by = eval(subsub_feats[,i])])
    if (k > 1){
      subgroup_name <- apply( subgroup_results_i[,1:k], 1, paste, collapse=",")
    }else{
      subgroup_name <- subgroup_results_i[,1]
    }
    d =subset(rbind.fill(data.frame(subgroup_results_i, group=subgroup_name)), N>thresh)
    subgroup_results_dbp<- rbind.fill(subgroup_results_dbp,d)
  }

}
print(subgroup_results_dbp)
nrow(subgroup_results_dbp)
print(under_thresh)

#For checking
#summary(lm(Systolic_Average ~ BMI, data = subset(data_clean,Sex.x =="Female" & Age_cat == "(52,59]")))
#save(subgroup_results,file="~/BP_Paper/subgroup_results.RData")
summary(lm(Diastolic_Average ~ BMI, data = subset(data_clean,Household_Income != "<5k RMB"& Province == "Tibet" & tx == "No meds" & Ethnicity == "Tibet")))

# Figure 5 --------------------------------------------------------

#f5_a <- ggplot( subset(toplot, !is.na(tx)), aes(x=x,y=coeflm.BMI.by_.,ymin=lowerlm, ymax=upperlm, color=(tx) )) +  geom_point(size=0.9) + ylab('Coefs') + xlab('') + 
  #labs (colour = "Treatment")+ theme(axis.ticks.x=element_blank(), axis.text.x = element_blank())  + ggtitle('5a') + title_theme
#f5_b <- ggplot( subset(toplot, !is.na(Sex.x)), aes(x=x,y=coeflm.BMI.by_.,ymin=lowerlm, ymax=upperlm, color=(Sex.x) )) +  geom_point(size=0.9) + ylab('Coefs') + xlab('')+ 
  #labs(colour="Sex") + theme(axis.ticks.x=element_blank(), axis.text.x = element_blank())+ ggtitle('5b') + title_theme
#f5_c <- ggplot( subset(toplot, !is.na(Province) & tx == "No meds"), aes(x=x,y=coeflm.BMI.by_.,ymin=lowerlm, ymax=upperlm, color=factor(Province == "Tibet" ) )) +  geom_point(size=0.9) + ylab('Coefs') + xlab('')+ 
  #labs(colour="Tibet") + theme(axis.ticks.x=element_blank(), axis.text.x = element_blank())+ ggtitle('5b') + title_theme
#grid.arrange(f5_a, f5_c, f5_b, ncol=1) # 750 368


toplot <- cbind(x=1:nrow(subgroup_results), subgroup_results[order(-subgroup_results$N),])
f6_a <- ggplot(toplot, aes(x=coeflm.BMI.by_.)) + geom_histogram(bins=40) + 
  xlab("Regression Coefficient") +ylab('Number of Subgroups') + ggtitle('(a)') + title_theme
f6_b <- ggplot(subset(toplot,!is.na(toplot$tx)), aes(x =coeflm.BMI.by_. , fill=factor(tx))) + xlab("Regression Coefficients")  + 
  geom_density( position="identity", alpha=0.5 ) +labs (fill="Treatment")+ ggtitle('(b)') + title_theme

toplot <- cbind(x=1:nrow(subgroup_results_dbp), subgroup_results_dbp[order(-subgroup_results_dbp$N),])
f7_a <- ggplot(toplot, aes(x=coeflm.BMI.by_.)) + geom_histogram(bins=40) +ylab('Number of Subgroups')+ xlab("Regression Coefficients") + 
  ggtitle('(c)') + title_theme
f7_b <- ggplot(subset(toplot,!is.na(toplot$tx)), aes(x =coeflm.BMI.by_. , fill=factor(tx))) + xlab("Subgroup Coefficients")  + 
  geom_density( position="identity", alpha=0.5 ) +labs (fill="Treatment")+ ggtitle('(d)') + title_theme
g<-grid.arrange(f6_a, f6_b,f7_a,f7_b,ncol=2)

ggsave(g,filename = "figure3.pdf",width = 8,height=6)



toplot <- cbind(x=1:nrow(subgroup_results), subgroup_results[order(-subgroup_results$N),])
f8_a <- ggplot(toplot, aes(x=spear)) + geom_histogram(bins=40) + 
  xlab("Spearman Correlation") +ylab('Number of Subgroups') + ggtitle('(a)') + title_theme
f8_b <- ggplot(subset(toplot,!is.na(toplot$tx)), aes(x =spear , fill=factor(tx))) + xlab("Spearman Correlation")  + 
  geom_density( position="identity", alpha=0.5 ) +labs (fill="Treatment")+ ggtitle('(b)') + title_theme

toplot <- cbind(x=1:nrow(subgroup_results_dbp), subgroup_results_dbp[order(-subgroup_results_dbp$N),])
f8_c <- ggplot(toplot, aes(x=spear)) + geom_histogram(bins=40)+ylab('Number of Subgroups') + xlab("Spearman Correlation") + 
  ggtitle('(c)') + title_theme
f8_d <- ggplot(subset(toplot,!is.na(toplot$tx)), aes(x =spear , fill=factor(tx))) + xlab("Spearman Correlation")  + 
  geom_density( position="identity", alpha=0.5 ) +labs (fill="Treatment")+ ggtitle('(d)') + title_theme
g<-grid.arrange(f8_a, f8_b,f8_c,f8_d,ncol=2)

ggsave(g,filename = "eFigure3.pdf",width = 8,height=6)





# Quantiles --------------------------------------------------------
quantile(toplot$coeflm.BMI.by_.[toplot$tx=="No meds"])
#quantile(subset(toplot252,727 and 61,842 )$coeflm.BMI.by_.,probs = c(.025,.975))
quantile(subset(toplot, !is.na(tx) & tx == "No meds" )$coeflm.BMI.by_.,probs = c(.025,.975))
quantile(subset(toplot, !is.na(tx) & tx == "On meds" )$coeflm.BMI.by_.,probs = c(.025,.975))
sd(subset(toplot, !is.na(tx) & tx == "No meds" )$coeflm.BMI.by_.)

# Figure 6 Matching ----------------------------------------------------------------

Systolic_Average_Cut <- cut(data_clean$Systolic_Average,breaks = 100)

#Simple density estimators
p_density<- table(Systolic_Average_Cut[data_clean$tx == "No meds"])/sum(data_clean$tx == "No meds")
q_density<- table(Systolic_Average_Cut[data_clean$tx == "On meds"])/sum(data_clean$tx == "On meds")

#Evaluate density of both distribution at each point of group 1 
p <- p_density[Systolic_Average_Cut[data_clean$tx == "No meds"]]
q <- q_density[Systolic_Average_Cut[data_clean$tx=="No meds"]]

matching_idx <- sample(which(data_clean$tx == "No meds"),sum(data_clean$tx=="On meds"), FALSE, q/p )

data_clean$untreated_matching <- NA
data_clean$untreated_matching[data_clean$tx == "On meds"] = "On meds"
data_clean$untreated_matching[matching_idx] = "No meds (matching)"
data_clean$untreated_matching <- factor(data_clean$untreated_matching)

legend_pos<-  theme(legend.position = "bottom")
f_unmatched <- ggplot(subset(data_clean), aes(x =Systolic_Average , fill=tx))  + geom_density( position="identity", alpha=0.5 )  +
  labs(fill="",x="SBP") + legend_pos + ggtitle('(a)') + title_theme
f_smooth_unmatched <- ggplot( subset(data_clean, BMI < 35 & BMI > 13.5 ), aes(y=Systolic_Average, x=BMI, colour=(tx))  )  +
  theme(legend.title=element_blank()) + ylab("SBP") +xlab ("BMI") + geom_smooth()  + legend_pos+ ggtitle('(b)') + title_theme
f_matched <- ggplot(subset(data_clean,!is.na(untreated_matching) ), aes(x =Systolic_Average , fill=untreated_matching))  +
  geom_density( position="identity", alpha=0.5 ) + labs(fill="", x = "SBP")+ legend_pos+ggtitle('(c)') + title_theme
f_smooth_matched <- ggplot( subset(data_clean, BMI < 35 & BMI > 13.5 & !is.na(untreated_matching)), aes(y=Systolic_Average, x=BMI, colour=(untreated_matching))  )  + 
  theme(legend.title=element_blank()) + ylab("SBP") +xlab ("BMI") + geom_smooth()  + legend_pos+ggtitle('(d)') + title_theme
g<-grid.arrange(f_unmatched,   f_smooth_unmatched, f_matched, f_smooth_matched, ncol=2)

ggsave(g,filename = "eFigure4.pdf",width = 8,height=6)

f2_nomeds<- ggplot(subset(data_clean, BMI < 35 & BMI > 13.5 & untreated_matching == "No meds (matching)"  ), aes(BMI, Systolic_Average)) + stat_binhex(bins=50)  + geom_smooth()+  scale_fill_gradient(low="#ffe5e5", high = "#ff0000")+ ggtitle('5d') + title_theme
f2_meds<- ggplot(subset(data_clean, BMI < 35 & BMI > 13.5 & untreated_matching == "On meds"  ), aes(BMI, Systolic_Average)) + stat_binhex(bins=50)  + geom_smooth()+  scale_fill_gradient(low="#ffe5e5", high = "#ff0000")
grid.arrange(f2_nomeds,  f2_meds, ncol=2)

summary(lm(Systolic_Average ~ BMI, data = subset(data_clean, untreated_matching == "On meds")))
summary(lm(Systolic_Average ~ BMI, data = subset(data_clean, untreated_matching == "No meds (matching)")))
summary(lm(Systolic_Average ~ BMI, data = subset(data_clean, tx == "No meds")))



# Figure 7 Projecting --------------------------------------------------------------
setwd("~/bmi-bp-paper/")




predicted_2025s <- read.csv('BMI_draws/predicted_2025s.csv')
predicted_2025s$weights[1:4] <- 0
predicted_2025s$weights[14:15] <- 0
predicted_2025s$weights <- predicted_2025s$weights/sum(predicted_2025s$weights)


predictors_htn_3 <- c("BMI", 'Education', 'Household_Income', 
                    'Hukou', 'Occupation', 'Marital_Status', 'Smoking_Currently', "Hx_Stroke", "Province");
predicted_sbps <- data.frame("Male Predicted 2025 SBP"=rep(0,nrow(predicted_2025s)), "Female Predicted 2025 SBP"=rep(0,nrow(predicted_2025s)));
current_sbps <- data.frame("Male MPP SBP"=rep(0,nrow(predicted_2025s)), "Female MPP SBP"=rep(0,nrow(predicted_2025s)));
rownames(current_sbps) <-rownames(predicted_2025s)
current_bmis_mpp <- data.frame("Male MPP BMI"=rep(0,nrow(predicted_2025s)), "Female MPP BMI"=rep(0,nrow(predicted_2025s)));

rownames(predicted_2025s) <- predicted_2025s$X
rownames(predicted_sbps) <- predicted_2025s$X

bmi_coefs <- data.frame("Male Coef"=rep(0,nrow(predicted_2025s)), "Female Coef"=rep(0,nrow(predicted_2025s)));
mf <- c("male", "female")

IschemicStrokeRR<- rep(0,nrow(predicted_2025s))
IschemicStrokeRR[5:13] <- c(2.05,2.05, 1.83,1.83, 1.63,1.63, 1.44, 1.44, 1.28)

IHDRR<- rep(0,nrow(predicted_2025s))
IHDRR[5:13] <- c(1.68, 1.68, 1.56,1.56, 1.45,1.45, 1.33, 1.33, 1.26)

names(IschemicStrokeRR) <- rownames(current_sbps)
names(IHDRR) <- rownames(IHDRR)

#num_strokes_per_100k<- rep(0,nrow(predicted_2025s))
#num_strokes_per_100k[5:13] <- c(26.4, 139.1, 139.1, 528.5, 528.5, 908.6, 908.6, 1486.9, 1486.9 )
#names(num_strokes_per_100k) <- rownames(current_sbps)




num_census<- data.frame("male"=rep(0,nrow(predicted_2025s)), "female"=rep(0,nrow(predicted_2025s)))
rownames(num_census) <- rownames(current_sbps)
num_census[5:13,1] <- c(60391104,63608678,53776418,40363234,41082938,29834426,20748471,16403453,11278859)
num_census[5:13,2] <- c(57634855,61145286,51818135,38389937,40229536,28832856,20364811,16568944,12573274)

num_strokes<- data.frame("male"=rep(0,nrow(predicted_2025s)), "female"=rep(0,nrow(predicted_2025s)))
num_ihds<- data.frame("male"=rep(0,nrow(predicted_2025s)), "female"=rep(0,nrow(predicted_2025s)))

num_strokes_per_100k<- data.frame("male"=rep(0,nrow(predicted_2025s)), "female"=rep(0,nrow(predicted_2025s)))
num_strokes_per_100k[5:13,1] <- c(26.4, 139.1, 139.1, 528.5, 528.5, 908.6, 908.6, 1486.9, 1486.9 )/1E5
num_strokes_per_100k[5:13,2] <- c(18, 92.7, 92.7, 339.9,339.9, 738.0, 738.0, 1219.7, 1219.7)/1E5
rownames(num_strokes_per_100k) <- rownames(predicted_2025s)

load('/home/George/bmi-bp-paper/BMI_draws/fnfchd.rdata')
num_ihds_rates<- data.frame("male"=rep(0,nrow(predicted_2025s)), "female"=rep(0,nrow(predicted_2025s)))
num_ihds_rates[5:13,1] <- c(fnfchdCHN1$fnfchd[fnfchdCHN1$year==2025][2], fnfchdCHN1$fnfchd[fnfchdCHN1$year==2025][2:9])
num_ihds_rates[5:13,2] <- c(fnfchdCHN2$fnfchd[fnfchdCHN2$year==2025][2],fnfchdCHN2$fnfchd[fnfchdCHN2$year==2025][2:9])

#num_ihds<- data.frame("male"=rep(0,nrow(predicted_2025s)), "female"=rep(0,nrow(predicted_2025s)))
#num_ihds[5:13,1] <- c(506,0,0, 63982,0,0,0, 248536,0)*3
#num_ihds[5:13,2] <- c(253,0,0, 20963,0,0,0, 106465,0)*3
##rownames(num_ihds) <- rownames(predicted_2025s)

PAFs <- data.frame("Male PAF"=rep(0,nrow(predicted_2025s)), "Female PAF"=rep(0,nrow(predicted_2025s)))
IHDPAFs <- data.frame("Male IHD PAF"=rep(0,nrow(predicted_2025s)), "Female IHD PAF"=rep(0,nrow(predicted_2025s)))
rownames(PAFs) <- rownames(predicted_2025s)
rownames(IHDPAFs) <- rownames(predicted_2025s)

for (i in 5:13){
  for (j in 1:2){
   minAge = 35 + (i -5)*5 
   maxAge = 35 + (i-4)*5  
   if (j==1) subgroup_sex = "Male" else subgroup_sex = "Female"
   data_sub_sub <- subset(data_clean, Age < maxAge & Age >= minAge & Sex.x== subgroup_sex)
   outlm_adj <- lm (Systolic_Average ~ ., data = data_sub_sub[,c("Systolic_Average", predictors_htn_3)])
   bmi_coef <- outlm_adj$coefficients["BMI"]
   increase_bmi <- predicted_2025s[i,mf[j]]- current_bmis[i,j]
   current_bmis_mpp[i,j] <- mean(data_sub_sub$BMI)
   predicted_sbps[i,j] <- mean(data_sub_sub$Systolic_Average)  + increase_bmi*bmi_coef
   current_sbps[i,j] <- mean(data_sub_sub$Systolic_Average)  
   bmi_coefs[i,j] <- bmi_coef
   
   #Compute PAF
   increase_sbp <- predicted_sbps[i,j] - current_sbps[i,j] 
   sbp_sd_subgroup <- sd(data_sub_sub$Systolic_Average)
   
   PAFs[i,j]<- PAFpnorm(ExpMn=predicted_sbps[i,j], ExpSD=sbp_sd_subgroup, LnRR=log(IschemicStrokeRR[i])/10, CFmean=current_sbps[i,j], CFsd=sbp_sd_subgroup, Min=70, Max=220, By=0.1)
   IHDPAFs[i,j]<- PAFpnorm(ExpMn=predicted_sbps[i,j], ExpSD=sbp_sd_subgroup, LnRR=log(IHDRR[i])/10, CFmean=current_sbps[i,j], CFsd=sbp_sd_subgroup, Min=70, Max=220, By=0.1)
   
   num_strokes[i,j] <- round(PAFs[i,j]*num_strokes_per_100k[i,j]*num_census[i,j])
   num_ihds[i,j] <- round(IHDPAFs[i,j]*num_ihds_rates[i,j]*num_census[i,j])
   
   print(sprintf("%d,%d, %s, Predicted SBP: %f, Current SBP: %f, sd: %f, RR: %f, PAF: %f, Num Strokes: %f", minAge,maxAge, mf[j], predicted_sbps[i,j], current_sbps[i,j],sbp_sd_subgroup, IschemicStrokeRR[i],PAFs[i,j], num_strokes[i,j]))
  }
}

#IHDPAFs_aves <- IHDPAFs;
#IHDPAFs_aves[5:7,] <- data.frame(t(matrix(rep(colMeans(IHDPAFs[5:7,]),3),ncol=3)))
#IHDPAFs_aves[8:11,] <- data.frame(t(matrix(rep(colMeans(IHDPAFs[8:11,]),4),ncol=4)))
#IHDPAFs_aves[12:13,] <- data.frame(t(matrix(rep(colMeans(IHDPAFs[12:13,]),2),ncol=2)))

paf_result_male <- round(cbind("BMI Increase"=predicted_2025s[,"male"] - current_bmis[,1] ,"Coef"=bmi_coefs[,1],"SBP Increase"=predicted_sbps[,1] - current_sbps[,1],"Ischemic Stroke RR"=IschemicStrokeRR, "PAF"=PAFs[,1], "Total_Strokes"=round(num_strokes_per_100k[,1]*num_census[,1]), "Attributable_Strokes"=num_strokes[,1] ),2)
rownames(paf_result_male) <- rownames(predicted_2025s)
#write.csv(paf_result_male[5:13,],"BMI_draws/paf_result_male.csv")

paf_result_female <- round(cbind("BMI Increase"=predicted_2025s[,"female"] - current_bmis[,2] ,"Coef"=bmi_coefs[,2],"SBP Increase"=predicted_sbps[,2] - current_sbps[,2],"Ischemic Stroke RR"=IschemicStrokeRR, "PAF"=PAFs[,2], "Total_Strokes"=round(num_strokes_per_100k[,2]*num_census[,2]), "Attributable_Strokes"=num_strokes[,2] ),2)
rownames(paf_result_female) <- rownames(predicted_2025s)
#write.csv(paf_result_female[5:13,],"BMI_draws/paf_result_female.csv")
#write.csv(paf_result_male[5:13,],"BMI_draws/paf_result_male.csv")

paf_ihd_result_male <- round(cbind("BMI Increase"=predicted_2025s[,"male"] - current_bmis[,1] ,"Coef"=bmi_coefs[,1],"SBP Increase"=predicted_sbps[,1] - current_sbps[,1],"IHD RR"=IHDRR, "PAF"=IHDPAFs[,1], 
                                          "Total_IHD"=round(num_ihds_rates[,1]*num_census[,1]), "Attributable_IHD"=round(num_ihds[,1])) ,2)
rownames(paf_ihd_result_male) <- rownames(predicted_2025s)

paf_ihd_result_female <- round(cbind("BMI Increase"=predicted_2025s[,"female"] - current_bmis[,2] ,"Coef"=bmi_coefs[,2],"SBP Increase"=predicted_sbps[,2] - current_sbps[,2],"IHD RR"=IHDRR, "PAF"=IHDPAFs[,2], 
                                          "Total_IHD"=round(num_ihds_rates[,2]*num_census[,2]), "Attributable_IHD"=round(num_ihds[,2])) ,2)
rownames(paf_ihd_result_female) <- rownames(predicted_2025s)

write.csv(paf_ihd_result_female[5:13,],"BMI_draws/paf_ihd_result_female.csv")
write.csv(paf_ihd_result_male[5:13,],"BMI_draws/paf_ihd_result_male.csv")

colSums(paf_ihd_result_female,na.rm = TRUE) 
colSums(paf_ihd_result_male,na.rm=TRUE)

sum(paf_result_male[,"Attributable_Strokes"])
sum(paf_result_female[,"Attributable_Strokes"])



#View(round(cbind(current_bmis, predicted_2025s[,c("male","female")], current_sbps, bmi_coefs, predicted_sbps, IschemicStrokeRR, "PAF.Male"=PAFs[,1],"PAF.Female"=PAFs[,2] ),2))
#Check it
data_sub_sub <- subset(data_clean, Age < 45 & Age >= 40 & Sex.x== "Female")
predicted_test = mean(data_sub_sub$Systolic_Average) + (predicted_2025s$female[predicted_2025s$X == "40-44"]-current_bmis[predicted_2025s$X == "40-44",2])*lm(Systolic_Average ~ ., data = data_sub_sub[,c("Systolic_Average", predictors_htn_3)])$coefficients["BMI"]
PAFpnorm(ExpMn=predicted_test, ExpSD=sd(data_sub_sub$Systolic_Average), LnRR=log(1.68)/10, CFmean=mean(data_sub_sub$Systolic_Average), CFsd=sd(data_sub_sub$Systolic_Average), Min=70, Max=220, By=0.1)
PAFpnorm(ExpMn=predicted_test, ExpSD=sd(data_sub_sub$Systolic_Average), LnRR=log(2.05)/10, CFmean=mean(data_sub_sub$Systolic_Average), CFsd=sd(data_sub_sub$Systolic_Average), Min=70, Max=220, By=0.1)

 
#cbind("Male 2014"=current_bmis[5:13,1], "Female 2014"=current_bmis[5:13,2], "Male 2025"=predicted_2025s[5:13,1],
      #"Female 2025"=predicted_2025s[5:13,2], "Delta Male" = predicted_2025s[5:13,1] - current_bmis[5:13,1],"Delta Female"=predicted_2025s[5:13,2] - current_bmis[5:13,2] )
 
#current_bmi_male <- sum(current_bmis[5:13,1]*adj_weights)
#predicted_bmi_male <- sum(predicted_2025s[5:13,"male"]*adj_weights)
#current_bmi_female <- sum(current_bmis[5:13,2]*adj_weights)
#predicted_bmi_female <- sum(predicted_2025s[5:13,"female"]*adj_weights)
# 
# current_sbp_male <- sum(current_sbps[5:13,1]*adj_weights)
# current_sbp_female <- sum(current_sbps[5:13,2]*adj_weights)
# 
# predicted_sbp_2025_male <- sum(adj_weights*predicted_sbps$Male[5:13])
# predicted_sbp_2025_female <- sum(adj_weights*predicted_sbps$Female[5:13])
# 
# diff_sbps <- predicted_sbps - current_sbps
# 
# predicted_sbp_2025_male - current_sbp_male
# predicted_sbp_2025_female - current_sbp_female
# sprintf("Current SBP Male: %f, 2025 SBP Male: %f, Current SBP Female: %f, 2025 SBP Female: %f", current_sbp_male, predicted_sbp_2025_male, current_sbp_female, predicted_sbp_2025_female)
# 
# 
# sbp_sd_male <- sd(data_clean$Systolic_Average[data_clean$Sex.x=="Male"])
# sbp_sd_female <- sd(data_clean$Systolic_Average[data_clean$Sex.x=="Female"])
# 
# #toplot_current_male <- rnorm(1E4, mean=current_sbp_male, sd = sbp_sd_male)
# #toplot_2025_male <- rnorm(1E4, mean=predicted_sbp_2025_male, sd = sbp_sd_male)
# #toplot_male <- data.frame(label=as.factor(c(rep("2009", 1E4), rep("2025", 1E4))), sbp=c(toplot_current_male, toplot_2025_male) )
# #pnorm(140, mean=current_sbp_male, sd = sbp_sd_male,lower.tail = FALSE)
# #pnorm(140, mean=predicted_sbp_2025_male, sd = sbp_sd_male,lower.tail = FALSE)
# 
# #pnorm(140, mean=current_sbp_female, sd = sbp_sd_female,lower.tail = FALSE)
# #pnorm(140, mean=predicted_sbp_2025_female, sd = sbp_sd_female,lower.tail = FALSE)
# 
# #quantile(x = data_clean$Systolic_Average[data_clean$Sex.x=="Male"],0.6)
# #f8_a <- ggplot(toplot_male, aes(x=sbp, color=label)) + geom_density()

# Now, let's do some PAF --------------------------------------------------------
names(current_sbps[5:13,1])
predicted_2025s$weights[5:13]

mf <- c("male", "female")
for (i in 0:4){
  for (j in 1:2){
   minAge = 35 + (i)*10 
   maxAge = 35 + (i+1)*10 
   data_sub_sub <- subset(data_clean, Age < maxAge & Age >= minAge & Sex.x== subgroup_sex)
   print(c(minAge,maxAge, nrow(data_sub_sub)))
   outlm_adj <- lm (Systolic_Average ~ ., data = data_sub_sub[,c("Systolic_Average", predictors_htn_3)])
   bmi_coef <- outlm_adj$coefficients["BMI"]
   increase_bmi <- predicted_2025s[i+5,mf[j]]- current_bmis[i+5,j]
   print(sprintf("%s, %d,%d: %f",mf[i], minAge,maxAge, increase_bmi))
  }
}

names(IschemicStrokeRR) <- c("35-44","45-55", "55-64","65-74", "75-84")



# Exploring Tibet --------------------------------------------------------
ggplot( subset(toplot, !is.na(Province) & tx=="No meds") , aes(x=coeflm.BMI.by_.,fill=(factor(Province=="Tibet" )) )) +  geom_density()

#allcoef <-lm_subgroup(data_clean$Systolic_Average, data_clean$BMI) 
ggplot( subset(toplot, !is.na(tx)), aes(x=x,y=coeflm.BMI.by_.,ymin=lowerlm, ymax=upperlm, color=(tx) )) +  geom_point(size=0.9)
ggplot( subset(toplot, !is.na(Province) & tx=="No meds") , aes(x=x,y=coeflm.BMI.by_.,ymin=lowerlm, ymax=upperlm, color=(factor(Province=="Tibet" )) )) +  geom_point(size=0.7)

subset(subgroup_results,tx=="No meds" & coeflm.BMI.by_. <0.5)[,c("coeflm.BMI.by_.", "group")]

tibet_nomeds <- subset(subgroup_results,tx=="No meds" & Province=="Tibet")
tibet_nomeds[order(-tibet_nomeds$coeflm.BMI.by_.),c("group","N", "coeflm.BMI.by_." )]

table(subset(data_clean, Occupation =="Farmer" & Household_Income == "<5k RMB" & Province=="Tibet" & tx == "No meds")[,"Sex.x"])

ggplot(subset(data_clean,!is.na(tx) & !is.na(Province) & tx=="No meds" & Province == "Tibet" ), aes(x =Systolic_Average , fill=Household_Income))  + geom_density( position="identity", alpha=0.5 ) + labs(fill="", x = "SBP")+ legend_pos

ggplot( subset(data_clean, BMI < 35 & BMI > 13.5 & !is.na(Province) & tx=="No meds" & Province == "Tibet"  ), aes(y=Systolic_Average, x=BMI)  )  + theme(legend.title=element_blank()) + ylab("SBP") +xlab ("BMI") + stat_bin_hex()  
ggplot( subset(data_clean, BMI < 35 & BMI > 13.5 & !is.na(Province) & tx=="No meds" & Province == "Tibet"  ), aes(y=Systolic_Average, x=BMI)  )  + theme(legend.title=element_blank()) + ylab("SBP") +xlab ("BMI") + stat_bin_hex()  

ggplot( subset(data_clean, BMI < 35 & BMI > 13.5 &tx=="No meds" & Province == "Tibet" ), aes(y=Systolic_Average, x=BMI, colour=(Household_Income))  )  + theme(legend.title=element_blank()) + ylab("SBP") +xlab ("BMI") + geom_smooth()  + legend_pos
ggplot( subset(data_clean, BMI < 35 & BMI > 13.5 &tx=="No meds" ), aes(y=Systolic_Average, x=BMI, colour=(Household_Income))  )  + theme(legend.title=element_blank()) + ylab("SBP") +xlab ("BMI") + geom_smooth()  + legend_pos


ggplot( subset(data_clean, BMI < 35 & BMI > 13.5 &tx=="No meds" ), aes(y=Systolic_Average, x=BMI, colour=factor(Province == "Tibet"))  )  + theme(legend.title=element_blank()) + ylab("SBP") +xlab ("BMI") + geom_smooth()  + legend_pos

ggplot(subset(data_clean, BMI < 35 & BMI > 13.5 & tx=="No meds"), aes(x =BMI , fill=factor(Province == "Tibet")))  + geom_density( position="identity", alpha=0.5 )
ggplot(subset(data_clean, BMI < 35 & BMI > 13.5 & tx=="No meds"), aes(x =Systolic_Average , fill=factor(Province == "Tibet")))  + geom_density( position="identity", alpha=0.5 )

# Subgroup differences ----------------------------------------------------

# require(ggplot2)
# ggplot( toplot, aes(x=x,y=coeflm.BMI.by_.,ymin=lowerlm, ymax=upperlm ) ) +  geom_point() 
# #ggplot( toplot, aes(x=x,y=coeflm.BMI.by_.,ymin=lowerlm, ymax=upperlm, color=(untreated_matching) )) +  geom_point()
# 
# #ggplot( subset(toplot, is, aes(x=x,y=coeflm.BMI.by_.,ymin=lowerlm, ymax=upperlm, color=(tx) )) +  geom_point()
# ggplot( toplot, aes(x=x,y=coeflm.BMI.by_.,ymin=lowerlm, ymax=upperlm, color=N) ) +  geom_point() + scale_colour_gradientn(colours=rainbow(3))
# ggplot( toplot, aes(x=x,y=coeflm.BMI.by_.,ymin=lowerlm, ymax=upperlm, color=N) ) +  geom_point() + scale_colour_gradientn(colours=rainbow(3))+geom_errorbar()
# 
# 
# write.csv(subgroup_results[order(subgroup_results$coeflm),c("group", "N", "coeflm.BMI", "lowerlm", "upperlm")], file="k_4_subgroups.csv")
# 
# 
# ggplot(subset(toplot,!is.na(toplot$tx)), aes(x =coeflm.BMI.by_. , fill=factor(tx)))  + geom_density( position="identity", alpha=0.5 ) 
# ggplot(subset(toplot,!is.na(toplot$Province) & tx == "No meds"), aes(x =coeflm.BMI.by_. , fill=factor(Province == "Tibet")))  + geom_density( position="identity", alpha=0.5 ) 
# ggplot(subset(data_clean,!is.na(Province) ), aes(x =Systolic_Average , fill=factor(Province == "Tibet")))  + geom_density( position="identity", alpha=0.5 ) 
# ggplot(subset(data_clean,!is.na(Province) ), aes(x =BMI , fill=factor(Province == "Tibet")))  + geom_density( position="identity", alpha=0.5 ) 
# 
# 
# # Subgroup differences ----------------------------------------------------
# 
# difffeats <- c("Sex.x", "Age_cat",  "Ethnicity", "Hukou", "Occupation", "Education", "Household_Income",  "Marital_Status", "Province")
# #agg_by <- as.list(meddifs[,c()])
# agg_by <- as.list(meddifs[,c("Sex.x", "Age_cat", "Ethnicity", "Hukou")])
# meddifs <- subset(toplot,!is.na(tx))
# #meddifs_ <- aggregate(meddifs,by = agg_by, FUN = function(x) x)
# 
# library(data.table)
# setDT(meddifs)
# gcldiff <- function(coef,tx,by_I){
#   print(meddifs$tx[by_I])
#   print(coef)
#   print(tx)
#   print(by_I)
#   
#   if (length(by_I)==2) {
#         
#    coef 
#   }else{
#     NA;
#   }
# }
# meddifs_ <- data.frame(meddifs[,  as.list(gcldiff(coeflm.BMI.by_.,tx, .I)),by = difffeats])
# 
# # All combos l2 -----------------------------------------------------------
# 
# 
# lm_l2_subgroup <- function(SA, BMI,by_){
#   byvec <- rep(FALSE,length(SA))
#   byvec[by_] <- TRUE;
#   outlm <- lm (SA[by_] ~ BMI[by_])
#   lower <- confint.default(outlm)[2,1] 
#   upper <- confint.default(outlm)[2,2]
#   #l2val <- l2_dist_test(SA,BMI,byvec,50)
#   l2val <- l2_dist(SA,BMI,byvec)
#   c("coeflm"= outlm$coefficients[2], "lowerlm"= lower, "upperlm" =upper, "N"=length(by_), "l2val"=l2val);
# }
# 
# subgroup_feats2 <- c("Sex.x",  "Ethnicity", "Hukou", "Occupation",  "Household_Income",  "Province", "tx")
# data_sub2 <- cbind(data_sub, rand_test=factor(sample(1:3,size = nrow(data_sub),replace = TRUE)))
# ks <- 1:4;
# thresh <- 5E3
# subgroup_results <- data.frame();
# for (k in ks) {
#   subsub_feats <- combn(subgroup_feats2, k)
#   library(data.table)
#   data_subgroups <- data_sub2[,c("Systolic_Average", "BMI", subgroup_feats2)]
#   setDT(data_subgroups)
#   for (i in 1:ncol(subsub_feats)){
#     print(sprintf("%d/%d completed", i, ncol(subsub_feats)))
#     subgroup_results_i <- data.frame(data_subgroups[,  as.list(lm_l2_subgroup(data_sub2$Systolic_Average, data_sub$BMI,.I)),by = eval(subsub_feats[,i])])
#     if (k > 1){
#       subgroup_name <- apply( subgroup_results_i[,1:k], 1, paste, collapse=",")
#     }else{
#       subgroup_name <- subgroup_results_i[,1]
#     }
#     d =subset(rbind.fill(data.frame(subgroup_results_i, group=subgroup_name)), N>thresh)
#     subgroup_results<- rbind.fill(subgroup_results,d)
#   }
#   
# }
# print(subgroup_results[order(-subgroup_results$l2val),])
# nrow(subgroup_results)
# 
# 
# # L2 dists ----------------------------------------------------------------
# 
# print(subgroup_results[order(-subgroup_results$l2val)[1:100],c("group","N","coeflm.BMI.by_.",   "l2val")])
# subgroup <- data_sub$Ethnicity == "Han" & data_sub$Province == "Liaoning" & data_sub$Occupation == "Other"
# ft1<- ggplot(data_sub[subgroup,], aes(BMI, Systolic_Average)) + stat_binhex()  + geom_smooth()+  scale_fill_gradient(low="#ffe5e5", high = "#ff0000")
# ft2<- ggplot(data_sub[!subgroup,], aes(BMI, Systolic_Average)) + stat_binhex() + geom_smooth()+  scale_fill_gradient(low="#ffe5e5", high = "#ff0000")
# ft12 <- grid.arrange(ft1,ft2)
# 
# subgroup <-  data_sub$Province == "Liaoning" & data_sub$Occupation == "Other"
# ggplot(data_sub[subgroup,], aes(BMI, Systolic_Average)) + geom_point()   
# 
# subgroup <- data_sub$Sex.x == "Male" & data_sub$Province == "Henan"  
# subgroup <- data_sub$Ethnicity == "Han" & data_sub$Household_Income == "<5k RMB" & data_sub$Province == "InnerMongolia"  
# subgroup <- data_sub$Ethnicity == "Han" & data_sub$Hukou == "Unified" & data_sub$Province == "Henan"  
# subgroup <- data_sub$tx == "No meds" & data_sub$Hukou == "Urban" & data_sub$Province == "Shandong"  
# ggplot(data_sub[subgroup,], aes(BMI, Systolic_Average)) + geom_point()   
# ggplot(data_sub, aes(BMI, Systolic_Average, colour = factor(subgroup))) + geom_smooth()   

