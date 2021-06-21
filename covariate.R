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
                                          IS_SBP, IS_DBP,
                                         IS_smk,  IS_BMI,
                                         tx, IS_Ethnicity
))
names(data) <- c("Patient_Id", "Sex.x", "Age", "Household_Income", "Education", "Occupation", "Hukou",
                 "Province", "Marital_Status", 
                 "Systolic_Average", "Diastolic_Average",  
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

data_clean$tx <- factor(data_clean$tx)

levels(data_clean$Sex.x) <- c("Male", "Female")
#levels(data_clean$Occupation) <- c("Farmer", "Retire", "Unemployment", "Housework", "Unknown", "Refuse", "Worker", "Manager", "Admin", "Technician", "Serv. Ind.", "Bus. Owner", "Military", "Other")
levels(data_clean$Occupation) <- c("Farmer", "Workers", "Administrators", "Admin. Clerk", "Technician", "Business", "Bus. Owner", "Military", "Others", "Retire", "Unemployed", "Housework", "Unknown", "Refuse")
levels(data_clean$Marital_Status) <- c("Unmarried", "Married")
levels(data_clean$Hukou) <- c("Rural", "Urban", "Unified","No Hukou")

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



covariates <- c("Sex.x", "Age_cat2", "Household_Income","Occupation", "Ethnicity", "Marital_Status",  "Hukou","Province","Education","Smoking_Currently","tx")


plot_titles <- c(`Sex.x`="Sex",  `Hukou`="Hukou", `Occupation`="Occupation", `Education`="Education", `Household_Income`="Household Income",  `Marital_Status`="Marital Status", `Province`="Province" ,`tx`="BP Meds" , `Ethnicity`="Ethnicity", `Age_cat`="Age", `Age_cat2`="Age", `Smoking_Currently`="Currently Smoking") 



predictors_htn <- c("BMI", "Age",  'Sex.x', 'Education', 'Ethnicity','Household_Income', 'Hukou', 'Occupation', 'Marital_Status', 'Smoking_Currently', "Province","tx");


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



