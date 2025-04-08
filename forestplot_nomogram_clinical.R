library(survivalROC)
library(timeROC)
library(survival)
library(glmnet)
library("survminer")
library(ggpubr)
library(dplyr)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(survivalsvm)
library(tibble)
library(BART)
library(linkET)
library(IOBR)

############################## This is the next step of survival analysis
#### Input the data

risk_scores<-readRDS('risk_score_list_for_clinical.rds')
tcga_info<-readRDS('tcga_tumor_info.rds')
#gse37_info<-readRDS('gse84437_sampleInfo.rds')
#gse33_info<-readRDS('gse84433_sampleInfo.rds')
gse59_info<-readRDS('gse15459_sampleInfo.rds')

gse26_info<-readRDS('gse84426_sampleInfo.rds')

#### subset tcga_info
tcga_info<-tcga_info[rownames(risk_scores$tcga),]

### update risk_score
# Assuming the dataframes are stored in a list called 'risk_scores'
# and the individual '_info' dataframes are in separate variables.

# List of info dataframes
info_df <- list(
  tcga_info = tcga_info,
 # gse37_info = gse37_info,
  #gse33_info = gse33_info,
  gse59_info = gse59_info,
  gse26_info = gse26_info
)

# Get the names to map them correctly
risk_names <- names(risk_scores)

# Merge each dataframe in 'risk_scores' with the corresponding '_info' dataframe by row names
for (name in risk_names) {
  # Perform the merge using the rownames as keys
  merged_df <- merge(
    risk_scores[[name]], 
    info_df[[paste0(name, "_info")]], 
    by = "row.names", 
    all.x = TRUE
  )
  
  # Restore the row names and drop the extra 'Row.names' column
  rownames(merged_df) <- merged_df$Row.names
  merged_df$Row.names <- NULL
  
  # Store the merged dataframe back in the list
  risk_scores[[name]] <- merged_df
}


###### Forest plot
data<-risk_scores$tcga
data$age_new<-ifelse(data$age<=60,"<=60",">60")
data$stage<-ifelse(data$stage %in% c('Stage I','Stage II'),'Stage I-II','Stage III-IV')
data$T_stage_new <- ifelse(data$T_stage %in% c('TX','T1','T2'), "Low", "High")
data$M_stage_new <- ifelse(data$M_stage %in% c('MX','M0'), "Low", "High")
data$N_stage_new <- ifelse(data$T_stage %in% c('NX','N0','N1'), "Low", "High")

# Convert categorical variables to numeric (ordinal encoding)
data$stage_numeric <- as.numeric(as.factor(data$stage))
data$T_stage_numeric <- as.numeric(as.factor(data$T_stage))
data$M_stage_numeric <- as.numeric(as.factor(data$M_stage))
data$N_stage_numeric <- as.numeric(as.factor(data$N_stage))
#### Adjust factor

data$age_new<-factor(data$age_new,levels = c("<=60",">60"))

data$Risk<-factor(data$Risk,levels = c("low risk", "high risk"))
data$stage<-factor(data$stage,levels = c('Stage I-II','Stage III-IV'))
data$T_stage<-factor(data$T_stage,levels = c('TX','T1','T2','T3','T4'))
data$M_stage<-factor(data$M_stage,levels = c('MX','M0','M1'))
data$N_stage<-factor(data$N_stage,levels = c('NX','N0','N1','N2','N3'))
data$sex<-factor(data$sex,levels = c("male", "female"))


### Univariate forest
fit_age <- coxph(Surv(time, status) ~ age, data = data)
fit_sex <- coxph(Surv(time, status) ~ sex, data = data)
fit_tumor_stage <- coxph(Surv(time, status) ~ stage, data = data)
fit_t_stage <- coxph(Surv(time, status) ~ T_stage, data = data)
fit_m_stage <- coxph(Surv(time, status) ~ M_stage, data = data)
fit_n_stage <- coxph(Surv(time, status) ~ N_stage, data = data)
fit_risk <- coxph(Surv(time, status) ~ Risk, data = data)

# List of all the fitted Cox models
fit_list <- list(
  age = fit_age,
  sex = fit_sex,
  tumor_stage = fit_tumor_stage,
  t_stage = fit_t_stage,
  m_stage = fit_m_stage,
  n_stage = fit_n_stage,
  risk = fit_risk
)

# Function to extract HR, 95% CI, and p-value from a Cox model
# Define the function to extract relevant Cox model information
extract_cox_info <- function(fit) {
  # Extract summary from the Cox model
  summary_fit <- summary(fit)
  
  # Extract only the first coefficient row (in case of multiple coefficients)
  coef_info <- summary_fit$coef[1, ]  
  ci_info <- summary_fit$conf.int[1, ]
  
  # Extract relevant values
  HR <- as.numeric(ci_info["exp(coef)"])            # Hazard Ratio
  lower_95CI <- as.numeric(ci_info["lower .95"])    # Lower 95% CI
  upper_95CI <- as.numeric(ci_info["upper .95"])    # Upper 95% CI
  p_value <- coef_info["Pr(>|z|)"]                  # P-value
  
  # Format p-value: "<0.0001" for very small values, otherwise 4 decimals
  p_value <- ifelse(p_value < 0.0001, "<0.0001", format(round(p_value, 4), nsmall = 4))
  
  # Return as a named vector
  return(c(HR = HR, lower95CI = lower_95CI, upper95CI = upper_95CI, pvalue = p_value))
}

# Example: List of Cox models to process
# fit_list <- list(fit_age, fit_sex, fit_t_stage, fit_n_stage, fit_risk)

# Apply the function to all Cox models and bind the results into a dataframe
result_df <- do.call(rbind, lapply(fit_list, extract_cox_info))

# Convert the result to a dataframe and ensure numeric columns
result_df <- data.frame(result_df, stringsAsFactors = FALSE)

# Ensure numeric conversion for HR and CI columns
result_df$HR <- as.numeric(result_df$HR)
result_df$lower95CI <- as.numeric(result_df$lower95CI)
result_df$upper95CI <- as.numeric(result_df$upper95CI)


# Ensure column names match the desired output
colnames(result_df) <- c( "HR", "lower95%CI", "upper95%CI", "pvalue")





library(forplo)
library(autoReg)


forplo(result_df[,1:3],
       add.columns = result_df[4],
       add.colnames = 'p-value',
       col= 'red',
       size = 2,
       ci.edge= T,
       char = 18,
       em = 'HR',
       shade.every=1,
       shade.col='gray40',
       font='Helvetica',
       leftbar.ticks=TRUE,
       right.bar=TRUE,
       rightbar.ticks=TRUE
)



###Multivariate forest plot


fit_multi_tcga<-coxph(Surv(time, status) ~ age + sex + stage+ T_stage_numeric + N_stage_numeric + M_stage_numeric+ Risk, data = data)
summary(fit_multi_tcga)

forplo(fit_multi_tcga,
       row.labels=c('Age',
                    'sex',
                    'stage',
                    'T stage',
                    'N stage',
                    'M stage',
                    'Risk'),
       col= 'red',
       size = 2,
       ci.edge= T,
       char = 18,  
       shade.every=1,
       
       shade.col='gray40',
       font='Helvetica',
       #title='Cox regression, sorted by HR, inverted HRs<1, scaled dots by effect size',
       right.bar=TRUE,
       rightbar.ticks=TRUE,
       leftbar.ticks=TRUE,
       
)

########### nomogram
library(rms)
library(nomogramEx)
library(regplot)


dd<-datadist(data)
options(datadist="dd")
options(na.action="na.delete")
summary(pbc$time)

# Adjust the labels for the variables
label(data$T_stage_numeric) <- "T stage"
label(data$N_stage_numeric) <- "N stage"
label(data$M_stage_numeric) <- "M stage"

coxpbc<-cph(formula = (Surv(time, status) ~ age + sex  +stage+T_stage + 
                         N_stage+M_stage + Risk), 
            data = data,x=T,y=T,surv = T,na.action=na.delete)  #,time.inc =2920
print(coxpbc)
surv<-Survival(coxpbc) 
surv1<-function(x) surv(12,x)
surv3<-function(x) surv(35,x)
surv5<-function(x) surv(60,x)

x<-nomogram(coxpbc,fun = list(surv1,surv3,surv5),lp=T,
            funlabel = c('1-year survival Probability', '3-year survival Probability','5-year survival Probability'),
            maxscale = 100,fun.at = c(0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))
pdf("nomogram_classical.pdf",width = 12,height = 10)
plot(x, lplabel="Linear Predictor",
     xfrac=.35,varname.label=TRUE, varname.label.sep="=", ia.space=.2, 
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, 
     cap.labels=FALSE,cex.var = 1.6,cex.axis = 1.05,lwd=5,
     label.every = 1,col.grid = gray(c(0.8, 0.95)))
dev.off()




##### Cach 2
library(nomogramEx)
mycol<-c("#A6CEE3","#1F78B4","#33adff","#2166AC")
names(mycol) = c("dencol","boxcocl","obscol","spkcol")
mycol<- as.list(mycol)
pbccox <- coxph(formula = Surv(time,status) ~  age + sex + +stage+T_stage + 
                  N_stage+M_stage + Risk , data = data)


regplot(pbccox,
        
        #observation=data[2,], #
        failtime = c(12,36,60), 
        #observation = data[20,],
        prfail = TRUE, 
        showP = T, 
        droplines = T,
        center = T,
        odds = F,
        #colors = mycol, 
        rank="sd",
        points = T,
        leftlabel=T,
        title = 'Nomogram',
        interval="confidence") #

#Draw a calibration curve for verification
f1<-cph(formula = Surv(time,status) ~  age + sex  +stage+T_stage_numeric + 
          N_stage_numeric+M_stage_numeric + Risk,data=data,x=T,y=T,surv = T,
        na.action=na.delete,time.inc = 12) 
f3<-cph(formula = Surv(time,status) ~  age + sex  +stage+T_stage_numeric + 
          N_stage_numeric+M_stage_numeric + Risk,data=data,x=T,y=T,surv = T,
        na.action=na.delete,time.inc = 36) 
f5<-cph(formula = Surv(time,status) ~  age + sex  +stage+T_stage_numeric + 
          N_stage_numeric+M_stage_numeric + Risk,data=data,x=T,y=T,surv = T,
        na.action=na.delete,time.inc = 60) 
#Parameter m=50 means repeated calculation of 50 samples in each group
cal1<-calibrate(f1, cmethod="KM", method="boot",u=12,m=50,B=1000) 
cal3<-calibrate(f3, cmethod="KM", method="boot",u=36,m=50,B=1000) 
cal5<-calibrate(f5, cmethod="KM", method="boot",u=60,m=50,B=1000) 


pdf("calibration_compare.pdf",width = 5,height = 5)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),subtitles = F,
     bty = "l",
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)


plot(cal3,lwd = 2,lty = 0,errbar.col = c("#B2182B"),subtitles = F,
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

plot(cal5,lwd = 2,lty = 0,errbar.col = c("orange"),subtitles = F,
     xlim = c(0,1),ylim= c(0,1),col = c("orange"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("orange"), pch = 16)


legend("bottomright", 
       legend = c("1-year", "3-year","5-year"), 
       col =c("#2166AC","#B2182B","orange"), 
       lwd = 2,
       cex = 1.5,
       bty = "n")
dev.off()

####### box plot compare risk group

## calculate signature score first (apply in TCGA-STAD)

sig_tme<-calculate_sig_score(pdata           = NULL,
                             eset            = mat_tumor,
                             signature       = signature_collection,
                             method          = "pca",
                             mini_gene_count = 2)

sig_tme <- t(column_to_rownames(sig_tme, var = "ID"))

###combine with sampleInfo

input_iobr <- combine_pd_eset(eset = sig_tme, pdata = risk_scores$tcga, scale = T)

res<- batch_surv(pdata    = input_iobr,
                 time     = "time", 
                 status   = "status", 
                 variable = colnames(input_iobr)[19:ncol(input_iobr)])


p2 <- sig_heatmap(input         = input_iobr, 
                  features      = res$ID[1:50],
                  group         = "Risk", 
                  palette_group = "jama", 
                  palette       = 6,
                  path          = "result" )


p3 <- sig_box(data           = input_iobr, 
              signature      = "TGFβ_myCAF",
              variable       = "Risk",
              jitter         = FALSE,
              cols           = NULL,
              palette        = "jama",
              
              show_pvalue    = T,
              angle_x_text   = 60, 
              hjust          = 1, 
              size_of_pvalue = 5, 
              size_of_font   = 8)
p3

##### get corr
p4<- get_cor(eset = sig_tme, id_eset=rownames(sig_tme), pdata = risk_scores$tcga, is.matrix = TRUE, var1 = "Fatty_Acid_Elongation", 
             var2 = "TGFβ_myCAF", subtype = "Risk", palette = "aaas")


##### 


###METHOD 1: CIBERSORT

cibersort<-deconvo_tme(eset = mat_tumor, method = "cibersort", arrays = FALSE, perm = 200 )
res_meta<-iobr_cor_plot(pdata_group           = input_iobr,
                        id1                   = "ID",
                        feature_data          = cibersort,
                        id2                   = "ID",
                        target                = NULL,
                        group                 = "Risk",
                        is_target_continuous  = F,
                        padj_cutoff           = 1,
                        index                 = 1,
                        category              = "signature",
                        signature_group       = (sig_group)[c(20)],
                        ProjectID             = "GBM-MCP",
                        
                        palette_box           = "jco",
                        palette_corplot       = "pheatmap",
                        palette_heatmap       = 2,
                        feature_limit         = 200,
                        character_limit       = 50,
                        show_heatmap_col_name = FALSE,
                        show_col              = FALSE,
                        show_plot             = TRUE)

###### box plot between risk group and clinical factor
data2<-select(data,RiskScore,subtype,KRAS,MSI,CIMP,TP53,sex,age_new,stage,T_stage,N_stage,M_stage)


#my_comparisons<-list( c("GASTRIC-CIMP", "GASTRIC-EBV-CIMP"), c("GASTRIC-EBV-CIMP", "OTHER"), c("GASTRIC-CIMP", "OTHER") )
# my_comparisons<-list( c("CIN", "EBV"), c("CIN", "GS"),c("CIN", "MSI"),
#                       c("EBV", "GS"), c("EBV", "MSI"),c("GS", "MSI"))

my_comparisons<-list( c("0", "1") )

data3<-data2[complete.cases(data2$TP53),]

data3<-data3[(data3$TP53 !='NA'),]

ggboxplot(data3, x = "TP53", y = "RiskScore",
          color = "TP53", palette =c("#00AFBB", "#E7B800", "#FC4E07","green4"),
          add = "jitter", shape = "TP53")+ stat_compare_means(comparisons = my_comparisons)+ theme_pubr()+
  
  theme(axis.text.x = element_text(size = 15,face = 'bold'),
        axis.text.y = element_text(size=14),
     
        strip.text = element_text(size = 10,face = 'bold'), # Enlarges the gene names in the facet labels
        legend.title = element_text(size = 10), # Enlarges the legend title
        legend.text = element_text(size = 15)) 
