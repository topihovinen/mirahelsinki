## This is an R code for statistical analysis and generating figures of MIRA Helsinki study ##
# ------------------------------------------------------------------------------------------ #

# -------------------------------0. PREFACE------------------------------------------------- #
{
  # Introduce needed libraries for R
  {
    library("ggplot2")
    library("digest")
    library("Hmisc")
    library("devtools")
    library("cairoDevice")
    library("dplyr")
    library("ggsignif")
    library("ggdendro")
    library("mvtnorm")
    library("dendextend")
    library("arm")
    library("svglite")
    library("ggsci")
    library("png")
    library("gridExtra")
    library("grid")
    library("xlsx")
    source("multiplot.R")
    source("permira.test.R") # Functions for permutation tests of our data,
  #                            provided in another supplementary
  }
  
  # Import data and choose the factor (group) order for figures and color palette
  {
    miradata <- read.csv("HUSLAB Data final.txt", sep='\t', dec=',', header=T)
    head(miradata)
    
    # This block discards all subjects with no blood sample
    {
      miradata <- miradata[!is.na(miradata$Bage),]
    }
    
    # Build additional information based on the imported data
    {
      miradata$BMI = miradata$Weight/((miradata$Height*0.01)^2)
      miradata$Cholestanol.per.Campesterol = miradata$Cholestanol / miradata$Campesterol
      miradata$Lathosterol.per.Sitosterol = miradata$Lathosterol / miradata$Sitosterol
      miradata$Lathosterol.per.Campesterol = miradata$Lathosterol / miradata$Campesterol
      miradata$TBAs = miradata$TCA + miradata$TCDCA + miradata$TDCA + miradata$TLCA + miradata$TUDCA #Taurine-conjugated bile acids
      miradata$GBAs = miradata$GCA + miradata$GCDCA + miradata$GDCA + miradata$GLCA + miradata$GUDCA #Glycine-conjugated bile acids
      miradata$PBAbalance = miradata$CDCA / miradata$CA #Unconjugated primary bile acid "balance"
      miradata$ConjuBAbalance = miradata$TBAs / miradata$GBAs #Conjugation balance of all bile acids
      miradata$CPBAs = miradata$TCDCA + miradata$GCDCA + miradata$TCA + miradata$GCA #Conjugated primary bile acids
      miradata$uCPBAs = miradata$CDCA + miradata$CA #Unconjugated primary bile acids
      miradata$PBAs = miradata$CPBAs + miradata$uCPBAs #Total primary bile acids
      miradata$CSBAs = miradata$TLCA + miradata$GLCA + miradata$TDCA + miradata$GDCA + miradata$TUDCA + miradata$GUDCA #Conjugated secondary bile acids
      miradata$uCSBAs = miradata$LCA + miradata$DCA + miradata$UDCA #Unconjugated secondary bile acids
      miradata$SBAs = miradata$CSBAs + miradata$uCSBAs #Total secondary bile acids
      miradata$TcPBAs = miradata$TCA + miradata$TCDCA #Taurine-conjugated primary bile acids
      miradata$GcPBAs = miradata$GCA + miradata$GCDCA #Glycine-conjugated primary bile acids
      miradata$TcSBAs = miradata$TDCA + miradata$TLCA + miradata$TUDCA #Taurine-conjugated secondary bile acids
      miradata$GcSBAs = miradata$GDCA + miradata$GLCA + miradata$GUDCA #Glycine-conjugated secondary bile acids
      miradata$TotalBAs = miradata$PBAs + miradata$SBAs + miradata$HDCA #Total bile acids
      miradata$BAconjugation = (miradata$CPBAs + miradata$CSBAs) / (miradata$uCPBAs + miradata$uCSBAs + miradata$CPBAs + miradata$CSBAs) #Conjugation proportion of all bile acids
      miradata$PBAconjugation = miradata$CPBAs / (miradata$CPBAs + miradata$uCPBAs) #Conjugation proportion of primary bile acids
      miradata$SBAconjugation = miradata$CSBAs / (miradata$CSBAs + miradata$uCSBAs) #Conjugation proportion of secondary bile acids
      miradata$ConjuPBAbalance = miradata$TcPBAs / miradata$GcPBAs #Conjugation balance of primary bile acids
      miradata$ConjuSBAbalance = miradata$TcSBAs / miradata$GcSBAs #Conjugation balance of secondary bile acids
      miradata$PBAperSBA = miradata$PBAs / miradata$SBAs #Primary vs. secondary bile acid balance
      miradata$TotCAs = miradata$CA + miradata$TCA + miradata$GCA #Total cholic acids
      miradata$TotCDCAs = miradata$CDCA + miradata$TCDCA + miradata$GCDCA #Total chenodeoxycholic acids
      miradata$TotLCAs = miradata$LCA + miradata$TLCA + miradata$GLCA #Total litocholic acids
      miradata$TotDCAs = miradata$DCA + miradata$TDCA + miradata$GDCA #Total deoxycholic acids
      miradata$TotUDCAs = miradata$UDCA + miradata$TUDCA + miradata$GUDCA #Total ursodeoxycholic acids
      miradata$Vit_A_status <- -((-15.277)*miradata$S.RBP - 7.013*miradata$S.Prealb..mg.l./55 + 0.367*miradata$S.CRP + 24.218) # Vitamin A status based on Talsma et al. 2015
      # Add group labels. Group 5 is used for the article, as traditionally vegetarians include also pesco-vegetarians
      miradata$Group3 <- factor(miradata$Group3, c("Control", "Vegetarian", "Pesco-vegetarian", "Vegan"))
      miradata$Group4 <- factor(miradata$Group4, c("Control", "Control (vegan in daycare)", "Vegetarian", "Pesco-vegetarian", "Vegan"))
      miradata$Group5 <- miradata$Group4
      miradata$Group5[miradata$Group5 == "Pesco-vegetarian"] = "Vegetarian"
      miradata$Group5 <- factor(miradata$Group5) # Drop empty levels from Group5
      miradata$Group6 <- miradata$Group5
      levels(miradata$Group6) <- c(levels(miradata$Group6), "Vegan since pregnancy")
      miradata$Group6[miradata$ID == 104] <- "Vegan since pregnancy"
      miradata$Group6[miradata$ID == 108] <- "Vegan since pregnancy"
      miradata$Group6[miradata$ID == 203] <- "Vegan since pregnancy"
      miradata$Group6[miradata$ID == 204] <- "Vegan since pregnancy"
      # Convert cholesterol absorption and synthesis biomarkers to ug/mg, original data was *100ug/mg
      miradata[,45:52] = miradata[,45:52]/100
    }
  }
}


# ---------------------1. GENERAL CHARACTERISTICS AND ANTHROPOMETRICS----------------------------- #

{
  # For height and weight analysis, we first create a data frame for height and weight
  # averages and corresponding 2SD curves of Finnish children (averages and 2SD data
  # based on Saari A, et al. Ann Med 2010.)
  {
    anthroFIN <- read.csv("thl_kasvukayrat.txt", sep = '\t', dec = ',', header = T)
    # As we don't have Finnish BMI curves from under two-year-olds, we'll have to make an extrapolation
    # Ages 2-4 in BMI curve form nice linear trend (all three parameters of the BMI SDS-function have a Pearson correlation coefficient > 0.99 between these ages), so we'll use this information for extrapolating.
    for(i in 100:1) {
      anthroFIN$BMI_mean_F[i] <- anthroFIN$BMI_mean_F[i+1] - (anthroFIN$BMI_mean_F[300] - anthroFIN$BMI_mean_F[101])/200
      anthroFIN$BMI_SD_F[i] <- anthroFIN$BMI_SD_F[i+1] - (anthroFIN$BMI_SD_F[300] - anthroFIN$BMI_SD_F[101])/200
      anthroFIN$BMI_nu_F[i] <- anthroFIN$BMI_nu_F[i+1] - (anthroFIN$BMI_nu_F[300] - anthroFIN$BMI_nu_F[101])/200
      anthroFIN$BMI_mean_M[i] <- anthroFIN$BMI_mean_M[i+1] - (anthroFIN$BMI_mean_M[300] - anthroFIN$BMI_mean_M[101])/200
      anthroFIN$BMI_SD_M[i] <- anthroFIN$BMI_SD_M[i+1] - (anthroFIN$BMI_SD_M[300] - anthroFIN$BMI_SD_M[101])/200
      anthroFIN$BMI_nu_M[i] <- anthroFIN$BMI_nu_M[i+1] - (anthroFIN$BMI_nu_M[300] - anthroFIN$BMI_nu_M[101])/200
    }
    # For figure printing, let's calculate "2 SD" limits of height and BMI for each age and form a dataframe
    anthroFIN$BMI_SD2u_F <- anthroFIN$BMI_mean_F * ((1 + 2*anthroFIN$BMI_nu_F*anthroFIN$BMI_SD_F)^(1/anthroFIN$BMI_nu_F))
    anthroFIN$BMI_SD2l_F <- anthroFIN$BMI_mean_F * ((1 - 2*anthroFIN$BMI_nu_F*anthroFIN$BMI_SD_F)^(1/anthroFIN$BMI_nu_F))
    anthroFIN$BMI_SD2u_M <- anthroFIN$BMI_mean_M * ((1 + 2*anthroFIN$BMI_nu_M*anthroFIN$BMI_SD_M)^(1/anthroFIN$BMI_nu_M))
    anthroFIN$BMI_SD2l_M <- anthroFIN$BMI_mean_M * ((1 - 2*anthroFIN$BMI_nu_M*anthroFIN$BMI_SD_M)^(1/anthroFIN$BMI_nu_M))
    anthroFIN$BMI_SD2u <- (anthroFIN$BMI_SD2u_F + anthroFIN$BMI_SD2u_M) / 2
    anthroFIN$BMI_SD2l <- (anthroFIN$BMI_SD2l_F + anthroFIN$BMI_SD2l_M) / 2
    anthroFIN$BMI_mean <- (anthroFIN$BMI_mean_F + anthroFIN$BMI_mean_M) / 2
    anthroFIN$Height_SD2l <- ((anthroFIN$Height_mean_F - 2*anthroFIN$Height_SD_F) + (anthroFIN$Height_mean_M - 2*anthroFIN$Height_SD_M)) / 2
    anthroFIN$Height_SD2u <- ((anthroFIN$Height_mean_F + 2*anthroFIN$Height_SD_F) + (anthroFIN$Height_mean_M + 2*anthroFIN$Height_SD_M)) / 2
    anthroFIN$Height_mean <- (anthroFIN$Height_mean_F + anthroFIN$Height_mean_M) / 2
    lengthmean <- data.frame(anthroFIN$Age, anthroFIN$Height_mean, anthroFIN$Height_SD2l, anthroFIN$Height_SD2u)
    colnames(lengthmean) <- c("age", "median", "ll", "ul")
    bmimean <- data.frame(anthroFIN$Age, anthroFIN$BMI_mean, anthroFIN$BMI_SD2l, anthroFIN$BMI_SD2u)
    colnames(bmimean) <- c("age", "median", "ll", "ul")
    # Calculate height and BMI z-values
    for(i in 1:length(miradata$ID)) {
      if(miradata$Sex[i] == "M") {
        miradata$Height_z[i] <- (miradata$Height[i] - anthroFIN$Height_mean_M[anthroFIN$Age == round(miradata$AntAge[i], 2)])/anthroFIN$Height_SD_M[anthroFIN$Age == round(miradata$AntAge[i], 2)]
        miradata$BMI_z[i] <- (((miradata$BMI[i]/anthroFIN$BMI_mean_M[anthroFIN$Age == round(miradata$AntAge[i], 2)])^(anthroFIN$BMI_nu_M[anthroFIN$Age == round(miradata$AntAge[i], 2)])) - 1) / (anthroFIN$BMI_nu_M[anthroFIN$Age == round(miradata$AntAge[i], 2)] * anthroFIN$BMI_SD_M[anthroFIN$Age == round(miradata$AntAge[i], 2)])
      } else if(miradata$Sex[i] == "N") {
        miradata$Height_z[i] <- (miradata$Height[i] - anthroFIN$Height_mean_F[anthroFIN$Age == round(miradata$AntAge[i], 2)])/anthroFIN$Height_SD_F[anthroFIN$Age == round(miradata$AntAge[i], 2)]
        miradata$BMI_z[i] <- (((miradata$BMI[i]/anthroFIN$BMI_mean_F[anthroFIN$Age == round(miradata$AntAge[i], 2)])^(anthroFIN$BMI_nu_F[anthroFIN$Age == round(miradata$AntAge[i], 2)])) - 1) / (anthroFIN$BMI_nu_F[anthroFIN$Age == round(miradata$AntAge[i], 2)] * anthroFIN$BMI_SD_F[anthroFIN$Age == round(miradata$AntAge[i], 2)])
      }
    }
  }
  # For MUAC, import data from WHO z-scores for upper arm circumference
  {
    muac_boys <- read.csv("WHO_MUAC_z-scores_boys.txt", sep = '\t', dec = '.', header = T)
    muac_boys_ext <- read.csv("WHO_MUAC_z-scores_boys_extension.txt", sep = '\t', dec = ',', header = T)
    muac_boys <- rbind(muac_boys, muac_boys_ext) # Combine WHO 0-5 years old MUAC with the extension provided by Mramba et al. 2017
    muac_girls <- read.csv("WHO_MUAC_z-scores_girls.txt", sep = '\t', dec = '.', header = T)
    muac_girls_ext <- read.csv("WHO_MUAC_z-scores_girls_extension.txt", sep = '\t', dec = ',', header = T)
    muac_girls <- rbind(muac_girls, muac_girls_ext)
    muac_boys[,1] <- muac_boys[,1] / 365
    muac_girls[,1] <- muac_girls[,1] / 365
    colnames(muac_boys)[1] <- "Age"
    colnames(muac_girls)[1] <- "Age"
    miradata$MUAC_z <- rep(0, length(miradata$ID))
    for(i in 1:length(miradata$ID)) {
      if(is.na(miradata$MUAC[i])) {
        miradata$MUAC_z[i] <- NA
      } else {
        if(miradata$Sex[i] == "M") {
          age <- muac_boys$Age[which.min(abs(miradata$AntAge[i] - muac_boys$Age))] #match the age with MUAC-table ages
          for(j in 3:9) {
            if(muac_boys[muac_boys$Age == age, j] < miradata$MUAC[i]) {
              sd_limit_lower <- muac_boys[muac_boys$Age == age, j]
              sd_limit_lower_id <- j - 6 # This is meant to generate e.g. -4 if z-score of MUAC is between -4 and -3
            } else {
              sd_limit_upper <- muac_boys[muac_boys$Age == age, j+1]
              sd_limit_upper_id <- j - 6 # This is meant to generate e.g. -3 if z-score of MUAC is between -4 and -3
              break
            }
          }
          miradata$MUAC_z[i] <- sd_limit_lower_id + ((miradata$MUAC[i] - sd_limit_lower) / (sd_limit_upper - sd_limit_lower))
        } else if(miradata$Sex[i] == "N") {
          age <- muac_girls$Age[which.min(abs(miradata$AntAge[i] - muac_girls$Age))] #match the age with MUAC-table ages
          for(j in 3:9) {
            if(muac_girls[muac_girls$Age == age, j] < miradata$MUAC[i]) {
              sd_limit_lower <- muac_girls[muac_girls$Age == age, j]
              sd_limit_lower_id <- j - 6 # This is meant to generate e.g. -4 if z-score of MUAC is between -4 and -3
            } else {
              sd_limit_upper <- muac_girls[muac_girls$Age == age, j+1]
              sd_limit_upper_id <- j - 6 # This is meant to generate e.g. -3 if z-score of MUAC is between -4 and -3
              break
            }
          }
          miradata$MUAC_z[i] <- sd_limit_lower_id + ((miradata$MUAC[i] - sd_limit_lower) / (sd_limit_upper - sd_limit_lower))
        }
      }
    }
  }
  
  # Age vs. height plot
  {
    hplot = ggplot(data=miradata, aes(x = AntAge, y = Height, colour = Group5, shape = Group5, size = Group5)) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          labels = c("Control", "Vegetarian", "Vegan"),
                          values=c("#000000FF", "#999999FF", "#028910FF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan"),
                         guide = guide_legend(title.position = "right")) +
      scale_size_manual(name = "Dietary group", values=c(2,2,2),
                        labels = c("Control", "Vegetarian", "Vegan"),
                        guide = guide_legend(title.position = "right")) +
      geom_line(inherit.aes = FALSE, data = lengthmean, aes(x = age, y = median), colour = "black", show.legend = F, linetype = "dashed", size = 0.5) +
      annotate(geom = "text", label = "Finnish\nmedian", x = 7.6, y = 128, angle = 20, size = 4) +
      geom_ribbon(inherit.aes = FALSE, data = lengthmean,
                  aes(x = age, ymin = ll, ymax = ul),
                  alpha = 0.15, show.legend = F) +
      geom_point() +
      geom_point(data = subset(miradata, Group5 == "Vegan")) +
      geom_point(data = subset(miradata, Group5 == "Vegetarian")) +
      labs(x = "Age (years)", y = "Height (cm)") +
      scale_x_continuous(limits = c(1, 8), expand = c(0, 0)) +
      theme(legend.justification = c(0.95,0.05),
            legend.position = c(0.95,0.05),
            legend.text = element_text(size = 12),
            legend.text.align = 0,
            legend.title = element_blank(),
            legend.key.size = unit(0.75, "lines"),
            legend.background = element_rect(fill="#FFDC9199", colour="#676767", size = 0.3),
            legend.key = element_rect(fill=NA),
            axis.title = element_text(size = 12))
    print(hplot)
    ggsave(file = "manusc/figures/Age vs. height.pdf", plot = hplot, device = "pdf", width = 6, height = 4)
  }
  # Age vs. BMI plot
  {
    bplot = ggplot(data=miradata, aes(x = AntAge, y = BMI, colour = Group5, shape = Group5, size = Group5)) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          labels = c("Control", "Vegetarian", "Vegan"),
                          values=c("#000000FF", "#999999FF", "#028910FF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan"),
                         guide = guide_legend(title.position = "right")) +
      scale_size_manual(name = "Dietary group", values=c(2,2,2),
                        labels = c("Control", "Vegetarian", "Vegan"),
                        guide = guide_legend(title.position = "right")) +
      geom_smooth(inherit.aes = FALSE, data = bmimean, aes(x = age, y = median), colour = "black", show.legend = F, method="loess", se = FALSE, linetype = "dashed", size = 0.5) +
      annotate(geom = "text", label = "Finnish\nmedian", x = 7.6, y = 16.0, angle = 9, size = 4) +
      geom_ribbon(inherit.aes = FALSE, data = bmimean,
                  aes(x = age, ymin = ll, ymax = ul),
                  alpha = 0.15, show.legend = F) +
      geom_point() +
      geom_point(data = subset(miradata, Group5 == "Vegan")) +
      geom_point(data = subset(miradata, Group5 == "Vegetarian")) +
      labs(x = "Age (years)", y = "BMI (kg/m^2)") +
      scale_x_continuous(limits = c(1, 8), expand = c(0, 0)) +
      theme(legend.justification = c(0.95,0.05),
            legend.position = c(0.95,0.77),
            legend.text = element_text(size = 12),
            legend.text.align = 0,
            legend.title = element_blank(),
            legend.key.size = unit(0.75, "lines"),
            legend.background = element_rect(fill="#FFDC9199", colour="#676767", size = 0.3),
            legend.key = element_rect(fill=NA),
            axis.title = element_text(size = 12))
    print(bplot)
    ggsave(file = "manusc/figures/Age vs. BMI.pdf", plot = bplot, device = "pdf", width = 6, height = 4)
  }
}

# -------------------------2. SERUM BIOMARKER ANALYSES---------------------------- #
{
  # Create ad-hoc name and unit vectors to print smoothly
  {
    miranames = c(rep("", 13), "Vitamin D2+D3", "Vitamin D3", "Vitamin D2", "Zinc", "Transthyretin", "Triglycerides",
                  "LDL-cholesterol", "HDL-cholesterol", "Total cholesterol", "Total Chol per HDL ratio", 
                  "LDL per HDL ratio", "Leukocytes", "Erythrocytes", "Hemoglobin", "Hematocrit",
                  "MCV (Mean RBC volume)", "RDW (RBC distribution width)", "MCH (Mean RBC hemoglobin)",
                  "MCHC (Mean Hb conc. in RBC)", "Thrombocytes", "B12-TC2", "Folate (Erythrocytes)",
                  "Hematocrit", "Glucose (plasma)", "Ferritin", "TfR", "Body Iron Storages (BIS)", "RBP",
                  "CRP", "AGP", "Iodine (urine)", "Cholestanol", "Cholestenol", "Desmosterol",
                  "Lathosterol", "Campesterol", "Sitosterol", "Avenasterol", "Squalene",
                  "Total cholesterol (HPLC)", "Tauroursodeoxycholic acid", "Ursodeoxycholic acid",
                  "Glycochenodeoxycholic acid", "Taurocholic acid", "Hyodeoxycholic acid", "Glycocholic acid",
                  "Glycolithocholic acid", "Chenodeoxycholic acid", "Taurodeoxycholic acid",
                  "Glycoursodeoxycholic acid", "Deoxycholic acid", "Glycodeoxycholic acid", "Lithocholic acid",
                  "Taurolithocholic acid", "Taurochenodeoxycholic acid", "Cholic acid",
                  "7-alpha-hydroxy-4-cholesten-3-one", "Cholestanol per Campesterol ratio",
                  "Lathosterol per Sitosterol ratio", "Lathosterol per Campesterol ratio", "BMI",
                  "Taurine-conjugated bile acids", "Glycine-conjugated bile acids", "Primary bile acid balance (CDCA per CA)",
                  "Bile acid conjugation balance (tauro- per glyco-)", "Conjugated primary bile acids",
                  "Unconjugated primary bile acids", "Primary bile acids", "Conjugated secondary bile acids",
                  "Unconjugated secondary bile acids", "Secondary bile acids", "Taurine-conjugated primary bile acids",
                  "Glycine-conjugated primary bile acids", "Taurine-conjugated secondary bile acids",
                  "Glycine-conjugated secondary bile acids", "Total bile acids", "Bile acid conjugation ratio",
                  "Primary bile acid conjugation ratio", "Secondary bile acid conjugation ratio",
                  "PBA conjugation balance (tauro- per glyco-)", "SBA conjugation balance (tauro- per glyco-)",
                  "Primary per secondary bile acids", "Total cholic acids", "Total chenodeoxycholic acids",
                  "Total lithocholic acids", "Total deoxycholic acids", "Total ursodeoxycholic acids", "Vitamin A status")
    
    miraunits = c(rep("", 13), "nmol/l","nmol/l","nmol/l","µmol/l","mg/l","mmol/l","mmol/l","mmol/l","mmol/l","", "",
                  "*10^9/l","*10^12/l","g/l", "%", "fl","%","pg","g/l","*10^9/l","pmol/l","nmol/l","%",
                  "mmol/l", "µg/l", "mg/l", "mg/kg", "RE (µmol/l)", "mg/l", "g/l", "µg/l",
                  rep("µg/mg of total cholesterol",8), "mmol/l", rep("µmol/l",17), "", "", "", "", "µmol/l", "µmol/l",
                  "", "", rep("µmol/l", 11), "", "", "", "", "", "", rep("µmol/l", 5), "")
  }
  
  # Choose groups to be plotted (or plot all subjects by choosing data = miradata)
  {
    data = miradata #Choose all subjects
    # data = subset(mirablood, Group3 %in% c("Omnivore", "Vegan")) #Choose groups for figures
    # data = subset(mirablood, !(ID %in% c(104, 403, 901, 902))) # Choose IDs for figures
  }
    
  # Get p-values of desired metabolites (assessed in rows.to.measure) from
  # permira.test and adjust them for multiple comparisons (B-H)
  
  # CAUTION! THIS PART OF CODE WILL TAKE SOME TIME TO RUN.
  {
    # Calculate p-values with permutation tests (ad-hoc permira.test function)
    {
      rows.to.measure <- c(14:15,17:23,25:33,35,37:52,54:70,75:101) # Choose which metabolites are to be tested
      doms <- rep(NA,max(rows.to.measure))
      p.values <- rep(NA,max(rows.to.measure))
      for(i in rows.to.measure) {
        if(is.na(p.values[i])) {
          results <- permira.test(miradata,
                                  metabolite = colnames(miradata)[i],
                                  method = "pr.index",
                                  in.class = FALSE)
          doms[i] <- results[1] # difference of means
          p.values[i] <- results[2] # p-value
        }
      }
    }
    
    # Multiple testing adjustments with Benjamini-Hochberg,
    # according to which set the data belongs
    {
      # As the previous approach yields lists, we need to make it into vector for multiple comparison
      # adjustment to work
      p.values <- as.numeric(as.character(p.values)) #This replaces NULLs with NAs so that unlist -function does not drop non-calculated metabolites out
      p.valuesvec <- unlist(p.values)
      p.valuesvec[p.valuesvec == 0] <- NA # Replace all the zeros with NA
      # CAUTION: Check if there was any permutation test yielding p = 0 to avoid problems with previous line
      p.values.adjusted <- rep(NA,length(p.valuesvec))
      
      # First the targeted metabolites set
      set1 <- c(14:15,17:33,35,37,44,101) #Assign ID's of targeted metabolites that were tested
      p.adj.set1 <- p.adjust(p.valuesvec[set1],method = "BH")
      for(i in 1:length(p.adj.set1)) {
        p.values.adjusted[set1[i]] <- p.adj.set1[i]
      }
      
      # Second, iron and acute phase protein set
      set2 <- c(38:43)
      p.adj.set2 <- p.adjust(p.valuesvec[set2],method = "BH")
      for(i in 1:length(p.adj.set2)) {
        p.values.adjusted[set2[i]] <- p.adj.set2[i]
      }
      
      # Cholesterol metabolism
      set3 <- c(45:52)
      p.adj.set3 <- p.adjust(p.valuesvec[set3],method = "BH")
      for(i in 1:length(p.adj.set3)) {
        p.values.adjusted[set3[i]] <- p.adj.set3[i]
      }
      
      # Bile acids
      set4 <- c(54:70,75:100)
      p.adj.set4 <- p.adjust(p.valuesvec[set4],method = "BH")
      for(i in 1:length(p.adj.set4)) {
        p.values.adjusted[set4[i]] <- p.adj.set4[i]
      }
    }
  }
  
  # Form a table of all analyzed targeted metabolites
  {
    miradata_titles = c(rep("", 13),
                        "Vitamin D2 + D3 - nmol/l",
                        "Vitamin D3 - nmol/l",
                        "Vitamin D2 - nmol/l",
                        "Zinc - mol/l",
                        "Transthyretin - mg/l",
                        "Triglycerides - mmol/l",
                        "LDL Cholesterol - mmol/l",
                        "HDL Cholesterol - mmol/l",
                        "Total Cholesterol - mmol/l",
                        "Total / HDL Cholesterol",
                        "LDL / HDL Cholesterol",
                        "Leukocytes - E9/l",
                        "Erythrocytes - E12/l",
                        "Hemoglobin - g/l",
                        "Hematocrit - %",
                        "Mean corpuscular volume (MCV) - fl",
                        "Red cell distribution width (RDW) - %",
                        "Mean corpuscular hemoglobin (MCH) - pg",
                        "Mean corpuscular hemoglobin concentration (MCHC) - g/l",
                        "Trombocytes - E9/l",
                        "Vitamin B12 (Bound to TC2) - pmol/l",
                        "Folate (Vitamin B9) in erythrocytes - nmol/l",
                        "",
                        "Glucose - mmol/l",
                        "Ferritin - ug/l",
                        "Transferrin receptor (TfR) - mg/l",
                        "Body Iron Stores (BIS)",
                        "Retinol-binding protein (RBP) - umol/l",
                        "C-reactive protein (CRP) - mg/l",
                        "AGP - g/l",
                        "Iodine in urine - ug/l",
                        "Cholestanol - ug/mg of total cholesterol",
                        "Cholestenol - ug/mg of total cholesterol",
                        "Desmosterol - ug/mg of total cholesterol",
                        "Lathosterol - ug/mg of total cholesterol",
                        "Campesterol - ug/mg of total cholesterol",
                        "Sitosterol - ug/mg of total cholesterol",
                        "Avenasterol - ug/mg of total cholesterol",
                        "Squalene - ug/mg of total cholesterol",
                        "",
                        "Tauroursodeoxycholic acid - umol/l",
                        "Ursodeoxycholic acid - umol/l",
                        "Glycochenodeoxycholic acid - umol/l",
                        "Taurocholic acid - umol/l",
                        "Hyodeoxycholic acid - umol/l",
                        "Glycocholic acid - umol/l",
                        "Glycolithocholic acid - umol/l",
                        "Chenodeoxycholic acid - umol/l",
                        "Taurodeoxycholic acid - umol/l",
                        "Glycoursodeoxycholic acid - umol/l",
                        "Deoxycholic acid - umol/l",
                        "Glycodeoxycholic acid - umol/l",
                        "Lithocholic acid - umol/l",
                        "Taurolithocholic acid - umol/l",
                        "Taurochenodeoxycholic acid - umol/l",
                        "Cholic acid - umol/l",
                        "7-alpha-hydroxy-4-cholesten-3-one - umol/l",
                        "",
                        "Cholestanol / Campesterol",
                        "Lathosterol / Sitosterol",
                        "Lathosterol / Campesterol",
                        "Taurine-conjugated bile acids - umol/l",
                        "Glycine-conjugated bile acids - umol/l",
                        "Primary bile acid balance (CDCA per CA)",
                        "Bile acid conjugation balance (tauro- per glyco-)",
                        "Conjugated primary bile acids - umol/l",
                        "Unconjugated primary bile acids - umol/l",
                        "Primary bile acids - umol/l",
                        "Conjugated secondary bile acids - umol/l",
                        "Unconjugated secondary bile acids - umol/l",
                        "Secondary bile acids - umol/l",
                        "Taurine-conjugated primary bile acids - umol/l",
                        "Glycine-conjugated primary bile acids - umol/l",
                        "Taurine-conjugated secondary bile acids - umol/l",
                        "Glycine-conjugated secondary bile acids - umol/l",
                        "Total bile acids - umol/l",
                        "Bile acid conjugation ratio - % of all bile acids",
                        "Primary bile acid conjugation ratio - % of all primary bile acids",
                        "Secondary bile acid conjugation ratio - % of all secondary bile acids",
                        "PBA conjugation balance (tauro- per glyco-)",
                        "SBA conjugation balance (tauro- per glyco-)",
                        "Primary bile acids / secondary bile acids",
                        "Total cholic acids - umol/l",
                        "Total chenodeoxycholic acids - umol/l",
                        "Total lithocholic acids - umol/l",
                        "Total deoxycholic acids - umol/l",
                        "Total ursodeoxycholic acids - umol/l",
                        "Vitamin A status (logistic regression with RBP, TTR and CRP)")
    miradata_tissues <- c(rep("",13),
                          rep("Serum", 5),
                          rep("Plasma",6),
                          rep("EDTA blood", 9),
                          "Serum",
                          "Erythrocytes",
                          "EDTA blood",
                          "Plasma",
                          rep("Serum", 6),
                          "Urine",
                          rep("Serum", 26),
                          "",
                          rep("Serum", 29),
                          "")
    miradata_medians <- data.frame(miradata_titles)
    metabolites_to_table <- c(14:35,37:41,101,42:52,54:70,72:100)
    for(i in 1:length(miradata_medians$miradata_titles)) {
      if(!(i %in% metabolites_to_table) | !is.numeric(miradata[,i])) {
        miradata_medians$Nutrient[i] <- miradata_titles[i]
        miradata_medians$Tissue[i] <- miradata_tissues[i]
        miradata_medians$Omnivore[i] <- 0
        miradata_medians$Vegetarian[i] <- 0
        miradata_medians$Vegan[i] <- 0
        miradata_medians$p.value[i] <- 1
      } else {
        miradata_medians$Nutrient[i] <- miradata_titles[i]
        miradata_medians$Tissue[i] <- miradata_tissues[i]
        miradata_medians$Omnivore[i] <- paste(signif(median(miradata[miradata$Group5 == "Control",i], na.rm = T), digits = 3),
                                                      "\n[",
                                                      signif(min(miradata[miradata$Group5 == "Control",i], na.rm = T), digits = 3),
                                                      "-",
                                                      signif(max(miradata[miradata$Group5 == "Control",i], na.rm = T), digits = 3),
                                                      "]")
        miradata_medians$Vegetarian[i] <- paste(signif(median(miradata[miradata$Group5 == "Vegetarian",i], na.rm = T), digits = 3),
                                                        "\n[",
                                                        signif(min(miradata[miradata$Group5 == "Vegetarian",i], na.rm = T), digits = 3),
                                                        "-",
                                                        signif(max(miradata[miradata$Group5 == "Vegetarian",i], na.rm = T), digits = 3),
                                                        "]")
        miradata_medians$Vegan[i] <- paste(signif(median(miradata[miradata$Group5 == "Vegan",i], na.rm = T), digits = 3),
                                                   "\n[",
                                                   signif(min(miradata[miradata$Group5 == "Vegan",i], na.rm = T), digits = 3),
                                                   "-",
                                                   signif(max(miradata[miradata$Group5 == "Vegan",i], na.rm = T), digits = 3),
                                                   "]")
        miradata_medians$p.value[i] <- p.values.adjusted[i]
      }
    }
    miradata_medians <- miradata_medians[,2:7] #There was problems with renaming nutrient names in for loop, so a quickfix is to have two nutrient name -columns and remove the other one here
    colnames(miradata_medians) <- c("Measurement",
                                    "Tissue",
                                    paste("Omnivore\nN =", sum(miradata_ravinto$crv_yhd_5 == 5)),
                                    paste("Vegetarian\nN =", sum(miradata_ravinto$crv_yhd_5 == 2)),
                                    paste("Vegan\nN =", sum(miradata_ravinto$crv_yhd_5 == 1)),
                                    "p-value*\nOmnivore vs. Vegan")
    miradata_medians <- miradata_medians[miradata_medians$Measurement != "",] # Curate all non-titled, thus irrelevant rows (e.g. ID, sex,...)
    write.xlsx(miradata_medians, "bloodsample_data.xlsx")
  }
  
  # Print metabolites to a figure and form plot objects for final figures.
  ## (For technical reasons in R, you may need to name the "plot" objects manually after
  ## each loop. Generating plotXX -named objects automatically in the loop yields
  ## identical objects with parameters according to the latest value of i. PlotXX-objects
  ## are later used in the code for generating the publication plots)
  {
    for (i in c(14:15, 17:33, 35:52, 54:100)) { #All figures
    #for (i in c(14:15, 18, 20:23, 35, 41, 44:51, 61, 69, 78, 93, 94)) { #Figures in the publication
      plot <- ggplot(data = data, aes(x=data$Group5, y=data[,i])) +
        theme_light() +
        scale_colour_manual(name = "Dietary group",
                            values = c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                            guide = guide_legend(title.position = "right")) +
        scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                           labels = c("Control", "Vegetarian", "Vegan")) +
        scale_size_manual(name = "Dietary group", values=c(1,1,1),
                          labels = c("Control", "Vegetarian", "Vegan")) +
        scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
        geom_point(position = position_jitter(width = .1, height = 0),
                   aes(colour = Group6, shape = Group5, size = Group5)) +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
        stat_summary(geom = "errorbar", fun.y = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
        labs(y = miraunits[i], x = "", title = miranames[i]) +
        #labs(y = miraunits[14], x = "", title = "Total Vitamin 25(OH)D   ") +
        geom_signif(comparisons = list(c("Control", "Vegan")),
                    test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) + 
        # following geom_hlines are for printing possible limits for deficiency
        # ("dashed") or insufficiency ("dotted")
        #geom_hline(aes(yintercept=1.17, linetype="dotted", alpha = .5)) + #RBP insufficiency
        #geom_hline(aes(yintercept=0.83, linetype="dashed", alpha = .5)) + #RBP deficiency
        #geom_hline(aes(yintercept=50, linetype="dashed", alpha = .5)) + #VitD deficiency
        #geom_hline(aes(yintercept=30, linetype="dashed", alpha = .5)) + #VitB12 deficiency
        #geom_hline(aes(yintercept=0, linetype="dashed", alpha = .5)) + #Iron (BIS) and VitA deficiency
        coord_cartesian(ylim = c(min(data[,i], na.rm=TRUE)*0.95,
                                 max(data[,i], na.rm=TRUE)*1.2)) +
        theme(plot.title = element_text(hjust = 0.5, size = 10),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
              axis.text.y = element_text(size = 10),
              axis.title.y = element_text(size = 10),
              legend.position = "none") +
        guides(colour=F, linetype=F, size=F, shape=F)
      #print to file
      filename = paste("manusc/figures/Figure subplots/", miranames[i], ".png", sep="")
      ggsave(file=filename, plot=plot, width=2.5, height=3)
      message("Printing figure ", i) # Check up during the loop that something is happening
      Sys.sleep(0.01) # To make sure the progress message from previous line is printed each time.
    }
  }
}
  
# ----------------------3. UNTARGETED METABOLOMICS DATA------------------------- #

{
  # Import the data and create lists of text labels
  {
    mirametab <- read.csv("MIRA Metabolomics final.txt", header = T, sep = '\t', dec=',', stringsAsFactors = F)
    
    # AMINO ACID PLOT
    #Set column ID's and types for amino acids manually
    aminoacid_ID <- c(13,19,24,38,
                      42,44,48,61,
                      63,65,84,85,
                      86,92,106,126,
                      143,149,190)
    aminoacid_labels <- c("Gly", "Ala", "Ser", "Pro",
                          "Val#", "Thr", "Cys", "Leu & Ile",
                          "Asn", "Asp", "Gln", "Lys", 
                          "Glu", "Met", "His", "Phe",
                          "Arg", "Tyr", "Trp")
    type <- c("Non-essential", "Non-essential", "Non-essential", "Non-essential",
              "Essential", "Essential", "Non-essential", "Essential",
              "Non-essential", "Non-essential", "Non-essential", "Essential",
              "Non-essential", "Essential", "Essential", "Essential",
              "Non-essential", "Non-essential","Essential")
  }
  
  # Calculate differences of means and their 95% confidence interval lengths
  # assuming heteroscedasticity (Behrens-Fisher approach, originally from 1935).
  # Equation for degrees of freedom, v, is also called Welch-Satterthwaite equation.
  {
    AAdoms <- rep(0,length(aminoacid_ID))
    AAcilens <- rep(0, length(aminoacid_ID))
    for(i in 1:length(aminoacid_ID)) {
      n1 = length(mirametab$ID[mirametab$Group5=="Vegan"])
      n2 = length(mirametab$ID[mirametab$Group5=="Control"])
      m1 = mean(mirametab[mirametab$Group5=="Vegan",aminoacid_ID[i]])
      m2 = mean(mirametab[mirametab$Group5=="Control",aminoacid_ID[i]])
      # Means of peak intensities normalized to control group mean:
      norm_m1 = mean(mirametab[mirametab$Group5=="Vegan",aminoacid_ID[i]]/m2)
      norm_m2 = mean(mirametab[mirametab$Group5=="Control",aminoacid_ID[i]]/m2) #This should yield 1 every time.
      s1 = sd(mirametab[mirametab$Group5=="Vegan",aminoacid_ID[i]]/m2)
      s2 = sd(mirametab[mirametab$Group5=="Control",aminoacid_ID[i]]/m2)
      v <- (s1^2/n1+s2^2/n2)^2/(s1^4/(n1^2*(n1-1)) + s2^4/(n2^2*(n2-1)))
      #Calculate differences of normalized means (normalized to control group mean)
      AAdoms[i] <- norm_m1 - norm_m2
      AAcilens[i] <- qt(0.975,v)*sqrt(s1^2/n1 + s2^2/n2)
    }
    AAdata <- data.frame(AAdoms, AAcilens, aminoacid_labels, type)
    colnames(AAdata) <- c("dom", "cilen", "label", "type")
  }
  
  # Make the plot
  {
    aminoplot <- ggplot(data = AAdata) +
      theme_light() +
      scale_fill_manual(name = "",
                          values=c("#BC3C29FF", "#E18727FF")) +
      geom_bar(aes(x=label, y=dom, fill=type), stat="identity") +
      geom_errorbar(aes(x=label, ymin=dom-cilen, ymax=dom+cilen),
                    width=0.4, colour="darkgray", size=0.5, alpha = 1) +
      facet_wrap(~type, scales = "free_x") +
      labs(y = "Average normalized peak intensity\nin vegans compared to controls", x = "", title = "") +
      theme(plot.title = element_text(hjust = 0.5, size = 10),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 8),
            axis.title.y = element_text(size = 10),
            legend.position = "none") +
      guides(colour=F, linetype=F, size=F, shape=F)
    ggsave(file="manusc/figures/aminoacids.png",
           plot = aminoplot, 
           width = 5, 
           height = 3, 
           device = "png")
  }
}

# -----------------------4. BUILD AND PRINT THE FIGURES-------------------------- #

# Figure 1
{
  # Create titles
  df <- data.frame()
  titleplot <- ggplot(df) +
    theme_void() +
    geom_point() +
    xlim(0, 10) +
    ylim(0.5, 1.5) +
    annotate("label", x = 5, y = 1,
             label = "A                Age vs. height of study subjects                   ", size = 6, fill = "#FFDC9199")
  titleplot2 <- ggplot(df) +
    theme_void() +
    geom_point() +
    xlim(0, 10) +
    ylim(0.5, 1.5) +
    annotate("label", x = 5, y = 1,
             label = "B                Age vs. BMI of study subjects                  ", size = 6, fill = "#FFDC9199")
  # Assign the layout matrix (how the plots are arranged)
  ployout = matrix(c(1,1,1,1,2,2,2,2,
                     3,3,3,3,4,4,4,4,
                     3,3,3,3,4,4,4,4,
                     3,3,3,3,4,4,4,4,
                     3,3,3,3,4,4,4,4,
                     3,3,3,3,4,4,4,4,
                     3,3,3,3,4,4,4,4,
                     3,3,3,3,4,4,4,4), ncol = 8, byrow = T)
  # Bring the plot together and print it to pdf
  ggsave(file="manusc/figures/Figure 1.pdf",
         plot=multiplot(titleplot, titleplot2,
                        hplot, bplot,
                        layout = ployout), width = 12, height = 5, device = "pdf")
}

# Figure 2
{
  # Create titles
  df <- data.frame()
  titleplot <- ggplot(df) +
    theme_void() +
    geom_point() +
    xlim(0, 10) +
    ylim(0.5, 1.5) +
    annotate("label", x = 5, y = 1,
             label = "A                                                Vitamins and minerals                                                ",
             size = 6, fill = "#FFDC9199")
  titleplot2 <- ggplot(df) +
    theme_void() +
    geom_point() +
    xlim(0, 10) +
    ylim(0.5, 1.5) +
    annotate("label", x = 5, y = 1,
             label = "B                                               Cholesterol metabolism                                               ",
             size = 6, fill = "#FFDC9199")
  subtitleplot <- ggplot(df) +
    theme_void() +
    geom_point() +
    xlim(0, 10) +
    ylim(0.5, 1.5) +
    annotate("label", x = 5, y = 1,
             label = "                                    Cholesterol absorption                                  ",
             size = 6, fill = "#F6F3EE")
  subtitleplot2 <- ggplot(df) +
    theme_void() +
    geom_point() +
    xlim(0, 10) +
    ylim(0.5, 1.5) +
    annotate("label", x = 5, y = 1,
             label = "                                  Cholesterol biosynthesis                                 ",
             size = 6, fill = "#F6F3EE")
  # Create an empty plot if needed to adjust figures
  emptyplot <- ggplot(df) +
    theme_void() +
    geom_point() +
    xlim(0, 10) +
    ylim(0, 10)
  # For this figure, we need to import the cholesterol metabolism figure from png
  img <- readPNG("manusc/figures/Cholesterol metabolism 2.0.png")
  img <- rasterGrob(img, interpolate=TRUE)
  cholpathway <- ggplot(df) +
    theme_void() +
    annotation_custom(img, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    xlim(0, 10) +
    ylim(0, 10)
  
  # Assign the layout matrix (how the plots are arranged)
  ployout = matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                     2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,
                     2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,
                     2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,
                     2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,
                     7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,
                     7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,
                     7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,
                     7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,
                     12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,
                     13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,
                     13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,
                     13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,
                     13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,
                     13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,
                     15,15,15,16,16,16,16,16,16,16,16,16,16,16,16,
                     15,15,15,17,17,17,18,18,18,19,19,19,20,20,20,
                     15,15,15,17,17,17,18,18,18,19,19,19,20,20,20,
                     15,15,15,17,17,17,18,18,18,19,19,19,20,20,20,
                     15,15,15,17,17,17,18,18,18,19,19,19,20,20,20,
                     21,21,21,22,22,22,22,22,22,22,22,22,22,22,22,
                     21,21,21,23,23,23,24,24,24,25,25,25,26,26,26,
                     21,21,21,23,23,23,24,24,24,25,25,25,26,26,26,
                     21,21,21,23,23,23,24,24,24,25,25,25,26,26,26,
                     21,21,21,23,23,23,24,24,24,25,25,25,26,26,26), ncol = 15, byrow = T)
  # Bring the plot together and print it to pdf
  ggsave(file="manusc/figures/Figure 2.pdf",
         plot=multiplot(titleplot,
                        plot14 +
                          scale_colour_manual(name = "Dietary group",
                                                    values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                                    guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 10)),
                        plot15 +
                          scale_colour_manual(name = "Dietary group",
                                              values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                              guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 10)),
                        plot41 +
                          scale_colour_manual(name = "Dietary group",
                                              values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                              guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 10)),
                        plot18 +
                          scale_colour_manual(name = "Dietary group",
                                              values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                              guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 10)),
                        plot101 +
                          scale_colour_manual(name = "Dietary group",
                                              values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                              guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 10)),
                        plot35 +
                          scale_colour_manual(name = "Dietary group",
                                              values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                              guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 9),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 10)),
                        plot17 +
                          scale_colour_manual(name = "Dietary group",
                                              values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                              guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 10)),
                        plot38 +
                          scale_colour_manual(name = "Dietary group",
                                              values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                              guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 10)),
                        plot40 +
                          scale_colour_manual(name = "Dietary group",
                                              values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                              guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 9),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 10)),
                        plot44 +
                          scale_colour_manual(name = "Dietary group",
                                              values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                              guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 10)),
                        titleplot2,
                        plot22 +
                          scale_colour_manual(name = "Dietary group",
                                              values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                              guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 10)),
                        cholpathway,
                        plot20 +
                          scale_colour_manual(name = "Dietary group",
                                              values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                              guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 10)),
                        subtitleplot,
                        plot45 +
                          scale_colour_manual(name = "Dietary group",
                                              values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                              guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 8)),
                        plot49 +
                          scale_colour_manual(name = "Dietary group",
                                              values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                              guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 8)),
                        plot50 +
                          scale_colour_manual(name = "Dietary group",
                                              values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                              guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 8)),
                        plot51 +
                          scale_colour_manual(name = "Dietary group",
                                              values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                              guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 8)),
                        plot21 +
                          scale_colour_manual(name = "Dietary group",
                                              values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                              guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 10)),
                        subtitleplot2,
                        plot52 +
                          scale_colour_manual(name = "Dietary group",
                                              values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                              guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 8)),
                        plot47 +
                          scale_colour_manual(name = "Dietary group",
                                              values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                              guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 8)),
                        plot46 +
                          scale_colour_manual(name = "Dietary group",
                                              values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                              guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 8)),
                        plot48 +
                          scale_colour_manual(name = "Dietary group",
                                              values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                              guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
                          theme(plot.title = element_text(hjust = 0.5, size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                axis.title.y = element_text(size = 8)),
                        #emptyplot,
                        layout = ployout),
         width = 9,
         height = 14,
         device = "pdf")
}

# Figure 3
{
  # Create titles
  df <- data.frame()
  titleplot <- ggplot(df) +
    theme_void() +
    geom_point() +
    xlim(0, 10) +
    ylim(0.5, 1.5) +
    annotate("label", x = 5, y = 1,
             label = "A                                Pathway analysis of untargeted metabolomics                           ",
             size = 6, fill = "#FFDC9199")
  titleplot2 <- ggplot(df) +
    theme_void() +
    geom_point() +
    xlim(0, 10) +
    ylim(0.5, 1.5) +
    annotate("label", x = 5, y = 1,
             label = "B                                              Serum levels of bile acids                                             ",
             size = 6, fill = "#FFDC9199")
  titleplot3 <- ggplot(df) +
    theme_void() +
    geom_point() +
    xlim(0, 10) +
    ylim(0.5, 1.5) +
    annotate("label", x = 5, y = 1,
             label = "C                           Levels of amino acids in untargeted metabolomics                         ",
             size = 6, fill = "#FFDC9199")
  # Create an empty plot if needed to adjust figures
  emptyplot <- ggplot(df) +
    theme_void() +
    geom_point() +
    xlim(0, 10) +
    ylim(0, 10)
  
  # Assign the layout matrix (how the plots are arranged)
  ployout = matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                     2,2,2,2,2,2,2,12,3,3,3,3,3,3,3,
                     2,2,2,2,2,2,2,12,3,3,3,3,3,3,3,
                     2,2,2,2,2,2,2,12,3,3,3,3,3,3,3,
                     2,2,2,2,2,2,2,12,3,3,3,3,3,3,3,
                     2,2,2,2,2,2,2,12,3,3,3,3,3,3,3,
                     4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
                     5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,
                     5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,
                     5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,
                     5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,
                     10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,
                     11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,
                     11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,
                     11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,
                     11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,
                     11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,
                     11,11,11,11,11,11,11,11,11,11,11,11,11,11,11), ncol = 15, byrow = T)
  # Bring the plot together and print it to pdf
  ggsave(file="manusc/figures/Figure 3.pdf",
         plot=multiplot(titleplot,
                        pathway_figure_pos, pathway_figure_neg,
                        titleplot2,
                        plot78 +
                          scale_colour_manual(name = "Dietary group",
                                                     values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                                     guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")),
                        plot93 +
                          scale_colour_manual(name = "Dietary group",
                                                     values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                                     guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")),
                        plot94 +
                          scale_colour_manual(name = "Dietary group",
                                                     values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                                     guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")),
                        plot61 +
                          scale_colour_manual(name = "Dietary group",
                                                     values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                                     guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")),
                        plot69 +
                          scale_colour_manual(name = "Dietary group",
                                                     values=c("#000000FF", "#999999FF", "#028910FF", "#028910FF"),
                                                     guide = guide_legend(title.position = "right")) + 
                          scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")),
                        titleplot3,
                        aminoplot,
                        layout = ployout), width = 9, height = 11, device = "pdf")
}