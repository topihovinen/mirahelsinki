##    This is an R code for statistical analysis and generating figures of MIRA Helsinki study    ##
# ------------------------------------------------------------------------------------------------ #
#             .-"''-.  _                                                                           #
# .'       `( \                                    __  __ ___ ____      _                          #
#         @/            ')   ,--,__,-"            |  \/  |_ _|  _ \    / \                         #
#         /        /      \ /     /   _/          | |\/| || || |_) |  / _ \                        #
#       __|           ,   |/         /            | |  | || ||  _ <  / ___ \                       #
#     .~  `\   / \ ,  |   /                       |_|  |_|___|_| \_\/_/   \_\                      #
#   .~      `\    `  /  _/   _/                   ---------------------------                      #
# .~          `\  ~~`__/    /                     |_| _  |  _  o __  |  o                          #
# ~             `--'/                             | |(/_ | _>  | | | |< |                          #
#              /   /    /                                                                          #
#             /  /'    /jgs                                                                        #
#                                                       http://github.com/topihovinen/mirahelsinki #
# ------------------------------------------------------------------------------------------------ #

# -------------------------------------------0. PREFACE------------------------------------------- #
{
  # Introduce needed libraries for R
  {
    library("ggplot2")
    library("ggrepel")
    library("digest")
    library("reshape2")
    library("Hmisc")
    library("devtools")
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
    library("extrafont")
    library("stringr")
    source("multiplot.R")
    source("permira.test.R") # Functions for permutation tests of our data,
  #                            provided in another supplementary
  }
  
  # Import data
  {
    # Import the "targeted" metabolomics data and choose the factor (group) order for figures and color palette
    {
      miradata <- read.csv("HUSLAB Data final.txt", sep='\t', dec=',', header=T)
      head(miradata)
      
      # This block discards all subjects with no blood sample
      {
        miradata <- miradata[!is.na(miradata$Bage),]
      }
      
      # Build additional information that can be yielded from the imported data
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
        {
          mirakrea <- read.csv("HUSLAB Data Krea.txt", sep='\t', dec=',', header=T)
          miradata$U.Krea <- rep(NA, length(miradata[,1]))
          for(i in 1:length(miradata[,1])) {
            if(miradata$ID[i] %in% mirakrea$ID) {
              miradata$U.Krea[i] <- mirakrea$U.Krea.mmol.l[mirakrea$ID==miradata$ID[i]]
            }
          }
          miradata$U.I.standard <- miradata$U.I/miradata$U.Krea # Urine iodine standardized to urine creatinine
          rm(mirakrea)
        }
      }
      # Add group labels. Group 5 is used for the article, as traditionally vegetarians include also pesco-vegetarians
      {
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
      }
      # Convert cholesterol absorption and synthesis biomarkers to ug/mg, original data was *100ug/mg
      {
        miradata[,45:52] = miradata[,45:52]/100
      }
    }
    
    # Import the "untargeted" metabolomics data for amino acid and fatty acid analysis
    {
      # Import and preprocess metabolomics dataset, include also a few metabolites from targeted metabolomics
      # for correlation analysis. Note that the original file includes data from both extraction methods for
      # metabolomics optimization and lipidomics optimization (see Article Appendix for details) and the dataset
      # is divided accordingly to miraLIPOmetab and miraMETAmetab further below
      {
        mirametab <- t(as.matrix(read.csv("MIRA Metabolomics new.txt", header = T, sep = '\t', dec=',', stringsAsFactors = F)))
        mirametab <- mirametab[-c(1:3),] #Remove IonMz and IonActive -rows
        for(i in seq(1,length(mirametab[,1]),2)) { # Calculate the means of each duplicate measurement
          mirametab[i,] = colMeans(matrix(c(mirametab[i,], mirametab[i+1,]), nrow = 2, byrow = TRUE))
        }
        mirametab <- mirametab[-(seq(2,length(mirametab[,1]),2)),]
        mirametab <- as.data.frame(mirametab, header = T)
        mirametab$ID <- str_sub(rownames(mirametab),-3)
        mirametab$Extraction.method <- str_sub(rownames(mirametab),14,17)
        rownames(mirametab) <- c()
        mirametab <- mirametab[mirametab$ID %in% miradata$ID,]
        mirametab$Group5 <- rep(as.factor("Control"), length(mirametab$ID))
        levels(mirametab$Group5) <- c("Control", "Vegetarian", "Vegan")
        for(i in 1:length(mirametab$ID)) {
          if(mirametab$ID[i] %in% miradata$ID) {
            mirametab$Group5[i] <- factor(miradata$Group5[miradata$ID == mirametab$ID[i]])
          }
        }
        #Transfer wanted data from HUSLAB to metabolomics matrix
        {
          mirametab$Total.Cholesterol = rep(0,length(mirametab[,1]))
          mirametab$Folate.Eryt = rep(0,length(mirametab[,1]))
          mirametab$Triglycerides = rep(0,length(mirametab[,1]))
          mirametab$Transthyretin = rep(0,length(mirametab[,1]))
          for(i in 1:length(mirametab[,1])) {
            mirametab$Total.Cholesterol[i] <- miradata$fP.Kol..mmol.l.[miradata$ID==mirametab$ID[i]]
            mirametab$Folate.Eryt[i] <- miradata$fE.Folaat..nmol.l.[miradata$ID==mirametab$ID[i]]
            mirametab$Triglycerides[i] <- miradata$fP.Trigly..mmol.l.[miradata$ID==mirametab$ID[i]]
            mirametab$Transthyretin[i] <- miradata$S.Prealb..mg.l.[miradata$ID==mirametab$ID[i]]
          }
        }
      }
      
      # Create a variable for metabolite names (as they were provided from the equipment in a separate dataframe)
      {
        mirametab_names <- read.csv("MIRA Metabolomics new_metanames.txt", header = T, sep = '\t', dec = ',', stringsAsFactors = F)
        mirametab_names <- data.frame(ionIdx = mirametab_names$ionIdx, ionMz = mirametab_names$ionMz, name = mirametab_names$label..bona.fide., hmdb = mirametab_names$HMDB.ids, stringsAsFactors = F)
      }
    }
  }
}

# -------------------------1. GENERAL CHARACTERISTICS AND ANTHROPOMETRICS------------------------- #
{
  # For height and weight analysis, we first create a data frame for height and weight
  # averages and corresponding 2SD curves of Finnish children (averages and 2SD data
  # based on Saari A, et al. Ann Med 2010.) Approximate and discrete scale is represented
  # in this public code due to copyright reasons of the full original data from Finnish
  # Institute for Health and Welfare (THL)
  {
    lengthmean = data.frame(data = matrix(c(1,72.25,76.75,81,
                                            2,82.75,88.5,94.25,
                                            3,89.75,96.5,103.25,
                                            4,96.25,104,111.5,
                                            5,102.5,111,119.5,
                                            6,108.75,118.25,127.25,
                                            7,114.25,124.25,134.25,
                                            8,119.75,130.25,140.75),
                                          ncol=4, nrow=8, byrow = T))
    colnames(lengthmean) <- c("age","ll","median","ul") #Rename the columns for clarity
    bmimean <- data.frame(data = matrix(c(1,14.55,16.25,17.5,
                                          2,14.4,16.1,17.4,
                                          3,14.2,15.95,17.25,
                                          4,14.0,15.8,17.1,
                                          5,13.8,15.6,17.1,
                                          6,13.7,15.65,17.25,
                                          7,13.75,15.8,17.7,
                                          8,13.8,16.15,18.3),
                                        ncol=4, nrow=8, byrow = T))
    colnames(bmimean) <- c("age","ll","median","ul")
  }
  # For MUAC, import data from WHO z-scores for upper arm circumference
  {
    muac_boys <- read.csv("WHO_MUAC_z-scores_boys.txt", sep = '\t', dec = '.', header = T)
    muac_boys_ext <- read.csv("WHO_MUAC_z-scores_boys_extension.txt", sep = '\t', dec = ',', header = T)
    muac_boys <- rbind(muac_boys, muac_boys_ext) # Combine WHO 0-5 years old MUAC with the extension provided by Mramba et al. 2017
    muac_girls <- read.csv("WHO_MUAC_z-scores_girls.txt", sep = '\t', dec = '.', header = T)
    muac_girls_ext <- read.csv("WHO_MUAC_z-scores_girls_extension.txt", sep = '\t', dec = ',', header = T)
    muac_girls <- rbind(muac_girls, muac_girls_ext)
    {
      rm(muac_boys_ext)
      rm(muac_girls_ext)
    }
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
    # Include MUAC_z also in mirametab -dataset (for correlation analysis of MUAC and essential amino acids)
    mirametab$MUAC_z = rep(0,length(mirametab[,1]))
    for(i in 1:length(mirametab[,1])) {
      mirametab$MUAC_z[i] <- miradata$MUAC_z[miradata$ID==mirametab$ID[i]]
    }
    # Divide metabolomics and lipidomics extractions to different data frames now that all the needed
    # data is included in mirametab
    miraMETAmetab <- subset(mirametab, Extraction.method == "META")
    miraLIPOmetab <- subset(mirametab, Extraction.method == "LIPO")
    rm(mirametab) #The original mirametab-frame is not needed after this
  }
  
  # Age vs. height plot
  {
    hplot = ggplot(data=miradata, aes(x = AntAge, y = Height, colour = Group5, shape = Group5, size = Group5)) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          labels = c("Omnivore", "Vegetarian", "Vegan"),
                          values=c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Omnivore", "Vegetarian", "Vegan"),
                         guide = guide_legend(title.position = "right")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Omnivore", "Vegetarian", "Vegan"),
                        guide = guide_legend(title.position = "right")) +
      geom_line(inherit.aes = FALSE, data = lengthmean, aes(x = age, y = median), colour = "black", show.legend = F, linetype = "dashed", size = 0.5) +
      annotate(geom = "text", label = "Finnish\nmedian", x = 7.6, y = 128, angle = 20, size = 2) +
      geom_ribbon(inherit.aes = FALSE, data = lengthmean,
                  aes(x = age, ymin = ll, ymax = ul),
                  alpha = 0.15, show.legend = F) +
      geom_point() +
      geom_point(data = subset(miradata, Group5 == "Vegan")) +
      geom_point(data = subset(miradata, Group5 == "Vegetarian")) +
      labs(x = "Age (years)", y = "Height (cm)") +
      scale_x_continuous(limits = c(1, 8), expand = c(0, 0)) +
      theme(text = element_text(family = "Helvetica"),
            legend.position = "top",
            legend.text = element_text(size = 8),
            legend.text.align = 0,
            legend.title = element_blank(),
            legend.direction = "horizontal",
            legend.key.size = unit(0.5, "lines"),
            legend.spacing.x = unit(0, "in"),
            legend.key = element_rect(fill=NA),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 8, colour = "black"),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank())
    print(hplot)
    ggsave(file = "manusc/figures/Lancet/Age vs. height.pdf", plot = hplot, device = "pdf", width = 6, height = 4)
  }
  
  # Age vs. BMI plot
  {
    bplot = ggplot(data=miradata, aes(x = AntAge, y = BMI, colour = Group5, shape = Group5, size = Group5)) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          labels = c("Omnivore", "Vegetarian", "Vegan"),
                          values=c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Omnivore", "Vegetarian", "Vegan"),
                         guide = guide_legend(title.position = "right")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Omnivore", "Vegetarian", "Vegan"),
                        guide = guide_legend(title.position = "right")) +
      geom_smooth(inherit.aes = FALSE, data = bmimean, aes(x = age, y = median), colour = "black", show.legend = F, method="loess", se = FALSE, linetype = "dashed", size = 0.5) +
      annotate(geom = "text", label = "Finnish\nmedian", x = 7.6, y = 16.0, angle = 9, size = 2) +
      geom_ribbon(inherit.aes = FALSE, data = bmimean,
                  aes(x = age, ymin = ll, ymax = ul),
                  alpha = 0.15, show.legend = F) +
      geom_point() +
      geom_point(data = subset(miradata, Group5 == "Vegan")) +
      geom_point(data = subset(miradata, Group5 == "Vegetarian")) +
      labs(x = "Age (years)", y = expression(paste(BMI (kg/m^2)))) +
      scale_x_continuous(limits = c(1, 8), expand = c(0, 0)) +
      theme(text = element_text(family = "Helvetica"),
            legend.position = "top",
            legend.text = element_text(size = 8),
            legend.text.align = 0,
            legend.title = element_blank(),
            legend.direction = "horizontal",
            legend.key.size = unit(0.5, "lines"),
            legend.spacing.x = unit(0, "in"),
            legend.key = element_rect(fill=NA),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 8, colour = "black"),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank())
    print(bplot)
    ggsave(file = "manusc/figures/Lancet/Age vs. BMI.pdf", plot = bplot, device = "pdf", width = 6, height = 4)
  }
}

# -----------------------------------2. HIERARCHICAL CLUSTERING----------------------------------- #
{
  clustdata <- c(14:15,17:22,25:33,35,37:39,41:43,45:52,54:70) # Choose all originally independently measured 49 metabolites
  miradata_clust <- miradata[,c(1,104,clustdata)] # Column 1 is ID and 104 is groups for labels and colors of the plot
  # Create color scale for plotting
  {
    miradata_clust$col <- as.character(miradata_clust$Group5)
    miradata_clust$col[miradata_clust$col == "Control"] <- "#4D6B73FF"
    miradata_clust$col[miradata_clust$col == "Vegetarian"] <- "#009FC3FF"
    miradata_clust$col[miradata_clust$col == "Vegan"] <- "#6CB43FFF"
  }
  # Calculate distance object and plot it using hclust
  {
    clust_dist <- miradata_clust[,c(3:(length(clustdata) + 2))] %>%
      scale %>%
      dist(method = "euclidean") %>%
      hclust(method = "ward.D2") %>%
      as.dendrogram %>%
      set("leaves_pch", 19) %>%
      set("leaves_cex", 2) %>%
      set("branches_lwd", 0.5) %>%
      set("branches_k_color", k = 2, value = c("#4D6B73FF","#009FC3FF")) %>%
      set("branches_k_color", k = 3, value = c("#4D6B73FF","#009FC3FF","#6CB43FFF"))
    clust_dist <- set(clust_dist, what = "leaves_col", value = miradata_clust$col[order(order(as.numeric(labels(clust_dist))))])
    #labels(clust_dist) <- as.character(miradata_clust$Group5[order(order(as.numeric(labels(clust_dist))))])
    labels(clust_dist) <- rep("",length(miradata_clust[,1]))
    # Plot using ggplot
    # Legend is tough to create so let's do it manually to img-variable
    img <- readPNG("manusc/figures/clustlegend.png")
    img <- rasterGrob(img, interpolate=TRUE)
    clustplot <- ggplot(data = as.ggdend(clust_dist)) + 
      theme_light() +
      labs(x = "", y = "Dissimilarity (arbitrary units)") +
      annotate("rect", xmin = 0.5, xmax = 30.5, ymin = -2, ymax = 1, fill = "#4D6B73FF", alpha = 0.2) +
      annotate("text", x = 15.5, y = -1, label = "A") +
      annotate("rect", xmin = 30.51, xmax = 35.5, ymin = -2, ymax = 1, fill = "#009FC3FF", alpha = 0.2) +
      annotate("text", x = 33, y = -1, label = "B") +
      annotate("rect", xmin = 35.51, xmax = 40.5, ymin = -2, ymax = 1, fill = "#6CB43FFF", alpha = 0.2) +
      annotate("text", x = 38, y = -1, label = "C") +
      annotation_custom(img, xmin = 13, xmax = 29, ymin = 20.2, ymax = 21.8) +
      theme(text = element_text(colour = "black", size = 8, family = "Helvetica"),
            plot.title = element_text(size = 8),
            axis.title.y = element_text(size = 8, colour = "black"),
            axis.ticks.y = element_line(colour = "black"),
            axis.ticks.x = element_blank(),
            axis.line.y = element_line(colour = "black"),
            axis.line.x = element_blank(),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.text.x = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank())
    
    clust_unt_dist <- miradata_unt_clust[,1:872] %>%
      scale %>%
      dist(method = "euclidean") %>%
      hclust(method = "ward.D2") %>%
      as.dendrogram %>%
      set("leaves_pch", 19) %>%
      set("leaves_cex", 2) %>%
      set("branches_lwd", 0.5) %>%
      set("branches_k_color", k = 2, value = c("#4D6B73FF","#4D6B73FF")) %>%
      set("branches_k_color", k = 4, value = c("#4D6B73FF","#4D6B73FF","#4D6B73FF","#4D6B73FF"))
    clust_unt_dist <- set(clust_unt_dist, what = "leaves_col", value = miradata_unt_clust$col[order(order(as.numeric(labels(clust_unt_dist))))])
    labels(clust_unt_dist) <- rep("",length(miradata_unt_clust[,1]))
    clust_unt_plot <- ggplot(data = as.ggdend(clust_unt_dist)) + 
      theme_light() +
      labs(x = "", y = "Dissimilarity (arbitrary units)") +
      annotate("rect", xmin = 0.5, xmax = 14.5, ymin = -5, ymax = 66, fill = "#4D6B73FF", alpha = 0.2) +
      annotate("text", x = 7.5, y = -2, label = "A") +
      annotate("rect", xmin = 14.51, xmax = 23.5, ymin = -5, ymax = 66, fill = "#6CB43FFF", alpha = 0.2) +
      annotate("text", x = 19, y = -2, label = "B") +
      annotate("text", x = 24.5, y = -2, label = "C") +
      annotate("rect", xmin = 25.51, xmax = 40.5, ymin = -5, ymax = 66, fill = "#4D6B73FF", alpha = 0.2) +
      annotate("text", x = 32.5, y = -2, label = "D") +
      annotation_custom(img, xmin = 11, xmax = 29, ymin = 77, ymax = 80.5) +
      geom_hline(aes(yintercept=66, linetype="dotted", alpha = .5)) +
      theme(text = element_text(colour = "black", size = 8, family = "Helvetica"),
            plot.title = element_text(size = 8),
            axis.title.y = element_text(size = 8, colour = "black"),
            axis.ticks.y = element_line(colour = "black"),
            axis.ticks.x = element_blank(),
            axis.line.y = element_line(colour = "black"),
            axis.line.x = element_blank(),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.text.x = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none")
    
    ggsave(file="manusc/figures/Figure EV1.pdf",
           plot=clust_unt_plot,
           width = 6, height = 4, dpi = 1200, units = "in", device = "pdf")
  }
}

# -----------------------------------3. SERUM BIOMARKER ANALYSES---------------------------------- #
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
                  "Total lithocholic acids", "Total deoxycholic acids", "Total ursodeoxycholic acids", "Vitamin A status",
                  "Urine creatinine", "Iodine in urine (adjusted to creatinine)", rep("", 5))
    
    miraunits = c(rep("", 13), "nmol/l","nmol/l","nmol/l","µmol/l","mg/l","mmol/l","mmol/l","mmol/l","mmol/l","", "",
                  "*10^9/l","*10^12/l","g/l", "%", "fl","%","pg","g/l","*10^9/l","pmol/l","nmol/l","%",
                  "mmol/l", "µg/l", "mg/l", "mg/kg", "RE (µmol/l)", "mg/l", "g/l", "µg/l",
                  rep("µg/mg of TC",8), "mmol/l", rep("µmol/l",17), "", "", "", "", "µmol/l", "µmol/l",
                  "", "", rep("µmol/l", 11), "", "", "", "", "", "", rep("µmol/l", 5), "", "mmol/l",
                  "ug/l per mmol/l creatinine", rep("",5))
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
    # Calculate p-values with permutation tests (ad-hoc permira.test function), comparing probability indexes and all possible pairs
    {
      rows.to.measure <- c(14:15,17:23,25:33,35,37:52,54:70,75:101,103) # Choose which metabolites are to be tested
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
      set1 <- c(14:15,17:33,35,37,101,103) #Assign ID's of targeted metabolites that were tested
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
    #for (i in c(14:15, 17:33, 35:52, 54:100)) { #All figures
    #for (i in c(14:15, 18, 20:23, 35, 41, 44:51, 61, 69, 78, 93, 94)) { #Figures in the publication
      plot103 <- ggplot(data = data, aes(x=Group5, y=data[,103])) +
        theme_light() +
        scale_colour_manual(name = "Dietary group",
                            values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF", "#6CB43FFF"),
                            guide = guide_legend(title.position = "right")) +
        scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                           labels = c("Control", "Vegetarian", "Vegan")) +
        scale_size_manual(name = "Dietary group", values=c(1,1,1),
                          labels = c("Control", "Vegetarian", "Vegan")) +
        scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
        geom_point(position = position_jitter(width = .1, height = 0),
                   aes(colour = Group5, shape = Group5, size = Group5)) +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
        stat_summary(geom = "errorbar", fun.y = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
        #labs(y = miraunits[39], x = "", title = miranames[39]) +
        labs(y = miraunits[103], x = "", title = "Iodine/creatinine\n(urine)") +
        geom_signif(comparisons = list(c("Control", "Vegan")),
                    test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) + 
        # following geom_hlines are for printing possible limits for deficiency
        # ("dashed") or insufficiency ("dotted")
        #geom_hline(aes(yintercept=1.17, linetype="dotted", alpha = .5)) + #RBP insufficiency
        #geom_hline(aes(yintercept=0.83, linetype="dashed", alpha = .5)) + #RBP deficiency
        #geom_hline(aes(yintercept=50, linetype="dashed", alpha = .5)) + #VitD deficiency
        #geom_hline(aes(yintercept=30, linetype="dashed", alpha = .5)) + #VitB12 deficiency
        #geom_hline(aes(yintercept=0, linetype="dashed", alpha = .5)) + #Iron (BIS) and VitA deficiency
        coord_cartesian(ylim = c(min(data[,103], na.rm=TRUE)*0.95,
                                 max(data[,103], na.rm=TRUE)*1.2)) +
        theme(plot.title = element_text(hjust = 0.5, size = 10),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
              axis.text.y = element_text(size = 10, colour = "black"),
              axis.title.y = element_text(size = 10),
              axis.ticks = element_line(colour = "black"),
              axis.line = element_line(colour = "black"),
              panel.border = element_blank(),
              legend.position = "none") +
        guides(colour=F, linetype=F, size=F, shape=F)
      #print to file
      filename = paste("manusc/figures/Figure subplots/BH adjusted/", miranames[103], ".png", sep="")
      ggsave(file=filename, plot=plot, width=2.5, height=3)
      message("Printing figure ", i) # Check up during the loop that something is happening
      Sys.sleep(0.01) # To make sure the progress message from previous line is printed each time.
    #}
  }
}

# -------------------------------4. THREE-GROUP STATISTICAL ANALYSES------------------------------ #
{
  # Kruskal-Wallis for anthropometrics
  {
    height_krusk <- kruskal.test(miradata$Height_z, g = miradata$Group5)
    bmi_krusk <- kruskal.test(miradata$BMI_z, g = miradata$Group5)
    muac_krusk <- kruskal.test(miradata$MUAC_z, g = miradata$Group5)
    krusk.wall.data.anthro <- data.frame(variable = c("Height z-score", "BMI SDS-score", "MUAC z-score"),
                                         p = c(height_krusk[[3]], bmi_krusk[[3]], muac_krusk[[3]]))
    rm(height_krusk, bmi_krusk, muac_krusk)
  }
  # Kruskal-Wallis for targeted metabolomics
  {
    alldata <- c(41,18,101,35,14,15,17,38,39,103,19,22,20,21,23,37,25:33,42,43,45,49,50,51,46,47,48,52,74,70,69,59,57,61,56,68,64,65,62,66,60,67,55,63,54,58,93,94,78) # Select all targeted independent metabolites in the wanted order
    krusk.wall.results <- sapply(miradata[,alldata], kruskal.test, g = miradata$Group5)
    krusk.wall.stats <- as.vector(unlist(krusk.wall.results[seq(1,length(krusk.wall.results),5)]))
    krusk.wall.p <- as.vector(unlist(krusk.wall.results[seq(3,length(krusk.wall.results),5)]))
    krusk.wall.p.adj <- rep(NA, length(krusk.wall.p))
    krusk.wall.p.adj[which(alldata %in% c(set1))] <- p.adjust(krusk.wall.p[which(alldata %in% c(set1))], method = "BH")
    krusk.wall.p.adj[which(alldata %in% c(set2))] <- p.adjust(krusk.wall.p[which(alldata %in% c(set2))], method = "BH")
    krusk.wall.p.adj[which(alldata %in% c(set3))] <- p.adjust(krusk.wall.p[which(alldata %in% c(set3))], method = "BH")
    krusk.wall.p.adj[which(alldata %in% c(set4))] <- p.adjust(krusk.wall.p[which(alldata %in% c(set4))], method = "BH")
    # Save the data
    krusk.wall.data.met <- data.frame(metabolite = colnames(miradata)[alldata],
                                      p = krusk.wall.p,
                                      p.adj = krusk.wall.p.adj)
    rm(krusk.wall.results, krusk.wall.stats, krusk.wall.p, krusk.wall.p.adj) # Remove intermediate variables
  }
  # Kruskal-Wallis for nutrient intake data
  {
    alldata_nutrients <- c(24,90,88,91,92,93,94,58,59,109,60,110,61,32,35,37,95,97,40,98,44,99,56,100,46,101,47,102,57,103,39,111,42,112,43,96,41,45,48,104,50,52,105,49,106,51,108,53,54,107,55) # Select daily intake columns of all nutrients in intake data in the wanted order
    krusk.wall.results.nutri <- sapply(miradata_ravinto[,alldata_nutrients], kruskal.test, g = miradata_ravinto$Group5)
    krusk.wall.stats.nutri <- as.vector(unlist(krusk.wall.results.nutri[seq(1,length(krusk.wall.results.nutri),5)]))
    krusk.wall.p.nutri <- as.vector(unlist(krusk.wall.results.nutri[seq(3,length(krusk.wall.results.nutri),5)]))
    krusk.wall.p.adj.nutri <- p.adjust(krusk.wall.p.nutri, method = "BH")
    # Print the data
    krusk.wall.data.nutri <- data.frame(nutrient = colnames(miradata_ravinto)[alldata_nutrients],
                                        p = krusk.wall.p.nutri,
                                        p.adj = krusk.wall.p.adj.nutri)
    rm(krusk.wall.results.nutri, krusk.wall.stats.nutri, krusk.wall.p.nutri, krusk.wall.p.adj.nutri) # Remove intermediate variables
  }
}

# ---------------------------------5. UNTARGETED METABOLOMICS DATA-------------------------------- #
{
  ## ------------ Amino acid analysis -------------- ##
  {
    # Find amino acids
    aminonames_orig <- c("Glycine", "Alanine", "Serine", "Proline", "Valine; Betaine", "Threonine", "Cys", "(Iso)Leucine", "Asparagine", "Aspartate", "Glutamine", "Lysine", "Glutamate", "Methionine", "Histidine", "Phenylalanine", "Arginine", "Tyrosine", "Tryptophan")
    aminoacid_ID <- rep(0,length(aminonames_orig))
    for(i in 1:length(aminonames_orig)) {
      aminoacid_ID[i] <- mirametab_names$ionIdx[which(mirametab_names$name == aminonames_orig[i])[1]]
    }
    
    # Calculate differences of means and their 95% confidence interval lengths
    # assuming heteroscedasticity (Behrens-Fisher approach, originally from 1935).
    # Equation for degrees of freedom, v, is also called Welch-Satterthwaite equation.
    {
      AAdoms <- rep(0,length(aminoacid_ID))
      AAcilens <- rep(0, length(aminoacid_ID))
      AAfc <- rep(0,length(aminoacid_ID))
      for(i in 1:length(aminoacid_ID)) {
        n1 = length(miraMETAmetab$ID[miraMETAmetab$Group5=="Vegan"])
        n2 = length(miraMETAmetab$ID[miraMETAmetab$Group5=="Control"])
        m1 = mean(miraMETAmetab[miraMETAmetab$Group5=="Vegan",aminoacid_ID[i]])
        m2 = mean(miraMETAmetab[miraMETAmetab$Group5=="Control",aminoacid_ID[i]])
        # Means of peak intensities normalized to control group mean:
        norm_m1 = mean(miraMETAmetab[miraMETAmetab$Group5=="Vegan",aminoacid_ID[i]]/m2)
        norm_m2 = mean(miraMETAmetab[miraMETAmetab$Group5=="Control",aminoacid_ID[i]]/m2) #This should yield 1 every time.
        s1 = sd(miraMETAmetab[miraMETAmetab$Group5=="Vegan",aminoacid_ID[i]]/m2)
        s2 = sd(miraMETAmetab[miraMETAmetab$Group5=="Control",aminoacid_ID[i]]/m2)
        v <- (s1^2/n1+s2^2/n2)^2/(s1^4/(n1^2*(n1-1)) + s2^4/(n2^2*(n2-1)))
        #Calculate differences of normalized means (normalized to control group mean)
        AAdoms[i] <- norm_m1 - norm_m2
        AAcilens[i] <- qt(0.975,v)*sqrt(s1^2/n1 + s2^2/n2)
        AAfc[i] <- norm_m1
      }
      # Provide proper labelling for the plot
      {
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
      # This combined dataframe will be used for plotting
      AAdata <- data.frame(dom = AAdoms,
                           cilen = AAcilens,
                           fc = AAfc,
                           label = aminoacid_labels,
                           type = type)
    }
    
    # Amino acid plot
    {
      aminoplot <- ggplot(data = AAdata) +
        theme_light() +
        scale_fill_manual(name = "",
                            values=c("#B30437FF", "#9170B4FF")) +
        #non-log2 version:
        geom_bar(aes(x=label, y=dom, fill=type), stat="identity") +
        geom_errorbar(aes(x=label, ymin=dom-cilen, ymax=dom+cilen),
                     width=0.4, colour="black", size=0.5, alpha = 1) +
        facet_wrap(~type, scales = "free_x") +
        labs(y = "Average normalized peak intensity\nin vegans compared to controls", x = "", title = "") +
        # #log2 version:
        # geom_bar(aes(x=label, y=log2(fc), fill=type), stat="identity") +
        # geom_errorbar(aes(x=label, ymin=log2(fc-cilen), ymax=log2(fc+cilen)),
        #               width=0.4, colour="black", size=0.5, alpha = 1) +
        # facet_wrap(~type, scales = "free_x") +
        # labs(y = "log2 FC in vegans compared to controls", x = "", title = "") +
        theme(text = element_text(size = 8, family = "Helvetica"),
              plot.title = element_text(hjust = 0.5, size = 8),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
              axis.text.y = element_text(size = 8, colour = "black"),
              axis.title.y = element_text(size = 8),
              legend.position = "none",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              strip.background = element_blank(),
              strip.text = element_text(size = 8, colour = "black"),
              axis.ticks = element_line(colour = "black"),
              axis.line = element_line(colour = "black")) +
        guides(colour=F, linetype=F, size=F, shape=F)
      ggsave(file="manusc/figures/Med/aminoacids.png",
             plot = aminoplot, 
             width = 5, 
             height = 3, 
             device = "png")
    }
  }
  
  ## ------------ Fatty acid analysis -------------- ##
  {
    # To better analyze untargeted metabolomics, normalize the data by dividing each measurement by the mean of omnivore group for the corresponding metabolite
    {
      miraLIPO_norm <- miraLIPOmetab
      for(i in 1:length(miraLIPO_norm$ID)) {
        if(miraLIPO_norm$Group5[i] == "Vegetarian") {
          miraLIPO_norm$ID2[i] <- paste("Vegetarian_", order(order(as.numeric(miraLIPO_norm$ID)))[i], sep = "")
        } else if(miraLIPO_norm$Group5[i] == "Vegan") {
          miraLIPO_norm$ID2[i] <- paste("Vegan_", order(order(as.numeric(miraLIPO_norm$ID)))[i], sep = "")
        } else {
          miraLIPO_norm$ID2[i] <- paste("Omnivore_", order(order(as.numeric(miraLIPO_norm$ID)))[i], sep = "")
        }
        for(j in 1:872)
          miraLIPO_norm[i,j] <- miraLIPO_norm[i,j] / mean(miraLIPOmetab[(miraLIPOmetab$Group5 == "Control"),j])
      }
    }
    # For combined analysis (and plotting) of triglycerides, it's convenient to combine them into one tall dataframe
    {
      miraLIPO_norm_MCFA <- melt(miraLIPO_norm[,c("Group5", "V699","V746","V759")], id.vars="Group5")
      miraLIPO_norm_VLCFA <- melt(miraLIPO_norm[,c("Group5", "V834","V845","V858","V861","V855","V869")], id.vars="Group5")
    }
    # This is a generic code for boxplots used in the manuscript. The labels and metabolite indexes should be changed manually.
    # To plot a single metabolite, use miraLIPO_norm, and to plot combined triglycerides you may use the corresponding dataframes above
    fattyacid_plot <- ggplot(data = miraLIPO_norm_VLCFA, aes(x=Group5,
                                                             y=log2(miraLIPO_norm_VLCFA[,3]))) + 
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right"),
                          labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_fill_manual(name = "Dietary group",
                        values = c("#4D6B7380", "#009FC380", "#6CB43F80"),
                        guide = guide_legend(title.position = "right"),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(position = position_jitter(width = .1, height = 0), aes(color = Group5,
                                                                         shape = Group5,
                                                                         size = Group5)) + 
      geom_boxplot(aes(fill = Group5,
                       color = Group5),
                   outlier.alpha = 0,
                   coef = 10,
                   lwd = .25,
                   width = .75) +
      labs(y = "log2 FC to\nomnivore mean", x = "LCFA/VLCFA triglycerides combined\nC50:0-o, C53:4, C58:0-o, C58:4, C58:9-o, C62:10", title = "") +
      #geom_signif(comparisons = list(c("Control", "DC-Vegan")), map_signif_level = T, margin_top = 0.05) +
      #geom_signif(comparisons = list(c("Control", "Vegan")), map_signif_level = T, margin_top = 0.15) + 
      #coord_cartesian(ylim = c(min(miraLIPO_norm[,i], na.rm=TRUE)*0.95, max(miraLIPO_norm[,i], na.rm=TRUE)*1.15)) +
      coord_cartesian(ylim = c(-3,2.5)) +
      theme(plot.title = element_text(hjust = 0, size = 15),
            #axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line = element_line(colour = "black"),
            #axis.text.y = element_text(size = 7, colour = "black"),
            #axis.title.y = element_text(size = 7),
            #axis.ticks.y = element_line(colour = "black"),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.title.x = element_text(size = 7),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none",
            #legend.position = "right",
            legend.background = element_blank(),
            legend.box.background = element_rect(colour = "black"),
            legend.text = element_text(size = 7),
            legend.title = element_blank()) +
      guides(linetype=F, size=F)
  }
}

# ---------------------------------6. BUILD AND PRINT THE FIGURES--------------------------------- #
{
  # Figure 1
  {
    # Create titles
    df <- data.frame()
    titleplot <- ggplot(df) +
      theme_void() +
      geom_point() +
      xlim(0, 10) +
      ylim(0.5, 1.5) +
      annotate("text", x = 5, y = 1,
               label = 'bold("A                                             ")',
               size = 5, parse = TRUE, family = "Helvetica")
    titleplot2 <- ggplot(df) +
      theme_void() +
      geom_point() +
      xlim(0, 10) +
      ylim(0.5, 1.5) +
      annotate("text", x = 5, y = 1,
               label = 'bold("B                                           ")',
               size = 5, parse = TRUE, family = "Helvetica")
    titleplot3 <- ggplot(df) +
      theme_void() +
      geom_point() +
      xlim(0, 10) +
      ylim(0.5, 1.5) +
      annotate("text", x = 5, y = 1,
               label = 'bold("C                                                                                                              ")',
               size = 5, parse = TRUE, family = "Helvetica")
    # Assign the layout matrix (how the plots are arranged)
    ployout = matrix(c(1,1,1,1,1,2,2,2,2,2,
                       3,3,3,3,3,4,4,4,4,4,
                       3,3,3,3,3,4,4,4,4,4,
                       3,3,3,3,3,4,4,4,4,4,
                       3,3,3,3,3,4,4,4,4,4,
                       3,3,3,3,3,4,4,4,4,4,
                       3,3,3,3,3,4,4,4,4,4,
                       3,3,3,3,3,4,4,4,4,4,
                       3,3,3,3,3,4,4,4,4,4,
                       3,3,3,3,3,4,4,4,4,4,
                       5,5,5,5,5,5,5,5,5,5,
                       6,6,7,7,8,8,9,9,10,10,
                       6,6,7,7,8,8,9,9,10,10,
                       6,6,7,7,8,8,9,9,10,10,
                       6,6,7,7,8,8,9,9,10,10,
                       6,6,7,7,8,8,9,9,10,10,
                       6,6,7,7,8,8,9,9,10,10,
                       11,11,12,12,13,13,14,14,15,15,
                       11,11,12,12,13,13,14,14,15,15,
                       11,11,12,12,13,13,14,14,15,15,
                       11,11,12,12,13,13,14,14,15,15,
                       11,11,12,12,13,13,14,14,15,15,
                       11,11,12,12,13,13,14,14,15,15), ncol = 10, byrow = T)
    # Bring the plot together and print it to pdf
    ggsave(file="manusc/figures/Figure 1.pdf",
           plot=multiplot(titleplot, titleplot2,
                          hplot, bplot,
                          titleplot3,
                          plotnene, plotnpro, plotnfat, plotnsafa, plotnchol,
                          plotnvita, plotnvitd, plotnfol, plotnzn, plotnfe,
                          layout = ployout), width = 7, height = 8, dpi = 1200, unit = "in", device = "pdf")
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
      annotate("text", x = 5, y = 1,
               label = 'bold("A                                                                                                              ")',
               size = 5, parse = TRUE, family = "Helvetica")
    titleplot2 <- ggplot(df) +
      theme_void() +
      geom_point() +
      xlim(0, 10) +
      ylim(0.5, 1.5) +
      annotate("text", x = 5, y = 1,
               label = 'bold("B                                                                                                              ")',
               size = 5, parse = TRUE, family = "Helvetica")
    subtitleplot <- ggplot(df) +
      theme_void() +
      geom_point() +
      xlim(0, 10) +
      ylim(0.5, 1.5) +
      annotate("text", x = 5, y = 1.1,
               label = "                                    Cholesterol absorption                                  ",
               size = 4, family = "Helvetica") +
      geom_hline(yintercept = 0.7)
    subtitleplot2 <- ggplot(df) +
      theme_void() +
      geom_point() +
      xlim(0, 10) +
      ylim(0.5, 1.5) +
      annotate("text", x = 5, y = 1.1,
               label = "                                  Cholesterol biosynthesis                                 ",
               size = 4, family = "Helvetica") +
      geom_hline(yintercept = 0.7)
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
                       rep(c(2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,
                       2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,
                       2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,
                       2,2,2,3,3,3,4,4,4,5,5,5,6,6,6),2),
                       rep(c(7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,
                       7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,
                       7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,
                       7,7,7,8,8,8,9,9,9,10,10,10,11,11,11),2),
                       12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,
                       rep(c(13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,
                       13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,
                       13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,
                       13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,
                       13,13,13,14,14,14,14,14,14,14,14,14,14,14,14),2),
                       15,15,15,16,16,16,16,16,16,16,16,16,16,16,16,
                       rep(c(15,15,15,17,17,17,18,18,18,19,19,19,20,20,20,
                       15,15,15,17,17,17,18,18,18,19,19,19,20,20,20,
                       15,15,15,17,17,17,18,18,18,19,19,19,20,20,20,
                       15,15,15,17,17,17,18,18,18,19,19,19,20,20,20),2),
                       21,21,21,22,22,22,22,22,22,22,22,22,22,22,22,
                       rep(c(21,21,21,23,23,23,24,24,24,25,25,25,26,26,26,
                       21,21,21,23,23,23,24,24,24,25,25,25,26,26,26,
                       21,21,21,23,23,23,24,24,24,25,25,25,26,26,26,
                       21,21,21,23,23,23,24,24,24,25,25,25,26,26,26),2)), ncol = 15, byrow = T)
    # Bring the plot together and print it to pdf
    ggsave(file="manusc/figures/Figure 2.pdf",
           plot=multiplot(titleplot,
                          plot14, plot15, plot41, plot18, plot101,
                          plot35, plot17, plot38, plot39, plot103,
                          titleplot2,
                          plot22, cholpathway,
                          plot20, subtitleplot,
                          plot45, plot49, plot50, plot51,
                          plot21, subtitleplot2,
                          plot52, plot47,
                          plot46,
                          plot48,
                          layout = ployout),
           width = 6.4,
           height = 9,
           dpi = 1200,
           unit = "in",
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
      annotate("text", x = 5, y = 1,
               label = 'bold("A                                                                                                              ")',
               size = 5, parse = TRUE, family = "Helvetica")
    titleplot2 <- ggplot(df) +
      theme_void() +
      geom_point() +
      xlim(0, 10) +
      ylim(0.5, 1.5) +
      annotate("text", x = 5, y = 1,
               label = 'bold("B                                                                                                              ")',
               size = 5, parse = TRUE, family = "Helvetica")
    titleplot3 <- ggplot(df) +
      theme_void() +
      geom_point() +
      xlim(0, 10) +
      ylim(0.5, 1.5) +
      annotate("text", x = 5, y = 1,
               label = 'bold("C                                                                                                              ")',
               size = 5, parse = TRUE, family = "Helvetica")
    # Create an empty plot if needed to adjust figures
    emptyplot <- ggplot(df) +
      theme_void() +
      geom_point() +
      xlim(0, 10) +
      ylim(0, 10)
    
    # For this plot, import png-files containing readily edited figures of pathway analysis
    #img <- readPNG("manusc/figures/pathwayenrichment.png")
    #img <- rasterGrob(img, interpolate=TRUE)
    #penric <- ggplot(df) +
    #  theme_void() +
    #  annotation_custom(img, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    #  xlim(0, 10) +
    #  ylim(0, 10)
    
    # Assign the layout matrix (how the plots are arranged)
    ployout = matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                       2,2,2,2,2,2,2,12,3,3,3,3,3,3,3,
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
                       11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,
                       11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,
                       11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,
                       11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,
                       11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,
                       11,11,11,11,11,11,11,11,11,11,11,12,12,12,12), ncol = 15, byrow = T)
    # Bring the plot together and print it to pdf
    ggsave(file="manusc/figures/Figure 3.pdf",
           plot=multiplot(titleplot,
                          pathway_figure_pos, pathway_figure_neg,
                          titleplot2,
                          plot78, plot93, plot94, plot61, plot69,
                          titleplot3,
                          aminoplot,
                          emptyplot,
                          layout = ployout), width = 7, height = 8, dpi = 1200, units = "in", device = "pdf")
  }
  
  # Figure 4
  {
    # Create titles
    df <- data.frame()
    emptyplot <- ggplot(df) +
      theme_void() +
      geom_point() +
      xlim(0, 10) +
      ylim(0, 10)
    # Assign the layout matrix (how the plots are arranged)
    ployout = matrix(c(1,1,1,1,2,2,2,3,3,3,4,4,5,5,5,5,
                       1,1,1,1,2,2,2,6,6,6,7,7,8,8,9,9,
                       12,12,12,10,10,10,10,10,12,11,11,11,11,12,12,12,
                       12,12,12,10,10,10,10,10,12,11,11,11,11,12,12,12), ncol = 16, byrow = T)
    # Bring the plot together and print it to pdf
    ggsave(file="manusc/figures/Figure 4.pdf",
           plot=multiplot(plot_ala, plot_dha, plot_carnitine183, plot_carnitine184, plot_carnitine204,
                          plot_lysopc160, plot_lysopc181, plot_lysope160, plot_lysope181,
                          plot_mcfa, plot_vlcfa,
                          emptyplot,
                          layout = ployout),
           width = 7, height = 5, dpi = 1200, units = "in", device = "pdf")
  }
}

# -------------------------------A. EXTRACT SOURCE DATA OF FIGURES-------------------------------- #
{
  # During the study, there were very few vegan children in Helsinki daycares. To protect them from
  # identification, the data is represented separately for each figure without IDs
  
  # Figure 1
  {
    src_fig1a <- data.frame(group5 = miradata$Group5,
                            age = miradata$Bage,
                            height = miradata$Height)
    src_fig1a <- src_fig1a[order(src_fig1a$height),]
    row.names(src_fig1a) <- NULL
    src_fig1b <- data.frame(group5 = miradata$Group5,
                            age = miradata$Bage,
                            bmi = miradata$BMI)
    src_fig1b <- src_fig1b[order(src_fig1b$bmi),]
    row.names(src_fig1b) <- NULL
    #fig1c_cols <- c(24,90,91,92,32,
    #                111,112,101,107,106)
    src_fig1c1 <- data.frame(group5 = miradata_ravinto$Group5,
                             energy = miradata_ravinto[,24])
    src_fig1c1 <- src_fig1c1[order(src_fig1c1$energy),]
    row.names(src_fig1c1) <- NULL
    src_fig1c2 <- data.frame(group5 = miradata_ravinto$Group5,
                             protein = miradata_ravinto[,90])
    src_fig1c2 <- src_fig1c2[order(src_fig1c2$protein),]
    row.names(src_fig1c2) <- NULL
    src_fig1c3 <- data.frame(group5 = miradata_ravinto$Group5,
                             fat = miradata_ravinto[,91])
    src_fig1c3 <- src_fig1c3[order(src_fig1c3$fat),]
    row.names(src_fig1c3) <- NULL
    src_fig1c4 <- data.frame(group5 = miradata_ravinto$Group5,
                             safa = miradata_ravinto[,92])
    src_fig1c4 <- src_fig1c4[order(src_fig1c4$safa),]
    row.names(src_fig1c4) <- NULL
    src_fig1c5 <- data.frame(group5 = miradata_ravinto$Group5,
                             cholesterol = miradata_ravinto[,32])
    src_fig1c5 <- src_fig1c5[order(src_fig1c5$cholesterol),]
    row.names(src_fig1c5) <- NULL
    src_fig1c6 <- data.frame(group5 = miradata_ravinto$Group5,
                             vit_a = miradata_ravinto[,111])
    src_fig1c6 <- src_fig1c6[order(src_fig1c6$vit_a),]
    row.names(src_fig1c6) <- NULL
    src_fig1c7 <- data.frame(group5 = miradata_ravinto$Group5,
                             vit_d = miradata_ravinto[,112])
    src_fig1c7 <- src_fig1c7[order(src_fig1c7$vit_d),]
    row.names(src_fig1c7) <- NULL
    src_fig1c8 <- data.frame(group5 = miradata_ravinto$Group5,
                             folate = miradata_ravinto[,101])
    src_fig1c8 <- src_fig1c8[order(src_fig1c8$folate),]
    row.names(src_fig1c8) <- NULL
    src_fig1c9 <- data.frame(group5 = miradata_ravinto$Group5,
                             zinc = miradata_ravinto[,107])
    src_fig1c9 <- src_fig1c9[order(src_fig1c9$zinc),]
    row.names(src_fig1c9) <- NULL
    src_fig1c10 <- data.frame(group5 = miradata_ravinto$Group5,
                             iron = miradata_ravinto[,106])
    src_fig1c10 <- src_fig1c10[order(src_fig1c10$iron),]
    row.names(src_fig1c10) <- NULL
    
    src_fig1 <- list(src_fig1a,
                     src_fig1b,
                     src_fig1c1, src_fig1c2, src_fig1c3, src_fig1c4, src_fig1c5,
                     src_fig1c6, src_fig1c7, src_fig1c8, src_fig1c9, src_fig1c10)
    sheetnames <- c("1A Height", "1B Weight",
                    "1C Energy", "1C Protein", "1C Fat", "1C SAFA", "1C Cholesterol",
                    "1C Vitamin A", "1C Vitamin D", "1C Folate", "1C Zinc", "1C Iron")
    for(i in 1:length(src_fig1)) {
      write.xlsx(src_fig1[[i]],
                 file = "manusc/Open data/EMM-2020-13492_SourceDataForFigure1A-C.xlsx",
                 sheetName = sheetnames[i],
                 append = TRUE)
    }
  }
  
  # Figure 2
  {
    # 2A columns: 14, 15, 41, 18, 101, 35, 17, 38, 39, 103
    src_fig2a1 <- data.frame(group5 = miradata$Group5,
                             vit_d_tot = miradata[,14])
    src_fig2a1 <- src_fig2a1[order(src_fig2a1$vit_d_tot),]
    row.names(src_fig2a1) <- NULL
    src_fig2a2 <- data.frame(group5 = miradata$Group5,
                             vit_d3 = miradata[,15])
    src_fig2a2 <- src_fig2a2[order(src_fig2a2$vit_d3),]
    row.names(src_fig2a2) <- NULL
    src_fig2a3 <- data.frame(group5 = miradata$Group5,
                             rbp = miradata[,41])
    src_fig2a3 <- src_fig2a3[order(src_fig2a3$rbp),]
    row.names(src_fig2a3) <- NULL
    src_fig2a4 <- data.frame(group5 = miradata$Group5,
                             transthyretin = miradata[,18])
    src_fig2a4 <- src_fig2a4[order(src_fig2a4$transthyretin),]
    row.names(src_fig2a4) <- NULL
    src_fig2a5 <- data.frame(group5 = miradata$Group5,
                             vit_a_status = miradata[,101])
    src_fig2a5 <- src_fig2a5[order(src_fig2a5$vit_a_status),]
    row.names(src_fig2a5) <- NULL
    src_fig2a6 <- data.frame(group5 = miradata$Group5,
                             folate = miradata[,35])
    src_fig2a6 <- src_fig2a6[order(src_fig2a6$folate),]
    row.names(src_fig2a6) <- NULL
    src_fig2a7 <- data.frame(group5 = miradata$Group5,
                             zinc = miradata[,17])
    src_fig2a7 <- src_fig2a7[order(src_fig2a7$zinc),]
    row.names(src_fig2a7) <- NULL
    src_fig2a8 <- data.frame(group5 = miradata$Group5,
                             ferritin = miradata[,38])
    src_fig2a8 <- src_fig2a8[order(src_fig2a8$ferritin),]
    row.names(src_fig2a8) <- NULL
    src_fig2a9 <- data.frame(group5 = miradata$Group5,
                             tfr = miradata[,39])
    src_fig2a9 <- src_fig2a9[order(src_fig2a9$tfr),]
    row.names(src_fig2a9) <- NULL
    src_fig2a10 <- data.frame(group5 = miradata$Group5,
                             iodine_to_crea = miradata[,103])
    src_fig2a10 <- src_fig2a10[order(src_fig2a10$iodine_to_crea),]
    row.names(src_fig2a10) <- NULL
    # Fig 2B cholesterols: 22, 20, 21
    src_fig2b1 <- data.frame(group5 = miradata$Group5,
                             chol_tot = miradata[,22])
    src_fig2b1 <- src_fig2b1[order(src_fig2b1$chol_tot),]
    row.names(src_fig2b1) <- NULL
    src_fig2b2 <- data.frame(group5 = miradata$Group5,
                             chol_ldl = miradata[,20])
    src_fig2b2 <- src_fig2b2[order(src_fig2b2$chol_ldl),]
    row.names(src_fig2b2) <- NULL
    src_fig2b3 <- data.frame(group5 = miradata$Group5,
                             chol_hdl = miradata[,21])
    src_fig2b3 <- src_fig2b3[order(src_fig2b3$chol_hdl),]
    row.names(src_fig2b3) <- NULL
    # 2B cholesterol metabolism: 45, 49, 50, 51, 52, 47, 46, 48
    src_fig2b4 <- data.frame(group5 = miradata$Group5,
                             cholestanol = miradata[,45])
    src_fig2b4 <- src_fig2b4[order(src_fig2b4$cholestanol),]
    row.names(src_fig2b4) <- NULL
    src_fig2b5 <- data.frame(group5 = miradata$Group5,
                             campesterol = miradata[,49])
    src_fig2b5 <- src_fig2b5[order(src_fig2b5$campesterol),]
    row.names(src_fig2b5) <- NULL
    src_fig2b6 <- data.frame(group5 = miradata$Group5,
                             sitosterol = miradata[,50])
    src_fig2b6 <- src_fig2b6[order(src_fig2b6$sitosterol),]
    row.names(src_fig2b6) <- NULL
    src_fig2b7 <- data.frame(group5 = miradata$Group5,
                             avenasterol = miradata[,51])
    src_fig2b7 <- src_fig2b7[order(src_fig2b7$avenasterol),]
    row.names(src_fig2b7) <- NULL
    src_fig2b8 <- data.frame(group5 = miradata$Group5,
                             squalene = miradata[,52])
    src_fig2b8 <- src_fig2b8[order(src_fig2b8$squalene),]
    row.names(src_fig2b8) <- NULL
    src_fig2b9 <- data.frame(group5 = miradata$Group5,
                             desmosterol = miradata[,47])
    src_fig2b9 <- src_fig2b9[order(src_fig2b9$desmosterol),]
    row.names(src_fig2b9) <- NULL
    src_fig2b10 <- data.frame(group5 = miradata$Group5,
                              cholestenol = miradata[,46])
    src_fig2b10 <- src_fig2b10[order(src_fig2b10$cholestenol),]
    row.names(src_fig2b10) <- NULL
    src_fig2b11 <- data.frame(group5 = miradata$Group5,
                              lathosterol = miradata[,48])
    src_fig2b11 <- src_fig2b11[order(src_fig2b11$lathosterol),]
    row.names(src_fig2b11) <- NULL

    src_fig2 <- list(src_fig2a1, src_fig2a2, src_fig2a3, src_fig2a4, src_fig2a5,
                     src_fig2a6, src_fig2a7, src_fig2a8, src_fig2a9, src_fig2a10,
                     src_fig2b1, src_fig2b2, src_fig2b3,
                     src_fig2b4, src_fig2b5, src_fig2b6, src_fig2b7,
                     src_fig2b8, src_fig2b9, src_fig2b10, src_fig2b11)
    sheetnames <- c("2A Total Vitamin D", "2A Vitamin D3", "2A RBP", "2A Transthyretin", "2A Vitamin A status",
                    "2A Folate", "2A Zinc", "2A Ferritin", "2A TfR", "2A Iodine to creatinine",
                    "2B Total cholesterol", "2B LDL cholesterol", "2B HDL cholesterol",
                    "2B Cholestanol", "2B Campesterol", "2B Sitosterol", "2B Avenasterol",
                    "2B Squalene", "2B Desmosterol", "2B Cholestenol", "2B Lathosterol")
    for(i in 1:length(src_fig2)) {
      write.xlsx(src_fig2[[i]],
                 file = "manusc/Open data/EMM-2020-13492_SourceDataForFigure2A-B.xlsx",
                 sheetName = sheetnames[i],
                 append = TRUE)
    }
  }
  
  # Figure 3
  {
    # 2A columns: 14, 15, 41, 18, 101, 35, 17, 38, 39, 103
    src_fig3a <- data.frame(pathway = pa_results.forfigs$pathway,
                            p.value_pos = pa_results.forfigs$p.value_pos,
                            p.value_neg = pa_results.forfigs$p.value_neg,
                            enrichmentscore = pa_results.forfigs$es)
    # Fig 3b: 78, 93, 94, 61, 69
    src_fig3b1 <- data.frame(group5 = miradata$Group5,
                             ba_conjug_tot = miradata[,78])
    src_fig3b1 <- src_fig3b1[order(src_fig3b1$ba_conjug_tot),]
    row.names(src_fig3b1) <- NULL
    src_fig3b2 <- data.frame(group5 = miradata$Group5,
                             ba_conjug_prim = miradata[,93])
    src_fig3b2 <- src_fig3b2[order(src_fig3b2$ba_conjug_prim),]
    row.names(src_fig3b2) <- NULL
    src_fig3b3 <- data.frame(group5 = miradata$Group5,
                             ba_conjug_sec = miradata[,94])
    src_fig3b3 <- src_fig3b3[order(src_fig3b3$ba_conjug_sec),]
    row.names(src_fig3b3) <- NULL
    src_fig3b4 <- data.frame(group5 = miradata$Group5,
                             cdca = miradata[,61])
    src_fig3b4 <- src_fig3b4[order(src_fig3b4$cdca),]
    row.names(src_fig3b4) <- NULL
    src_fig3b5 <- data.frame(group5 = miradata$Group5,
                             ca = miradata[,69])
    src_fig3b5 <- src_fig3b5[order(src_fig3b5$ca),]
    row.names(src_fig3b5) <- NULL
    src_fig3c <- AAdata
    colnames(src_fig3c) <- c("diff_means", "conf_int_length", "aminoacid", "type")
    
    src_fig3 <- list(src_fig3a,
                     src_fig3b1, src_fig3b2, src_fig3b3, src_fig3b4, src_fig3b5,
                     src_fig3c)
    sheetnames <- c("3A Pathway Analysis",
                    "3B Bile acid conjugation (all BAs)",
                    "3B Bile acid conjugation (primary BAs)",
                    "3B Bile acid conjugation (sec. BAs)",
                    "3B Chenodeoxycholic acid",
                    "3B Cholic acid",
                    "3C Amino acids")
    for(i in 1:length(src_fig3)) {
      write.xlsx(src_fig3[[i]],
                 file = "manusc/Open data/EMM-2020-13492_SourceDataForFigure3A-C.xlsx",
                 sheetName = sheetnames[i],
                 append = TRUE)
    }
  }
  
  # Figure 4
  {
    # 4A: 416, 500
    src_fig4a1 <- data.frame(group5 = miraLIPO_norm$Group5,
                             ala = miraLIPO_norm[,416])
    src_fig4a1 <- src_fig4a1[order(src_fig4a1$ala),]
    row.names(src_fig4a1) <- NULL
    src_fig4a2 <- data.frame(group5 = miraLIPO_norm$Group5,
                             dha = miraLIPO_norm[,500])
    src_fig4a2 <- src_fig4a2[order(src_fig4a2$dha),]
    row.names(src_fig4a2) <- NULL
    # 4B: 615, 613, 657
    src_fig4b1 <- data.frame(group5 = miraLIPO_norm$Group5,
                             carnitine183 = miraLIPO_norm[,615])
    src_fig4b1 <- src_fig4b1[order(src_fig4b1$carnitine183),]
    row.names(src_fig4b1) <- NULL
    src_fig4b2 <- data.frame(group5 = miraLIPO_norm$Group5,
                             carnitine184 = miraLIPO_norm[,613])
    src_fig4b2 <- src_fig4b2[order(src_fig4b2$carnitine184),]
    row.names(src_fig4b2) <- NULL
    src_fig4b3 <- data.frame(group5 = miraLIPO_norm$Group5,
                             carnitine204 = miraLIPO_norm[,657])
    src_fig4b3 <- src_fig4b3[order(src_fig4b3$carnitine204),]
    row.names(src_fig4b3) <- NULL
    # 4C: 698, 724, 654, 686
    src_fig4c1 <- data.frame(group5 = miraLIPO_norm$Group5,
                             lysopc160 = miraLIPO_norm[,698])
    src_fig4c1 <- src_fig4c1[order(src_fig4c1$lysopc160),]
    row.names(src_fig4c1) <- NULL
    src_fig4c2 <- data.frame(group5 = miraLIPO_norm$Group5,
                             lysopc181 = miraLIPO_norm[,724])
    src_fig4c2 <- src_fig4c2[order(src_fig4c2$lysopc181),]
    row.names(src_fig4c2) <- NULL
    src_fig4c3 <- data.frame(group5 = miraLIPO_norm$Group5,
                             lysope160 = miraLIPO_norm[,654])
    src_fig4c3 <- src_fig4c3[order(src_fig4c3$lysope160),]
    row.names(src_fig4c3) <- NULL
    src_fig4c4 <- data.frame(group5 = miraLIPO_norm$Group5,
                             lysope181 = miraLIPO_norm[,686])
    src_fig4c4 <- src_fig4c4[order(src_fig4c4$lysope181),]
    row.names(src_fig4c4) <- NULL
    # 4D: MCFA
    src_fig4d1 <- miraLIPO_norm_MCFA
    src_fig4d1 <- src_fig4d1[order(src_fig4d1$variable, src_fig4d1$value),]
    row.names(src_fig4d1) <- NULL
    src_fig4d1$variable <- as.character(src_fig4d1$variable)
    src_fig4d1$variable[src_fig4d1$variable == "V699"] <- "tg_c260"
    src_fig4d1$variable[src_fig4d1$variable == "V746"] <- "tg_c310"
    src_fig4d1$variable[src_fig4d1$variable == "V759"] <- "tg_c330"
    # 4D: (V)LCFA
    src_fig4d2 <- miraLIPO_norm_VLCFA
    src_fig4d2 <- src_fig4d2[order(src_fig4d2$variable, src_fig4d2$value),]
    row.names(src_fig4d2) <- NULL
    src_fig4d2$variable <- as.character(src_fig4d2$variable)
    src_fig4d2$variable[src_fig4d2$variable == "V834"] <- "tg_c500-o"
    src_fig4d2$variable[src_fig4d2$variable == "V845"] <- "tg_c534"
    src_fig4d2$variable[src_fig4d2$variable == "V858"] <- "tg_c580-o"
    src_fig4d2$variable[src_fig4d2$variable == "V861"] <- "tg_c584"
    src_fig4d2$variable[src_fig4d2$variable == "V855"] <- "tg_c589-o"
    src_fig4d2$variable[src_fig4d2$variable == "V869"] <- "tg_c6210"
    
    src_fig4 <- list(src_fig4a1, src_fig4a2,
                     src_fig4b1, src_fig4b2, src_fig4b3,
                     src_fig4c1, src_fig4c2, src_fig4c3, src_fig4c4,
                     src_fig4d1, src_fig4d2)
    sheetnames <- c("4A Alpha-linolenic acid", "4A Docosahexaenic acid",
                    "4B Carnitine 18-3", "4B Carnitine 18-4", "4B Carnitine 20-4",
                    "4C LysoPC 16-0", "4C LysoPC 18-1", "4C LysoPE 16-0", "4C LysoPE 18-1",
                    "4D MCFA Triglycerides", "4D (V)LCFA Triglycerides")
    for(i in 1:length(src_fig4)) {
      write.xlsx(src_fig4[[i]],
                 file = "manusc/Open data/EMM-2020-13492_SourceDataForFigure4A-D.xlsx",
                 sheetName = sheetnames[i],
                 append = TRUE)
    }
  }
  
  # Figure EV1
  {
    src_figev1 <- as.data.frame(t(as.matrix(dist(scale(miradata_clust[,c(3:(length(clustdata) + 2))]), method = "euclidean"))))
    src_figev1$group5 <- miradata_clust$Group5
    row.names(src_figev1) <- NULL
    colnames(src_figev1) <- c(1:40, "group5")
    write.xlsx(src_figev1,
               file = "manusc/Open data/EMM-2020-13492_SourceDataForFigureEV1.xlsx",
               sheetName = "Eucl. distance matrix")
  }
}

# -----------------------------B. INDIVIDUAL CODE FOR EACH SUBFIGURE------------------------------ #
{
  # There was a struggle with looping over id's as if you save e.g. plot14 with parts referring to i:th
  # column, you will need to have i = 14 each time you would like to use the saved plot anywhere.
  
  # Nutrient intake figures
  {
    # For these figures, we need to use intake data frame instead of "miradata"
    data_rav = miradata_ravinto
    data = data_rav
    # plotnene (Energy intake, column 24)
    {
      plotnene <- ggplot(data = data_rav, aes(x=Group5, y=data[,24])) +
        theme_light() +
        scale_colour_manual(name = "Dietary group",
                            values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF", "#6CB43FFF"),
                            guide = guide_legend(title.position = "right")) +
        scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                           labels = c("Control", "Vegetarian", "Vegan")) +
        scale_size_manual(name = "Dietary group", values=c(1,1,1),
                          labels = c("Control", "Vegetarian", "Vegan")) +
        scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
        geom_point(position = position_jitter(width = .1, height = 0),
                   aes(colour = Group5, shape = Group5, size = Group5)) +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
        stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
        labs(y = "kJ/d", x = "", title = "Energy") +
        #geom_signif(comparisons = list(c("Control", "Vegan")),
        #            test = 'permira.test.gg.rav', map_signif_level = T, margin_top = 0.1, textsize = 3) +
        coord_cartesian(ylim = c(min(data[,24], na.rm=TRUE)*0.95,
                                 max(data[,24], na.rm=TRUE)*1.2)) +
        theme(text = element_text(family = "Helvetica"),
              plot.title = element_text(hjust = 0.5, size = 8),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
              axis.text.y = element_text(size = 8, colour = "black"),
              axis.title.y = element_text(size = 8),
              axis.ticks = element_line(colour = "black"),
              axis.line = element_line(colour = "black"),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              legend.position = "none") +
        guides(colour=F, linetype=F, size=F, shape=F)
    }
    # plotnpro (Protein intake, column 90)
    {
      plotnpro <- ggplot(data = data_rav, aes(x=Group5, y=data[,90])) +
        theme_light() +
        scale_colour_manual(name = "Dietary group",
                            values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF", "#6CB43FFF"),
                            guide = guide_legend(title.position = "right")) +
        scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                           labels = c("Control", "Vegetarian", "Vegan")) +
        scale_size_manual(name = "Dietary group", values=c(1,1,1),
                          labels = c("Control", "Vegetarian", "Vegan")) +
        scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
        geom_point(position = position_jitter(width = .1, height = 0),
                   aes(colour = Group5, shape = Group5, size = Group5)) +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
        stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
        labs(y = "% of daily energy", x = "", title = "Protein") +
        geom_signif(comparisons = list(c("Control", "Vegan")),
                    test = 'permira.test.gg.rav', map_signif_level = T, margin_top = 0.1, textsize = 3) +
        coord_cartesian(ylim = c(min(data[,90], na.rm=TRUE)*0.95,
                                 max(data[,90], na.rm=TRUE)*1.2)) +
        theme(text = element_text(family = "Helvetica"),
              plot.title = element_text(hjust = 0.5, size = 8),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
              axis.text.y = element_text(size = 8, colour = "black"),
              axis.title.y = element_text(size = 8),
              axis.ticks = element_line(colour = "black"),
              axis.line = element_line(colour = "black"),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              legend.position = "none") +
        guides(colour=F, linetype=F, size=F, shape=F)
    }
    # plotnfat (Fat intake, column 91)
    {
      plotnfat <- ggplot(data = data_rav, aes(x=Group5, y=data[,91])) +
        theme_light() +
        scale_colour_manual(name = "Dietary group",
                            values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF", "#6CB43FFF"),
                            guide = guide_legend(title.position = "right")) +
        scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                           labels = c("Control", "Vegetarian", "Vegan")) +
        scale_size_manual(name = "Dietary group", values=c(1,1,1),
                          labels = c("Control", "Vegetarian", "Vegan")) +
        scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
        geom_point(position = position_jitter(width = .1, height = 0),
                   aes(colour = Group5, shape = Group5, size = Group5)) +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
        stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
        labs(y = "% of daily energy", x = "", title = "Fat") +
        #geom_signif(comparisons = list(c("Control", "Vegan")),
        #            test = 'permira.test.gg.rav', map_signif_level = T, margin_top = 0.1, textsize = 3) +
        coord_cartesian(ylim = c(min(data[,91], na.rm=TRUE)*0.95,
                                 max(data[,91], na.rm=TRUE)*1.2)) +
        theme(text = element_text(family = "Helvetica"),
              plot.title = element_text(hjust = 0.5, size = 8),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
              axis.text.y = element_text(size = 8, colour = "black"),
              axis.title.y = element_text(size = 8),
              axis.ticks = element_line(colour = "black"),
              axis.line = element_line(colour = "black"),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              legend.position = "none") +
        guides(colour=F, linetype=F, size=F, shape=F)
    }
    # plotnsafa (Saturated fatty acid intake, column 92)
    {
      plotnsafa <- ggplot(data = data_rav, aes(x=Group5, y=data[,92])) +
        theme_light() +
        scale_colour_manual(name = "Dietary group",
                            values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF", "#6CB43FFF"),
                            guide = guide_legend(title.position = "right")) +
        scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                           labels = c("Control", "Vegetarian", "Vegan")) +
        scale_size_manual(name = "Dietary group", values=c(1,1,1),
                          labels = c("Control", "Vegetarian", "Vegan")) +
        scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
        geom_point(position = position_jitter(width = .1, height = 0),
                   aes(colour = Group5, shape = Group5, size = Group5)) +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
        stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
        labs(y = "% of daily energy", x = "", title = "Saturated\nfatty acids") +
        geom_signif(comparisons = list(c("Control", "Vegan")),
                    test = 'permira.test.gg.rav', map_signif_level = T, margin_top = 0.1, textsize = 3) +
        geom_segment(x = 1, xend = 2, y = 20, yend = 20) +
        geom_segment(x = 1, xend = 1, y = 19.5, yend = 20) +
        geom_segment(x = 2, xend = 2, y = 19.5, yend = 20) +
        annotate("text", x = 1.5, y = 20.5, label = "*", size = 3) +
        coord_cartesian(ylim = c(min(data[,92], na.rm=TRUE)*0.95,
                                 max(data[,92], na.rm=TRUE)*1.3)) +
        theme(text = element_text(family = "Helvetica"),
              plot.title = element_text(hjust = 0.5, size = 8),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
              axis.text.y = element_text(size = 8, colour = "black"),
              axis.title.y = element_text(size = 8),
              axis.ticks = element_line(colour = "black"),
              axis.line = element_line(colour = "black"),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              legend.position = "none") +
        guides(colour=F, linetype=F, size=F, shape=F)
    }
    # plotnchol (Saturated fatty acid intake, column 32)
    {
      plotnchol <- ggplot(data = data_rav, aes(x=Group5, y=data[,32])) +
        theme_light() +
        scale_colour_manual(name = "Dietary group",
                            values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF", "#6CB43FFF"),
                            guide = guide_legend(title.position = "right")) +
        scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                           labels = c("Control", "Vegetarian", "Vegan")) +
        scale_size_manual(name = "Dietary group", values=c(1,1,1),
                          labels = c("Control", "Vegetarian", "Vegan")) +
        scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
        geom_point(position = position_jitter(width = .1, height = 0),
                   aes(colour = Group5, shape = Group5, size = Group5)) +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
        stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
        labs(y = "g/d", x = "", title = "Cholesterol") +
        geom_signif(comparisons = list(c("Control", "Vegan")),
                    test = 'permira.test.gg.rav', map_signif_level = T, margin_top = 0.1, textsize = 3) +
        geom_segment(x = 1, xend = 2, y = 290, yend = 290) +
        geom_segment(x = 1, xend = 1, y = 280, yend = 290) +
        geom_segment(x = 2, xend = 2, y = 280, yend = 290) +
        annotate("text", x = 1.5, y = 295, label = "*", size = 3) +
        geom_segment(x = 2, xend = 3, y = 290, yend = 290) +
        geom_segment(x = 2, xend = 2, y = 280, yend = 290) +
        geom_segment(x = 3, xend = 3, y = 280, yend = 290) +
        annotate("text", x = 2.5, y = 295, label = "*", size = 3) +
        coord_cartesian(ylim = c(min(data[,32], na.rm=TRUE)*0.95,
                                 max(data[,32], na.rm=TRUE)*1.25)) +
        theme(text = element_text(family = "Helvetica"),
              plot.title = element_text(hjust = 0.5, size = 8),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
              axis.text.y = element_text(size = 8, colour = "black"),
              axis.title.y = element_text(size = 8),
              axis.ticks = element_line(colour = "black"),
              axis.line = element_line(colour = "black"),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              legend.position = "none") +
        guides(colour=F, linetype=F, size=F, shape=F)
    }
    # plotnvita (Vitamin A intake, column 111)
    {
      plotnvita <- ggplot(data = data_rav, aes(x=Group5, y=data[,111])) +
        theme_light() +
        scale_colour_manual(name = "Dietary group",
                            values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF", "#6CB43FFF"),
                            guide = guide_legend(title.position = "right")) +
        scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                           labels = c("Control", "Vegetarian", "Vegan")) +
        scale_size_manual(name = "Dietary group", values=c(1,1,1),
                          labels = c("Control", "Vegetarian", "Vegan")) +
        scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
        geom_point(position = position_jitter(width = .1, height = 0),
                   aes(colour = Group5, shape = Group5, size = Group5)) +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
        stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
        labs(y = expression(paste(mu, "g (RAE)/d")), x = "", title = "Vitamin A") +
        #geom_signif(comparisons = list(c("Control", "Vegan")),
        #            test = 'permira.test.gg.rav', map_signif_level = T, margin_top = 0.1, textsize = 3) +
        coord_cartesian(ylim = c(min(data[,111], na.rm=TRUE)*0.95,
                                 max(data[,111], na.rm=TRUE)*1.2)) +
        theme(text = element_text(family = "Helvetica"),
              plot.title = element_text(hjust = 0.5, size = 8),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
              axis.text.y = element_text(size = 8, colour = "black"),
              axis.title.y = element_text(size = 8),
              axis.ticks = element_line(colour = "black"),
              axis.line = element_line(colour = "black"),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              legend.position = "none") +
        guides(colour=F, linetype=F, size=F, shape=F)
    }
    # plotnvitd (Vitamin D intake, column 112)
    {
      plotnvitd <- ggplot(data = data_rav, aes(x=Group5, y=data[,112])) +
        theme_light() +
        scale_colour_manual(name = "Dietary group",
                            values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF", "#6CB43FFF"),
                            guide = guide_legend(title.position = "right")) +
        scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                           labels = c("Control", "Vegetarian", "Vegan")) +
        scale_size_manual(name = "Dietary group", values=c(1,1,1),
                          labels = c("Control", "Vegetarian", "Vegan")) +
        scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
        geom_point(position = position_jitter(width = .1, height = 0),
                   aes(colour = Group5, shape = Group5, size = Group5)) +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
        stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
        labs(y = expression(paste(mu, "g/d")), x = "", title = "Vitamin D") +
        #geom_signif(comparisons = list(c("Control", "Vegan")),
        #            test = 'permira.test.gg.rav', map_signif_level = T, margin_top = 0.1, textsize = 3) +
        coord_cartesian(ylim = c(min(data[,112], na.rm=TRUE)*0.95,
                                 max(data[,112], na.rm=TRUE)*1.2)) +
        theme(text = element_text(family = "Helvetica"),
              plot.title = element_text(hjust = 0.5, size = 8),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
              axis.text.y = element_text(size = 8, colour = "black"),
              axis.title.y = element_text(size = 8),
              axis.ticks = element_line(colour = "black"),
              axis.line = element_line(colour = "black"),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              legend.position = "none") +
        guides(colour=F, linetype=F, size=F, shape=F)
    }
    # plotnfol (Folate intake, column 107)
    {
      plotnfol <- ggplot(data = data_rav, aes(x=Group5, y=data[,101])) +
        theme_light() +
        scale_colour_manual(name = "Dietary group",
                            values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF", "#6CB43FFF"),
                            guide = guide_legend(title.position = "right")) +
        scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                           labels = c("Control", "Vegetarian", "Vegan")) +
        scale_size_manual(name = "Dietary group", values=c(1,1,1),
                          labels = c("Control", "Vegetarian", "Vegan")) +
        scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
        geom_point(position = position_jitter(width = .1, height = 0),
                   aes(colour = Group5, shape = Group5, size = Group5)) +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
        stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
        labs(y = expression(paste(mu, "g/d")), x = "", title = "Folate") +
        geom_signif(comparisons = list(c("Control", "Vegan")),
                    test = 'permira.test.gg.rav', map_signif_level = T, margin_top = 0.1, textsize = 3) +
        geom_segment(x = 1, xend = 2, y = 600, yend = 600) +
        geom_segment(x = 1, xend = 1, y = 580, yend = 600) +
        geom_segment(x = 2, xend = 2, y = 580, yend = 600) +
        annotate("text", x = 1.5, y = 610, label = "**", size = 3) +
        geom_segment(x = 2, xend = 3, y = 600, yend = 600) +
        geom_segment(x = 2, xend = 2, y = 580, yend = 600) +
        geom_segment(x = 3, xend = 3, y = 580, yend = 600) +
        annotate("text", x = 2.5, y = 610, label = "*", size = 3) +
        coord_cartesian(ylim = c(min(data[,101], na.rm=TRUE)*0.95,
                                 max(data[,101], na.rm=TRUE)*1.25)) +
        theme(text = element_text(family = "Helvetica"),
              plot.title = element_text(hjust = 0.5, size = 8),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
              axis.text.y = element_text(size = 8, colour = "black"),
              axis.title.y = element_text(size = 8),
              axis.ticks = element_line(colour = "black"),
              axis.line = element_line(colour = "black"),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              legend.position = "none") +
        guides(colour=F, linetype=F, size=F, shape=F)
    }
    # plotnzn (Zinc intake, column 107)
    {
      plotnzn <- ggplot(data = data_rav, aes(x=Group5, y=data[,107])) +
        theme_light() +
        scale_colour_manual(name = "Dietary group",
                            values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF", "#6CB43FFF"),
                            guide = guide_legend(title.position = "right")) +
        scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                           labels = c("Control", "Vegetarian", "Vegan")) +
        scale_size_manual(name = "Dietary group", values=c(1,1,1),
                          labels = c("Control", "Vegetarian", "Vegan")) +
        scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
        geom_point(position = position_jitter(width = .1, height = 0),
                   aes(colour = Group5, shape = Group5, size = Group5)) +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
        stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
        labs(y = "mg/d", x = "", title = "Zinc") +
        geom_signif(comparisons = list(c("Control", "Vegan")),
                    test = 'permira.test.gg.rav', map_signif_level = T, margin_top = 0.1, textsize = 3) +
        coord_cartesian(ylim = c(min(data[,107], na.rm=TRUE)*0.95,
                                 max(data[,107], na.rm=TRUE)*1.2)) +
        theme(text = element_text(family = "Helvetica"),
              plot.title = element_text(hjust = 0.5, size = 8),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
              axis.text.y = element_text(size = 8, colour = "black"),
              axis.title.y = element_text(size = 8),
              axis.ticks = element_line(colour = "black"),
              axis.line = element_line(colour = "black"),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              legend.position = "none") +
        guides(colour=F, linetype=F, size=F, shape=F)
    }
    # plotnfe (Iron intake, column 106)
    {
      plotnfe <- ggplot(data = data_rav, aes(x=Group5, y=data[,106])) +
        theme_light() +
        scale_colour_manual(name = "Dietary group",
                            values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF", "#6CB43FFF"),
                            guide = guide_legend(title.position = "right")) +
        scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                           labels = c("Control", "Vegetarian", "Vegan")) +
        scale_size_manual(name = "Dietary group", values=c(1,1,1),
                          labels = c("Control", "Vegetarian", "Vegan")) +
        scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
        geom_point(position = position_jitter(width = .1, height = 0),
                   aes(colour = Group5, shape = Group5, size = Group5)) +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
        stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
        labs(y = "mg/d", x = "", title = "Iron") +
        geom_signif(comparisons = list(c("Control", "Vegan")),
                    test = 'permira.test.gg.rav', map_signif_level = T, margin_top = 0.1, textsize = 3) +
        geom_segment(x = 1, xend = 2, y = 20, yend = 20) +
        geom_segment(x = 1, xend = 1, y = 19.5, yend = 20) +
        geom_segment(x = 2, xend = 2, y = 19.5, yend = 20) +
        annotate("text", x = 1.5, y = 20.5, label = "**", size = 3) +
        coord_cartesian(ylim = c(min(data[,106], na.rm=TRUE)*0.95,
                                 max(data[,106], na.rm=TRUE)*1.25)) +
        theme(text = element_text(family = "Helvetica"),
              plot.title = element_text(hjust = 0.5, size = 8),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
              axis.text.y = element_text(size = 8, colour = "black"),
              axis.title.y = element_text(size = 8),
              axis.ticks = element_line(colour = "black"),
              axis.line = element_line(colour = "black"),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              legend.position = "none") +
        guides(colour=F, linetype=F, size=F, shape=F)
    }
  }
  
  # plot14 (Total vitamin D)
  {
    plot14 <- ggplot(data = data, aes(x=Group5, y=data[,14])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[14], x = "", title = miranames[14]) +
      geom_signif(comparisons = list(c("Control", "Vegan")),
                  test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      geom_hline(aes(yintercept=50, linetype="dashed", alpha = .5)) +
      coord_cartesian(ylim = c(min(data[,14], na.rm=TRUE)*0.95,
                               max(data[,14], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot15 (Vitamin D3)
  {
    plot15 <- ggplot(data = data, aes(x=Group5, y=data[,15])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[15], x = "", title = miranames[15]) +
      geom_signif(comparisons = list(c("Control", "Vegan")),
                  test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      geom_hline(aes(yintercept=50, linetype="dashed", alpha = .5)) +
      coord_cartesian(ylim = c(min(data[,15], na.rm=TRUE)*0.95,
                               max(data[,15], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot41 (RBP)
  {
    plot41 <- ggplot(data = data, aes(x=Group5, y=data[,41])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[41], x = "", title = miranames[41]) +
      geom_signif(comparisons = list(c("Control", "Vegan")),
                  test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      geom_hline(aes(yintercept=1.17, linetype="dotted", alpha = .5)) +
      geom_hline(aes(yintercept=0.83, linetype="dashed", alpha = .5)) +
      coord_cartesian(ylim = c(min(data[,41], na.rm=TRUE)*0.95,
                               max(data[,41], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot18 (Transthyretin)
  {
    plot18 <- ggplot(data = data, aes(x=Group5, y=data[,18])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[18], x = "", title = miranames[18]) +
      geom_signif(comparisons = list(c("Control", "Vegan")),
                  test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      coord_cartesian(ylim = c(min(data[,18], na.rm=TRUE)*0.95,
                               max(data[,18], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot101 (Vitamin A Status)
  {
    plot101 <- ggplot(data = data, aes(x=Group5, y=data[,101])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = "Arbitrary units", x = "", title = miranames[101]) +
      geom_signif(comparisons = list(c("Control", "Vegan")),
                  test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      geom_hline(aes(yintercept=0, linetype="dashed", alpha = .5)) +
      coord_cartesian(ylim = c(min(data[,101], na.rm=TRUE)*0.95,
                               max(data[,101], na.rm=TRUE)*1.3)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot35 (Folate)
  {
    plot35 <- ggplot(data = data, aes(x=Group5, y=data[,35])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[35], x = "", title = "Folate (Erythrocytes)         ") +
      geom_signif(comparisons = list(c("Control", "Vegan")),
                  test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      coord_cartesian(ylim = c(min(data[,35], na.rm=TRUE)*0.95,
                               max(data[,35], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot17 (Zinc)
  {
    plot17 <- ggplot(data = data, aes(x=Group5, y=data[,17])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[17], x = "", title = miranames[17]) +
      #geom_signif(comparisons = list(c("Control", "Vegan")),
      #            test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      geom_segment(x = 1, xend = 2, y = 26, yend = 26) +
      geom_segment(x = 1, xend = 1, y = 25.5, yend = 26) +
      geom_segment(x = 2, xend = 2, y = 25.5, yend = 26) +
      annotate("text", x = 1.5, y = 26.5, label = "*", size = 3) +
      coord_cartesian(ylim = c(min(data[,17], na.rm=TRUE)*0.95,
                               max(data[,17], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot38 (Ferritin)
  {
    plot38 <- ggplot(data = data, aes(x=Group5, y=data[,38])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[38], x = "", title = miranames[38]) +
      #geom_signif(comparisons = list(c("Control", "Vegan")),
      #            test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      coord_cartesian(ylim = c(min(data[,38], na.rm=TRUE)*0.95,
                               max(data[,38], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot39 (TfR)
  {
    plot39 <- ggplot(data = data, aes(x=Group5, y=data[,39])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[39], x = "", title = miranames[39]) +
      #geom_signif(comparisons = list(c("Control", "Vegan")),
      #            test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      coord_cartesian(ylim = c(min(data[,39], na.rm=TRUE)*0.95,
                               max(data[,39], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot103 (Iodine to creatinine)
  {
    plot103 <- ggplot(data = data, aes(x=Group5, y=data[,103])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[103], x = "", title = "Iodine (Urine)") +
      #geom_signif(comparisons = list(c("Control", "Vegan")),
      #            test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) + 
      coord_cartesian(ylim = c(min(data[,103], na.rm=TRUE)*0.95,
                               max(data[,103], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot22 (Total cholesterol)
  {
    plot22 <- ggplot(data = data, aes(x=Group5, y=data[,22])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[22], x = "", title = miranames[22]) +
      geom_signif(comparisons = list(c("Control", "Vegan")),
                  test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      coord_cartesian(ylim = c(min(data[,22], na.rm=TRUE)*0.95,
                               max(data[,22], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none") +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot20 (LDL cholesterol)
  {
    plot20 <- ggplot(data = data, aes(x=Group5, y=data[,20])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[20], x = "", title = miranames[20]) +
      geom_signif(comparisons = list(c("Control", "Vegan")),
                  test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      coord_cartesian(ylim = c(min(data[,20], na.rm=TRUE)*0.95,
                               max(data[,20], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none") +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot21 (HDL cholesterol)
  {
    plot21 <- ggplot(data = data, aes(x=Group5, y=data[,21])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[21], x = "", title = miranames[21]) +
      geom_signif(comparisons = list(c("Control", "Vegan")),
                  test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      coord_cartesian(ylim = c(min(data[,21], na.rm=TRUE)*0.95,
                               max(data[,21], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none") +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot45 (Cholestanol) 45, 49, 50, 51, 52, 47, 46, 48,
  {
    plot45 <- ggplot(data = data, aes(x=Group5, y=data[,45])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[45], x = "", title = miranames[45]) +
      geom_signif(comparisons = list(c("Control", "Vegan")),
                  test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      geom_segment(x = 2, xend = 3, y = 2.25, yend = 2.25) +
      geom_segment(x = 2, xend = 2, y = 2.21, yend = 2.25) +
      geom_segment(x = 3, xend = 3, y = 2.21, yend = 2.25) +
      annotate("text", x = 2.5, y = 2.29, label = "*", size = 3) +
      coord_cartesian(ylim = c(min(data[,45], na.rm=TRUE)*0.95,
                               max(data[,45], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot49 (Campesterol)
  {
    plot49 <- ggplot(data = data, aes(x=Group5, y=data[,49])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[49], x = "", title = miranames[49]) +
      geom_signif(comparisons = list(c("Control", "Vegan")),
                  test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      coord_cartesian(ylim = c(min(data[,49], na.rm=TRUE)*0.95,
                               max(data[,49], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot50 (Sitosterol)
  {
    plot50 <- ggplot(data = data, aes(x=Group5, y=data[,50])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[50], x = "", title = miranames[50]) +
      geom_signif(comparisons = list(c("Control", "Vegan")),
                  test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      geom_segment(x = 1, xend = 2, y = 3.53, yend = 3.53) +
      geom_segment(x = 1, xend = 1, y = 3.45, yend = 3.53) +
      geom_segment(x = 2, xend = 2, y = 3.45, yend = 3.53) +
      annotate("text", x = 1.5, y = 3.6, label = "*", size = 3) +
      coord_cartesian(ylim = c(min(data[,50], na.rm=TRUE)*0.95,
                               max(data[,50], na.rm=TRUE)*1.25)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot51 (Avenasterol)
  {
    plot51 <- ggplot(data = data, aes(x=Group5, y=data[,51])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[51], x = "", title = miranames[51]) +
      geom_signif(comparisons = list(c("Control", "Vegan")),
                  test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      geom_segment(x = 1, xend = 2, y = 1.09, yend = 1.09) +
      geom_segment(x = 1, xend = 1, y = 1.06, yend = 1.09) +
      geom_segment(x = 2, xend = 2, y = 1.06, yend = 1.09) +
      annotate("text", x = 1.5, y = 1.12, label = "*", size = 3) +
      coord_cartesian(ylim = c(min(data[,51], na.rm=TRUE)*0.95,
                               max(data[,51], na.rm=TRUE)*1.3)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot52 (Squalene)
  {
    plot52 <- ggplot(data = data, aes(x=Group5, y=data[,52])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[52], x = "", title = miranames[52]) +
      #geom_signif(comparisons = list(c("Control", "Vegan")),
      #            test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      coord_cartesian(ylim = c(min(data[,52], na.rm=TRUE)*0.95,
                               max(data[,52], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot47 (Desmosterol)
  {
    plot47 <- ggplot(data = data, aes(x=Group5, y=data[,47])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[47], x = "", title = miranames[47]) +
      #geom_signif(comparisons = list(c("Control", "Vegan")),
      #            test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      geom_segment(x = 2, xend = 3, y = 1.4, yend = 1.4) +
      geom_segment(x = 2, xend = 2, y = 1.38, yend = 1.4) +
      geom_segment(x = 3, xend = 3, y = 1.38, yend = 1.4) +
      annotate("text", x = 2.5, y = 1.42, label = "*", size = 3) +
      coord_cartesian(ylim = c(min(data[,47], na.rm=TRUE)*0.95,
                               max(data[,47], na.rm=TRUE)*1.15)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot46 (Cholestenol)
  {
    plot46 <- ggplot(data = data, aes(x=Group5, y=data[,46])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[46], x = "", title = miranames[46]) +
      #geom_signif(comparisons = list(c("Control", "Vegan")),
      #            test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      geom_segment(x = 1, xend = 2, y = 0.5, yend = 0.5) +
      geom_segment(x = 1, xend = 1, y = 0.49, yend = 0.5) +
      geom_segment(x = 2, xend = 2, y = 0.49, yend = 0.5) +
      annotate("text", x = 1.5, y = 0.51, label = "*", size = 3) +
      coord_cartesian(ylim = c(min(data[,46], na.rm=TRUE)*0.95,
                               max(data[,46], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot48 (Lathosterol)
  {
    plot48 <- ggplot(data = data, aes(x=Group5, y=data[,48])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[48], x = "", title = miranames[48]) +
      #geom_signif(comparisons = list(c("Control", "Vegan")),
      #            test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      geom_segment(x = 1, xend = 2, y = 3.1, yend = 3.1) +
      geom_segment(x = 1, xend = 1, y = 3.03, yend = 3.1) +
      geom_segment(x = 2, xend = 2, y = 3.03, yend = 3.1) +
      annotate("text", x = 1.5, y = 3.15, label = "*", size = 3) +
      coord_cartesian(ylim = c(min(data[,48], na.rm=TRUE)*0.95,
                               max(data[,48], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot78 (Total BA conjugation tauro/glyco)
  {
    plot78 <- ggplot(data = data, aes(x=Group5, y=data[,78])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[78], x = "", title = "Total BA conjugation\n(tauro-/glyco)") +
      geom_signif(comparisons = list(c("Control", "Vegan")),
                  test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      coord_cartesian(ylim = c(min(data[,78], na.rm=TRUE)*0.95,
                               max(data[,78], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot93 (Primary BA conjugation tauro/glyco)
  {
    plot93 <- ggplot(data = data, aes(x=Group5, y=data[,93])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[93], x = "", title = "Prim. BA conjugation\n(tauro-/glyco)") +
      geom_signif(comparisons = list(c("Control", "Vegan")),
                  test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      coord_cartesian(ylim = c(min(data[,93], na.rm=TRUE)*0.95,
                               max(data[,93], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot94 (Secondary BA conjugation tauro/glyco)
  {
    plot94 <- ggplot(data = data, aes(x=Group5, y=data[,94])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[94], x = "", title = "Sec. BA conjugation\n(tauro-/glyco)") +
      #geom_signif(comparisons = list(c("Control", "Vegan")),
      #            test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      coord_cartesian(ylim = c(min(data[,94], na.rm=TRUE)*0.95,
                               max(data[,94], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot61 (Chenodeoxycholic acid)
  {
    plot61 <- ggplot(data = data, aes(x=Group5, y=data[,61])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[61], x = "", title = miranames[61]) +
      geom_signif(comparisons = list(c("Control", "Vegan")),
                  test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      coord_cartesian(ylim = c(min(data[,61], na.rm=TRUE)*0.95,
                               max(data[,61], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot69 (Cholic acid)
  {
    plot69 <- ggplot(data = data, aes(x=Group5, y=data[,69])) +
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Control", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Control", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("OMN", "VGTR", "VGN")) +
      geom_point(position = position_jitter(width = .1, height = 0),
                 aes(colour = Group5, shape = Group5, size = Group5)) +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = .15) +
      stat_summary(geom = "errorbar", fun = mean, aes(ymin = ..y.., ymax = ..y..), width = .3) +
      labs(y = miraunits[69], x = "", title = miranames[69]) +
      geom_signif(comparisons = list(c("Control", "Vegan")),
                  test = 'permira.test.gg', map_signif_level = T, margin_top = 0.1, textsize = 3) +
      coord_cartesian(ylim = c(min(data[,69], na.rm=TRUE)*0.95,
                               max(data[,69], na.rm=TRUE)*1.2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks = element_line(colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = "none",
            plot.margin = margin(0.1,0.3,0,0.3, unit = "lines")) +
      guides(colour=F, linetype=F, size=F, shape=F)
  }
  # plot_ala (Alpha-Linolenic Acid (ion ID 416))
  {
    plot_ala <- ggplot(data = miraLIPO_norm, aes(x=Group5,
                                                 y=log2(miraLIPO_norm[,416]))) + 
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right"),
                          labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_fill_manual(name = "Dietary group",
                        values = c("#4D6B7380", "#009FC380", "#6CB43F80"),
                        guide = guide_legend(title.position = "right"),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(position = position_jitter(width = .1, height = 0), aes(color = Group5,
                                                                         shape = Group5,
                                                                         size = Group5)) + 
      geom_boxplot(aes(fill = Group5,
                       color = Group5),
                   outlier.alpha = 0,
                   coef = 10,
                   lwd = .25,
                   width = .75) +
      labs(y = "log2 FC to\nomnivore mean", x = "Alpha-linolenic acid\n(C18:3)", title = "A") +
      #geom_segment(x = 1, xend = 3, y = 1.75, yend = 1.75) +
      #geom_segment(x = 1, xend = 1, y = 1.7, yend = 1.75) +
      #geom_segment(x = 3, xend = 3, y = 1.7, yend = 1.75) +
      #annotate("text", x = 2, y = 1.9, label = "NS.", size = 3) +
      coord_cartesian(ylim = c(-2,2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = -0.3, size = 15, face = "bold"),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = 8),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks.y = element_line(colour = "black"),
            # axis.text.y = element_blank(),
            # axis.title.y = element_blank(),
            # axis.ticks.y = element_blank(),
            # axis.line.y = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none",
            #legend.position = "right",
            legend.background = element_blank(),
            legend.box.background = element_rect(colour = "black"),
            legend.text = element_text(size = 7),
            legend.title = element_blank()) +
      guides(linetype=F, size=F)
  }
  # plot_dha (Docosahexaenic acid (ion ID 500))
  {
    plot_dha <- ggplot(data = miraLIPO_norm, aes(x=Group5,
                                                 y=log2(miraLIPO_norm[,500]))) + 
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right"),
                          labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_fill_manual(name = "Dietary group",
                        values = c("#4D6B7380", "#009FC380", "#6CB43F80"),
                        guide = guide_legend(title.position = "right"),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(position = position_jitter(width = .1, height = 0), aes(color = Group5,
                                                                         shape = Group5,
                                                                         size = Group5)) + 
      geom_boxplot(aes(fill = Group5,
                       color = Group5),
                   outlier.alpha = 0,
                   coef = 10,
                   lwd = .25,
                   width = .75) +
      labs(y = "log2 FC to\nomnivore mean", x = "Docosahexaenic acid\n(C22:6)", title = "") +
      geom_segment(x = 1, xend = 3, y = 1.3, yend = 1.3) +
      geom_segment(x = 1, xend = 1, y = 1.25, yend = 1.3) +
      geom_segment(x = 3, xend = 3, y = 1.25, yend = 1.3) +
      annotate("text", x = 2, y = 1.4, label = "***", size = 3) +
      coord_cartesian(ylim = c(-2,2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = 8),
            axis.ticks.x = element_blank(),
            # axis.text.y = element_text(size = 8, colour = "black"),
            # axis.title.y = element_text(size = 8),
            # axis.ticks.y = element_line(colour = "black"),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none",
            #legend.position = "right",
            legend.background = element_blank(),
            legend.box.background = element_rect(colour = "black"),
            legend.text = element_text(size = 7),
            legend.title = element_blank()) +
      guides(linetype=F, size=F)
  }
  # plot_carnitine183 (Gamma-linolenyl carnitine, carnitine 18:3 (ion ID 615))
  {
    plot_carnitine183 <- ggplot(data = miraLIPO_norm, aes(x=Group5,
                                                          y=log2(miraLIPO_norm[,615]))) + 
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right"),
                          labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_fill_manual(name = "Dietary group",
                        values = c("#4D6B7380", "#009FC380", "#6CB43F80"),
                        guide = guide_legend(title.position = "right"),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(position = position_jitter(width = .1, height = 0), aes(color = Group5,
                                                                         shape = Group5,
                                                                         size = Group5)) + 
      geom_boxplot(aes(fill = Group5,
                       color = Group5),
                   outlier.alpha = 0,
                   coef = 10,
                   lwd = .25,
                   width = .75) +
      labs(y = "log2 FC to\nomnivore mean", x = "Carnitine 18:3", title = "B") +
      #geom_segment(x = 1, xend = 3, y = 0.65, yend = 0.65) +
      #geom_segment(x = 1, xend = 1, y = 0.6, yend = 0.65) +
      #geom_segment(x = 3, xend = 3, y = 0.6, yend = 0.65) +
      #annotate("text", x = 2, y = 0.9, label = "NS.", size = 3) +
      coord_cartesian(ylim = c(-0.75,1.5)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = -0.5, size = 15, face = "bold"),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = 8),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks.y = element_line(colour = "black"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none") +
      guides(linetype=F, size=F)
  }
  # plot_carnitine184 (Stearidonyl carnitine, carnitine 18:4 (ion ID 613))
  {
    plot_carnitine184 <- ggplot(data = miraLIPO_norm, aes(x=Group5,
                                                          y=log2(miraLIPO_norm[,613]))) + 
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right"),
                          labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_fill_manual(name = "Dietary group",
                        values = c("#4D6B7380", "#009FC380", "#6CB43F80"),
                        guide = guide_legend(title.position = "right"),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(position = position_jitter(width = .1, height = 0), aes(color = Group5,
                                                                         shape = Group5,
                                                                         size = Group5)) + 
      geom_boxplot(aes(fill = Group5,
                       color = Group5),
                   outlier.alpha = 0,
                   coef = 10,
                   lwd = .25,
                   width = .75) +
      labs(y = "log2 FC to\nomnivore mean", x = "Carnitine 18:4", title = "") +
      geom_segment(x = 1, xend = 3, y = 0.65, yend = 0.65) +
      geom_segment(x = 1, xend = 1, y = 0.6, yend = 0.65) +
      geom_segment(x = 3, xend = 3, y = 0.6, yend = 0.65) +
      annotate("text", x = 2, y = 0.75, label = "***", size = 3) +
      coord_cartesian(ylim = c(-0.75,1.5)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = 8),
            axis.ticks.x = element_blank(),
            # axis.text.y = element_text(size = 8, colour = "black"),
            # axis.title.y = element_text(size = 8),
            # axis.ticks.y = element_line(colour = "black"),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none",
            #legend.position = "right",
            legend.background = element_blank(),
            legend.box.background = element_rect(colour = "black"),
            legend.text = element_text(size = 7),
            legend.title = element_blank()) +
      guides(linetype=F, size=F)
  }
  # plot_carnitine204 (Arachidyl carnitine, carnitine 20:4 (ion ID 657))
  {
    plot_carnitine204 <- ggplot(data = miraLIPO_norm, aes(x=Group5,
                                                          y=log2(miraLIPO_norm[,657]))) + 
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right"),
                          labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_fill_manual(name = "Dietary group",
                        values = c("#4D6B7380", "#009FC380", "#6CB43F80"),
                        guide = guide_legend(title.position = "right"),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(position = position_jitter(width = .1, height = 0), aes(color = Group5,
                                                                         shape = Group5,
                                                                         size = Group5)) + 
      geom_boxplot(aes(fill = Group5,
                       color = Group5),
                   outlier.alpha = 0,
                   coef = 10,
                   lwd = .25,
                   width = .75) +
      labs(y = "log2 FC to\nomnivore mean", x = "Carnitine 20:4", title = "") +
      #geom_segment(x = 1, xend = 3, y = 1.25, yend = 1.25) +
      #geom_segment(x = 1, xend = 1, y = 1.2, yend = 1.25) +
      #geom_segment(x = 3, xend = 3, y = 1.2, yend = 1.25) +
      #annotate("text", x = 2, y = 1.45, label = "NS.", size = 3) +
      coord_cartesian(ylim = c(-0.75,1.5)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = 8),
            axis.ticks.x = element_blank(),
            # axis.text.y = element_text(size = 8, colour = "black"),
            # axis.title.y = element_text(size = 8),
            # axis.ticks.y = element_line(colour = "black"),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            #legend.position = "none",
            legend.position = "right",
            legend.background = element_blank(),
            legend.box.background = element_rect(colour = "black"),
            legend.text = element_text(size = 8),
            legend.box.margin = margin(0,0,0,0),
            legend.box.spacing = unit(0.05, "in"),
            legend.margin = margin(0.6,0.1,0.6,0.7),
            legend.key.height = unit(0.2, "in"),
            legend.spacing.x = unit(0.03, "in"),
            legend.title = element_blank()) +
      guides(linetype=F, size=F)
  }
  # plot_lysopc160 (2-Palmitoylglycerophosphocholine, LysoPC 16:0 (ion ID 698))
  {
    plot_lysopc160 <- ggplot(data = miraLIPO_norm, aes(x=Group5,
                                                       y=log2(miraLIPO_norm[,698]))) + 
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right"),
                          labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_fill_manual(name = "Dietary group",
                        values = c("#4D6B7380", "#009FC380", "#6CB43F80"),
                        guide = guide_legend(title.position = "right"),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(position = position_jitter(width = .1, height = 0), aes(color = Group5,
                                                                         shape = Group5,
                                                                         size = Group5)) + 
      geom_boxplot(aes(fill = Group5,
                       color = Group5),
                   outlier.alpha = 0,
                   coef = 10,
                   lwd = .25,
                   width = .75) +
      labs(y = "log2 FC to\nomnivore mean", x = "LysoPC 16:0", title = "C") +
      #geom_segment(x = 1, xend = 3, y = 1.25, yend = 1.25) +
      #geom_segment(x = 1, xend = 1, y = 1.2, yend = 1.25) +
      #geom_segment(x = 3, xend = 3, y = 1.2, yend = 1.25) +
      #annotate("text", x = 2, y = 1.5, label = "NS.", size = 3) +
      coord_cartesian(ylim = c(-1,2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = -0.3, size = 15, face = "bold"),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = 8),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks.y = element_line(colour = "black"),
            # axis.text.y = element_blank(),
            # axis.title.y = element_blank(),
            # axis.ticks.y = element_blank(),
            # axis.line.y = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none") +
      guides(linetype=F, size=F)
  }
  # plot_lysopc181 (2-Oleoylglycerophosphocholine, LysoPC 18:1 (ion ID 724))
  {
    plot_lysopc181 <- ggplot(data = miraLIPO_norm, aes(x=Group5,
                                                       y=log2(miraLIPO_norm[,724]))) + 
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right"),
                          labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_fill_manual(name = "Dietary group",
                        values = c("#4D6B7380", "#009FC380", "#6CB43F80"),
                        guide = guide_legend(title.position = "right"),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_segment(x = 1, xend = 3, y = 1.4, yend = 1.4) +
      geom_segment(x = 1, xend = 1, y = 1.3, yend = 1.4) +
      geom_segment(x = 3, xend = 3, y = 1.3, yend = 1.4) +
      annotate("text", x = 2, y = 1.5, label = "*", size = 3) +
      geom_point(position = position_jitter(width = .1, height = 0), aes(color = Group5,
                                                                         shape = Group5,
                                                                         size = Group5)) + 
      geom_boxplot(aes(fill = Group5,
                       color = Group5),
                   outlier.alpha = 0,
                   coef = 10,
                   lwd = .25,
                   width = .75) +
      labs(y = "log2 FC to\nomnivore mean", x = "LysoPC 18:1", title = "") +
      coord_cartesian(ylim = c(-1,2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = 8),
            axis.ticks.x = element_blank(),
            # axis.text.y = element_text(size = 8, colour = "black"),
            # axis.title.y = element_text(size = 8),
            # axis.ticks.y = element_line(colour = "black"),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none") +
      guides(linetype=F, size=F)
  }
  # plot_lysope160 (Lysophosphatidylethanolamine (LysoPE) 16:0 (ion ID 654))
  {
    plot_lysope160 <- ggplot(data = miraLIPO_norm, aes(x=Group5,
                                                       y=log2(miraLIPO_norm[,654]))) + 
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right"),
                          labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_fill_manual(name = "Dietary group",
                        values = c("#4D6B7380", "#009FC380", "#6CB43F80"),
                        guide = guide_legend(title.position = "right"),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(position = position_jitter(width = .1, height = 0), aes(color = Group5,
                                                                         shape = Group5,
                                                                         size = Group5)) + 
      geom_boxplot(aes(fill = Group5,
                       color = Group5),
                   outlier.alpha = 0,
                   coef = 10,
                   lwd = .25,
                   width = .75) +
      labs(y = "log2 FC to\nomnivore mean", x = "LysoPE 16:0", title = "") +
      #geom_segment(x = 1, xend = 3, y = 1.25, yend = 1.25) +
      #geom_segment(x = 1, xend = 1, y = 1.15, yend = 1.25) +
      #geom_segment(x = 3, xend = 3, y = 1.15, yend = 1.25) +
      #annotate("text", x = 2, y = 1.5, label = "NS.", size = 3) +
      coord_cartesian(ylim = c(-1,2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = 8),
            axis.ticks.x = element_blank(),
            # axis.text.y = element_text(size = 8, colour = "black"),
            # axis.title.y = element_text(size = 8),
            # axis.ticks.y = element_line(colour = "black"),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none") +
      guides(linetype=F, size=F)
  }
  # plot_lysope181 (Lysophosphatidylethanolamine (LysoPE) 18:1 (ion ID 686))
  {
    plot_lysope181 <- ggplot(data = miraLIPO_norm, aes(x=Group5,
                                                       y=log2(miraLIPO_norm[,686]))) + 
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right"),
                          labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_fill_manual(name = "Dietary group",
                        values = c("#4D6B7380", "#009FC380", "#6CB43F80"),
                        guide = guide_legend(title.position = "right"),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      #geom_segment(x = 1, xend = 3, y = 1.45, yend = 1.45) +
      #geom_segment(x = 1, xend = 1, y = 1.35, yend = 1.45) +
      #geom_segment(x = 3, xend = 3, y = 1.35, yend = 1.45) +
      #annotate("text", x = 2, y = 1.7, label = "NS.", size = 3) +
      geom_point(position = position_jitter(width = .1, height = 0), aes(color = Group5,
                                                                         shape = Group5,
                                                                         size = Group5)) + 
      geom_boxplot(aes(fill = Group5,
                       color = Group5),
                   outlier.alpha = 0,
                   coef = 10,
                   lwd = .25,
                   width = .75) +
      labs(y = "log2 FC to\nomnivore mean", x = "LysoPE 18:1", title = "") +
      coord_cartesian(ylim = c(-1,2)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = 8),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none") +
      guides(linetype=F, size=F)
  }
  # plot_mcfa (MCFA triglycerides combined)
  {
    plot_mcfa <- ggplot(data = miraLIPO_norm_MCFA, aes(x=Group5,
                                                       y=log2(miraLIPO_norm_MCFA[,3]))) + 
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right"),
                          labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_fill_manual(name = "Dietary group",
                        values = c("#4D6B7380", "#009FC380", "#6CB43F80"),
                        guide = guide_legend(title.position = "right"),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_segment(x = 1, xend = 3, y = 2.25, yend = 2.25) +
      geom_segment(x = 1, xend = 1, y = 2.15, yend = 2.25) +
      geom_segment(x = 3, xend = 3, y = 2.15, yend = 2.25) +
      annotate("text", x = 2, y = 2.5, label = "NA.^{1}", parse = TRUE, size = 3) +
      geom_point(position = position_jitter(width = .1, height = 0), aes(color = Group5,
                                                                         shape = Group5,
                                                                         size = Group5)) + 
      geom_boxplot(aes(fill = Group5,
                       color = Group5),
                   outlier.alpha = 0,
                   coef = 10,
                   lwd = .25,
                   width = .75) +
      labs(y = "log2 FC to\nomnivore mean", x = "MCFA triglycerides combined     \nC26:0, C31:0, C33:0     ", title = "D") +
      coord_cartesian(ylim = c(-3,2.5)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = -0.1, size = 15, face = "bold"),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = 8),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks.y = element_line(colour = "black"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none") +
      guides(linetype=F, size=F)
  }
  # plot_vlcfa (LCFA/VLCFA triglycerides combined)
  {
    plot_vlcfa <- ggplot(data = miraLIPO_norm_VLCFA, aes(x=Group5,
                                                         y=log2(miraLIPO_norm_VLCFA[,3]))) + 
      theme_light() +
      scale_colour_manual(name = "Dietary group",
                          values = c("#4D6B73FF", "#009FC3FF", "#6CB43FFF"),
                          guide = guide_legend(title.position = "right"),
                          labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_fill_manual(name = "Dietary group",
                        values = c("#4D6B7380", "#009FC380", "#6CB43F80"),
                        guide = guide_legend(title.position = "right"),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_shape_manual(name = "Dietary group", values=c(15,19,17),
                         labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_size_manual(name = "Dietary group", values=c(1,1,1),
                        labels = c("Omnivore", "Vegetarian", "Vegan")) +
      scale_x_discrete(labels = c("Omnivore", "Vegetarian", "Vegan")) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(position = position_jitter(width = .1, height = 0), aes(color = Group5,
                                                                         shape = Group5,
                                                                         size = Group5)) + 
      geom_boxplot(aes(fill = Group5,
                       color = Group5),
                   outlier.alpha = 0,
                   coef = 10,
                   lwd = .25,
                   width = .75) +
      labs(y = "log2 FC to\nomnivore mean", x = "     LCFA/VLCFA triglycerides combined\n     C50:0-o, C53:4, C58:0-o, C58:4, C58:9-o, C62:10",
           title = "") +
      geom_segment(x = 1, xend = 3, y = 2.3, yend = 2.3) +
      geom_segment(x = 1, xend = 1, y = 2.2, yend = 2.3) +
      geom_segment(x = 3, xend = 3, y = 2.2, yend = 2.3) +
      annotate("text", x = 2, y = 2.5, label = "NA.^{1}", parse = TRUE, size = 3) +
      coord_cartesian(ylim = c(-3,2.5)) +
      theme(text = element_text(size = 8, family = "Helvetica"),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = 8),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none") +
      guides(linetype=F, size=F)
  }
}

## ---------------------------------------------NOTES-------------------------------------------- ##
{
  # These subplots are saved as the individual replotting may be time-consuming. They are also open access by request from
  # the author if someone wants to try this code
  save(plot_ala, plot_dha, plot_carnitine183, plot_carnitine184, plot_carnitine204,
       plot_lysopc160, plot_lysopc181, plot_lysope160, plot_lysope181,
       plot_mcfa, plot_vlcfa,
       hplot, bplot,
       plot14, plot15, plot17, plot18, plot20, plot21, plot22, plot35, plot38, plot39,
       plot41, plot45, plot46, plot47, plot48, plot49, plot50, plot51, plot52, plot61,
       plot69, plot78, plot93, plot94, plot101, plot103,
       cholpathway, aminoplot,
       pathway_figure_pos, pathway_figure_neg,
       file = "manusc/subplots_for_figures.RData")
}