library("xlsx")

# Import nutritional data
miradata_ravinto <- read.csv("muuttujat_analyysiin.txt", sep = '\t', dec = '.', header = T)
# For consistency and convenience, we include only subjects with blood sample available
miradata_ravinto <- miradata_ravinto[miradata_ravinto$ID %in% miradata$ID,]
# Make weight-, age-, sex- and diet group data consistent with miradata
miradata_ravinto$Sex <- miradata_ravinto$sukupuoli # For sex, we need to simply change the numeric "sukupuoli" variable to binary "M"/"N" for the ad-hoc permira.test to work
miradata_ravinto$Sex[miradata_ravinto$Sex == 1] <- "N"
miradata_ravinto$Sex[miradata_ravinto$Sex == 2] <- "M"
miradata_ravinto$Group5 <- miradata_ravinto$crv_yhd_5
miradata_ravinto$Group5[miradata_ravinto$Group5 == 1] <- "Vegan"
miradata_ravinto$Group5[miradata_ravinto$Group5 == 2] <- "Vegetarian"
miradata_ravinto$Group5[miradata_ravinto$Group5 == 5] <- "Control"
for(i in 1:length(miradata_ravinto$ID)) {
  miradata_ravinto$Weigth[i] <- miradata$Weight[miradata$ID == miradata_ravinto$ID[i]]
  miradata_ravinto$AntAge[i] <- miradata$AntAge[miradata$ID == miradata_ravinto$ID[i]]
}
# Choose variables for the table in their corresponding order
final_nutrients <- c("ENERJ_per_d",
                     "PROT_e_pros",
                     "CHO_e_pros",
                     "FAT_e_pros",
                     "FASAT_e_pros",
                     "FAMS_e_pros",
                     "FAPU_e_pros",
                     "FATRN_per_d",
                     "F18D2CN6_per_d",
                     "F18D2CN6_per_MJ",
                     "F18D3N3_per_d",
                     "F18D3N3_per_MJ",
                     "rs_F20D5N3_per_d",
                     "F20D5N3_per_d",
                     "saantil_EPA",
                     "F20D5N3_per_MJ",
                     "rs_F22D6N3_per_d",
                     "F22D6N3_per_d",
                     "saantil_DHA",
                     "F22D6N3_per_MJ",
                     "CHOL_per_d",
                     "CHOL_per_MJ",
                     "SUCS_per_d",
                     "FIBC_per_d",
                     "FIBC_per_MJ",
                     "rs_THIA_per_d",
                     "THIA_per_d",
                     "saantil_Tiamiini",
                     "THIA_per_MJ",
                     "rs_RIBF_per_d",
                     "RIBF_per_d",
                     "saantil_Riboflaviini",
                     "RIBF_per_MJ",
                     "rs_NIAEQ_per_d",
                     "NIAEQ_per_d",
                     "saantil_Niasiini",
                     "NIAEQ_per_MJ",
                     "rs_VITB6_per_d",
                     "VITB6_per_d",
                     "saantil_Pyridoksiini",
                     "VITB6_per_MJ",
                     "rs_FOL_per_d",
                     "FOL_per_d",
                     "saantil_Foolihappo",
                     "FOL_per_MJ",
                     "rs_VITB12_per_d",
                     "VITB12_per_d",
                     "saantil_B12vitamiini",
                     "VITB12_per_MJ",
                     "rs_VITC_per_d",
                     "VITC_per_d",
                     "saantil_Cvitamiini",
                     "VITC_per_MJ",
                     "rs_VITA_per_d",
                     "VITA_per_d",
                     "saantil_RAE_tot",
                     "VITA_per_MJ",
                     "rs_VITD_per_d",
                     "VITD_per_d",
                     "saantil_Dvit_tot",
                     "VITD_per_MJ",
                     "rs_VITE_per_d",
                     "VITE_per_d",
                     "saantil_Evitamiini",
                     "VITE_per_MJ",
                     "VITK_per_d",
                     "VITK_per_MJ",
                     "NAT_per_d",
                     "NAT_per_MJ",
                     "rs_CA_per_d",
                     "CA_per_d",
                     "saantil_Kalsium",
                     "CA_per_MJ",
                     "K_per_d",
                     "K_per_MJ",
                     "rs_MG_per_d",
                     "MG_per_d",
                     "saantil_Magnesium",
                     "MG_per_MJ",
                     "rs_FE_per_d",
                     "FE_per_d",
                     "saantil_Rauta",
                     "FE_per_MJ",
                     "JODI_per_d",
                     "rs_JODI_per_d",
                     "saantil_Jodi",
                     "JODI_per_MJ",
                     "P_per_d",
                     "P_per_MJ",
                     "rs_ZN_per_d",
                     "ZN_per_d",
                     "saantil_Sinkki",
                     "ZN_per_MJ")
final_nutrients_units <- c("kJ", #Energy
                           "%", #Protein
                           "%", #Carbohydrates
                           "%", #Fat
                           "%", #Sat. fat
                           "%", #Monounsaturated fatty acids
                           "%", #Polyunsaturated fatty acids
                           "g", #Trans fatty acids
                           "g", #Linolic acid
                           "g/MJ",
                           "g", #Alphalinolenic acid
                           "g/MJ",
                           "mg", #EPA
                           "mg",
                           "mg",
                           "mg/MJ",
                           "mg", #DHA
                           "mg",
                           "mg",
                           "mg/MJ",
                           "g", #Cholesterol
                           "g/MJ",
                           "g", #Sucrose
                           "g", #Fiber
                           "g/MJ",
                           "mg", #B1-vit
                           "mg",
                           "mg",
                           "mg/MJ",
                           "mg", #B2-vit
                           "mg",
                           "mg",
                           "mg/MJ",
                           "mg", #B3-vit
                           "mg",
                           "mg",
                           "mg/MJ",
                           "mg", #B6-vit
                           "mg",
                           "mg",
                           "mg/MJ",
                           "ug", #Folate
                           "ug",
                           "ug",
                           "ug/MJ",
                           "ug", #B12-vit
                           "ug",
                           "ug",
                           "ug/MJ",
                           "mg", #C-vit
                           "mg",
                           "mg",
                           "mg/MJ",
                           "ug", #A-vit
                           "ug",
                           "ug",
                           "ug/MJ",
                           "ug", #D-vit
                           "ug",
                           "ug",
                           "ug/MJ",
                           "mg", #E-vit
                           "mg",
                           "mg",
                           "mg/MJ",
                           "ug", #K-vit
                           "ug/MJ",
                           "mg", #Sodium
                           "mg/MJ",
                           "mg", #Calcium
                           "mg",
                           "mg",
                           "mg/MJ",
                           "mg", #Potassium
                           "mg/MJ",
                           "mg", #Magnesium
                           "mg",
                           "mg",
                           "mg/MJ",
                           "mg", #Iron
                           "mg",
                           "mg",
                           "mg/MJ",
                           "ug", #Iodine
                           "ug",
                           "ug",
                           "ug/MJ",
                           "mg", #Phosphorus
                           "mg/MJ",
                           "mg", #Zinc
                           "mg",
                           "mg",
                           "mg/MJ")
final_nutrients_names <- c("Daily energy",
                           "Protein",
                           "Carbohydrates",
                           "Fat",
                           "Saturated fatty acids",
                           "Monounsaturated fatty acids",
                           "Polyunsaturated fatty acids",
                           "Trans fatty acids",
                           "Linoleic acid",
                           "\t per daily energy intake",
                           "Alpha-linolenic acid",
                           "\t per daily energy intake",
                           "EPA (20:5 n-3)",
                           "\t of which from food",
                           "\t of which from supplements",
                           "\t per daily energy intake",
                           "DHA (22:6 n-6)",
                           "\t of which from food",
                           "\t of which from supplements",
                           "\t per daily energy intake",
                           "Cholesterol",
                           "\t per daily energy intake",
                           "Sucrose",
                           "Fiber",
                           "\t per daily energy intake",
                           "Thiamine (B1)",
                           "\t of which from food",
                           "\t of which from supplements",
                           "\t per daily energy intake",
                           "Riboflavin (B2)",
                           "\t of which from food",
                           "\t of which from supplements",
                           "\t per daily energy intake",
                           "Niacin equivalents (B3)",
                           "\t of which from food",
                           "\t of which from supplements",
                           "\t per daily energy intake",
                           "Pyridoxine (B6)",
                           "\t of which from food",
                           "\t of which from supplements",
                           "\t per daily energy intake",
                           "Folate (B9)",
                           "\t of which from food",
                           "\t of which from supplements",
                           "\t per daily energy intake",
                           "Cobalamine (B12)",
                           "\t of which from food",
                           "\t of which from supplements",
                           "\t per daily energy intake",
                           "Vitamin C",
                           "\t of which from food",
                           "\t of which from supplements",
                           "\t per daily energy intake",
                           "Vitamin A**",
                           "\t of which from food",
                           "\t of which from supplements",
                           "\t per daily energy intake",
                           "Vitamin D",
                           "\t of which from food",
                           "\t of which from supplements",
                           "\t per daily energy intake",
                           "Vitamin E",
                           "\t of which from food",
                           "\t of which from supplements",
                           "\t per daily energy intake",
                           "Vitamin K",
                           "\t per daily energy intake",
                           "Sodium",
                           "\t per daily energy intake",
                           "Calcium",
                           "\t of which from food",
                           "\t of which from supplements",
                           "\t per daily energy intake",
                           "Potassium",
                           "\t per daily energy intake",
                           "Magnesium",
                           "\t of which from food",
                           "\t of which from supplements",
                           "\t per daily energy intake",
                           "Iron",
                           "\t of which from food",
                           "\t of which from supplements",
                           "\t per daily energy intake",
                           "Iodine",
                           "\t of which from food",
                           "\t of which from supplements",
                           "\t per daily energy intake",
                           "Phosphorus",
                           "\t per daily energy intake",
                           "Zinc",
                           "\t of which from food",
                           "\t of which from supplements",
                           "\t per daily energy intake")

# Run permutation tests on chosen nutrients
{
  #The first row is for starting comparison from zero
  #miradata_ravinto_pvalues <- data.frame(Nutrient = character(),
  #                                       p.value = double())
  for(i in 1:length(final_nutrients)) {
    nutriname <- final_nutrients[i]
    if(!(nutriname %in% miradata_ravinto_pvalues$Nutrient)) {
      new_row <- data.frame(Nutrient = nutriname,
                            p.value = permira.test(data = miradata_ravinto,
                                                   metabolite = nutriname,
                                                   comparison = c("Vegan", "Control"),
                                                   ageclass = c(4),
                                                   n = 47500,
                                                   method = "pr.index",
                                                   in.class = TRUE)$p.value,
                            p.value.corr = 0)
      miradata_ravinto_pvalues <- rbind(miradata_ravinto_pvalues,
                                        new_row)
    }
  }
}
# Benjamini-Hochberg correction for multiple testing
miradata_ravinto_pvalues$p.value.corr <- p.adjust(miradata_ravinto_pvalues$p.value, method = "BH")

# Combine data into single table in .xlsx form
{
  miradata_ravinto_medians <- data.frame(final_nutrients_names)
  for(i in 1:length(miradata_ravinto_medians$final_nutrients_names)) {
    nutri <- colnames(miradata_ravinto) == final_nutrients[i]
    miradata_ravinto_medians$Nutrient[i] <- paste(final_nutrients_names[i], "-", final_nutrients_units[i])
    miradata_ravinto_medians$Omnivore[i] <- paste(signif(median(miradata_ravinto[miradata_ravinto$crv_yhd_5 == 5,nutri], na.rm = T), digits = 3),
                                                  "\n[",
                                                  signif(min(miradata_ravinto[miradata_ravinto$crv_yhd_5 == 5,nutri], na.rm = T), digits = 3),
                                                  "-",
                                                  signif(max(miradata_ravinto[miradata_ravinto$crv_yhd_5 == 5,nutri], na.rm = T), digits = 3),
                                                  "]")
    miradata_ravinto_medians$Vegetarian[i] <- paste(signif(median(miradata_ravinto[miradata_ravinto$crv_yhd_5 == 2,nutri], na.rm = T), digits = 3),
                                                    "\n[",
                                                    signif(min(miradata_ravinto[miradata_ravinto$crv_yhd_5 == 2,nutri], na.rm = T), digits = 3),
                                                    "-",
                                                    signif(max(miradata_ravinto[miradata_ravinto$crv_yhd_5 == 2,nutri], na.rm = T), digits = 3),
                                                    "]")
    miradata_ravinto_medians$Vegan[i] <- paste(signif(median(miradata_ravinto[miradata_ravinto$crv_yhd_5 == 1,nutri], na.rm = T), digits = 3),
                                               "\n[",
                                               signif(min(miradata_ravinto[miradata_ravinto$crv_yhd_5 == 1,nutri], na.rm = T), digits = 3),
                                               "-",
                                               signif(max(miradata_ravinto[miradata_ravinto$crv_yhd_5 == 1,nutri], na.rm = T), digits = 3),
                                               "]")
    if(sum(miradata_ravinto_pvalues$Nutrient == final_nutrients[i])) {
      miradata_ravinto_medians$pvalues[i] <- miradata_ravinto_pvalues$p.value.corr[miradata_ravinto_pvalues$Nutrient == final_nutrients[i]]
    } else {
      miradata_ravinto_medians$pvalues[i] <- NA
    }
  }
  miradata_ravinto_medians <- miradata_ravinto_medians[,2:6] #There was problems with renaming nutrient names in for loop, so a quickfix is to have two nutrient name -columns and remove the other one here
  colnames(miradata_ravinto_medians) <- c("Nutrients",
                                          paste("Omnivore\nN =", sum(miradata_ravinto$crv_yhd_5 == 5)),
                                          paste("Vegetarian\nN =", sum(miradata_ravinto$crv_yhd_5 == 2)),
                                          paste("Vegan\nN =", sum(miradata_ravinto$crv_yhd_5 == 1)),
                                          "p-value*\nOmnivore vs. Vegan")
  write.xlsx(miradata_ravinto_medians, "nutrient_data.xlsx")
}
