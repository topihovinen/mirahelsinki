## ---------------------------- Pathway analysis R-Tool -------------------------- ##
library("stringr")
library("ggplot2")
library("ggrepel")
library("colorspace")
library("xlsx")
library("shiny")

# Define basic settings for the algorithm
{
  datafile = "pathwaydata.txt"	# Write your data filename with quote marks (has to be a csv/txt-file!)
  rank.method = "fc"			# This does not change anything at this point, keep it as it is
  pathwaysfile = "hmdbPtw_v3.0.mo"	# Write your pathwaylist filename with quote marks
  enrichment.type = "ions"		# Possible enrichment types are "ions" and "metabolites"
  fisher.alternative = "greater"		# Possible alternative hypotheses are "two.sided", "greater" and "less"
  print.figures = TRUE			# Do you want the script to print figures?
  table.output = TRUE			# Do you want the script to write excel-table of results?
  fc_cutoff = 0.25			# Fold change -cutoff for pathway analysis
  p_cutoff = 0.05				# p-value -cutoff for pathway analysis
  magnitude.measure = "es"		# In the figures, what is the measure of magnitude of change in pathway. Possible values "mean.change", "mean.effect" and "es".
}

{
  # Test if input parameters are valid
  try(if(!rank.method %in% c("fc", "p")) stop("Invalid rank.method"))
  try(if(!enrichment.type %in% c("metabolites", "ions", "ions_2")) stop("Invalid enrichment.type"))
  try(if(!fisher.alternative %in% c("two.sided", "greater", "less")) stop("Invalid fisher.alternative"))
  
  message("---------- PATHWAY ANALYSIS R-TOOL v0.1 ----------\n")
  message("Parameters OK.\n")
  # Read data
  {
    message("Reading data.....\n")
    diffdata <- read.csv(datafile, sep = '\t', dec = ',', header = T)
    diffdata <- data.frame(ionIdx = diffdata$ionIdx,
                           log2.FC = diffdata$log2.FC.,
                           p.value = diffdata$p.value,
                           q.value.Storey = diffdata$q.value..Storey.,
                           hmdb_id = as.character(diffdata$hmdb_id),
                           stringsAsFactors = FALSE)
    message("Finished reading data.\n\tIons found: ", length(diffdata$ionIdx))
  }
  # Read and arrange pathways
  {
    message("Reading and arranging pathways.....")
    pathways <- read.csv(pathwaysfile, header = F)
    pwy_idx <- grep(">#", pathways[,1])
    pwy_names <- as.character(pathways[pwy_idx,1])
    pwy_names <- substr(pwy_names, 3, nchar(pwy_names))
    pwy_n <- length(pwy_names)
    pwy_longest <- max(diff(pwy_idx)) - 1
    pathways_curated <- data.frame(matrix(ncol = pwy_n, nrow = pwy_longest))
    colnames(pathways_curated) <- pwy_names
    for(i in setdiff(1:length(pathways[,1]), pwy_idx)) {
      pwy_col <- sum(pwy_idx < i)
      pwy_row <- i - max(pwy_idx[pwy_idx < i])
      pathways_curated[pwy_row, pwy_col] <- substr(as.character(pathways[i,]), 5, nchar(as.character(pathways[i,])))
    }
    message("Finished arranging pathways.\n\tPathways found: ", pwy_n, "\n\tLongest pathway length: ", pwy_longest)
  }
  # Divide HMDB-codes to different lines
  {
    message("Dividing ions to metabolites.....\n")
    ions_n <- length(diffdata$hmdb_id)
    for(i in 1:ions_n) {
      ids <- strsplit(diffdata$hmdb_id[i], "; ")[[1]]
      met_n <- length(ids)
      diffdata$hmdb_id[i] <- ids[1]
      for(j in 2:met_n) {
        diffdata <- rbind(diffdata,
                          data.frame(ionIdx = diffdata$ionIdx[i],
                                     log2.FC = diffdata$log2.FC[i],
                                     p.value = diffdata$p.value[i],
                                     q.value.Storey = diffdata$q.value.Storey[i],
                                     hmdb_id = ids[j]))
      }
      # Make a timer feature so the user won't worry too much
      if(i%%(round(ions_n/10)) == 0) {
        Sys.sleep(0.1)
        message("\r", round(i/ions_n*100, 0), "%")
        flush.console()
      }
    }
    metabolites_n <- length(diffdata$hmdb_id)
    message("Finished dividing ions.\n\tMetabolites found: ", metabolites_n)
    message("Found ", sum(is.na(diffdata$hmdb_id)), " metabolites with not assigned HMDB ID. Removing them.")
    diffdata <- diffdata[!is.na(diffdata$hmdb_id),] # Remove metabolites with unknown HMDB ID as they can not be used in analysis using HMDB ID's
  }
  # For each pathway, calculate how many metabolites have a change and conduct Fisher exact test
  {
    # Find the subset of data that is changed
    pa_results <- data.frame(pathway = as.character(),
                             n.changed_pos = as.numeric(), #count sig. positively changed
                             n.changed_neg = as.numeric(), #count sig. negatively changed
                             n.found = as.numeric(), #count found metabolites of the pathway in data
                             n.tot = as.numeric(), #total amount of metabolites in pathway
                             pc.changed_pos = as.numeric(),
                             pc.changed_neg = as.numeric(),
                             p.value_pos = as.numeric(),
                             p.value_neg = as.numeric(),
                             p.value_both = as.numeric(),
                             stringsAsFactors = FALSE)
    fc_rankorder <- rank(diffdata$log2.FC[1:ions_n], na.last = TRUE)
    if(enrichment.type == "metabolites") {
      message("Finding significant changes.....")
      diffdata_pos <- diffdata[diffdata$log2.FC > fc_cutoff & diffdata$p.value < p_cutoff,]
      n.sig_changes_pos <- length(diffdata_pos$ionIdx)
      diffdata_neg <- diffdata[diffdata$log2.FC < -fc_cutoff & diffdata$p.value < p_cutoff,]
      n.sig_changes_neg <- length(diffdata_neg$ionIdx)
      message("\tSignificant positive changes found: ", n.sig_changes_pos,
              "\n\tSignificant negative changes found: ", n.sig_changes_neg)
      message("Running tests for pathways.....")
      for(i in 1:length(colnames(pathways_curated))) {
        new_pathway <- colnames(pathways_curated)[i]
        new_n.changed_pos <- sum(na.omit(pathways_curated[,i]) %in% diffdata_pos$hmdb_id)
        new_n.changed_neg <- sum(na.omit(pathways_curated[,i]) %in% diffdata_neg$hmdb_id)
        new_n.found <- sum(na.omit(pathways_curated[,i]) %in% diffdata$hmdb_id)
        new_n.tot <- length(na.omit(pathways_curated[,i]))
        new_pc.changed_pos <- new_n.changed_pos / new_n.found
        new_pc.changed_neg <- new_n.changed_neg / new_n.found
        new_p.value_pos <- fisher.test(matrix(c(new_n.changed_pos,
                                                (new_n.found - new_n.changed_pos),
                                                n.sig_changes_pos,
                                                (metabolites_n - n.sig_changes_pos)),
                                              nrow = 2),
                                       alternative = fisher.alternative)$p.value
        new_p.value_neg <- fisher.test(matrix(c(new_n.changed_neg,
                                                (new_n.found - new_n.changed_neg),
                                                n.sig_changes_neg,
                                                (metabolites_n - n.sig_changes_neg)),
                                              nrow = 2),
                                       alternative = fisher.alternative)$p.value
        new_p.value_both <- fisher.test(matrix(c((new_n.changed_neg + new_n.changed_pos),
                                                 (new_n.found - new_n.changed_neg - new_n.changed_pos),
                                                 (n.sig_changes_neg + n.sig_changes_pos),
                                                 (metabolites_n - n.sig_changes_neg - n.sig_changes_pos)),
                                               nrow = 2),
                                        alternative = fisher.alternative)$p.value
        pa_results <- rbind(pa_results,
                            data.frame(pathway = new_pathway,
                                       n.changed_pos = new_n.changed_pos,
                                       n.changed_neg = new_n.changed_neg,
                                       n.found = new_n.found,
                                       n.tot = new_n.tot,
                                       pc.changed_pos = new_pc.changed_pos,
                                       pc.changed_neg = new_pc.changed_neg,
                                       p.value_pos = new_p.value_pos,
                                       p.value_neg = new_p.value_neg,
                                       p.value_both = new_p.value_both))
      }
    }
    else if(enrichment.type == "ions") {
      # Format some new variables to result dataframe
      {
        pa_results$mean.changed_pos <- as.numeric()
        pa_results$mean.changed_neg <- as.numeric()
        pa_results$mean.effect <- as.numeric()
        pa_results$mean.effect_abs <- as.numeric()
        if(magnitude.measure == "es") {
          pa_results$es <- as.numeric()
        }
      }
      message("Finding significant changes.....")
      diffdata_pos <- diffdata[diffdata$log2.FC > fc_cutoff & diffdata$p.value < p_cutoff,]
      n.sig_changes_pos <- length(diffdata_pos$ionIdx) - sum(duplicated(diffdata_pos$ionIdx))
      diffdata_neg <- diffdata[diffdata$log2.FC < -fc_cutoff & diffdata$p.value < p_cutoff,]
      n.sig_changes_neg <- length(diffdata_neg$ionIdx) - sum(duplicated(diffdata_neg$ionIdx))
      message("\tSignificant positive changes found: ", n.sig_changes_pos,
              "\n\tSignificant negative changes found: ", n.sig_changes_neg)
      message("Running tests for pathways.....")
      for(i in 1:length(colnames(pathways_curated))) {
        new_pathway <- colnames(pathways_curated)[i]
        new_n.changed_pos <- sum(na.omit(pathways_curated[,i]) %in% diffdata_pos$hmdb_id)
        new_n.changed_neg <- sum(na.omit(pathways_curated[,i]) %in% diffdata_neg$hmdb_id)
        new_n.found <- sum(na.omit(pathways_curated[,i]) %in% diffdata$hmdb_id)
        # Curate n's by removing metabolites that occur under same ion peak
        {
          changed_pos <- which(na.omit(diffdata_pos$hmdb_id) %in% pathways_curated[,i])
          changed_neg <- which(na.omit(diffdata_neg$hmdb_id) %in% pathways_curated[,i])
          found_id <- which(na.omit(diffdata$hmdb_id) %in% pathways_curated[,i])
          found_id_curated <- found_id[!duplicated(diffdata$ionIdx[found_id])] # Remove all duplicated occurrences of the same ion index
          new_n.changed_pos <- new_n.changed_pos - sum(duplicated(diffdata_pos$ionIdx[changed_pos]))
          new_n.changed_neg <- new_n.changed_neg - sum(duplicated(diffdata_neg$ionIdx[changed_neg]))
          new_n.found <- new_n.found - sum(duplicated(diffdata$ionIdx[found_id]))
        }
        # Calculate average fold change for e.g. plotting
        {
          new_mean.changed_pos <- mean(diffdata_pos$log2.FC[changed_pos])
          if(length(found_id_curated) > 1) {
            new_mean.effect <- mean(diffdata$log2.FC[found_id_curated] / sd(diffdata$log2.FC[found_id_curated], na.rm = T))
            new_mean.effect_abs <- mean(abs(diffdata$log2.FC[found_id_curated]) / sd(diffdata$log2.FC[found_id_curated], na.rm = T))
          } else if(length(found_id_curated) <= 1) {
            new_mean.effect <- NA
            new_mean.effect_abs <- NA
          }
          new_mean.changed_neg <- mean(diffdata_neg$log2.FC[changed_neg])
        }
        # Calculate values for the final result sheet
        new_n.tot <- length(na.omit(pathways_curated[,i]))
        new_pc.changed_pos <- new_n.changed_pos / new_n.found
        new_pc.changed_neg <- new_n.changed_neg / new_n.found
        new_p.value_pos <- fisher.test(matrix(c(new_n.changed_pos,
                                                (new_n.found - new_n.changed_pos),
                                                n.sig_changes_pos,
                                                (ions_n - n.sig_changes_pos)),
                                              nrow = 2),
                                       alternative = fisher.alternative)$p.value
        new_p.value_neg <- fisher.test(matrix(c(new_n.changed_neg,
                                                (new_n.found - new_n.changed_neg),
                                                n.sig_changes_neg,
                                                (ions_n - n.sig_changes_neg)),
                                              nrow = 2),
                                       alternative = fisher.alternative)$p.value
        new_p.value_both <- fisher.test(matrix(c((new_n.changed_neg + new_n.changed_pos),
                                                 (new_n.found - new_n.changed_neg - new_n.changed_pos),
                                                 (n.sig_changes_neg + n.sig_changes_pos),
                                                 (ions_n - n.sig_changes_neg - n.sig_changes_pos)),
                                               nrow = 2),
                                        alternative = fisher.alternative)$p.value
        # Include also Enrichment Score calculation
        if(magnitude.measure == "es") {
          # Calculate ES_calc, that has the size and direction of each step in enrichment score calculation. Then by cumulative sum conduct the actual ES to ES_result
          pathwayhit_idx <- diffdata$ionIdx[found_id_curated]
          diffdata_sorted <- diffdata[order(diffdata$log2.FC),] #Sort diffdata based on fold change. This does not affect reducing metabolite-level data to ion-level data, as metabolites under the same ion peak will sort next to each other
          if(enrichment.type == "ions") {
            diffdata_sorted <- diffdata_sorted[!duplicated(diffdata_sorted$ionIdx),] # Remove duplicates to measure only ions
            ions_n_sorted <- length(diffdata_sorted$ionIdx) #Refresh ions_n if there was a difference after removing duplicated rows
            pathwayhit_idx_sorted <- which(diffdata_sorted$ionIdx %in% pathwayhit_idx)
            ES_calc <- rep(0, ions_n_sorted) #Create empty vector for calculating step sizes
            ES_calc[pathwayhit_idx_sorted] <- abs(diffdata_sorted$log2.FC[pathwayhit_idx_sorted]) # Set ES_calc step heights for "hits" based on the knowledge of their ionIdx's and fold changes
            average_stepdown <- sum(abs(diffdata_sorted$log2.FC[pathwayhit_idx_sorted])) / (ions_n_sorted - length(pathwayhit_idx_sorted))
            ES_calc[-pathwayhit_idx_sorted] <- -average_stepdown
            ES_result <- cumsum(ES_calc) / max(c((new_n.changed_neg + new_n.changed_pos),1))
            ES <- max(abs(ES_result))
            if(ES > max(pa_results$es)) {
              ES_max_idx <- logical(length = ions_n_sorted) #Create a vector that would have TRUE in the index corresponding maximal ES of all pathways
              ES_max_idx[pathwayhit_idx_sorted] <- TRUE
              ES_max_function <- data.frame(fES = ES_result,
                                            x = c(1:ions_n_sorted),
                                            idx = ES_max_idx) # Save ES function for plotting from the highest ES metabolite
              ES_max_name <- new_pathway
            }
          }
          # else {
          #  ES_calc[pathwayhit_idx] <- diffdata$log2.FC[pathwayhit_idx]
          #  average_stepdown <- sum(diffdata$log2.FC[pathwayhit_idx]) / (ions_n - length(found_id_curated))
          #  ES_calc[-pathwayhit_idx] <- -average_stepdown
          #  ES_result <- cumsum(ES_calc)
          #  ES <- max(abs(ES_result))
          #}
          #if(ES > max(pa_results$es)) {
          #  ES_max_idx <- logical(length = ions_n) #Create a vector that would have TRUE in the index corresponding maximal ES of all pathways
          #  ES_max_idx[pathwayhit_idx] <- TRUE
          #  ES_max_function <- data.frame(fES = ES_result,
          #                                x = c(1:ions_n),
          #                                idx = ES_max_idx) # Save ES function for plotting from the highest ES metabolite
          #  ES_max_name <- new_pathway
          #}
          pa_results <- rbind(pa_results,
                              data.frame(pathway = new_pathway,
                                         n.changed_pos = new_n.changed_pos,
                                         n.changed_neg = new_n.changed_neg,
                                         n.found = new_n.found,
                                         n.tot = new_n.tot,
                                         pc.changed_pos = new_pc.changed_pos,
                                         pc.changed_neg = new_pc.changed_neg,
                                         p.value_pos = new_p.value_pos,
                                         p.value_neg = new_p.value_neg,
                                         p.value_both = new_p.value_both,
                                         mean.changed_pos = new_mean.changed_pos,
                                         mean.changed_neg = new_mean.changed_neg,
                                         mean.effect = new_mean.effect,
                                         mean.effect_abs = new_mean.effect_abs,
                                         es = ES))
        }
        else {
          pa_results <- rbind(pa_results,
                              data.frame(pathway = new_pathway,
                                         n.changed_pos = new_n.changed_pos,
                                         n.changed_neg = new_n.changed_neg,
                                         n.found = new_n.found,
                                         n.tot = new_n.tot,
                                         pc.changed_pos = new_pc.changed_pos,
                                         pc.changed_neg = new_pc.changed_neg,
                                         p.value_pos = new_p.value_pos,
                                         p.value_neg = new_p.value_neg,
                                         p.value_both = new_p.value_both,
                                         mean.changed_pos = new_mean.changed_pos,
                                         mean.changed_neg = new_mean.changed_neg,
                                         mean.effect = new_mean.effect,
                                         mean.effect_abs = new_mean.effect_abs))
        }
      }
    }
    message("Pathway testing done.\nMost changed pathways\n\tPositive change: ",
            pa_results$pathway[which.min(pa_results$p.value_pos)],
            " (p-value ", pa_results$p.value_pos[which.min(pa_results$p.value_pos)],
            ")\n\tNegative change: ",
            pa_results$pathway[which.min(pa_results$p.value_neg)],
            " (p-value ", pa_results$p.value_neg[which.min(pa_results$p.value_neg)],
            ")\n\tAny change: ",
            pa_results$pathway[which.min(pa_results$p.value_both)],
            " (p-value ", pa_results$p.value_both[which.min(pa_results$p.value_both)], ")")
    # For figures and tables, correct some names of pathways as the import gets messed up with ','s in pathway names
    levels(pa_results$pathway)[levels(pa_results$pathway) == "Valine"] <- "Valine, Leucine and Isoleucine degradation"
  }
  # Print figures
  if(print.figures) {
  # Find 5 most extreme values
    {
        extreme5_idx_pos <- data.frame(pathwayname = rep(NA, 5),
                                       dist = rep(0,5),
                                       stringsAsFactors = FALSE)
        extreme5_idx_neg <- data.frame(pathwayname = rep(NA, 5),
                                       dist = rep(0,5),
                                       stringsAsFactors = FALSE)
        if(magnitude.measure == "mean.effect") {
          for(i in 1:length(pa_results$pathway)) {
            dist_pos <- (((-log10(pa_results$p.value_pos[i])) / max(-log10(pa_results$p.value_pos), na.rm = T))^2 + (pa_results$mean.effect[i] / max(pa_results$mean.effect, na.rm = T))^2)
            if(dist_pos > min(extreme5_idx_pos$dist) & !is.na(dist_pos)) {
              extreme5_idx_pos$pathwayname[which.min(extreme5_idx_pos$dist)] <- pa_results$pathway[i]
              extreme5_idx_pos$dist[which.min(extreme5_idx_pos$dist)] <- dist_pos
            }
            dist_neg <- (((-log10(pa_results$p.value_neg[i])) / max(-log10(pa_results$p.value_neg), na.rm = T))^2 + (pa_results$mean.effect[i] / max(pa_results$mean.effect, na.rm = T))^2)
            if(dist_neg > min(extreme5_idx_neg$dist) & !is.na(dist_neg)) {
              extreme5_idx_neg$pathwayname[which.min(extreme5_idx_neg$dist)] <- pa_results$pathway[i]
              extreme5_idx_neg$dist[which.min(extreme5_idx_neg$dist)] <- dist_neg
            }
          }
        }
        else if(magnitude.measure == "mean.change") {
          for(i in 1:length(pa_results$pathway)) {
            dist_pos <- (((-log10(pa_results$p.value_pos[i])) / max(-log10(pa_results$p.value_pos), na.rm = T))^2 + (pa_results$mean.changed_pos[i] / max(pa_results$mean.changed_pos, na.rm = T))^2)
            if(dist_pos > min(extreme5_idx_pos$dist) & !is.na(dist_pos)) {
              extreme5_idx_pos$pathwayname[which.min(extreme5_idx_pos$dist)] <- pa_results$pathway[i]
              extreme5_idx_pos$dist[which.min(extreme5_idx_pos$dist)] <- dist_pos
            }
            dist_neg <- (((-log10(pa_results$p.value_neg[i])) / max(-log10(pa_results$p.value_neg), na.rm = T))^2 + (pa_results$mean.changed_neg[i] / max(pa_results$mean.changed_neg, na.rm = T))^2)
            if(dist_neg > min(extreme5_idx_neg$dist) & !is.na(dist_neg)) {
              extreme5_idx_neg$pathwayname[which.min(extreme5_idx_neg$dist)] <- pa_results$pathway[i]
              extreme5_idx_neg$dist[which.min(extreme5_idx_neg$dist)] <- dist_neg
            }
          }
        }
        else if(magnitude.measure == "es") {
          for(i in 1:length(pa_results$pathway)) {
            dist_pos <- (((-log10(pa_results$p.value_pos[i])) / max(-log10(pa_results$p.value_pos), na.rm = T))^2 + (pa_results$es[i] / max(pa_results$es, na.rm = T))^2)
            if(dist_pos > min(extreme5_idx_pos$dist) & !is.na(dist_pos)) {
              extreme5_idx_pos$pathwayname[which.min(extreme5_idx_pos$dist)] <- pa_results$pathway[i]
              extreme5_idx_pos$dist[which.min(extreme5_idx_pos$dist)] <- dist_pos
            }
            dist_neg <- (((-log10(pa_results$p.value_neg[i])) / max(-log10(pa_results$p.value_neg), na.rm = T))^2 + (pa_results$es[i] / max(pa_results$es, na.rm = T))^2)
            if(dist_neg > min(extreme5_idx_neg$dist) & !is.na(dist_neg)) {
              extreme5_idx_neg$pathwayname[which.min(extreme5_idx_neg$dist)] <- pa_results$pathway[i]
              extreme5_idx_neg$dist[which.min(extreme5_idx_neg$dist)] <- dist_neg
            }
          }
        }
      }
    message("Printing figures.....")
    if(magnitude.measure == "mean.effect") {
        pathway_figure_pos <- ggplot(data = pa_results,
                                     aes(x = mean.effect,
                                         y = -log10(p.value_pos))) +
          scale_color_gradient(low = "#000000", high = "#ff0000") +
          theme_light() +
          geom_point(aes(colour = -log10(p.value_pos),
                         size = mean.effect)) +
          geom_text(data = pa_results[extreme5_idx_pos$pathwayname,], 
                    aes(label = pathway, x = mean.effect, y = -log10(p.value_pos)), size = 3, hjust = -.05) +
          guides(size = "none",
                 colour = "none") +
          labs(x = "Mean effect size",
               y = "-log10 p-value of Fisher test")
        pathway_figure_neg <- ggplot(data = pa_results,
                                     aes(x = mean.effect,
                                         y = -log10(p.value_neg))) +
          scale_color_gradient(low = "#000000", high = "#ff0000") +
          theme_light() +
          geom_point(aes(colour = -log10(p.value_neg),
                         size = mean.effect)) +
          geom_text(data = pa_results[extreme5_idx_neg$pathwayname,], 
                    aes(label = pathway, x = mean.effect, y = -log10(p.value_neg)), size = 3, hjust = -.05) +
          guides(size = "none",
                 colour = "none") +
          labs(x = "Mean effect size",
               y = "-log10 p-value of Fisher test")
        ggsave(file = "Pathway analysis (positive change).pdf", plot = pathway_figure_pos, device = "pdf", width = 6, height = 4)
        ggsave(file = "Pathway analysis (positive change).png", plot = pathway_figure_pos, device = "png", width = 6, height = 4)
        ggsave(file = "Pathway analysis (negative change).pdf", plot = pathway_figure_neg, device = "pdf", width = 6, height = 4)
        ggsave(file = "Pathway analysis (negative change).png", plot = pathway_figure_neg, device = "png", width = 6, height = 4)
      }
    else if(magnitude.measure == "mean.change") {
        pathway_figure_pos <- ggplot(data = pa_results,
                                     aes(x = mean.changed_pos,
                                         y = -log10(p.value_pos))) +
          scale_color_gradient(low = "#000000", high = "#ff0000") +
          theme_light() +
          geom_point(aes(colour = -log10(p.value_pos),
                         size = mean.changed_pos)) +
          geom_text(data = pa_results[extreme5_idx_pos$pathwayname,], 
                    aes(label = pathway, x = mean.changed_pos, y = -log10(p.value_pos)), size = 3, hjust = -.05) +
          guides(size = "none",
                 colour = "none") +
          labs(x = "Mean of significant positive changes",
               y = "-log10 p-value of Fisher test")
        pathway_figure_neg <- ggplot(data = pa_results,
                                     aes(x = mean.changed_neg,
                                         y = -log10(p.value_neg))) +
          scale_color_gradient(low = "#000000", high = "#ff0000") +
          theme_light() +
          geom_point(aes(colour = -log10(p.value_neg),
                         size = mean.changed_neg)) +
          geom_text(data = pa_results[extreme5_idx_neg$pathwayname,], 
                    aes(label = pathway, x = mean.changed_neg, y = -log10(p.value_neg)), size = 3, hjust = -.05) +
          guides(size = "none",
                 colour = "none") +
          labs(x = "Mean of significant negative changes",
               y = "-log10 p-value of Fisher test")
        ggsave(file = "results/Pathway analysis (positive change).pdf", plot = pathway_figure_pos, device = "pdf", width = 6, height = 4)
        ggsave(file = "results/Pathway analysis (positive change).png", plot = pathway_figure_pos, device = "png", width = 6, height = 4)
        ggsave(file = "results/Pathway analysis (negative change).pdf", plot = pathway_figure_neg, device = "pdf", width = 6, height = 4)
        ggsave(file = "results/Pathway analysis (negative change).png", plot = pathway_figure_neg, device = "png", width = 6, height = 4)
      }
    else if(magnitude.measure == "es") {
      pathway_figure_pos <- ggplot(data = pa_results,
                                   aes(x = es,
                                       y = -log10(p.value_pos))) +
        scale_color_gradient(low = "#000000", high = "#ff0000") +
        theme_light() +
        geom_point(aes(colour = -log10(p.value_pos),
                       size = es)) +
        geom_text_repel(data = pa_results[extreme5_idx_pos$pathwayname[extreme5_idx_pos$dist %in% sort(extreme5_idx_pos$dist)[3:5]],], # A very ugly implementation to actually use only three most remarkably changed pathways in labeling
                  aes(label = pathway, x = es, y = -log10(p.value_pos)),
                  size = 3,
                  hjust = -.05,
                  box.padding = 1.5,
                  nudge_x = 0.2,
                  nudge_y = 0.2) +
        guides(size = "none",
               colour = "none") +
        labs(x = "Enrichment score",
             y = "-log10 p-value of Fisher test",
             title = "Higher metabolite levels in vegans than in omnivores (+)") +
        theme(plot.title = element_text(size = 11))
      pathway_figure_neg <- ggplot(data = pa_results,
                                   aes(x = es,
                                       y = -log10(p.value_neg))) +
        scale_color_gradient(low = "#000000", high = "#ff0000") +
        theme_light() +
        geom_point(aes(colour = -log10(p.value_neg),
                       size = es)) +
        geom_text_repel(data = pa_results[extreme5_idx_neg$pathwayname,], 
                  aes(label = pathway, x = es, y = -log10(p.value_neg)),
                  size = 3,
                  hjust = -.05,
                  box.padding = 1.5) +
        guides(size = "none",
               colour = "none") +
        labs(x = "Enrichment score",
             y = "-log10 p-value of Fisher test",
             title = "Lower metabolite levels in vegans than in omnivores (-)") +
        theme(plot.title = element_text(size = 11))
      es_figure1 <- ggplot(data = ES_max_function,
                           aes(x = x,
                               y = fES)) +
        theme_minimal() +
        theme(panel.grid = element_blank(),
              axis.line = element_line()) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
        geom_line(color = "red") +
        xlim(0, max(ES_max_function$x)) +
        labs(y = "Enrichment Score ES(x)",
             x = "<< Negative changes | Positive changes >>",
             title = paste(" Enrichment score function for\n", ES_max_name))
    }
    ggsave(file = "results/Pathway analysis (positive change).pdf", plot = pathway_figure_pos, device = "pdf", width = 6, height = 4)
    ggsave(file = "results/Pathway analysis (positive change).png", plot = pathway_figure_pos, device = "png", width = 6, height = 4)
    ggsave(file = "results/Pathway analysis (negative change).pdf", plot = pathway_figure_neg, device = "pdf", width = 6, height = 4)
    ggsave(file = "results/Pathway analysis (negative change).png", plot = pathway_figure_neg, device = "png", width = 6, height = 4)
    ggsave(file = "results/ES plot for highest ES.pdf", plot = es_figure1, device = "pdf", width = 6, height = 4)
    ggsave(file = "results/ES plot for highest ES.png", plot = es_figure1, device = "png", width = 6, height = 4)
  }
  # Form output tables
  if(table.output) {
    if(magnitude.measure == "es") {
      colnames(pa_results) <- c("Pathway name",
                                "Sig. positive changes",
                                "Sig. negative changes",
                                "Found ions/metabolites",
                                "Total ions/metabolites",
                                "Positive change %",
                                "Negative change %",
                                "p-value (positive changes)",
                                "p-value (negative changes)",
                                "p-value (all changes)",
                                "Mean of positive sig. changes",
                                "Mean of negative sig. changes",
                                "Mean effect size of all changes",
                                "Mean absolute effect size of all changes",
                                "Enrichment score")
    }
    else {
      colnames(pa_results) <- c("Pathway name",
                                "Sig. positive changes",
                                "Sig. negative changes",
                                "Found ions/metabolites",
                                "Total ions/metabolites",
                                "Positive change %",
                                "Negative change %",
                                "p-value (positive changes)",
                                "p-value (negative changes)",
                                "p-value (all changes)",
                                "Mean of positive sig. changes",
                                "Mean of negative sig. changes",
                                "Mean effect size of all changes",
                                "Mean absolute effect size of all changes")
    }
    write.xlsx(pa_results, "results/Pathway Analysis results.xlsx")
  }
}
