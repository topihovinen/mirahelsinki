# --------------------FUNCTIONS FOR PERMUTATION TESTS IN MIRA STUDY----------------------- #

# ~~~~~~~~permira.test~~~~~~~~ #
# The actual permutation test function. Output is the observed difference/test value of
# method of choice and significance level from permutation test.
{
  permira.test <- function(data,
                               metabolite = "fP.Kol..mmol.l.",
                               comparison = c("Vegan", "Control"),
                               ageclass = c(4), #how to divide ages to strata
                               n = 47500,
                               method = "diff.mean",
                               in.class = TRUE) {
    # Check if input data is properly formatted
    # Methods for permutation test supported by this function would be
    # difference of means (diff.mean), difference of weighed class means (diff.mean.c),
    # probability index (pr.index) and hodges-lehmann-sen estimator (hls)
    proper.methods <- c("diff.mean", "diff.mean.c", "pr.index", "hls")
    {
      try(if(length(comparison) != 2)
        stop("You need to define exactly 2 groups for comparison."))
      try(if(!(metabolite %in% colnames(data)))
        stop("The metabolite is not found in your data."))
      try(if(!(method %in% proper.methods))
        stop("You have chosen an invalid method ", method))
    }
    
    # Start of the actual function
    # Import and trim the appropriate data rows to variable 'data' and print stuff for user
    {
      message("----------Permutation test for MIRA dataset----------\n")
      data = data[(data$Group5 == comparison[1] | data$Group5 == comparison[2]),]
      data$Group5 <- factor(data$Group5)
      message("Analyzing: ", metabolite, " between ", comparison[1], " and ", comparison[2])
      message("Missing data points: ", sum(is.na(data[,metabolite])))
      data <- data[!is.na(data[,metabolite]),] #Clear all NA samples from the data
    }
    
    # Prepare required variables
    {
      # Create age classes according to users wishes
      data$agegroup <- cut(data$AntAge, c(-Inf,ageclass,Inf),labels = 1:(length(ageclass)+1))
      
      # Create empty vector for differences of means to be calculated.
      # Amount of permutations is defined by parameter n for the function
      message("Number of permutations: ", n, "\n")
      nulldist <- rep(0,n)
    }
    
    # Calculate the p-values according to chosen method
    {
      message("Method: ", method)
      # Difference of means
      if(method == "diff.mean") {
        # Run the permutations and for each permutation, calculate difference of mean
        # for the null distribution.
        message("CAUTION: This permutation test algorithm (with comparison of differences of mean) has not been validated yet.\nUsing comparison of probability indexes is warmly recommended.")
        for(i in 1:n) {
          data$permgroup <- data$Group5
          for(j in 1:3) {
            for(k in c("N","M")) {
              data$permgroup[data$agegroup == j & data$Sex == k] <- sample(data$permgroup[data$agegroup == j & data$Sex == k])
            }
          }
          nulldist[i] <- mean(data[data$permgroup == "Vegan",metabolite], na.rm = T) - mean(data[data$permgroup == "Control",metabolite], na.rm = T)
          if(i %% 500 == 0) {
            message("Permutation round: ", i)
          }
        }
        
        # Calculate the real difference in means
        {
          test.value <- mean(data[data$Group5 == "Vegan",metabolite], na.rm = T) - mean(data[data$Group5 == "Control",metabolite], na.rm = T)
          # P value equals the percentage of values in null distribution as far or further from 0 as the test value.
          p.value <- (sum(nulldist >= abs(test.value))/n) + (sum(nulldist <= -abs(test.value))/n)
          message("Difference of means: ", test.value, "\tp = ", p.value)
        }
      }
      # Difference of weighed class means
      if(method == "diff.mean.c") {
        message("Method for calculating differences of weighed class means is not ready yet.\nReturning NA's")
        return(data.frame(NA,NA))
      }
      # Probability index
      if(method == "pr.index") {
        # Run the permutations and for each permutation, calculate probability index
        # for the null distribution.
        for(i in 1:n) {
          data$permgroup <- data$Group5
          if(in.class) {
            pr.index = 0
            pairamount = 0
            for(j in 1:(length(ageclass)+1)) {
              for(k in c("N","M")) {
                data$permgroup[data$agegroup == j & data$Sex == k] <- sample(data$permgroup[data$agegroup == j & data$Sex == k])
                if(in.class) { #make comparisons only in group
                  for(l in data$ID[data$permgroup == comparison[1] &
                                   data$agegroup == j &
                                   data$Sex == k]) {
                    for(m in data$ID[data$permgroup == comparison[2] &
                                     data$agegroup == j &
                                     data$Sex == k]) {
                      pairamount = pairamount + 1
                      if(data[data$ID == l, metabolite] > data[data$ID == m, metabolite]) {
                        pr.index = pr.index + 1
                      }
                    }
                  }
                }
              }
            }
          }
          if(in.class) {
            pr.index = pr.index / pairamount
          }
          # For probability index, make all possible pairings of (X_i, Y_i) and calculate
          # how high proportion of pairs has X_i > Y_i. This part is applied, if
          # comparison is done between all samples, not just the ones in the same
          # age/sex class.
          if(!in.class) {
            for(j in 1:(length(ageclass)+1)) {
              for(k in c("N","M")) {
                data$permgroup[data$agegroup == j & data$Sex == k] <- sample(data$permgroup[data$agegroup == j & data$Sex == k])
              }
            }
            pairamount <- sum(data$permgroup == comparison[1]) * sum(data$permgroup == comparison[2])
            pr.index = 0
            for(l in data$ID[data$permgroup == comparison[1]]) {
              for(m in data$ID[data$permgroup == comparison[2]]) {
                if(data[data$ID == l,metabolite] > data[data$ID == m,metabolite]) {
                  pr.index = pr.index + 1/pairamount
                }
              }
            }
          }
          nulldist[i] <- pr.index
          if(i %% 500 == 0) {
            message("Permutation round: ", i)
          }
        }
        
        # Calculate the real probability index to test.value
        {
          test.value = 0
          # Comparison with all pairs
          if(!in.class) {
            for(j in data$ID[data$Group5 == comparison[1]]) {
              for(k in data$ID[data$Group5 == comparison[2]]) {
                if(data[data$ID == j,metabolite] > data[data$ID == k,metabolite]) {
                  test.value = test.value + 1/pairamount
                }
              }
            }
          }
          # Comparison with only class-matched pairs
          if(in.class) {
            pairamount = 0
            for(j in 1:(length(ageclass)+1)) {
              for(k in c("N","M")) {
                for(l in data$ID[data$Group5 == comparison[1] &
                                 data$agegroup == j &
                                 data$Sex == k]) {
                  for(m in data$ID[data$Group5 == comparison[2] &
                                   data$agegroup == j &
                                   data$Sex == k]) {
                    pairamount = pairamount + 1
                    if(data[data$ID == l, metabolite] > data[data$ID == m, metabolite]) {
                      test.value = test.value + 1
                    }
                  }
                }
              }
            }
            test.value = test.value / pairamount
          }
          # Final test value will be (p-0.5)^2, p-value for two-sided test can
          # be assessed from the positive tail. Let us also modify the null distribution
          # accordingly (as it now includes only probability indexes)
          test.value = (test.value - 0.5)^2
          nulldist = (nulldist - 0.5)^2
          p.value = (sum(nulldist >= test.value))/n
          
          ## P value equals the percentage of values in null distribution as far or further from 0 as the test value.
          #if(test.value > 0.5) {
          #  p.value <- (sum(nulldist >= test.value))/n + (sum(nulldist <= 1-test.value))/n
          #}
          #if(test.value < 0.5) {
          #  p.value <- (sum(nulldist <= test.value))/n + (sum(nulldist >= 1-test.value))/n
          #}
          message("Test value ((p-0.5)^2): ", test.value, "\tp = ", p.value)
        }
      }
      # Hodges-Lehmann-Sen estimator
      if(method == "hls") {
        message("Method for calculating Hodges-Lehmann-Sen estimator is not ready yet.\nReturning NA's. Please use probability index instead.")
        return(data.frame(NA,NA))
      }
    }
    
    # Return a data frame with test value (actual difference in means) and p-value
    return(data.frame(test.value, p.value))
  }
}

# ~~~~~~permira.test.gg~~~~~~~ #
# An adapter code for ggplot's function geom_signif() for printing significance
# levels for the figures. Very quick and dirty solution, unfortunately, no elegance
# This is needed because geom_signif shows only two dataset vectors to the
# statistical test function, and we need information of the groups and ages, too,
# for the permutation test to function

# NOTE: This function requires precalculated vector of p-values in 'p.values' variable,
# where the index of the p-value must match the column index of the corresponding metabolite
# in the actual data file. For example, total cholesterol is found on the column 22 on the
# original data. Thus, p.values[22] is the p-value for total cholesterol.
{
  permira.test.gg <- function(data1, data2) {
    # Check if all data rows in given two datasets are included in our working data.
    # If not, try fixing the situation by importing the main data again to working data variable
    if(!(all(data %in% miradata)) & all(miradata %in% data)) {
      data = miradata
    }
    # After this check, match the given data columns (corresponding to metabolites) with
    # our working data matrix to find which metabolites geom_signif() wants to test.
    # If everything goes fine, entire column matches with one column and we use the
    # column index to match this with pre-calculated p.values vector indexes.
    # CAUTION! There's a considerable risk for mistakes in this, so please check the printed
    # figures and their significance levels.
    for(i in 1:length(colnames(data))) {
      if(all(data1 %in% data[,i])) {
        output <- data.frame(p.values.adjusted[[i]])
        # geom_signif() reads the p-value from "your_dataframe$p.value",
        # so we need to rename the column accordingly
        colnames(output) <- c("p.value")
      }
    }
    return(output)
  }
}