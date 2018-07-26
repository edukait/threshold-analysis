# This script analyzes data from the Angicart output TSV files at different
# intensity thresholds. It calculates the scaling exponents a and b using
# five different methods.

# First, set the working directory to the folder containing the .tsv files.
setwd("/Users/kaitlinylim/Documents/angicart/data/images/mouse_brain/tCar/Microsphere_inj/MCAO day 7 microspheres-pecam peri-infarct 2 001_vd_822_822_1000_561_sm")

# Also import all the functions needed for this analysis.
source("/Users/kaitlinylim/Documents/angicart/five_methods_R/five_methods_R/get_ratio.R")
source("/Users/kaitlinylim/Documents/angicart/five_methods_R/five_methods_R/get_conservation.R")
source("/Users/kaitlinylim/Documents/angicart/five_methods_R/five_methods_R/get_distribution_exp.R")
source("/Users/kaitlinylim/Documents/angicart/five_methods_R/five_methods_R/get_distal_tips.R")
source("/Users/kaitlinylim/Documents/angicart/five_methods_R/five_methods_R/get_hierarchical_exp.R")
source("/Users/kaitlinylim/Documents/angicart/five_methods_R/five_methods_R/get_regression_exp.R")
source("/Users/kaitlinylim/Documents/angicart/five_methods_R/five_methods_R/tsvreader.R")

library(ggplot2)

# Now we save all the file names of the .tsv files we will use.
fileNames <- Sys.glob("*.tsv")
# fileNames
nameLength <- length(fileNames)

# Global variables
a_avg_ratio_all <- c()
b_avg_ratio_all <- c()

a_avg_conserv <- c()
b_avg_conserv <- c()

a_dist_values <- c()
b_dist_values <- c()

a_reg_values <- c()
b_reg_values <- c()

a_ha_values <- c()
b_ha_values <- c()

# Need to import this in order to read and manipulate the .tsv files as needed.
library("data.table")

# NOTE: Need to add code that will save all the threshold values
chosen_th <- scan("thresholds.txt", numeric(), sep = " ") * 100
# chosen_th

# Now, we will loop through each .tsv file.
for (name in 1:nameLength) {
  # This reads one file at a time.
  treeDF <- tsvreader(fileNames[name])
  # This converts the radius values to diameter values.
  # ?read.table()
  treeDF[ , c(4,5)] <- treeDF[ , c(4,5)]*2
  names(treeDF)[c(4,5)] <- c("diameter_vol", "diameter")
  
  # Now, we begin using the five methods of analysis.
  ########## BRANCHING RATIOS ###########
  ## Calculating and adding branching ratios beta and gamma to the data table.  
  # We'll do this using the get_ratio() function.
  beta_vol <- get_ratio(treeDF, treeDF$diameter_vol)
  beta <- get_ratio(treeDF, treeDF$diameter)
  gamma <- get_ratio(treeDF, treeDF$length)
  
  # Redefine data table and add new columns.
  treeDF <- data.frame(treeDF, "beta_vol" = beta_vol, "beta" = beta, "gamma" = gamma)
  
  #### Calculating a and b for all beta and gamma values.
  ## Scaling exponent a from beta.
  # First remove some bad rows and/or outliers.
  if(length(which(treeDF$beta == Inf)) > 0){
    bad_rows <- which(treeDF$beta == Inf)
    beta_vec <- treeDF$beta[-bad_rows]
  }else{
    beta_vec <- treeDF$beta
  }
  # Calculate the scaling exponent a.
  a <- -log(mean(beta_vec, na.rm = T))/log(2)
  
  # Save this value to a global vector in order to make the graph.
  a_avg_ratio_all <- c(a_avg_ratio_all, a)
  
  ## Scaling exponent b from gamma.
  # First remove some bad rows and/or outliers.
  if(length(which(treeDF$gamma >= 10)) > 0){
    bad_rows <- which(treeDF$gamma >= 10)
    gamma_vec <- treeDF$gamma[-bad_rows]
  }else{
    gamma_vec <- treeDF$gamma
  }
  # Calculating the scaling exponent b.
  b <- -log(median(gamma_vec, na.rm = T))/log(2)
  
  # Save this value to a global vector in order to make the graph.
  b_avg_ratio_all <- c(b_avg_ratio_all, b)
  
  ########## CONSERVATION BASED SCALING EXPONENTS ###########
  treeDF <- get_conservation(treeDF)
  
  ## Now calculate mean a.
  # First remove bad rows and/or outliers.
  # Note that we are temporarily focusing on positive values of q, but the negative values can also be of interest.
  if(length(which(treeDF$q >= 0 & treeDF$q <= 1000))){
    bad_rows <- which(treeDF$q >= 0 & treeDF$q <= 1000)
    q_vec <- treeDF$q[which(treeDF$q >= 0 & treeDF$q <= 1000)]
  }else{
    q_vec <- treeDF$q
  }
  
  # Calculate a.
  a <- mean(q_vec)
  
  # Save the value to a global variable.
  a_avg_conserv <- c(a_avg_conserv, a)
  
  ## Now calculate mean b.
  # First remove bad rows and/or outliers.
  # Note that we are temporarily focusing on positive values of q, but the negative values can also be of interest.
  if(length(which(treeDF$s >= 0 & treeDF$s <= 1000))){
    bad_rows <- which(treeDF$s >= 0 & treeDF$s <= 1000)
    s_vec <- treeDF$s[which(treeDF$s >= 0 & treeDF$s <= 1000)]
  }else{
    s_vec <- treeDF$s
  }
  
  # Calculate b.
  b <- mean(s_vec)
  
  # Save the value to a global variable.
  b_avg_conserv <- c(b_avg_conserv, b)
  
  ########## DISTRIBUTION BASED SCALING EXPONENT ###########
  #### Get the a and b values for all beta and gamma.
  output <- get_distribution_exp(treeDF = treeDF, binning_type = 1, plot = FALSE)
  a_dist_values <- c(a_dist_values, output[1])
  b_dist_values <- c(b_dist_values, output[4])
  
  ########## REGRESSION BASED SCALING EXPONENT ###########
  treeDF <- get_distal_tips(treeDF)
  
  #### Get the a and b values for all beta and gamma.
  output <- get_regression_exp(treeDF = treeDF, plot = FALSE)
  a_reg_values <- c(a_reg_values, output[1])
  b_reg_values <- c(b_reg_values, output[4])
  
  ########## HIERARCHICAL AVERAGING BASED SCALING EXPONENT ###########
  # NOTE: For high threshold values, this may generate this error:
    # Error in `levels<-`(`*tmp*`, value = as.character(levels)) : 
    # factor level [2] is duplicated 
  
  #### Get the a and b values for all beta and gamma.
  # treeDF
  print(chosen_th[name])
  output <- get_hierarchical_exp(treeDF = treeDF, plot = FALSE, runs = 5000, seed = 1234)
  a_ha_values <- c(a_ha_values, output[1])
  b_ha_values <- c(b_ha_values, output[4])
}

# Now, using the global variables, we will graph our results.
## Graph for get_ratio
df_ratio <- data.frame(type = rep(c("all a", "all b"), each=length(a_avg_ratio_all)), 
                       thresholds = chosen_th[1:length(a_avg_ratio_all)],
                       vals = c(a_avg_ratio_all, b_avg_ratio_all))
# df_ratio
plot <- ggplot(df_ratio, aes(x=thresholds, y=vals, group=type)) + 
  geom_line(aes(linetype=type, color=type)) + 
  geom_point(aes(color=type)) + 
  scale_linetype_manual(values=c(rep("solid", 2))) +
  scale_color_manual(values=c("red", "blue")) + 
  ggtitle("Plot of Ratio-Based scaling exponents a and b vs. threshold value") +
  labs(x="Threshold Value", y="Scaling exponents a and b") +
  theme(plot.title = element_text(size=11, hjust=0.5, face="bold")) +
  theme(legend.title=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  theme(panel.background = element_rect(fill = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size=0.25, linetype="solid", colour="black"),
        panel.grid.minor = element_line(size=0.25, linetype="solid", colour="black"))

ggsave(filename="ratio_plot.pdf", plot=plot)

## Graph for get_conservation
df_conserv <- data.frame(type=rep(c("all a", "all b"), each=length(a_avg_conserv)), 
                         thresholds = chosen_th[1:length(a_avg_conserv)],
                         vals = c(a_avg_conserv, b_avg_conserv))
# df_conserv
plot <- ggplot(df_conserv, aes(x=thresholds, y=vals, group=type)) + 
  geom_line(aes(linetype=type, color=type)) + 
  geom_point(aes(color=type)) + 
  scale_linetype_manual(values=c(rep("solid", 2))) +
  scale_color_manual(values=c("red", "blue")) + 
  ggtitle("Plot of Conservation-Based scaling exponents a and b vs. threshold value") +
  labs(x="Threshold Value", y="Scaling exponents a and b") +
  theme(plot.title = element_text(size=10, hjust=0.5, face="bold")) +
  theme(legend.title=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  theme(panel.background = element_rect(fill = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size=0.25, linetype="solid", colour="black"),
        panel.grid.minor = element_line(size=0.25, linetype="solid", colour="black"))
ggsave(filename="conserv_plot.pdf", plot=plot)

## Graph for get_distribution_exp
df_dist <- data.frame(type=rep(c("all a", "all b"), each=length(a_dist_values)), 
                      thresholds = chosen_th[1:length(a_dist_values)],
                      vals = c(a_dist_values, b_dist_values)) 

plot <- ggplot(df_dist, aes(x=thresholds, y=vals, group=type)) + 
  geom_line(aes(linetype=type, color=type)) + 
  geom_point(aes(color=type)) + 
  scale_linetype_manual(values=c(rep("solid", 2))) +
  scale_color_manual(values=c("red", "blue")) + 
  ggtitle("Plot of Distribution-Based scaling exponents a and b vs. threshold value") +
  labs(x="Threshold Value", y="Scaling exponents a and b") +
  theme(plot.title = element_text(size=10, hjust=0.5, face="bold")) +
  theme(legend.title=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  theme(panel.background = element_rect(fill = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size=0.25, linetype="solid", colour="black"),
        panel.grid.minor = element_line(size=0.25, linetype="solid", colour="black"))

ggsave(filename="dist_plot.pdf", plot=plot)

## Graph for get_regression_exp

df_reg <- data.frame(type=rep(c("all a", "all b"), each=length(a_reg_values)),
                     thresholds = chosen_th[1:length(a_reg_values)],
                     vals = c(a_reg_values, b_reg_values))
# dr_reg
plot <- ggplot(df_reg, aes(x=thresholds, y=vals, group=type)) + 
  geom_line(aes(linetype=type, color=type)) + 
  geom_point(aes(color=type)) + 
  scale_linetype_manual(values=c(rep("solid", 2))) +
  scale_color_manual(values=c("red", "blue")) + 
  ggtitle("Plot of Regression-Based scaling exponents a and b vs. threshold value") +
  labs(x="Threshold Value", y="Scaling exponents a and b") +
  theme(plot.title = element_text(size=10, hjust=0.5, face="bold")) +
  theme(legend.title=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  theme(panel.background = element_rect(fill = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size=0.25, linetype="solid", colour="black"),
        panel.grid.minor = element_line(size=0.25, linetype="solid", colour="black"))

ggsave(filename="reg_plot.pdf", plot=plot)

## Graph for get_hierarchical_exp
df_ha <- data.frame(type=rep(c("all a", "all b"), each=length(a_ha_values)),
                    thresholds = chosen_th[1:length(a_ha_values)],
                    vals = c(a_ha_values, b_ha_values))
# df_ha
plot <- ggplot(df_ha, aes(x=thresholds, y=vals, group=type)) + 
  geom_line(aes(linetype=type, color=type)) + 
  geom_point(aes(color=type)) + 
  scale_linetype_manual(values=c(rep("solid", 2))) +
  scale_color_manual(values=c("red", "blue")) + 
  ggtitle("Plot of Hierarchical Averaging scaling exponents a and b vs. threshold value") +
  labs(x="Threshold Value", y="Scaling exponents a and b") +
  theme(plot.title = element_text(size=9, hjust=0.5, face="bold")) +
  theme(legend.title=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  theme(panel.background = element_rect(fill = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size=0.25, linetype="solid", colour="black"),
        panel.grid.minor = element_line(size=0.25, linetype="solid", colour="black"))

ggsave(filename="ha_plot.pdf", plot=plot)

