# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggpmisc)
library(forcats)
library(RColorBrewer)

# Function to convert zero values to NA
convert_zeros_to_NA <- function(data) {
  data[data == 0] <- NA
  return(data)
}

# Function to perform and summarize linear regression
perform_regression <- function(response, predictor) {
  model <- lm(log(response) ~ log(predictor))
  return(summary(model))
}

# Read and preprocess data
mydata <- read.csv("mmv.csv")
targetindex <- read.csv("targetindex.csv")[,-1] # Remove first column

# Merge data
mydata <- merge(mydata, targetindex)

# Read and preprocess essentiality data
mydata_essentiality <- read.csv("science_mmv.csv")
names(mydata_essentiality)[names(mydata_essentiality) == "Gene_ID"] <- "Predicted_target_protein_ID"
mydata$Predicted_target_protein_ID <- substr(mydata$Predicted_target_protein_ID, 1, nchar(mydata$Predicted_target_protein_ID) - 2)
mydata <- merge(mydata, mydata_essentiality, all.x = TRUE, all.y = FALSE)

# Force specific variables to be numeric and convert zeros to NA
numeric_vars <- c("BLAST_BitScore", "LiverStage_48h_EC50", "BloodStage_48h_EC50", "BloodStage_72h_EC50", 
                  "HEK_CytotoxicConcentration50", "targetindex", "MIS", "MFS", 
                  "BLAST_PercentageIdentity", "TargetDruggabilityIndex", "ConSurf_SimilarityPercentage")
mydata[numeric_vars] <- lapply(mydata[numeric_vars], as.numeric)
mydata[numeric_vars] <- lapply(mydata[numeric_vars], convert_zeros_to_NA)

mydata$GeneEssentiality[mydata$GeneEssentiality==""] <- NA

# Perform regression analyses
regression_results <- list()
for (var in c("BloodStage_48h_EC50", "BloodStage_72h_EC50", "LiverStage_48h_EC50", "HEK_CytotoxicConcentration50")) {
  for (response_var in c("BLAST_PercentageIdentity", "BLAST_BitScore", "ConSurf_SimilarityPercentage", "targetindex", "MIS", "MFS")) {
    regression_results[[paste(response_var, var, sep = "_")]] <- perform_regression(mydata[[response_var]], mydata[[var]])
  }
}

print(regression_results)

# Create and save plots

# Reverse the factor levels for 'Essentiality'
mydata$GeneEssentiality <- factor(mydata$GeneEssentiality, 
                              levels = rev(c("Essential", "Uncertain","Slow", "No changes", "Sterile",  "Dispensable")))

p1 <- ggplot(data = na.omit(mydata), aes(x = GeneEssentiality, y = BloodStage_48h_EC50, fill = GeneEssentiality)) +
  geom_boxplot(outlier.shape = NA) +  
  theme_bw() +
  scale_fill_brewer(palette = "Blues") +
  ylab("EC50 (uM)") +
  ggtitle("A). Asexual blood stage assay EC50 (48 hrs) by essentiality of predicted protein targets")

print(p1)


p2 <- ggplot(data = na.omit(mydata), aes(x = GeneEssentiality, y = BloodStage_72h_EC50, fill = GeneEssentiality)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  scale_fill_brewer(palette = "Blues") +  # Use a blue palette with reversed order
  ylab("EC50 (uM)") +
  ggtitle("B). Asexual blood stage assay EC50 (72 hrs) by essentiality of predicted protein targets") 

print(p2)

#plot 3

# Calculate the interquartile range (IQR)
iqr <- IQR(mydata$BloodStage_48h_EC50, na.rm = TRUE)
# Calculate the first and third quartiles
q1 <- quantile(mydata$BloodStage_48h_EC50, 0.25, na.rm = TRUE)
q3 <- quantile(mydata$BloodStage_48h_EC50, 0.75, na.rm = TRUE)

# Define the y-limits with some space above the upper limit
lower_limit <- max(0, q1 - 1.5 * iqr)
upper_limit <- q3 + 1.5 * iqr

# Add a buffer to the upper limit for extra space
buffer <- (upper_limit - lower_limit) * 0.01 # 1% buffer
upper_limit_adjusted <- upper_limit + buffer

p3 <- ggplot(data = na.omit(mydata), aes(x = factor(TargetDruggabilityIndex), y = BloodStage_48h_EC50, fill = TargetDruggabilityIndex)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme_bw() +
  ylab("EC50 (uM)") +
  xlab("Druggability Index") +
  ggtitle("C). Asexual blood stage assay EC50 (48 hrs) by Druggability index") +
  coord_cartesian(ylim = c(lower_limit, upper_limit_adjusted)) # Adjust the y-axis limits with buffer

print(p3)

#plot 4

# Calculate the interquartile range (IQR)
iqr <- IQR(mydata$BloodStage_72h_EC50, na.rm = TRUE)
# Calculate the first and third quartiles
q1 <- quantile(mydata$BloodStage_72h_EC50, 0.25, na.rm = TRUE)
q3 <- quantile(mydata$BloodStage_72h_EC50, 0.75, na.rm = TRUE)

# Define the y-limits with some space above the upper limit
lower_limit <- max(0, q1 - 1.5 * iqr)
upper_limit <- q3 + 0.7 * iqr

# Add a buffer to the upper limit for extra space
buffer <- (upper_limit - lower_limit) * 0.001 # 20% buffer
upper_limit_adjusted <- upper_limit + buffer

p4 <- ggplot(data = na.omit(mydata), aes(x = factor(TargetDruggabilityIndex), y = BloodStage_72h_EC50, fill = TargetDruggabilityIndex)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme_bw() +
  ylab("EC50 (uM)") +
  xlab("Druggability Index") +
  ggtitle("D). Asexual blood stage assay EC50 (72 hrs) by Druggability index") +
  coord_cartesian(ylim = c(lower_limit, upper_limit_adjusted)) # Adjust the y-axis limits with buffer

print(p4)


# First, remove legends from individual plots
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")

# Next, extract the legends from a plot with a legend
legend1 <- get_legend(p1 + theme(legend.position = "right"))
legend2 <- get_legend(p3 + theme(legend.position = "right"))

# Finally, arrange all the components
combined_plot <- ggarrange(
  p1, p2, legend1, 
  p3, p4, legend2,
  ncol = 3, nrow = 2,
  widths = c(2, 2, 0.2), # Adjust the relative widths of plots and legends
  heights = c(1, 1)      # Adjust if necessary
)

# Print the combined plot
print(combined_plot)

ggsave("Fig 3.pdf", combined_plot, height = 11, width = 16)

#Supplementary Fig 2

# Read the dataset
mydata = read.csv("mmv.csv")

# Convert Druggability to numeric and replace zeros with NA for regression analysis
mydata$TargetDruggabilityIndex <- as.numeric(mydata$TargetDruggabilityIndex)
mydata$TargetDruggabilityIndex[which(mydata$TargetDruggabilityIndex == 0)] <- NA

# Removing duplicate Predicted_target_protein_ID for each Druggability group
distinct_data <- mydata %>%
  distinct(Predicted_target_protein_ID, TargetDruggabilityIndex, .keep_all = TRUE)

# Creating a frequency table for Druggability
Druggability_table <- as.data.frame(table(distinct_data$TargetDruggabilityIndex))
names(Druggability_table) <- c("DruggabilityIndex", "Frequency")

# Ensure consistent ordering of DruggabilityIndex
Druggability_table$DruggabilityIndex <- factor(Druggability_table$DruggabilityIndex,
                                               levels = sort(unique(distinct_data$TargetDruggabilityIndex)))

# Add a numeric version of DruggabilityIndex for the fill aesthetic in the plot
Druggability_table$DruggabilityNum <- as.numeric(as.character(Druggability_table$DruggabilityIndex))

# Plotting Druggability Frequency
plot_dg <- ggplot(Druggability_table, aes(x = DruggabilityIndex, y = Frequency)) +
  geom_bar(stat = "identity", aes(fill = DruggabilityNum), color = "black") +
  geom_text(aes(label = Frequency), vjust = -0.3, size = 3.5) +
  scale_fill_gradientn(colors = brewer.pal(9, "Blues")) +
  theme_bw() +
  labs(title = "Frequency of Druggability indices for all Pf proteins with data",
       x = "Druggability index",
       y = "Frequency") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

# Print the Druggability plot
print(plot_dg)

# Processing Essentiality data
mydata$GeneEssentiality[which(mydata$GeneEssentiality == 0)] <- NA
mydata$GeneEssentiality <- fct_collapse(mydata$GeneEssentiality, 
                                    Slow = c("slow", "Slow"), 
                                    Essential = c("Essential", "Embryonic lethal", "Larval arrest"))

# Setting factor levels for GeneEssentiality
mydata$GeneEssentiality <- factor(mydata$GeneEssentiality, 
                              levels = rev(c("Essential", "Unfit", "Uncertain", "Slow", "No changes", "Sterile", "Dispensable")))

# Removing duplicates for GeneEssentiality data
distinct_data <- mydata %>%
  distinct(Predicted_target_protein_ID, GeneEssentiality, .keep_all = TRUE)

# Creating a frequency table for GeneEssentiality
essentiality_table <- as.data.frame(table(distinct_data$GeneEssentiality))
names(essentiality_table) <- c("Essentiality", "Frequency")

# Plotting Essentiality Frequency
plot_es <- ggplot(essentiality_table, aes(x = Essentiality, y = Frequency, fill = Essentiality)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = Frequency), vjust = -0.3, size = 3.5) +
  scale_fill_brewer(palette = "Blues") +
  theme_bw() +
  labs(title = "Frequency of essentiality categories for all Pf proteins with essentiality data",
       x = "Essentiality",
       y = "Frequency") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

# Print the Essentiality plot
print(plot_es)

# Combining Druggability and Essentiality plots
combined_plot <- ggarrange(plot_es, plot_dg, ncol = 1, nrow = 2, common.legend = FALSE, legend = "right")

# Print the combined plot
print(combined_plot)

# Save the combined plot to a file
ggsave("Supplementary Fig 2.pdf", combined_plot, height = 9, width = 9)

# Correlation plots

# Reading the dataset for analysis
mydata <- read.csv("mmv_143.csv")
targetindex <- read.csv("targetindex.csv")
names(targetindex)[names(targetindex) == "Compound.Name"] <- "CompoundName"

# Merging target index data with the main dataset
mydata <- merge(mydata, targetindex, all.x = TRUE, by = "CompoundName")

# Reading essentiality data
mydata_essentiality <- read.csv("science_mmv.csv")
names(mydata_essentiality)[names(mydata_essentiality) == "Gene_ID"] <- "Predicted_target_protein_ID"
mydata$Predicted_target_protein_ID <- substr(mydata$Predicted_target_protein_ID, 1, nchar(mydata$Predicted_target_protein_ID) - 2)

# Merging essentiality data with the main dataset
mydata <- merge(mydata, mydata_essentiality, all.x = TRUE, all.y = FALSE)

# Converting various columns to numeric and replacing zeros with NAs for regression analysis
numeric_vars <- c("BLAST_BitScore", "LiverStage_48h_EC50", "BloodStage_48h_EC50", "BloodStage_72h_EC50", "MIS", "MFS", "BLAST_PercentageIdentity", "TargetDruggabilityIndex", "ConSurf_SimilarityPercentage")
mydata[numeric_vars] <- lapply(mydata[numeric_vars], as.numeric)
mydata[numeric_vars] <- lapply(mydata[numeric_vars], function(x) replace(x, x == 0, NA))

# Conducting Shapiro-Wilk tests for normality on various variables
shapiro_tests <- lapply(mydata[numeric_vars], shapiro.test)

# Define a formula for regression analysis
my.formula <- y ~ x

# Function to create regression plots
create_regression_plot <- function(data, x_var, y_var, title, y_label, x_label) {
  p <- ggplot(data, aes_string(x = paste("log(", x_var, ")", sep = ""), y = paste("log(", y_var, ")", sep = ""))) +
    geom_smooth(method = "lm", se = TRUE, formula = my.formula) +
    stat_poly_eq(formula = my.formula, label.y = 2, label.x = "right",
                 eq.with.lhs = "italic(hat(y))~`=`~",
                 aes(label = paste(..eq.label.., sep = "*plain(\",\")~")),
                 parse = TRUE) +
    stat_fit_glance(method = 'lm',
                    method.args = list(formula = my.formula),
                    geom = 'text',
                    aes(label = paste("P-value = ", signif(..p.value.., digits = 2), sep = "")),
                    label.x = 'right', label.y = .7, size = 3.5) +
    stat_poly_eq(aes(label = paste(..rr.label..)), 
                 label.x = "right", label.y = .9,
                 formula = my.formula, parse = TRUE, size = 3.5) +
    ggtitle(title) +
    labs(y = y_label, x = x_label) +
    geom_point() +
    theme(axis.line = element_line(colour = 'black', size = 0.5, linetype = 'solid'))
  return(p)
}

# Creating various regression plots
pPI_B48 <- create_regression_plot(mydata, "BLAST_PercentageIdentity", "BloodStage_48h_EC50", 
                                  "BLAST percentage identity vs EC50 at 48 hours (asexual blood stage)", 
                                  "EC50 (uM) \n (log scale)", "BLAST percentage identity (%) \n (log scale)")

pPI_B72 <- create_regression_plot(mydata, "BLAST_PercentageIdentity", "BloodStage_72h_EC50", 
                                  "BLAST percentage identity vs EC50 at 72 hours (asexual blood stage)", 
                                  "EC50 (uM) \n (log scale)", "BLAST percentage identity (%) \n (log scale)")

pBS_B48 <- create_regression_plot(mydata, "BLAST_BitScore", "BloodStage_48h_EC50", 
                                  "BLAST bit score vs EC50 at 48 hours (asexual blood stage)", 
                                  "EC50 (uM) \n (log scale)", "BLAST bit score \n (log scale)")

pBS_B72 <- create_regression_plot(mydata, "BLAST_BitScore", "BloodStage_72h_EC50", 
                                  "BLAST bit score vs EC50 at 72 hours (asexual blood stage)", 
                                  "EC50 (uM) \n (log scale)", "BLAST bit score \n (log scale)")

pCS_B48 <- create_regression_plot(mydata, "ConSurf_SimilarityPercentage", "BloodStage_48h_EC50", 
                                  "Conserved amino acid similarity vs EC50 at 48 hours (asexual blood stage)", 
                                  "EC50 (uM) \n (log scale)", "Conserved amino acid similarity \n (log scale)")

pCS_B72 <- create_regression_plot(mydata, "ConSurf_SimilarityPercentage", "BloodStage_72h_EC50", 
                                  "Conserved amino acid similarity vs EC50 at 72 hours (asexual blood stage)", 
                                  "EC50 (uM) \n (log scale)", "Conserved amino acid similarity \n (log scale)")

# Arranging all plots in a grid
all_plots <- ggarrange(pPI_B48, pPI_B72, pBS_B48, pBS_B72, pCS_B48, pCS_B72, ncol = 2, nrow = 3)

print(all_plots)
# Saving the combined regression plots to a file
ggsave("Fig 2.pdf", all_plots, height = 14, width = 12)

#                                       ---END----
  