# Load required libraries
library(ggplot2)
library(dplyr)

# Complete exposure category mapping
exposure_category_mapping <- c(
  # Brain-Related
  "Chronotype" = "Brain-Related",
  "Fluid Intelligence Score" = "Brain-Related",
  "Neuroticism" = "Brain-Related",
  "Schizophrenia" = "Brain-Related",
  "Sleep Duration" = "Brain-Related",
  "Sleeplessness or Insomnia" = "Brain-Related",
  
  # Cardiovascular
  "Atrial Fibrillation" = "Cardiovascular",
  "Coronary Artery Disease" = "Cardiovascular",
  "Self-Reported Hypertension" = "Cardiovascular",
  
  # Dermatological
  "Basal Cell Carcinoma" = "Dermatological",
  "Hair or Balding Pattern 2" = "Dermatological",
  "Hair or Balding Pattern 3" = "Dermatological",
  "Hair or Balding Pattern 4" = "Dermatological",
  "Self-Reported Eczema or Dermatitis" = "Dermatological",
  "Self-Reported Gout" = "Dermatological",
  "Self-Reported Psoriasis" = "Dermatological",
  "Self-Reported Rheumatoid Arthritis" = "Dermatological",
  
  # Gastrointestinal
  "Crohn's Disease" = "Gastrointestinal",
  "Inflammatory Bowel Disease" = "Gastrointestinal",
  "Ulcerative Colitis" = "Gastrointestinal",
  
  # Immunological
  "Diagnosed Asthma" = "Immunological",
  "Eosinophil Count" = "Immunological",
  "Granulocyte Count" = "Immunological",
  "Hayfever or Allergic Rhinitis" = "Immunological",
  "High Light Scatter Reticulocyte Count" = "Immunological",
  "Lymphocyte Count" = "Immunological",
  "Monocyte Count" = "Immunological",
  "Myeloid White Cell Count" = "Immunological",
  "Neutrophil Count" = "Immunological",
  "Platelet Count" = "Immunological",
  "Red Blood Cell Count" = "Immunological",
  "Reticulocyte Count" = "Immunological",
  "Self-Reported Ankylosing Spondylitis" = "Immunological",
  "Self-Reported Asthma" = "Immunological",
  "Self-Reported Multiple Sclerosis" = "Immunological",
  "Sum of Basophil and Neutrophil Counts" = "Immunological",
  "Sum of Eosinophil and Basophil Counts" = "Immunological",
  "Sum of Neutrophil and Eosinophil Counts" = "Immunological",
  "Systemic Lupus Erythematosus" = "Immunological",
  "White Blood Cell Count" = "Immunological",
  
  # Metabolic
  "Birth Weight" = "Metabolic",
  "Body Fat Percentage" = "Metabolic",
  "ER Negative Breast Cancer" = "Metabolic",
  "IDL Triglycerides" = "Metabolic",
  "LDL Cholesterol" = "Metabolic",
  "Overall Breast Cancer" = "Metabolic",
  "Self-Reported High Cholesterol" = "Metabolic",
  "Self-Reported Hyperthyroidism" = "Metabolic",
  "Self-Reported Hypothyroidism" = "Metabolic",
  "Self-Reported Type 1 Diabetes" = "Metabolic",
  "Type 2 Diabetes" = "Metabolic"
)

# Numeric mapping for exposures
exposure_numeric_mapping <- c(
  # Chronotype and Cognitive Traits
  "Chronotype" = 1,
  "Fluid Intelligence Score" = 2,
  "Neuroticism" = 3,
  "Schizophrenia" = 4,
  "Sleep Duration" = 5,
  "Sleeplessness or Insomnia" = 6,
  
  # Cardiovascular Conditions
  "Atrial Fibrillation" = 7,
  "Coronary Artery Disease" = 8,
  "Self-Reported Hypertension" = 9,
  
  # Dermatological Traits
  "Basal Cell Carcinoma" = 10,
  "Hair or Balding Pattern 2" = 11,
  "Hair or Balding Pattern 3" = 12,
  "Hair or Balding Pattern 4" = 13,
  "Self-Reported Eczema or Dermatitis" = 14,
  "Self-Reported Gout" = 15,
  "Self-Reported Psoriasis" = 16,
  "Self-Reported Rheumatoid Arthritis" = 17,
  
  # Gastrointestinal Conditions
  "Crohn's Disease" = 18,
  "Inflammatory Bowel Disease" = 19,
  "Ulcerative Colitis" = 20,
  
  # Immunological/Blood Traits
  "Diagnosed Asthma" = 21,
  "Eosinophil Count" = 22,
  "Granulocyte Count" = 23,
  "Hayfever or Allergic Rhinitis" = 24,
  "High Light Scatter Reticulocyte Count" = 25,
  "Lymphocyte Count" = 26,
  "Monocyte Count" = 27,
  "Myeloid White Cell Count" = 28,
  "Neutrophil Count" = 29,
  "Platelet Count" = 30,
  "Red Blood Cell Count" = 31,
  "Reticulocyte Count" = 32,
  "Self-Reported Ankylosing Spondylitis" = 33,
  "Self-Reported Asthma" = 34,
  "Self-Reported Multiple Sclerosis" = 35,
  "Sum of Basophil and Neutrophil Counts" = 36,
  "Sum of Eosinophil and Basophil Counts" = 37,
  "Sum of Neutrophil and Eosinophil Counts" = 38,
  "Systemic Lupus Erythematosus" = 39,
  "White Blood Cell Count" = 40,
  
  # Metabolic Traits
  "Birth Weight" = 41,
  "Body Fat Percentage" = 42,
  "ER Negative Breast Cancer" = 43,
  "IDL Triglycerides" = 44,
  "LDL Cholesterol" = 45,
  "Overall Breast Cancer" = 46,
  "Self-Reported High Cholesterol" = 47,
  "Self-Reported Hyperthyroidism" = 48,
  "Self-Reported Hypothyroidism" = 49,
  "Self-Reported Type 1 Diabetes" = 50,
  "Type 2 Diabetes" = 51
)

# Add outcome and type columns for each dataset
result_amyloid <- read.csv("~/Dropbox/sexbiasedmr/data/adpathology/int2MR_AD_result/AD_result_amyloid.assoc.linear.gz.csv", sep="")
result_amyloid <- result_amyloid %>%
  mutate(category = exposure_category_mapping[exposure])
result_gpath <- read.csv("~/Dropbox/sexbiasedmr/data/adpathology/int2MR_AD_result/AD_result_gpath.assoc.linear.gz.csv", sep="")
result_gpath <- result_gpath %>%
  mutate(category = exposure_category_mapping[exposure])
result_tangles <- read.csv("~/Dropbox/sexbiasedmr/data/adpathology/int2MR_AD_result/AD_result_tangles.assoc.linear.gz.csv", sep="")
result_tangles <- result_tangles %>%
  mutate(category = exposure_category_mapping[exposure])
result_dlbany <- read.csv("~/Dropbox/sexbiasedmr/data/adpathology/int2MR_AD_result/AD_result_dlbany.assoc.logistic.gz.csv", sep="")
result_dlbany <- result_dlbany %>%
  mutate(category = exposure_category_mapping[exposure])
result_hspath_typ <- read.csv("~/Dropbox/sexbiasedmr/data/adpathology/int2MR_AD_result/AD_result_hspath_typ.assoc.logistic.gz.csv", sep="")
result_hspath_typ <- result_hspath_typ %>%
  mutate(category = exposure_category_mapping[exposure])
result_tdp_st4_binary <- read.csv("~/Dropbox/sexbiasedmr/data/adpathology/int2MR_AD_result/AD_result_tdp_st4_binary.assoc.logistic.gz.csv", sep="")
result_tdp_st4_binary <- result_tdp_st4_binary %>%
  mutate(category = exposure_category_mapping[exposure])

result_amyloid$outcome <- "amyloid"
result_gpath$outcome <- "gpath"
result_tangles$outcome <- "tangles"
result_amyloid$type    <- "AD pathology"
result_gpath$type      <- "AD pathology"
result_tangles$type    <- "AD pathology"
result_dlbany$outcome  <- "dlbany"
result_hspath_typ$outcome <- "hspath_typ"
result_tdp_st4_binary$outcome <- "tdp_st4_binary"
result_dlbany$type     <- "Non-AD pathology"
result_hspath_typ$type <- "Non-AD pathology"
result_tdp_st4_binary$type <- "Non-AD pathology"

# Combine all results into one dataframe
result_heatmap <- rbind(
  result_amyloid, result_gpath, result_tangles,
  result_dlbany, result_hspath_typ, result_tdp_st4_binary
)

# Compute Z-score columns
result_heatmap$Zscore                 <- result_heatmap$est_beta_int   / result_heatmap$se_beta_int
result_heatmap$Zscore_beta           <- result_heatmap$est_beta       / result_heatmap$se_beta
result_heatmap$Zscore_beta_plus_int  <- result_heatmap$est_beta_plus_int / result_heatmap$se_beta_plus_int

# Assign categories and numeric order to exposures
result_heatmap$category   <- exposure_category_mapping[result_heatmap$exposure]
result_heatmap$categoryNum<- exposure_numeric_mapping[result_heatmap$exposure]
result_heatmap <- na.omit(result_heatmap)

# Define colors for each category
category_colors <- c(
  "Brain-Related"   = "cyan",
  "Cardiovascular"  = "purple",
  "Dermatological"  = "orange",
  "Gastrointestinal"= "blue",
  "Immunological"   = "red",
  "Metabolic"       = "green"
)

# Order exposures and assign colors
result_heatmap <- result_heatmap %>% arrange(categoryNum)
result_heatmap$exposure_colors <- category_colors[result_heatmap$category]
result_heatmap$exposure <- factor(
  result_heatmap$exposure,
  levels = unique(result_heatmap$exposure[order(result_heatmap$categoryNum)])
)
exposure_colors <- category_colors[result_heatmap$category]

# Calculate category boundaries and centers for vertical separators
category_counts <- table(exposure_category_mapping)
category_boundaries <- cumsum(category_counts) + 0.5
category_centers    <- cumsum(category_counts) - category_counts/2 + 0.5

# Specify outcome order
result_heatmap$outcome <- factor(
  result_heatmap$outcome,
  levels = c("tangles", "gpath", "amyloid", "tdp_st4_binary", "hspath_typ", "dlbany")
)

# Determine boundaries between AD and non-AD pathologies
ad_pathology     <- c("tangles", "gpath", "amyloid")
non_ad_pathology <- c("tdp_st4_binary", "hspath_typ", "dlbany")
boundary_y       <- which(levels(result_heatmap$outcome) == non_ad_pathology[1]) - 0.5

# x-coordinate for right-side bracket
x_bracket <- max(as.numeric(result_heatmap$exposure)) + 1
x_label   <- x_bracket + 0.5  # shift label right

# Create the heatmap
p <- ggplot(result_heatmap, aes(x = exposure, y = outcome, fill = Zscore)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "Z-score"
  ) +
  geom_vline(xintercept = category_boundaries, linetype = "dashed") +
  geom_hline(yintercept = boundary_y, linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.text.x   = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5, margin = margin(t = 5)),
    axis.text.y   = element_text(size = 6),
    legend.position = "top",
    legend.title   = element_text(size = 8),
    legend.text    = element_text(size = 6),
    panel.grid     = element_blank(),
    plot.margin    = margin(t = 50, r = 120, b = 100, l = 20)
  ) +
  labs(x = "Exposure", y = "Outcome", title = "Heatmap of Z Score") +
  coord_cartesian(clip = "off")

# Add right-side bracket and label for AD pathology
ad_start <- which(levels(result_heatmap$outcome) == ad_pathology[1])
ad_end   <- which(levels(result_heatmap$outcome) == ad_pathology[3])
p <- p +
  annotate("segment", x = x_bracket, xend = x_bracket, y = ad_start, yend = ad_end) +
  annotate("segment", x = x_bracket, xend = x_bracket - 0.5, y = ad_start, yend = ad_start) +
  annotate("segment", x = x_bracket, xend = x_bracket - 0.5, y = ad_end,   yend = ad_end)   +
  annotate(
    "text", x = x_label, y = mean(c(ad_start, ad_end)),
    label = "AD pathology", size = 2.25, fontface = "bold", angle = 270, hjust = 0.5
  )

# Add right-side bracket and label for Non-AD pathology
non_ad_start <- which(levels(result_heatmap$outcome) == non_ad_pathology[1])
non_ad_end   <- which(levels(result_heatmap$outcome) == non_ad_pathology[3])
p <- p +
  annotate("segment", x = x_bracket, xend = x_bracket, y = non_ad_start, yend = non_ad_end) +
  annotate("segment", x = x_bracket, xend = x_bracket - 0.5, y = non_ad_start, yend = non_ad_start) +
  annotate("segment", x = x_bracket, xend = x_bracket - 0.5, y = non_ad_end,   yend = non_ad_end)   +
  annotate(
    "text", x = x_label, y = mean(c(non_ad_start, non_ad_end)),
    label = "Non-AD pathology", size = 2.25, fontface = "bold", angle = 270, hjust = 0.5
  )

# Add bottom brackets and labels for each category
category_positions <- split(result_heatmap$exposure, result_heatmap$category)
y_bracket <- -0.5
y_label_cat <- y_bracket - 0.7

for (cat in names(category_positions)) {
  exposures <- category_positions[[cat]]
  x1 <- min(which(levels(result_heatmap$exposure) %in% exposures))
  x2 <- max(which(levels(result_heatmap$exposure) %in% exposures))
  
  p <- p +
    annotate("segment", x = x1, xend = x2, y = y_bracket, yend = y_bracket) +
    annotate("segment", x = x1, xend = x1, y = y_bracket, yend = y_bracket + 0.5) +
    annotate("segment", x = x2, xend = x2, y = y_bracket, yend = y_bracket + 0.5) +
    annotate(
      "text", x = (x1 + x2)/2, y = y_label_cat,
      label = cat, size = 2.5, fontface = "italic", vjust = 1
    )
}

# Display the plot
print(p)