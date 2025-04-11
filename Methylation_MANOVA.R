# Load necessary libraries
options(java.parameters = "-Xmx8g")  # For Java memory, if using rJava

if (!requireNamespace("irlba", quietly = TRUE)) install.packages("irlba")
if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel")
library(irlba)
library(parallel)
library(dplyr)
library(data.table)
library(readr)
library(tidyr)
library(limma)
library(ggplot2)
library(e1071)
library(plotly)

setwd("~/Desktop/DATA/004 MSC Characterisation/004.1 MSC sequencing Project/")
vmr_threshold <- 0.1  # Set to 0 for all CpGs
fast_pca <- FALSE #partial SVD calculations - calculates only the first 2x PCs.

#Data selector
file_identifier <- "fullepith"
age_group_id <- c("paediatric") 
tissue_id <- c("epithelium") 
gut_segment_id <- c("TI") # TI, SC
treatment_id <- c("none") #none
diagnosis_id <- c("control") # "IBDu","control", "UC", "CD", "IBDu". "control", "", IBDu "control",

#File paths
file1_path <- "~/Desktop/DATA/004 MSC Characterisation/004.1 MSC sequencing Project/raw/sample.csv"
file2_path <- "~/Desktop/DATA/004 MSC Characterisation/004.1 MSC sequencing Project/raw_methyl_set1o.csv"



######### Actual Code Below ########

# Read CSV files
if (!exists("data1")){
  data1 <- fread(file1_path)
  print("loading reference sheet")
} else {
  print("reference already loaded")
}
if (!exists("data2")){
  data2 <- fread(file2_path)
  print("loading raw data")
} else {
  print("raw data already loaded")
}

# Filter rows in data1 based on the conditions
rows_to_remove <- data1 %>%
  filter(tissue %in% tissue_id,
         age_group %in% age_group_id,
         gut_segment %in% gut_segment_id,
         treatment %in% treatment_id,
         diagnosis %in% diagnosis_id)

# Find matching columns in data2
matching_columns <- intersect(rows_to_remove$array_id, colnames(data2))

rows_to_remove <- rows_to_remove %>%
  filter(array_id %in% matching_columns)
# Create metadata rows
extracted_tissue <- c("tissue", rows_to_remove$tissue)
extracted_age_group <- c("age_group", rows_to_remove$age_group)
extracted_gut_segment <- c("gut_segment", rows_to_remove$gut_segment)
extracted_treatment <- c("treatment", rows_to_remove$treatment)
extracted_diagnosis <- c("diagnosis", rows_to_remove$diagnosis)
extracted_patient_id <- c("diagnosis", rows_to_remove$patient_id)

##


# Subset data2 based on matching columns
subset_data2 <- data2[, ..matching_columns]
subset_data2 <- data2[, c(names(data2)[1], names(subset_data2)), with = FALSE]


# Create a data.table from the existing data
# Convert extracted variables into a matrix first
extracted_data <- rbind(extracted_tissue, 
                        extracted_gut_segment, 
                        extracted_age_group, 
                        extracted_treatment,
                        extracted_diagnosis,
                        extracted_patient_id)
rownames(extracted_data) <- NULL
extracted_data <- as.data.frame(extracted_data)
colnames(extracted_data) <- colnames(subset_data2)
# Convert to a data.table and combine with the subset_data2


#VMR filter 
calculate_vmr_optimized <- function(x) {
  # Remove NA values once for the entire vector
  x_clean <- x[!is.na(x)]
  
  # Calculate mean and variance in one pass
  n <- length(x_clean)
  mean_x <- sum(x_clean) / n
  var_x <- sum((x_clean - mean_x)^2) / (n - 1)
  
  # VMR calculation
  vmr <- var_x / mean_x
  
  # Handle division by zero
  if (is.nan(vmr) || is.infinite(vmr)) {
    return(NA)
  } else {
    return(vmr)
  }
}

# For applying to multiple rows efficiently
calculate_vmr_rows <- function(mat) {
  # Apply the optimized VMR function to each row
  apply(mat, 1, calculate_vmr_optimized)
}

# Usage in your script:
if (vmr_threshold > 0) {
  vmrs <- calculate_vmr_rows(subset_data2[, -1])
  filtered_data <- subset_data2[vmrs > vmr_threshold & !is.na(vmrs), ]
} else {
  filtered_data <- subset_data2
}





#new file2
epic_path <- "EPIC-8v2-0_A1.csv"
# Function to process the CSV files
if (exists("epicv2.0_database") == FALSE){
  epicv2.0_database <- read.csv(epic_path, stringsAsFactors = FALSE)
}
# Check if the "UCSC_RefGene_Name" column exists in the second file
if (!"UCSC_RefGene_Name" %in% colnames(epicv2.0_database)) {
  stop("Column 'UCSC_RefGene_Name' not found in the second CSV file.")
}

matching_rows <- filtered_data[[1]]

##Gene checker - just make sure that the correct genes are selected. 
result <- epicv2.0_database %>%
  # Apply sub() to both columns before filtering
  filter(sub("_.*", "", epicv2.0_database[[1]]) %in% sub("_.*", "", matching_rows)) %>%
  dplyr::select(Name, UCSC_RefGene_Name, Regulatory_Feature_Group, UCSC_RefGene_Group, CHR)
result
### isolate the relevant columns and run PCA with FastPCA.R
top_cpgs <- matching_rows






############## 

#########
#trying to incorporate tissue identity
tissue_id2 <- extracted_data[1,] %>% unlist()
tissue_id2 <- tissue_id2[-1]
tissue_id2 <- as.factor(tissue_id2)

age_group_id2 <- extracted_data[3,] %>% unlist()
age_group_id2 <- age_group_id2[-1]
age_group_id2 <- as.factor(age_group_id2)

treatment_id2 <- extracted_data[4,] %>% unlist()
treatment_id2 <- treatment_id2[-1]
treatment_id2 <- as.factor(treatment_id2)

diagnosis_id2 <- extracted_data[5,] %>% unlist()
diagnosis_id2 <- diagnosis_id2[-1]
diagnosis_id2 <- as.factor(diagnosis_id2)


# Separate segment_id from the gene expression data
segment_id <- extracted_data[2,] %>% unlist()
segment_id <- segment_id[-1]
segment_id <- as.factor(segment_id)

### 

# conventional v semi-automated filter
patient_id2 <- extracted_data[6,] %>% unlist()
threshold <- 873
patient_id2 <- gsub("T", "", patient_id2)
patient_id2 <- ifelse(patient_id2 > threshold, "Semi-Automated", "Conventional")
patient_id2 <- patient_id2[-1]
patient_id2 <- as.factor(patient_id2)

expression_data <- filtered_data
cpgs <- expression_data[[1]]


expression_data_numeric <- as.data.frame(sapply(expression_data, as.numeric))
expression_data_numeric <- expression_data_numeric[-1]
rownames(expression_data_numeric) <- cpgs
#modify this to make it applicable to all variables
if (length(levels(factor(diagnosis_id2))) >= 2 ){
  design <- model.matrix(~ 0+ factor(patient_id2))
  print("Factoring in Segment & Diagnosis source")
} else {
  design <- model.matrix(~ 0+ factor(patient_id2))
  print("Factoring in Segment only")
}


fit <- lmFit(expression_data_numeric, design)


fit <- eBayes(fit)  # Apply empirical Bayes smoothing
# Get the results
results <- topTable(fit, adjust = "BH", number = Inf)

# Function to perform optimized PCA
optimized_pca <- function(data, n_components = 2, scale. = TRUE) {
  # Determine the number of cores to use (leave one core free)
  n_cores <- max(1, parallel::detectCores() - 1)
  
  # Ensure data is a matrix
  data <- as.matrix(data)
  
  # Transpose the data if necessary (samples should be rows)
  if (ncol(data) < nrow(data)) {
    data <- t(data)
  }
  
  # Handle scaling
  if (scale.) {
    # Calculate column means and standard deviations
    col_means <- colMeans(data, na.rm = TRUE)
    col_sds <- apply(data, 2, sd, na.rm = TRUE)
    
    # Replace zero standard deviations with 1 to avoid division by zero
    col_sds[col_sds == 0] <- 1
    
    # Center and scale the data
    data <- scale(data, center = col_means, scale = col_sds)
  }
  
  # Perform truncated SVD (equivalent to PCA) using multiple cores
  tryCatch({
    svd_result <- irlba(data, nv = min(n_components, ncol(data) - 1), 
                        nu = min(n_components, ncol(data) - 1),
                        center = FALSE, scale = FALSE, 
                        fastpath = TRUE, work = n_cores)
    
    # Calculate explained variance
    var_explained <- (svd_result$d^2 / sum(svd_result$d^2)) * 100
    
    # Prepare the result in a format similar to prcomp
    result <- list(
      x = data %*% svd_result$v,  # PC scores
      rotation = svd_result$v,    # Loadings
      sdev = svd_result$d / sqrt(max(1, nrow(data) - 1)),  # Standard deviations
      var_explained = var_explained,
      center = if(scale.) col_means else FALSE,
      scale = if(scale.) col_sds else FALSE
    )
    
    class(result) <- "prcomp"  # Set class for compatibility
    
    return(result)
  }, error = function(e) {
    message("Error in PCA computation: ", e$message)
    return(NULL)
  })
}


############## PCA
if (length(rownames(results)) > 900000) { 
    pca_var <- 900000
    top_cpgs <- rownames(results)[order(results$P.Value)][1:pca_var]
    paste0("PCA on ", pca_var, " CpGs")
  } else {pca_var <- length(rownames(results))
  top_cpgs <- rownames(expression_data_numeric)
  print("Running PCA on all CpGs")}


# Messing around with custom numbers
#pca_var <- 4
#msc_top_cpgs <- top_cpgs
#top_cpgs <- ("a", "b," "c")
# # Adjust the number as needed


# Messing around with custom numbers
#top_cpgs <- tail(top_cpgs, length(top_cpgs)-2)

# Subset the expression data to include only top CpGs
top_expression_data <- expression_data_numeric[top_cpgs, ]
####
#if looking for custom cpgs
top_expression_data <- na.omit(top_expression_data)
# Transpose the data for PCA (samples as rows)
pca_data <- t(top_expression_data)
# Perform PCA
if (fast_pca == TRUE)
{
  pca_result <- optimized_pca(pca_data, n_components = 2, scale = TRUE)
} else {
  pca_result <- prcomp(pca_data, scale. = TRUE)
  
}

var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
# Extract principal components
pc1 <- pca_result$x[,1]
pc2 <- pca_result$x[,2]
pca_df <- data.frame(PC1 = pc1, PC2 = pc2, Dissociation = patient_id2, Age = age_group_id2, Diagnosis = diagnosis_id2, Sample = colnames(top_expression_data))

pca_name <- paste0("PCA of ", length(top_cpgs), " top CpGs gated by VMR > ", vmr_threshold)



# Create the ggplot object (similar to before, but we'll add a tooltip)
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Dissociation, text = Sample)) + # fill = Diagnosis,
  geom_point(size = 1) +
  theme_minimal() +
  labs(title = pca_name,
       x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 2), "%)")) +
  theme(legend.position = "right") + 
  scale_fill_manual(values = c("control" = "cyan", "UC" = "yellow", "IBDu" = "orange", "CD" = "darkred"))

manova_res <- manova(cbind(pca_result$x[,1],  pca_result$x[,2],  pca_result$x[,3],  pca_result$x[,4]) ~ Dissociation, data = pca_df)
summary(manova_res)
pca_final <- data.frame(PC1 = pc1, PC2 = pc2, PC3 = pca_result$x[,3], PC4 = pca_result$x[,4], PC5 = pca_result$x[,5], PC6 = pca_result$x[,6], Dissociation = patient_id2, Age = age_group_id2, Diagnosis = diagnosis_id2, Sample = colnames(top_expression_data))


# Convert to interactive plotly object 
interactive_pca_plot <- ggplotly(pca_plot, tooltip = "text")
print(interactive_pca_plot)
htmlwidgets::saveWidget(interactive_pca_plot, "interactive_pca_plot.html")