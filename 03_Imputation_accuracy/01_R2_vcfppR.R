# Load necessary libraries
library(vcfppR)

# Parse command-line arguments for Tool and Coverage
args <- commandArgs(trailingOnly = TRUE)
tool <- args[1]  # First argument: Tool (Glimpse2, Quilt, or Stitch)
coverage <- args[2]  # Second argument: Coverage (0.5x, 1x, or 3x)

# Define file path for truth/standard set
truthvcf <- "True_data.bcf"

# Define file paths for imputed VCFs for each tool and coverage level
# Adjust based on your analysis
file_paths <- list(
  Glimpse2 = list(
    "0.5x" = "0.5x_imputed.Glimpse2.bcf",
    "1x"   = "1x_imputed.Glimpse2.bcf",
    "3x"   = "3x_imputed.Glimpse2.bcf"
  ),
  Quilt = list(
    "0.5x" = "0.5x_imputed.Quilt.bcf",
    "1x"   = "1x_imputed.Quilt.bcf",
    "3x"   = "3x_imputed.Quilt.bcf"
  ),
  Stitch_Ref = list(
    "0.5x" = "0.5x_imputed.Stitch_Ref.bcf",
    "1x"   = "1x_imputed.Stitch_Ref.bcf",
    "3x"   = "3x_imputed.Stitch_Ref.bcf"
  ),
  Stitch = list(
    "0.5x" = "0.5x_imputed.Stitch.bcf",
    "1x"   = "1x_imputed.Stitch.bcf",
    "3x"   = "3x_imputed.Stitch.bcf"
  )
)

# Get the current imputed VCF file path
imputed_vcf <- file_paths[[tool]][[coverage]]

# Define MAF bins
maf_bins <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)

# Perform comparison for the current tool and coverage
res <- vcfcomp(test = imputed_vcf, truth = truthvcf,
               stats = "r2", bins = maf_bins, formats = c("DS","GT"))

# Convert the results to a dataframe
r2_data <- as.data.frame(res$r2)
r2_data$MAF <- rownames(res$r2)

# Add columns for tool and coverage
r2_data$Tool <- tool
r2_data$Coverage <- coverage

# Convert the MAF column to factor and assign new labels
r2_data$MAF <- factor(r2_data$MAF,
                      levels = c("[0,0.05]", "(0.05,0.1]", "(0.1,0.15]", "(0.15,0.2]",
                                 "(0.2,0.25]", "(0.25,0.3]", "(0.3,0.35]", "(0.35,0.4]",
                                 "(0.4,0.45]", "(0.45,0.5]"),
                      labels = c("0-0.05", "0.05-0.1", "0.1-0.15", "0.15-0.2",
                                 "0.2-0.25", "0.25-0.3", "0.3-0.35", "0.35-0.4",
                                 "0.4-0.45", "0.45-0.5"))

# Save the results to a CSV
output_file <- paste0("r2_results_", tool, "_", coverage, ".csv")
write.csv(r2_data, output_file, row.names = FALSE)
