# Get all R script for data processing
PROCESSING_SCRIPTS := $(wildcard 1_processing/*.R)
# Get all R script for data modelling
MODELLING_SCRIPTS := $(wildcard 2_modelling/*.R)
# Get all R script for generating the figure
FIGURE_SCRIPTS := $(wildcard 3_make_figures/*.R)

# Data should be stored in this folder
DATAFOLDER := data/

# Define the default target
# The data analysis in the paper consist of 3 main steps:
# 1. Process the data using the minimal data processing workflow
# 2. Model the data using scplainer
# 3. Make figures from the modelled data
all: processing modelling figures clean

# Rule to run the data processing
processing:
	@echo "Processing the data"
	@if [ ! -f $(DATAFOLDER) ]; then \
		mkdir $(DATAFOLDER); \
		echo "Folder '$(DATAFOLDER)' created."; \
	fi
	@for file in $(PROCESSING_SCRIPTS); do \
		echo "Working on $$file"; \
		R CMD BATCH $$file; \
	done

# Rule to run the data modelling
modelling:
	@echo "Modelling the data"
	@for file in $(MODELLING_SCRIPTS); do \
		echo "Working on $$file"; \
		R CMD BATCH $$file; \
	done

# Rule to generate the figure
figures:
	@echo "Generating figures"
	@for file in $(FIGURE_SCRIPTS); do \
		echo "Working on $$file"; \
		R CMD BATCH $$file; \
	done

# Remove any useless artefacts
# - Rplots.pdf: generated everytime a plot is generated
# - .Rhistory & .RData: generated when closing an R session
# -*$.Rout: log file generated by R CMD BATCH
# - cured_data.tsv: generated by HarmonizR
clean:
	@rm -f Rplots.pdf .Rhistory .RData *.Rout  cured_data.tsv
