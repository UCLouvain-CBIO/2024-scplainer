# Get all R script for data processing
PROCESSING_SCRIPTS := $(wildcard ./1_processing/*.R)
# Get all R script for data modelling
MODELLING_SCRIPTS := $(wildcard 2_modelling/*.R)
# Get all R script for generating the figure
FIGURE_SCRIPTS := $(wildcard 3_make_figures/*.R)

# Define the default target
# The data analysis in the paper consist of 3 main steps:
# 1. Process the data using the minimal data processing workflow
# 2. Model the data using scplainer
# 3. Make figures from the modelled data
all: processing modelling figures clean

# Rule to run the data processing
processing:
	@echo "Processing the data"
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
clean:
	@rm Rplots.pdf .Rhistory *.Rout
