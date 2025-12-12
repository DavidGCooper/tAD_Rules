# tAD_Rules
Data and analyses of the grammar rules that define transcriptional activation domains using high throughput in vivo growth experiments paired with next gen sequencing.
The following is a description of the contents of this repository.

Data Processing and Analysis Files:
tAD_Rules_Processing.rmd
- contains code used to process sequencing read counts into functionality scores for each sequence
tAD_Rules_Analysis.rmd
- contains code used to analyze sets of sequences and to create figures
RulesFeatureFunctions.R
- function for calculating score for sequence grammar rules 
run5ML-Function.R
- function for training and testing machine learning models

Directories:
Figures
- figures generated in R and included in associated manuscript
Regressions
- lasso regression models
RulesDataFrames
- output tables with sequence grammar rules scores for designed, natural, and random tAD sequences
RulesSTARprocessed
- feature (sequence) count output files from using STAR to map sequencing reads to the designed library mock-genome file

Additional files consist of intial and intermediate input as well as final output files generated while processing the sequence data.
These files are described and/or referenced within the main tAD_Rules_Processing.rmd file.
