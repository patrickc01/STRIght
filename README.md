# STRight
Human and Mouse STR analysis using NGS

## Installation and Requirements
STRight requires Python3 and the pandas package

## Usage
To use STRight.py to analyze human STRs:
1. Copy STRight.py and STRight_key.csv to the folder containing your joined STR fastq files
2. Run STRight.py

To use STRight.py to analyze mouse STRs:
1. Copy STRight_mouse.py and STRight_mouse_key.csv to the folder containing your joined STR fastq files
2. Run STRight_mouse.py


## Output
Two files are created in the run folder:  
1.  STR_results_df.csv : contains the fastq file names and the STR found, total reads and top 4 STR lengths and fractions found in each file.
2.  STR_analysis.txt : contains a summary of analyzed fastq files along with the top STR sequences that were found in each fastq file.
