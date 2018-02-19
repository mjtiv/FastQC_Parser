# FastQC_Parser
Parses base quality scores from fastqc for an overall summary report

Note: Code was Developed in the Abasht Laboratory at the University of Delaware under the 
supervision of Dr. Behnam Abasht website: "http://canr.udel.edu/faculty/behnam-abasht/"

This program quickly walks through a directories of FastQC results and extracts the average quality scores 
from the "fastqc_data.txt" file. The "fastqc_data.txt" file contains 
data tables of all the information used to visualize the html version of the results.

Note: This code works with older FastQC output where two folders are produced (zipped and unzipped data). If
using a newer version of FastQC slight modification of code will be needed to extract zipped data. 

The reason this script was needed was to try and understand overall global behavior of samples and try to
quantify the visualization of the data seen in the html files.

#Required Files for Running Script
1. FASTQC_Parameter_File.txt (open up file and change directory information---self-explanatory)
2. FASTQC_Parsing.py

#Running Program
Open up the program and run in python, output will put in the same directory as program

#Output Files from Program
Various FASTQC outputs from the statistical analysis performed by FastQC (mean, median, lower quartile, upper quartile
tenth percentile, ninetieth percentile and summary report file). In total seven files are outputed. A directory is created by the program called "Log_File_Directory," which contains all the log files for all the statistical tests performed by the program (6 files).

References

FastQC - https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
