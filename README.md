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
1. FASTQC_Sample_Final_Report.txt
-contains the average of all the lanes quality scores (breaks up reads by forward and reverse)
2. FASTQC_Sample_Log_File.txt
-log file of all the averages from each lane (use for de-bugging of code ---if issues occur)


References

FastQC - https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
