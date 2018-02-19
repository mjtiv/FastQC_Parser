#!/usr/bin/env python3.6

# Modules used in the code
import numpy as np
import os


# function opens up a fastqc_data.txt file and begins to
# parse apart the results for the qc scores
def get_fastqc_data():

    # opening fastqc results
    fastqc_results = open('fastqc_data.txt', 'r')

    # list to store Per Base Quality Scores
    data = []

    # creating on/off switches for reading data
    read_data = 'stop'

    for line in fastqc_results:
        if line.startswith('#Base'):
            read_data = 'start'
            pass

        elif line.startswith(">>END_MODULE"):
            read_data = 'stop'

        # once the next module is reached kill the loop
        elif line.startswith('>>Per sequence quality scores'):
            break

        elif read_data == 'start':
            data.append(line)

    fastqc_results.close()
    return data


def get_stats_of_data(fastqc_results, position):

    # Get average from fastqc
    quality_scores = []

    for data in fastqc_results:
        data = data.split("\t")
        quality_scores.append(float(data[position]))

    first_quarter_index = int(0.25*len(quality_scores))
    last_quarter_index = int(0.75*len(quality_scores))

    # Getting Stats from Quality Scores
    # print(quality_scores)
    avg_qual_score = np.mean(quality_scores)

    # middle blocks (double)
    first_block = quality_scores[:first_quarter_index]
    avg_first_block = np.mean(first_block)

    middle_blocks = quality_scores[first_quarter_index:last_quarter_index]
    avg_middle_blocks = np.mean(middle_blocks)

    # last block of data
    last_block = quality_scores[last_quarter_index:]
    avg_last_block = np.mean(last_block)

    return {'avg_qual_score': avg_qual_score, 'avg_first_block': avg_first_block,
            'avg_middle_blocks': avg_middle_blocks, 'avg_last_block': avg_last_block}


# reads in the parameter file to see what directories need to be parsed
# for analysis
def read_parameter_file():

    # Opens the parameter file to get all the required inputs for the rest of the code
    parameter_file = open('FASTQC_Parameter_File.txt', 'r')  # Open the file in python

    # List to store directories
    directories = []

    # creating on/off switches for reading parameter
    read_data = 'stop'

    for line in parameter_file:
        if line.startswith('DIRECTORIES'):
            read_data = 'start'

        elif read_data == 'start':
            line = line.rstrip()
            directories.append(line)
        else:
            pass

    # Close the initial file
    parameter_file.close()

    return directories


# reads in the specific
def parse_and_print_results(report_file, log_file, list_of_directories, position):

    report_file.write("Sample\tR1_1st_Block_Avg\tR1_Middle_Blocks_Avg\tR1_Final_Block_Avg\tR1_Avg\t"
                      "R2_1st_Block_Avg\tR2_Middle_Blocks_Avg\tR2_Final_Block_Avg\tR2_Block_Avg\t"
                      "Overall_Sample_Avg\n")

    log_file.write("Sample\tR1_1st_Block_Avg\tR1_Middle_Blocks_Avg\tR1_Final_Block_Avg\tR1_Avg\t"
                   "R2_1st_Block_Avg\tR2_Middle_Blocks_Avg\tR2_Final_Block_Avg\tR2_Block_Avg\t"
                   "Overall_Sample_Avg\n")

    for i in range(len(list_of_directories)):

        directory_pathway = list_of_directories[i]

        directory_list = []
        found_directories = os.listdir(list_of_directories[i])
        for found_directory in found_directories:
            if os.path.isdir(found_directory):
                if found_directory.startswith('Log'):
                    pass
                else:
                    directory_list.append(found_directory)
                    # print(found_directory)

        # Loop searches through the directories and retrieves all the data from the various sub-directories
        # Sample specific lane runs from FASTQC
        for directory in directory_list:

            # Parsing the directory name to get the sample name
            sample_name = directory

            report_file.write(sample_name+"\t")
            log_file.write(sample_name+"\t")

            # R1 lanes from sample
            sample_r1_avg_first_block = []
            sample_r1_avg_middle_blocks = []
            sample_r1_avg_final_block = []
            sample_r1_avg = []

            # R2 lanes from sample
            sample_r2_avg_first_block = []
            sample_r2_avg_middle_blocks = []
            sample_r2_avg_final_block = []
            sample_r2_avg = []

            # Overall Average
            sample_avg = []

            # Changes the directory of the program to go inside the sample directory
            new_directory = directory_pathway + "/" + directory + "/"
            os.chdir(new_directory)

            # In directory of ALL lane runs for a sample
            fastqc_directories = os.listdir(new_directory)

            for fastqc_directory in fastqc_directories:

                if fastqc_directory.endswith(".zip"):
                    pass

                elif fastqc_directory.endswith("_fastqc"):
                    final_directory_step = new_directory + fastqc_directory + "/"

                    # Change to the final directory step (where data is at)
                    os.chdir(final_directory_step)

                    if 'R1' in fastqc_directory:
                        # Parsing through each samples fastqc lane results
                        fastqc_results = get_fastqc_data()

                        # Getting the lane specific stats
                        stats_on_lane = get_stats_of_data(fastqc_results, position)
                        lane_avg_qual_score = stats_on_lane['avg_qual_score']
                        lane_avg_first_block = stats_on_lane['avg_first_block']
                        lane_avg_middle_blocks = stats_on_lane['avg_middle_blocks']
                        lane_avg_last_block = stats_on_lane['avg_last_block']

                        # Adding to the sample list (creating averages of averages)
                        sample_r1_avg_first_block.append(lane_avg_first_block)
                        sample_r1_avg_middle_blocks.append(lane_avg_middle_blocks)
                        sample_r1_avg_final_block.append(lane_avg_last_block)
                        sample_r1_avg.append(lane_avg_qual_score)

                        # Overall Average
                        sample_avg.append(lane_avg_qual_score)

                    elif 'R2' in fastqc_directory:
                        # Parsing through each samples fastqc lane results
                        fastqc_results = get_fastqc_data()

                        # Getting the lane specific stats
                        stats_on_lane = get_stats_of_data(fastqc_results, position)
                        lane_avg_qual_score = stats_on_lane['avg_qual_score']
                        lane_avg_first_block = stats_on_lane['avg_first_block']
                        lane_avg_middle_blocks = stats_on_lane['avg_middle_blocks']
                        lane_avg_last_block = stats_on_lane['avg_last_block']

                        # Adding to the sample list (creating averages of averages)
                        sample_r2_avg_first_block.append(lane_avg_first_block)
                        sample_r2_avg_middle_blocks.append(lane_avg_middle_blocks)
                        sample_r2_avg_final_block.append(lane_avg_last_block)
                        sample_r2_avg.append(lane_avg_qual_score)

                        # Overall Average
                        sample_avg.append(lane_avg_qual_score)

                    else:
                        print("Error with parsing R1 and R2")

                # Tells the program to skip other not important files
                else:
                    pass

            # Convert lists per lane to strings for log file
            sample_avg_str_list = (','.join(map(str, sample_avg)))

            # R1 lists per lane to strings for log file
            sample_r1_avg_first_block_str_list = (','.join(map(str, sample_r1_avg_first_block)))
            sample_r1_avg_middle_blocks_str_list = (','.join(map(str, sample_r1_avg_middle_blocks)))
            sample_r1_avg_final_block_str_list = (','.join(map(str, sample_r1_avg_final_block)))
            sample_r1_avg_str_list = (','.join(map(str, sample_r1_avg)))

            # R2 lists per lane to strings for log file
            sample_r2_avg_str_list = (','.join(map(str, sample_r2_avg)))
            sample_r2_avg_first_block_str_list = (','.join(map(str, sample_r2_avg_first_block)))
            sample_r2_avg_middle_blocks_str_list = (','.join(map(str, sample_r2_avg_middle_blocks)))
            sample_r2_avg_final_block_str_list = (','.join(map(str, sample_r2_avg_final_block)))

            # Writing R1 Lane Info to Log File
            log_file.write(sample_r1_avg_first_block_str_list+"\t")
            log_file.write(sample_r1_avg_middle_blocks_str_list + "\t")
            log_file.write(sample_r1_avg_final_block_str_list + "\t")
            log_file.write(sample_r1_avg_str_list + "\t")

            # Writing R2 Lane Info to Log File
            log_file.write(sample_r2_avg_first_block_str_list + "\t")
            log_file.write(sample_r2_avg_middle_blocks_str_list + "\t")
            log_file.write(sample_r2_avg_final_block_str_list + "\t")
            log_file.write(sample_r2_avg_str_list + "\t")

            # Writing Overall Average to file
            log_file.write(sample_avg_str_list + "\n")

            # Writing sample results to a file
            report_file.write(str(np.mean(sample_r1_avg_first_block))+'\t')
            report_file.write(str(np.mean(sample_r1_avg_middle_blocks)) + '\t')
            report_file.write(str(np.mean(sample_r1_avg_final_block)) + '\t')
            report_file.write(str(np.mean(sample_r1_avg)) + '\t')

            # Writing sample results to a file
            report_file.write(str(np.mean(sample_r2_avg_first_block)) + '\t')
            report_file.write(str(np.mean(sample_r2_avg_middle_blocks)) + '\t')
            report_file.write(str(np.mean(sample_r2_avg_final_block)) + '\t')
            report_file.write(str(np.mean(sample_r2_avg)) + '\t')

            # Writing overall sample average to file
            report_file.write(str(np.mean(sample_avg)) + '\n')

    return()


# Performing Summary Report Analysis (walks through directories
# again to perform the analysis
def parse_for_summary_results(list_of_directories, fastqc_summary_report):

    for i in range(len(list_of_directories)):

        directory_pathway = list_of_directories[i]

        directory_list = []
        found_directories = os.listdir(list_of_directories[i])
        for found_directory in found_directories:
            if os.path.isdir(found_directory):
                if found_directory.startswith('Log'):
                    pass
                else:
                    directory_list.append(found_directory)
                    # print(found_directory)

        # Loop searches through the directories and retrieves all the data from the various sub-directories
        # Sample specific lane runs from FASTQC
        for directory in directory_list:

            # Changes the directory of the program to go inside the sample directory
            new_directory = directory_pathway + "/" + directory + "/"
            os.chdir(new_directory)

            # In directory of ALL lane runs for a sample
            fastqc_directories = os.listdir(new_directory)

            for fastqc_directory in fastqc_directories:

                if fastqc_directory.endswith(".zip"):
                    pass

                elif fastqc_directory.endswith("_fastqc"):
                    final_directory_step = new_directory + fastqc_directory + "/"

                    # Change to the final directory step (where data is at)
                    os.chdir(final_directory_step)

                    if 'R1' in fastqc_directory:
                        # Parsing through each samples fastqc lane results
                        fastqc_summary_report = fastqc_summary(fastqc_directory, fastqc_summary_report)

                    elif 'R2' in fastqc_directory:
                        # Parsing through each samples fastqc lane results
                        fastqc_summary_report = fastqc_summary(fastqc_directory, fastqc_summary_report)

                    else:
                        print("Error with parsing R1 and R2")

                # Tells the program to skip other not important files
                else:
                    pass

    fastqc_summary_report.close()

    return()


# prints the final summary reports data to file (simple line
# loop parsing)
def fastqc_summary(directory_name, fastqc_summary_report):
    # opening fastqc results
    fastqc_results = open('summary.txt', 'r')
    fastqc_summary_report.write(str(directory_name)+"\t")
    for line in fastqc_results:
        data = line.split("\t")
        fastqc_summary_report.write(str(data[0])+"\t")

    pathway_of_file = os.getcwd()
    fastqc_summary_report.write(str(pathway_of_file)+"\t")

    fastqc_summary_report.write("\n")

    return fastqc_summary_report


def main():

    print("There is nothing either good or bad, but thinking makes it so.")
    print("Hamlet - Hamlet (William Shakespeare)")
    print("")

    print("FastQC Parsing Program is Starting to Run")
    print("")

    # Get current working directory
    home_directory = os.getcwd()

    # Get the list of directories (provided by user)
    list_of_directories = read_parameter_file()

    # See if directory exists otherwise make it
    verdict = os.path.exists('Log_File_Directory')
    if str(verdict) == 'False':
        os.makedirs('Log_File_Directory')
    else:
        print("Log Directory already exists")
        print("")

    ########################################################################################
    # Mean Report
    # Create Analysis files for Mean Reports of FastQC Analysis
    print("Starting run of mean report of FastQC data")
    print("")
    # Change to home working directory
    os.chdir(home_directory)
    mean_report_file = open("FASTQC_Mean_Report.txt", "w")
    mean_log_file = open("Log_File_Directory/FASTQC_Mean_Report_LOG.txt", "w")
    position = 1
    parse_and_print_results(mean_report_file, mean_log_file, list_of_directories, position)
    # Closing final files
    mean_report_file.close()
    mean_log_file.close()

    ########################################################################################
    # Median Report
    # Create Analysis files for Median Reports of FastQC Analysis
    print("Starting run of median report of FastQC data")
    print("")
    # Change back to home working directory
    os.chdir(home_directory)
    median_report_file = open("FASTQC_Median_Report.txt", "w")
    median_log_file = open("Log_File_Directory/FASTQC_Median_Report_LOG.txt", "w")
    position = 2
    parse_and_print_results(median_report_file, median_log_file, list_of_directories, position)
    # Closing final files
    median_report_file.close()
    median_log_file.close()

    ########################################################################################
    # Lower Quartile Report
    # Create Analysis files for Lower Quartile Reports of analysis
    print("Starting run of Lower Quartile report of FastQC data")
    print("")
    # Change back to home working directory
    os.chdir(home_directory)
    lower_quartile_report_file = open("FASTQC_Lower_Quartile_Report.txt", "w")
    lower_quartile_log_file = open("Log_File_Directory/FASTQC_Lower_Quartile_LOG.txt", "w")
    position = 3
    parse_and_print_results(lower_quartile_report_file, lower_quartile_log_file, list_of_directories, position)
    # Closing final files
    lower_quartile_report_file.close()
    lower_quartile_log_file.close()

    ########################################################################################
    # Upper Quartile Report
    # Create Analysis files for Upper Quartile Reports of analysis
    print("Starting run of Upper Quartile report of FastQC data")
    print("")
    # Change back to home working directory
    os.chdir(home_directory)
    upper_quartile_report_file = open("FASTQC_Upper_Quartile_Report.txt", "w")
    upper_quartile_log_file = open("Log_File_Directory/FASTQC_Upper_Quartile_LOG.txt", "w")
    position = 4
    parse_and_print_results(upper_quartile_report_file, upper_quartile_log_file, list_of_directories, position)
    # Closing final files
    upper_quartile_report_file.close()
    upper_quartile_log_file.close()

    ########################################################################################
    # 10th Percentile Report
    # Create Analysis files for Upper Quartile Reports of analysis
    print("Starting run of 10th Percentile report of FastQC data")
    print("")
    # Change back to home working directory
    os.chdir(home_directory)
    tenth_percentile_report_file = open("FASTQC_Tenth_Percentile_Report.txt", "w")
    tenth_percentile_log_file = open("Log_File_Directory/FASTQC_Tenth_Percentile_LOG.txt", "w")
    position = 5
    parse_and_print_results(tenth_percentile_report_file, tenth_percentile_log_file, list_of_directories, position)
    # Closing final files
    tenth_percentile_report_file.close()
    tenth_percentile_log_file.close()

    ########################################################################################
    # 90th Percentile Report
    # Create Analysis files for Upper Quartile Reports of analysis
    print("Starting run of 90th Percentile report of FastQC data")
    print("")
    # Change back to home working directory
    os.chdir(home_directory)
    ninetieth_percentile_report_file = open("FASTQC_Ninetieth_Percentile_Report.txt", "w")
    ninetieth_percentile_log_file = open("Log_File_Directory/FASTQC_Ninetieth_Percentile_LOG.txt", "w")
    position = 6
    parse_and_print_results(ninetieth_percentile_report_file, ninetieth_percentile_log_file,
                            list_of_directories, position)
    # Closing final files
    ninetieth_percentile_report_file.close()
    ninetieth_percentile_log_file.close()

    ########################################################################################
    # Create Summary Report of All Data (Pass/Warning/Fail)
    print("Starting run of summary.txt files from FastQC")
    print("")
    # Change back to working directory
    os.chdir(home_directory)
    fastqc_summary_report = open("FASTQC_Summary_Report.txt", 'w')
    fastqc_summary_report.write("Title\tBasic_Stats\tPer_base_sequence_quality\t"
                                "Per_sequence_quality_scores\tPer_base_sequence_content\t"
                                "Per_base_GC_content\tPer_sequence_GC_content\t"
                                "Per_base_N_content\tSequence_Length_Distribution\t"
                                "Sequence_Duplication Levels\tOverrepresented_sequences\t"
                                "Kmer_Content\tPathway_of_File\n")
    parse_for_summary_results(list_of_directories, fastqc_summary_report)
    fastqc_summary_report.close()

    print("Done Running of Parsing Program")


main()
