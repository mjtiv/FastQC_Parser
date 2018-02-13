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


def get_stats_of_data(fastqc_results):

    # Get average from fastqc
    quality_scores = []

    for data in fastqc_results:
        data = data.split("\t")
        quality_scores.append(float(data[1]))

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


def main():

    # Create new files for reports of analysis
    final_report_file = open("FASTQC_Sample_Final_Report.txt", "w")
    sample_log_file = open("FASTQC_Sample_Log_File.txt", "w")

    final_report_file.write("Sample\tR1_1st_Block_Avg\tR1_Middle_Blocks_Avg\tR1_Final_Block_Avg\tR1_Avg\t"
                            "R2_1st_Block_Avg\tR2_Middle_Blocks_Avg\tR2_Final_Block_Avg\tR2_Block_Avg\t"
                            "Overall_Sample_Avg\n")

    sample_log_file.write("Sample\tR1_1st_Block_Avg\tR1_Middle_Blocks_Avg\tR1_Final_Block_Avg\tR1_Avg\t"
                          "R2_1st_Block_Avg\tR2_Middle_Blocks_Avg\tR2_Final_Block_Avg\tR2_Block_Avg\t"
                          "Overall_Sample_Avg\n")

    # Get the list of directories (provided by user)
    list_of_directories = read_parameter_file()

    for i in range(len(list_of_directories)):

        directory_pathway = list_of_directories[i]

        # How parse the directory found on stackoverflow link (only gets directories)
        # Code from comment 1. http://stackoverflow.com/questions/7781545/how-to-get-all-folder
        # -only-in-a-given-path-in-python
        directory_list = [d for d in os.listdir(list_of_directories[i]) if os.path.isdir(
            os.path.join(list_of_directories[i],
                         d))]

        # Loop searches through the directories and retrieves all the data from the various sub-directories
        # Sample specific lane runs from FASTQC
        for directory in directory_list:

            # Parsing the directory name to get the sample name
            sample_name = directory

            final_report_file.write(sample_name+"\t")
            sample_log_file.write(sample_name+"\t")

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
                        stats_on_lane = get_stats_of_data(fastqc_results)
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
                        stats_on_lane = get_stats_of_data(fastqc_results)
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
            sample_log_file.write(sample_r1_avg_first_block_str_list+"\t")
            sample_log_file.write(sample_r1_avg_middle_blocks_str_list + "\t")
            sample_log_file.write(sample_r1_avg_final_block_str_list + "\t")
            sample_log_file.write(sample_r1_avg_str_list + "\t")

            # Writing R2 Lane Info to Log File
            sample_log_file.write(sample_r2_avg_first_block_str_list + "\t")
            sample_log_file.write(sample_r2_avg_middle_blocks_str_list + "\t")
            sample_log_file.write(sample_r2_avg_final_block_str_list + "\t")
            sample_log_file.write(sample_r2_avg_str_list + "\t")

            # Writing Overall Average to file
            sample_log_file.write(sample_avg_str_list + "\n")

            # Writing sample results to a file
            final_report_file.write(str(np.mean(sample_r1_avg_first_block))+'\t')
            final_report_file.write(str(np.mean(sample_r1_avg_middle_blocks)) + '\t')
            final_report_file.write(str(np.mean(sample_r1_avg_final_block)) + '\t')
            final_report_file.write(str(np.mean(sample_r1_avg)) + '\t')

            # Writing sample results to a file
            final_report_file.write(str(np.mean(sample_r2_avg_first_block)) + '\t')
            final_report_file.write(str(np.mean(sample_r2_avg_middle_blocks)) + '\t')
            final_report_file.write(str(np.mean(sample_r2_avg_final_block)) + '\t')
            final_report_file.write(str(np.mean(sample_r2_avg)) + '\t')

            # Writing overall sample average to file
            final_report_file.write(str(np.mean(sample_avg)) + '\n')

    # Closing final files
    final_report_file.close()
    sample_log_file.close()


main()
