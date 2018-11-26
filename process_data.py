# Example comment
import os  # <-- self explanatory
import pandas  # <-- self explanatory
import re  # <-- self explanatory
import csv  # <-- self explanatory
from sklearn.linear_model import LinearRegression
import numpy as np
import matplotlib.pyplot as plt  # Imports matplotlib and renames it plt 

INPUT_DIR = "data/YeastGenes"  # Sets the path of the Yeast Genome data
RNASEQ_DATA_PATH = "data/RNAseq_data.gff3"  # Sets the path of the RNA dataset
OUTPUT_DIR = "output"  # Sets the path for the output file with the processed data


def load_rnaseq_data(seqdata_path):  # Loads RNA dataset for processing
    rnaseq_data = dict()  # Sets the format of the dataset as dictionary
    parsed_file = pandas.read_csv(seqdata_path, comment="#", header=None, delimiter="\t")  # Parses the data 

    feature_format = r"ID=(.*);Name=(.*);log2_transcription_level=(.*);Note=(.*)"  # Sets the format of the relevant data from the set (Name and Transcription Level)
    for row_num, row in parsed_file.iterrows(): # Iterates through the lines of the dataset and pulls the 
        feature = row.values[-1]  # Finds last chunck of data (Tab delimited)
        regex_match = re.match(feature_format, feature)  # Finds the relevant data in the dataset
        orf_name_w_5utr = regex_match.group(1)  # Pulls the name of the gene
        orf_name = orf_name_w_5utr.replace("_5UTR", "")  # Removes the extra characters in the string
        transcription_level = regex_match.group(3)  # Pulls transcription level from the data
        rnaseq_data[orf_name] = transcription_level  # Associates transcription level with the gene name
    return rnaseq_data  # Outputs the name and transcription level


def calculate_gc_percent(sequence):
    count_gc = len(re.findall("[G|C]", sequence))  # Finds number of G's and C's in the sequence
    count_all = len(sequence)  # Finds overall length of sequence
    gc_percent = count_gc / count_all  # Expresses the GC content as # of G and C over the total length of the sequence
    return gc_percent


def update_file(input_dir, output_dir, rnaseq_info):  # Defines function to create updated files
    seq_data = list()  # Formats the data as a list
    dir_list = os.listdir(input_dir)  # Lists files in directory
    for file in dir_list:  # Iterates through the files in the input directory of Yeast genes
        path = os.path.join(input_dir, file)  # Creates path for each file 
        with open(path) as f:  # Opens each file in the set Yeast files
            text = f.read().splitlines()  # Grabs the text of each file, and splits along lines
        orf_name = file.replace(".txt", "")  # Removes the .txt from each file name

        if orf_name in rnaseq_info:  # Checks to see if the Orf name is in the RNA dataset
            trans_level = rnaseq_info[orf_name]  # If it is, the transcription level is pulled from the list
        else:  # If not...
            trans_level = "NA"  # Set the transcription level to "NA"
        text[0] += " " + orf_name + " " + trans_level + "\n"  # Rewrites the text of the yeast gene file to include the Orf name and transcription level

        output_path = os.path.join(output_dir, file)  # Defines the file path to output the updated yeast info to
        with open(output_path, 'w') as of:  # Writes the data to a new file
            of.writelines(text)  # Inserts the text to the file

        # GC Content calculations
        if len(text) == 2:  # If the data has 2 lines, continue. Otherwise, ignore it. (Basically just for YEL076C
            gc_percent = calculate_gc_percent(text[1])
        else:  # If there isnt data...
            gc_percent = 0  # %GC = 0

        seq_data.append((orf_name, trans_level, gc_percent))  # Adds name, transcription level and GC content
    return seq_data  # Returns data to main code for plotting later


def plot_transcription_gc(seq_data):  # Beginning of plotting
    trans_level_list = [float(x[1]) for x in seq_data if x[1] != "NA"]  # Finds transcription level and sets it to X data
    gc_pct_list = [float(x[2]) for x in seq_data if x[1] != "NA"]  # Finds GC content and sets it to Y data

    plot = plt.scatter(trans_level_list, gc_pct_list, s=3)  # Plots Y vs X
    plt.xlabel("Transcription Level")  # Sets X label
    plt.ylabel("GC%")  # Sets Y label
    plt.title("Transcription Level vs GC%")  # Sets Title

    trans_level_list_2d = np.array(trans_level_list).reshape(-1, 1)
    regression = LinearRegression()
    regression.fit(trans_level_list_2d, gc_pct_list)
    plt.plot(trans_level_list, regression.predict(trans_level_list_2d), color="green")

    plt.show()  # Displays plot


if __name__ == '__main__':  # Tells python to run this if this is the file is being run
    rnaseq_info = load_rnaseq_data(RNASEQ_DATA_PATH)  # Opens RNA data
    os.makedirs(OUTPUT_DIR, exist_ok=True)  # Creates output directory
    seq_data = update_file(INPUT_DIR, OUTPUT_DIR, rnaseq_info)  # Calls funciton to update files

    with open("sequence_info.csv", "w") as csv_file:  # Creates csv with data
        writer = csv.writer(csv_file)  # Creates CSV file
        writer.writerows(seq_data)  # Writes data to file

    plot_transcription_gc(seq_data)  # calls plot creation function


