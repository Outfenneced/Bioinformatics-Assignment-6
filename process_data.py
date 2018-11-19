import os
import pandas
import re
import csv
import matplotlib.pyplot as plt

INPUT_DIR = "data/YeastGenes"
RNASEQ_DATA_PATH = "data/RNAseq_data.gff3"
OUTPUT_DIR = "output"


def load_rnaseq_data(seqdata_path):
    rnaseq_data = dict()
    parsed_file = pandas.read_csv(seqdata_path, comment="#", header=None, delimiter="\t")

    feature_format = r"ID=(.*);Name=(.*);log2_transcription_level=(.*);Note=(.*)"
    for row_num, row in parsed_file.iterrows():
        feature = row.values[-1]
        regex_match = re.match(feature_format, feature)
        orf_name_w_5utr = regex_match.group(1)
        orf_name = orf_name_w_5utr.replace("_5UTR", "")
        transcription_level = regex_match.group(3)
        rnaseq_data[orf_name] = transcription_level
    return rnaseq_data


def update_file(input_dir, output_dir, rnaseq_info):
    seq_data = list()
    dir_list = os.listdir(input_dir)
    for file in dir_list:
        path = os.path.join(input_dir, file)
        with open(path) as f:
            text = f.read().splitlines()
        orf_name = file.replace(".txt", "")
        if orf_name in rnaseq_info:
            trans_level = rnaseq_info[orf_name]
        else:
            trans_level = "NA"
        text[0] += " " + orf_name + " " + trans_level + "\n"

        if len(text) == 2:
            count_gc = len(re.findall("[G|C]", text[1]))
            count_all = len(text[1])
            gc_percent = count_gc/count_all
        else:
            gc_percent = 0

        output_path = os.path.join(output_dir, file)
        with open(output_path, 'w') as of:
            of.writelines(text)

        seq_data.append((orf_name, trans_level, gc_percent))
    return seq_data


def plot_transcription_gc(seq_data):
    trans_level_list = [x[1] for x in seq_data if x[1] != "NA"]
    gc_pct_list = [x[2] for x in seq_data if x[1] != "NA"]

    plot = plt.scatter(trans_level_list, gc_pct_list, s=3)
    plt.xlabel("Transcription Level")
    plt.ylabel("GC%")
    plt.title("Transcription Level vs GC%")
    plt.show()


if __name__ == '__main__':
    rnaseq_info = load_rnaseq_data(RNASEQ_DATA_PATH)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    seq_data = update_file(INPUT_DIR, OUTPUT_DIR, rnaseq_info)

    with open("sequence_info.csv", "w") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerows(seq_data)

        plot_transcription_gc(seq_data)


