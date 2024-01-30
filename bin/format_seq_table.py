#!/usr/bin/env python

from argparse import ArgumentParser
import csv

def main():
    parser = ArgumentParser(description="Reformats context seqs table to fit easyquant format")
    parser.add_argument(
        "-i", 
        "--input_table", 
        dest="input_table",
        required=True,
        help="Specify input table (context_seqs.csv)",
    )
    parser.add_argument(
        "-o",
        "--output_table",
        dest="output_table",
        required=True,
        help="Specify output table",
    )

    args = parser.parse_args()

    # 5: FTID
    # 6: context_sequence_id
    # 18: context_sequence
    # 19: context_sequence_bp
    # 22: context_sequence_wt1
    # 23: context_sequence_wt2
    # 24: context_sequence_wt1_bp
    # 25: context_sequence_wt2_bp
    with open(args.input_table) as csvfile, open(args.output_table, "w") as outf:
        outf.write("name;sequence;position\n")
        reader = csv.DictReader(csvfile, delimiter=';')
        for row in reader:
            ftid = row['FTID']
            ctx_seq_id = row['context_sequence_id']
            ctx_seq_ft = row['context_sequence']
            ctx_seq_ft_bp = row['context_sequence_bp']
            ctx_seq_wt1 = row['context_sequence_wt1']
            ctx_seq_wt1_bp = row['context_sequence_wt1_bp']
            ctx_seq_wt2 = row['context_sequence_wt2']
            ctx_seq_wt2_bp = row['context_sequence_wt2_bp']
            
            outf.write("{}_{}_{}_ft;{};{}\n".format(ftid, ctx_seq_id, ctx_seq_ft_bp, ctx_seq_ft, ctx_seq_ft_bp))
            outf.write("{}_{}_{}_wt1;{};{}\n".format(ftid, ctx_seq_id, ctx_seq_wt1_bp, ctx_seq_wt1, ctx_seq_wt1_bp))
            outf.write("{}_{}_{}_wt2;{};{}\n".format(ftid, ctx_seq_id, ctx_seq_wt2_bp, ctx_seq_wt2, ctx_seq_wt2_bp))


if __name__ == "__main__":
    main()