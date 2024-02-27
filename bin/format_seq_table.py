#!/usr/bin/env python

from argparse import ArgumentParser
import csv
import sys

csv.field_size_limit(sys.maxsize)

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
        ctx_seq_dict = {}
        outf.write("name;sequence;position\n")
        reader = csv.DictReader(csvfile, delimiter=';')
        for row in reader:
            ctx_seq_id = row['context_sequence_id']
            
            if ctx_seq_id in ctx_seq_dict:
                print(row['context_sequence_wt1'] == ctx_seq_dict[ctx_seq_id]["wt1"])
            else:
                ctx_seq_dict[ctx_seq_id] = {
                    "ft": row['context_sequence'],
                    "ft_bp": row['context_sequence_bp'],
                    "wt1": row['context_sequence_wt1'],
                    "wt1_bp": row['context_sequence_wt1_bp'],
                    "wt2": row['context_sequence_wt2'],
                    "wt2_bp": row['context_sequence_wt2_bp']
                }
        for ctx_seq_id in ctx_seq_dict:
            values = ctx_seq_dict[ctx_seq_id]
            for seq_type in ("ft", "wt1", "wt2"):
                outf.write("{0}_{1}_{2};{3};{1}\n".format(
                    ctx_seq_id, 
                    values[seq_type + "_bp"], 
                    seq_type,
                    values[seq_type]
                ))


if __name__ == "__main__":
    main()