#!/usr/bin/env python3

"""
@author: BNT (URLA), TRON (PASO)
@version: 20190118
"""
import os
import os.path
import time

import easy_fuse
import pkg_resources
from logzero import logger
import pandas as pd
from easy_fuse.misc import queueing
from easy_fuse.misc.config import EasyFuseConfiguration


class FusionSummary(object):
    """Collect stats of the run and write them to file"""

    def __init__(
            self,
            input_fusions,
            input_fusion_context_seqs,
            input_requant_cpm,
            input_requant_counts,
            input_reads_stats,
            output_folder,
            config: EasyFuseConfiguration):

        assert os.path.exists(input_fusions) and os.path.isfile(input_fusions), \
            "Fusion file not found: {}".format(input_fusions)
        self.input_fusions = input_fusions
        assert os.path.exists(input_fusion_context_seqs) and os.path.isfile(input_fusion_context_seqs), \
            "Fusion context seq file not found: {}".format(input_fusion_context_seqs)
        self.input_fusion_context_seqs = input_fusion_context_seqs
        assert os.path.exists(input_requant_cpm) and os.path.isfile(input_requant_cpm), \
            "Requant CPM file not found: {}".format(input_requant_cpm)
        self.input_requant_cpm = input_requant_cpm
        assert os.path.exists(input_requant_counts) and os.path.isfile(input_requant_counts), \
            "Requant counts file not found: {}".format(input_requant_counts)
        self.input_requant_counts = input_requant_counts
        assert os.path.exists(input_reads_stats) and os.path.isfile(input_reads_stats), \
            "Read stats file not found: {}".format(input_reads_stats)
        self.input_reads_stats = input_reads_stats
        assert os.path.exists(output_folder) and os.path.isdir(output_folder), \
            "Output folder not found: {}".format(output_folder)
        self.output_folder = output_folder
        self.cfg = config

        # loads data
        with open(input_reads_stats, "r") as rfile:
            self.input_read_count = int(rfile.readline())
        fusion_tools = self.cfg["general"]["fusiontools"].split(",")
        requant_mode = self.cfg["general"]["requant_mode"]
        self.data = self.load_data(fusion_tools=fusion_tools, requant_mode=requant_mode)

        pd.options.mode.chained_assignment = None
        self.blacklist = []

    def load_data(self, fusion_tools, requant_mode):
        """Join the three data tables context_seq, detected_fusions and requantification"""
        # define path' to the context seq, detected fusion and re-quantification files

        context_data = pd.read_csv(self.input_fusion_context_seqs, sep=";")
        # append key for later join with requant data
        context_data["ftid_plus"] = (
                context_data["FTID"] + "_" + context_data["context_sequence_id"]
        )

        logger.info(
            "Appending normalized fusion counts from {} to the data table.".format(
                fusion_tools
            )
        )

        # read the original append CPM values from individual fusion prediction tools to context data
        detect_data = pd.read_csv(self.input_fusions, sep=";")
        context_data = self.append_tool_cnts_to_context_file(
            context_data, detect_data, fusion_tools
        )

        logger.info(
            "Appending (normalized) requantification counts from filtered and original mapping to the data table."
        )
        # perform subsequent joins on ftid_plus
        context_data.set_index("ftid_plus", inplace=True)

        # read and append normalized (cpm) and raw (counts) requantification data to context data
        # for mode in requant_mode:
        requant_cpm_data = pd.read_csv(self.input_requant_cpm, sep=";")
        context_data = context_data.join(
            requant_cpm_data.set_index("ftid_plus"), how="left"
        )
        requant_cnt_data = pd.read_csv(self.input_requant_counts, sep=";")
        context_data = context_data.join(
            requant_cnt_data.set_index("ftid_plus"),
            lsuffix="_{}".format(requant_mode),
            rsuffix="_cnt_{}".format(requant_mode),
            how="left",
        )

        return context_data.fillna(0)

    def run(self, model_predictions):
        """Execute individual methods"""

        # 0:type; 1:boundary; 2:frame; 3:pepseq; 4:counts; 5:blacklist; 6:prediction;
        # 7:exonError; 8:unfiltered; 9:allFilter; 10:allButPredFilter

        start_time = time.time()
        blacklist_file = self.cfg.get("general", "blacklist", fallback=None)
        self.load_blacklist(blacklist_file)

        # check whether output already exists
        cols_to_aggregate_on = ["BPID", "context_sequence_id", "FTID"]

        # in comparison to the previous implementation, pre-existing results will always be overwritten as it turned out that the
        # previous implementation is to error prone...
        joined_table_1b = self.data.groupby(
            by=cols_to_aggregate_on, sort=False, as_index=False
        )
        joined_table_1b.agg(self.custom_data_aggregation).to_csv(
            "{}_fusRank_1.csv".format(self.output_folder), sep=";", index=False
        )

        joined_table_1 = pd.read_csv("{}_fusRank_1.csv".format(self.output_folder), sep=";")

        for table in ["1", "2"]:
            summary_file = "{}_fusRank_{}.csv".format(self.output_folder, table)
            if model_predictions and os.path.exists(summary_file) and os.path.isfile(summary_file):
                model_path = pkg_resources.resource_filename(
                    easy_fuse.__name__,
                    "resources/model/Fusion_modeling_FFPE_train_v35.random_forest.model_full_data.EF_full.rds",
                )
                model_threshold = self.cfg["general"]["model_pred_threshold"]

                # append prediction scores based on pre-calculated model
                cmd_model = "{0} --fusion_summary {1} --model_file {2} --prediction_threshold {3} --output {4}".format(
                    pkg_resources.resource_filename(
                        easy_fuse.__name__, "resources/R/R_model_prediction.R"
                    ),
                    summary_file,
                    model_path,
                    model_threshold,
                    "{}.pred.csv".format(summary_file[:-4]),
                )

                queueing.submit(
                    "", cmd_model.split(" "), "", "", "", "", "", "", "", "", "", "none"
                )
                # re-read the table with prediction for filter counting
                # urla - note: there is probably a more elegant way using getattr/setattr but I'm not at all familiar with its pros/cons
                if table == "1":
                    joined_table_1 = pd.read_csv(
                        "{}_fusRank_1.pred.csv".format(self.output_folder), sep=";"
                    )

        joined_table_1.set_index(keys=cols_to_aggregate_on, drop=False, inplace=True)
        time_taken = time.time() - start_time
        logger.info(
            "Done. Processing time: {0:.2f} s".format(
                time_taken
            )
        )

        return self.count_records(joined_table_1, False, "F1N", model_predictions=model_predictions)


    def count_records(self, pd_data_frame, is_merged, name, model_predictions):
        """stub"""
        logger.info("Processing: {}".format(name))
        #        blacklist = ["HLA", "IG"]
        # list of sets to store (un)filtered fusion gene names
        # 0:unfiltered; 1:type; 2:boundary; 3:frame; 4:pepseq; 5:counts;
        # 6:blacklist; 7:prediction; 8:allFilter; 9:allButPredFilter
        # this list is initialized with empty sets which are either replaced with data or stay empty
        # this guarantees that filters are always at the same position in the list, even if they are empty
        fusion_gene_set_list = [set() for i in range(10)]
        # the dataframe may be empty if (a) no fusions were predicted or (b) no fusion remained after replicate merging
        if pd_data_frame.empty:
            return map(len, fusion_gene_set_list), set()
        df_rep = pd_data_frame
        # try to correct annotation issues
        ftid_dict = {}
        for i in df_rep.index:
            ftid = df_rep.loc[i, "FTID"]
            if ftid not in ftid_dict:
                ftid_dict[ftid] = [df_rep.loc[i, "Fusion_Gene"]]
            else:
                ftid_dict[ftid].append(df_rep.loc[i, "Fusion_Gene"])

        for ftid in ftid_dict:
            # check if multiple fusion_gene strings are associated with the ftid
            if len(ftid_dict[ftid]) > 1:
                # check if there is one of the multiple fusion_gene records where both partners correspond to the ftid
                for fusion_gene in ftid_dict[ftid]:
                    try:
                        left_fusion, right_fusion = fusion_gene.split("_")
                        if left_fusion in ftid and right_fusion in ftid:
                            ftid_dict[ftid] = fusion_gene
                            break
                    except ValueError:
                        logger.warning(
                            "Fusion gene is in an unexpected format (If you have gene names containing underscores, please contact the authors!)"
                        )
                # if none of the fusion_gene records fit to the ftid, take the shortest and print a warning
                if type(ftid_dict[ftid]) == list:
                    logger.warning(
                        'Unresolvable annotation bias detected between ftid "{0}" and fusion_gene records "{1}"'.format(
                            ftid, ftid_dict[ftid]
                        )
                    )
                    ftid_dict[ftid] = min(ftid_dict[ftid], key=len)
                    logger.warning(
                        "Chosing {} as fusion_gene record which might be wrong.".format(
                            ftid_dict[ftid]
                        )
                    )
            else:
                ftid_dict[ftid] = ftid_dict[ftid][0]

        # create sets of fusion gene names according to different fiters
        fusion_gene_set_list[1] = set(
            [
                ftid_dict[x]
                for x in df_rep[df_rep["type"].astype(str) != "cis_near"]["FTID"]
            ]
        )  # exclude cis near events as they are likely read through
        fusion_gene_set_list[2] = set(
            [
                ftid_dict[x]
                for x in df_rep[df_rep["exon_boundary"].astype(str) == "both"]["FTID"]
            ]
        )  # allow only fusions where the breakpoints are on exon boundaries in BOTH fusion partners
        fusion_gene_set_list[3] = set(
            [
                ftid_dict[x]
                for x in df_rep[df_rep["frame"].astype(str) != "no_frame"]["FTID"]
            ]
        )  # exclude no frame fusions
        # urla - note: I'm not 100% why str.isalpha is working everywhere I tested, but not on our dev servers... This way, it works, although I thought the pandas str method would this already...
        fusion_gene_set_list[4] = set(
            [
                ftid_dict[x]
                for x in df_rep[
                    df_rep["neo_peptide_sequence"].astype(str).str.isalpha()
                ]["FTID"]
            ]
        )  # the neo peptide sequence must be an alphabetic sequence of at least length 1
        if is_merged:
            junc_cnt_col = "ft_junc_cnt_org_1"
            if not junc_cnt_col in list(df_rep):
                junc_cnt_col = "ft_junc_cnt_best_1"
                if not junc_cnt_col in list(df_rep):
                    junc_cnt_col = "ft_junc_cnt_fltr_1"
            fusion_gene_set_list[5] = set(
                [
                    ftid_dict[x]
                    for x in df_rep[
                        df_rep[junc_cnt_col]
                        + df_rep[junc_cnt_col.replace("junc", "span")]
                        + df_rep[junc_cnt_col.replace("_1", "_2")]
                        + df_rep[
                            junc_cnt_col.replace("junc", "span").replace("_1", "_2")
                        ]
                        > 0
                    ]["FTID"]
                ]
            )  # at least 1 junction or spanning counts on the fusion transcript in either replicate
        else:
            junc_cnt_col = "ft_junc_cnt_org"
            if not junc_cnt_col in list(df_rep):
                junc_cnt_col = "ft_junc_cnt_best"
                if not junc_cnt_col in list(df_rep):
                    junc_cnt_col = "ft_junc_cnt_fltr"
            fusion_gene_set_list[5] = set(
                [
                    ftid_dict[x]
                    for x in df_rep[
                        df_rep[junc_cnt_col]
                        + df_rep[junc_cnt_col.replace("junc", "span")]
                        > 0
                    ]["FTID"]
                ]
            )  # at least 1 junction or spanning counts on the fusion transcript
        # both fusion partners are not allowed to start with anything mentioned in the "blacklisted" list (single gene/gene family exclusion)
        # and the fusion gene is not allowed to be identical to a string in the blacklist (gene pair exclusion)
        fusion_gene_set_list[6] = set(
            [
                ftid_dict[z]
                for z in df_rep[
                    [
                        all(
                            [
                                not (
                                    y.startswith(tuple(self.blacklist))
                                    or "_".join(x) in self.blacklist
                                )
                                for y in x
                            ]
                        )
                        for x in df_rep["Fusion_Gene"].str.split("_")
                    ]
                ]["FTID"]
            ]
        )
        if model_predictions:
            if is_merged:
                fusion_gene_set_list[7] = set(
                    [
                        ftid_dict[x]
                        for x in df_rep[
                            (df_rep["prediction_class_1"].astype(str) == "positive")
                            | (df_rep["prediction_class_2"].astype(str) == "positive")
                        ]["FTID"]
                    ]
                )
                # the random forrest model must have classified this at least once as "posititve" in either replicate
            else:
                fusion_gene_set_list[7] = set(
                    [
                        ftid_dict[x]
                        for x in df_rep[
                            df_rep["prediction_class"].astype(str) == "positive"
                        ]["FTID"]
                    ]
                )
                # the random forrest model must have classified this at least once as "posititve"
        fusion_gene_filter_wo_pred = set(ftid_dict.values())
        # intersect with all but predictions
        for i in range(1, 7):
            fusion_gene_filter_wo_pred = fusion_gene_filter_wo_pred.intersection(
                fusion_gene_set_list[i]
            )
        fusion_gene_filter_all = fusion_gene_filter_wo_pred.intersection(
            fusion_gene_set_list[7]
        )

        fusion_gene_set_list[0] = set(ftid_dict.values())
        fusion_gene_set_list[8] = fusion_gene_filter_all
        fusion_gene_set_list[9] = fusion_gene_filter_wo_pred

        if model_predictions:
            return map(len, fusion_gene_set_list), fusion_gene_filter_all
        return map(len, fusion_gene_set_list), fusion_gene_filter_wo_pred

    def append_tool_cnts_to_context_file(
        self, context_data, detected_fusions, fusion_tool_list
    ):
        """Add data to context seq table"""
        # Init new columns
        for tool in fusion_tool_list:
            context_data["{}_detected".format(tool.lower())] = 0
            context_data["{}_junc".format(tool.lower())] = 0
            context_data["{}_span".format(tool.lower())] = 0
        context_data["tool_count"] = 0
        context_data["tool_frac"] = 0.0

        # iterate over context seq data
        for i in context_data.index:
            # get tool prediction for each BPID
            tmp_slice = detected_fusions.loc[
                detected_fusions["BPID"] == context_data["BPID"][i]
            ]
            # check for each tool, whether a prediction was made by it and update the respective columns
            for tool in fusion_tool_list:
                # pre-set values to 0 or NA indicating no prediction with this tool
                # line wise initialization shouldn't be required as it was initialized with this values in the beginning
                for j in tmp_slice.index:
                    if tool == tmp_slice.loc[j, "Tool"]:
                        context_data.loc[i, "{}_detected".format(tool.lower())] = 1
                        context_data.loc[
                            i, "{}_junc".format(tool.lower())
                        ] = self.normalize_counts_cpm(
                            tmp_slice.loc[j, "Junction_Reads"]
                        )
                        context_data.loc[
                            i, "{}_span".format(tool.lower())
                        ] = self.normalize_counts_cpm(
                            tmp_slice.loc[j, "Spanning_Reads"]
                        )
                        context_data.loc[i, "tool_count"] += 1
            context_data.loc[i, "tool_frac"] = float(
                context_data.loc[i, "tool_count"]
            ) / float(len(fusion_tool_list))
        return context_data

    def load_blacklist(self, blacklist_file):
        """Read blacklist data from file"""
        if blacklist_file is not None:
            with open(blacklist_file, "r") as bl_file:
                for line in bl_file:
                    if line.startswith("#"):
                        continue
                    self.blacklist.append(line.strip())

    def normalize_counts_cpm(self, count):
        """Returns the CPM of a read count based on the number of original input reads"""
        if self.input_read_count == 0:
            return count
        return count / float(self.input_read_count) * 1000000.0

    @staticmethod
    def custom_data_aggregation(pd_series):
        """Removes duplicates from a pandas series fields and creates a joined string from the resulting list"""
        # urla - note: this is an adaptation from https://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-whilst-preserving-order
        #              in a way that it does not return a list but a string to print from the list
        seen = set()
        seen_add = seen.add
        return ";".join(
            str(e) for e in [x for x in pd_series if not (x in seen or seen_add(x))]
        )


def add_summarize_data_args(parser):

    parser.add_argument(
        "--input-fusions",
        dest="input_fusions",
        required=True,
        help="Path to input file with fusions",
    )
    parser.add_argument(
        "--input-fusion-context-seqs",
        dest="input_fusion_context_seqs",
        required=True,
        help="Path to input file with fusion context sequences",
    )
    parser.add_argument(
        "--input-requant-cpm",
        dest="input_requant_cpm",
        required=True,
        help="Path to input file with requant CPM",
    )
    parser.add_argument(
        "--input-requant-counts",
        dest="input_requant_counts",
        required=True,
        help="Path to input file with requant counts",
    )
    parser.add_argument(
        "--input-reads-stats",
        dest="input_reads_stats",
        required=True,
        help="Path to input file with reads stats",
    )
    parser.add_argument(
        "-c",
        "--config-file",
        dest="config_file",
        required=True,
        help="Specify alternative config file to use for your analysis",
    )
    parser.add_argument(
        "--model_predictions",
        default=False,
        action="store_true",
        help="Score predictions based on pre-train model",
    )
    parser.add_argument(
        "-o",
        "--output-folder",
        dest="output_folder",
        help="Specify the output folder for the summary files."
    )
    parser.set_defaults(func=summarize_data_command)


def summarize_data_command(args):
    config = EasyFuseConfiguration(args.config_file)
    summary = FusionSummary(
        input_fusions=args.input_fusions,
        input_fusion_context_seqs=args.input_fusion_context_seqs,
        input_requant_cpm=args.input_requant_cpm,
        input_requant_counts=args.input_requant_counts,
        input_reads_stats=args.input_reads_stats,
        output_folder=args.output_folder,
        config=config)
    summary.run(args.model_predictions)
