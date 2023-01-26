import argparse

from easy_fuse.misc.count_input_reads import add_count_reads_args
from easy_fuse.misc.qc_parser import add_qc_parser_args
from easy_fuse.read_selection import add_read_selection_args
from easy_fuse.requantify import add_requantify_args
from easy_fuse.summarize_data import add_summarize_data_args
from easy_fuse.fusionreadfilter import add_read_filter_args
from easy_fuse.fusiontoolparser import add_fusion_parse_args
from easy_fuse.fusionannotation import add_annotation_parser_args
from logzero import logger
import easy_fuse
from easy_fuse.processing import add_pipeline_parser_args
from easy_fuse.tool_wrapper.skewer_wrapper import add_skewer_args
from easy_fuse.tool_wrapper.soapfuse_wrapper import add_soapfuse_wrapper_args
from easy_fuse.tool_wrapper.star_custom_index import add_star_custom_index_args


epilog = "Copyright (c) 2023 TRON gGmbH (See LICENSE for licensing details)"


def easy_fuse_cli():
    # set up logger
    parser = argparse.ArgumentParser(
        description="EasyFuse v{}".format(easy_fuse.__version__),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=epilog,
    )

    subparsers = parser.add_subparsers(description="Commands")

    pipeline_parser = subparsers.add_parser(
        "pipeline", description="Runs the complete EasyFuse pipeline", epilog=epilog
    )
    add_pipeline_parser_args(pipeline_parser)

    annotation_parser = subparsers.add_parser(
        "annotation",
        description="Generate mapping stats for fusion detection",
        epilog=epilog,
    )
    add_annotation_parser_args(annotation_parser)

    fusion_parser = subparsers.add_parser(
        "fusion-parser",
        description="Parse predicted gene fusions from tools",
        epilog=epilog,
    )
    add_fusion_parse_args(fusion_parser)

    read_filter_parser = subparsers.add_parser(
        "read-filter",
        description="Filter reads adequate for fusion detection",
        epilog=epilog,
    )
    add_read_filter_args(read_filter_parser)

    summarize_data_parser = subparsers.add_parser(
        "summarize-data",
        description="Collect stats from fusion detection tools",
        epilog=epilog,
    )
    add_summarize_data_args(summarize_data_parser)

    qc_parser_parser = subparsers.add_parser(
        "qc-parser",
        description="Parses QC data",
        epilog=epilog,
    )
    add_qc_parser_args(qc_parser_parser)

    requantify_parser = subparsers.add_parser(
        "requantify",
        description="Requantifies based on junction/spanning reads",
        epilog=epilog,
    )
    add_requantify_args(requantify_parser)

    read_selection_parser = subparsers.add_parser(
        "read-selection",
        description="Read selection",
        epilog=epilog,
    )
    add_read_selection_args(read_selection_parser)

    skewer_wrapper_parser = subparsers.add_parser(
        "skewer-wrapper",
        description="Runs skewer",
        epilog=epilog,
    )
    add_skewer_args(skewer_wrapper_parser)

    soapfuse_wrapper_parser = subparsers.add_parser(
        "soapfuse-wrapper",
        description="Runs soapfuse",
        epilog=epilog,
    )
    add_soapfuse_wrapper_args(soapfuse_wrapper_parser)

    star_index_parser = subparsers.add_parser(
        "star-index",
        description="Runs a custom STAR index",
        epilog=epilog,
    )
    add_star_custom_index_args(star_index_parser)

    count_reads_parser = subparsers.add_parser(
        "count-reads",
        description="Read count",
        epilog=epilog,
    )
    add_count_reads_args(count_reads_parser)

    args = parser.parse_args()

    try:
        args.func(args)
    except AttributeError as e:
        logger.exception(e)
        parser.parse_args(["--help"])
