"""
Module for extracting information from the GFF DB file.
"""

# pylint: disable=E0401
import gffutils  # type: ignore

from bin.fusionannotation.src.breakpoint import Breakpoint
from bin.fusionannotation.src.cds import CDS
from bin.fusionannotation.src.exon import Exon
from bin.fusionannotation.src.feature_validation import filter_cds_by_exons
from bin.fusionannotation.src.transcript import Transcript

FEATURE_TYPE_EXON = "exon"
FEATURE_TYPE_CDS = "CDS"


class GffDbController:
    """Controller for GFF DB access."""
    def __init__(self, db_file):
        self.db = gffutils.FeatureDB(db_file)


    def get_features_overlapping_position(self, bp: Breakpoint, feature_type: str) -> list:
        """Get features overlapping the breakpoint.
        Feature types can be "exon" or "CDS".

        Args:
            bp (Breakpoint): Breakpoint to check for overlapping features
            feature_type (str): Either "exon" or "CDS"

        Returns:
            list: of feature objects
        """
        features = []
        for feature in self.db.region(
            region=f"{bp.chrom}:{bp.pos - 1}-{bp.pos + 1}",
            completely_within=False,
            featuretype=(feature_type),
        ):
            if feature_type == FEATURE_TYPE_EXON:
                features.append(
                    Exon(
                        feature.id,
                        feature.start,
                        feature.stop,
                        feature.attributes["Parent"][0]
                    )
                )
            elif feature_type == FEATURE_TYPE_CDS:
                features.append(
                    CDS(
                        feature.id,
                        feature.start,
                        feature.stop,
                        int(feature.frame),
                        feature.attributes["Parent"][0]
                    )
                )
        return features


    def get_exons_cds_overlapping_position(self, bp: Breakpoint) -> tuple:
        """Get two lists with breakpoint overlapping exons and cds from the database.

        Args:
            db (object): GFF database controller

        Returns:
            tuple: Exons and CDS overlapping the breakpoint
        """

        exons = self.get_features_overlapping_position(bp, "exon")
        cds = self.get_features_overlapping_position(bp, "CDS")
        cds_filt = filter_cds_by_exons(exons, cds)
        return (exons, cds_filt)


    def get_features_from_transcript(self, transcript_id: str, feature_type: str) -> list:
        """Get features from a transcript"""
        features = []
        for feature in self.db.children(
            transcript_id, featuretype=(feature_type), order_by="start", reverse=False
        ):
            if feature_type == FEATURE_TYPE_EXON:
                features.append(
                    Exon(feature.id, feature.start, feature.stop, "")
                )
            elif feature_type == FEATURE_TYPE_CDS:
                features.append(
                    CDS(feature.id, feature.start, feature.stop, int(feature.frame), "")
                )
        return features


    def get_parent_transcript(self, feature_id: str) -> Transcript:
        """Get transcript and gene parents from exon feature"""
        trans_id = ""
        trans_biotype = ""
        gene_id = ""
        gene_name = ""
        gene_biotype = ""
        description = ""
        for parent in self.db.parents(feature_id):
            if parent.featuretype == "transcript":
                trans_id = parent.id
                trans_biotype = parent.attributes["transcript_biotype"][0]
            if parent.featuretype == "gene":
                gene_id = parent.attributes.get("gene_id", "")
                gene_name = parent.attributes.get("gene_name", gene_id)[0]
                gene_biotype = parent.attributes["gene_biotype"][0]
                try:
                    description = parent.attributes["description"][0].replace(";", " ")
                except KeyError:
                    description = "NA"
        return Transcript(
            trans_id,
            trans_biotype,
            gene_name,
            gene_biotype,
            description
        )
