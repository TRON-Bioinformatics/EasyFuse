from configparser import NoOptionError, NoSectionError
from unittest import TestCase

import pkg_resources

import easy_fuse.tests
from easy_fuse.misc.config import EasyFuseConfiguration


class TestConfig(TestCase):
    def test_config_loading(self):
        config_file = pkg_resources.resource_filename(
            easy_fuse.tests.__name__, "resources/test.config.ini"
        )
        config = EasyFuseConfiguration(config_file=config_file)
        self.assertIsInstance(config, EasyFuseConfiguration)
        self.assertEqual(config.config_file, config_file)
        self.assertEqual(config.get("general", "pipeline_name"), "EasyFuse")
        self.assertEqual(config.get("general", "min_read_len_perc"), "0.75")
        self.assertEqual(config.get("general", "max_dist_proper_pair"), "200000")
        self.assertEqual(config.get("resources", "qc"), "6,10")
        self.assertEqual(config.get("commands", "fastqc"), "/path/to/fastqc_bin")
        self.assertEqual(
            config.get("references", "genome_fasta"),
            "/path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        )
        self.assertEqual(config.get("indices", "star"), "/path/to/STAR_idx/")
        self.assertEqual(
            config.get("other_files", "infusion_cfg"),
            "/path/to/infusion_index/infusion.cfg",
        )

    def test_non_existing_value(self):
        config_file = pkg_resources.resource_filename(
            easy_fuse.tests.__name__, "resources/test.config.ini"
        )
        config = EasyFuseConfiguration(config_file=config_file)
        with self.assertRaises(NoOptionError):
            config.get("general", "idontexist")

    def test_non_existing_section(self):
        config_file = pkg_resources.resource_filename(
            easy_fuse.tests.__name__, "resources/test.config.ini"
        )
        config = EasyFuseConfiguration(config_file=config_file)
        with self.assertRaises(NoSectionError):
            config.get("idontexist", "pipeline_name")

    def test_non_existing_config_file(self):
        config_file = pkg_resources.resource_filename(
            easy_fuse.tests.__name__, "resources/idontexist.config.ini"
        )
        with self.assertRaises(AssertionError):
            EasyFuseConfiguration(config_file)
