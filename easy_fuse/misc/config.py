import os
from configparser import ConfigParser


class EasyFuseConfiguration(ConfigParser):

    def __init__(self, config_file):
        super().__init__()
        self.config_file = os.path.abspath(config_file)
        assert os.path.exists(self.config_file), "Config file {} does not exist!".format(self.config_file)
        self.load_configuration()

    def load_configuration(self):
        self.read(self.config_file)

