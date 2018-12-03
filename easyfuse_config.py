"""
Reading/accessing "config.ini" file

@author: Tron (PASO), BNT (URLA)
@version: 20181126
"""

class Config(object):
    """Class to access all required parameter information about a run"""
    def __init__(self, infile):
        """inits the object and reads a file with state info"""
        self._file = infile
        self._categories, self._d, self._key_lists = self._read_config_file()

    def _read_config_file(self):
        """read key-value pairs from the state file self._file"""
        record_dict = {}
        key_lists = {}
        categories = []
        category = "NONE"
        with open(self._file) as cfg:
            for line in cfg:
                line = line.strip("\n")
                line = line.strip()
                if not (line.startswith("#") or line == ""):
                    if line.startswith("[") and line.endswith("]"):
                        category = line.strip("[").strip("]").strip()
                        categories.append(category)
                        record_dict[category] = {}
                        key_lists[category] = []
                    else:
                        words = line.split("=")
                        words = [x.strip() for x in words]
                        if len(words) < 2:
                            continue
                        temp_d = record_dict.get(category, {})
                        temp_d[words[0]] = words[1]
                        record_dict[category] = temp_d
                        temp_l = key_lists.get(category, [])
                        temp_l.append(words[0])
                        key_lists[category] = temp_l
        return categories, record_dict, key_lists

    def get(self, category, key):
        """get an dictionary entry named "key" from self._d"""
        return self._d[category][key]

    def get_path(self):
        """return the fullpath string of this config file"""
        return self._file

    def get_categories(self):
        """return a list of catgories"""
        return self._categories

    def get_keys(self, category):
        """return a list of keys belonging to a defined category"""
        return self._key_lists[category]
