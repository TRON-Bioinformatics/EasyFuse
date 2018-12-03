#!/usr/bin/env python

"""
Reading/accessing "config.ini" file

@author: Tron (PASO), BNT (URLA)
@version: 20181126
"""


import sys

class Config(object):
    """ Class to hold all state information about a run, 
    can be used to restart after a break point 
    """
    def __init__(self,infile):
        """ inits the object and reads a file with state info """
        self._file = infile
        self._categories, self._d, self._key_lists = self._read_config_file()
    
    def _read_config_file(self):
        """ read key-value pairs from the state file self._file """
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
                            continue;
                        temp_d = record_dict.get(category, {})
                        temp_d[words[0]] = words[1]
                        record_dict[category] = temp_d
                        temp_l = key_lists.get(category, [])
                        temp_l.append(words[0])
                        key_lists[category] = temp_l
        return categories, d, key_lists
 
  
    def get(self, category,key):
        """get an dictionary entry named "key" from self._d"""
        return self._d[category][key]
  
    def get_categories(self):
        return self._categories

    def get_keys(self,category):
        return self._key_lists[category]
