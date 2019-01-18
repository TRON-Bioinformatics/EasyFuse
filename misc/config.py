"""
Reading/accessing "config.ini" file

@author: Tron (PASO), BNT (URLA)
@version: 20190118
"""


import sys
import os.path
import misc.queue as Queueing

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
    
    def run_self_test(self):
        """Test if set commands work and provided file/dir path exist"""
        print("Starting self test of the config file...")
        count_errors = 0
        for category in self.get_categories():
            # ['general', 'resources', 'commands', 'references', 'indices', 'otherFiles', 'easyfuse_helper', 'liftover']
            if category == "general":
                print("Checking general settings...")
                if len(self.get(category, "tools")) == 0:
                    print("Error 99: No tools selected for running!")
                    count_errors += 1
                elif "Fetchdata" in self.get(category, "tools") and len(self.get(category, "fd_tools")) == 0:
                    print("Error 99: No tools selected for fetchdata run!")
                    count_errors += 1
            elif category == "resources":
                for ressource in self.get_keys(category):
                    try:
                        mem, cpu = self.get(category, ressource).split(",")
                        int(mem)
                        int(cpu)
                    except ValueError:
                        print("Error 99: Incorrect ressource definition for {}. Please provide comma separated integers for CPU and MEM limits!".format(ressource))
                        count_errors += 1
            elif category == "commands":
                print("Checking provided commands...")
                for command in self.get_keys(category):
                    try:
                        Queueing.submit_nonqueue_and_get_stdout(["which", self.get(category, command)])
                    except:
                        print("Error 99: The command assigned to \"{}\" could not be found in your PATH!".format(command))
                        count_errors += 1
            else:
                print("Checking provided path in \"{}\"".format(category))
                for path in self.get_keys(category):
                    if "bowtie" in path:
                        if not os.path.isfile("{}.rev.2.ebwt".format(self.get(category, path))):
                            print("Error 99: Bowtie indices not found. Please make sure that {} is the correct prefix for the index!".format(self.get(category, path)))
                            count_errors += 1
                    elif not os.path.isfile(self.get(category, path)) and not os.path.isdir(self.get(category, path)):
                        print("Error 99: File or directory \"{}\" does not exist. Please check definition of {}!".format(self.get(category, path), path))
                        count_errors += 1
        if count_errors > 0:
            print("Error 99: Identified {} problems with the config file => program execution aborted!".format(count_errors))
            sys.exit(99)
        else:
            print("All tests finished successfully!")
