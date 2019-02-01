#!/usr/bin/env python

"""
TODO: doc
"""

import re
import sys
import subprocess
from argparse import ArgumentParser

class VersCont(object):
    def __init__(self, dep_file):
        self.dep_file = dep_file
        self.dep_dict = {}
        self.py_dict = {}

    def load_dep_dict(self):
        """Reads the dependency file and stores its content in a dict"""
        with open(self.dep_file, "r") as deps:
            next(deps) # skip header
            for line in deps:
                line_splitter = line.rstrip("\r\n").split("\t")
                self.dep_dict[line_splitter[0]] = line_splitter[1:]

    def load_py_dict(self, with_pip_not_conda):
        """Load all installed python package version from a pip freeze call into a dict"""
        output = ""
        try:
            if with_pip_not_conda:
                output = subprocess.Popen(["pip", "freeze"], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            else:
                output = subprocess.Popen(["conda", "list"], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            (stdoutdata, stderrdata) = output.communicate()
            
            if stderrdata:
                raise IOError
            output = stdoutdata
        except IOError:
            print("Error 99: The was a problem calling \"pip freeze\". Is pip installed? Trying to get versions from conda...")
            self.load_py_dict(False)

        for py_line in stdoutdata.split("\n"):
            if py_line and not py_line.startswith("#"):
                if with_pip_not_conda:
                    self.py_dict[py_line.split("==")[0]] = py_line.split("==")[1]
                else:
                    self.py_dict[py_line.split()[0]] = py_line.split()[1]

    def get_and_print_tool_versions(self):
        """For each tool in the dep list, try to identify the installed version and compare against the tested"""
        # sanity check that we are not running on anything but linux
        if not "linux" in sys.platform:
            print("Error 1: Easyfuse is expected to be run under linux!")
            sys.exit(1)
        self.load_dep_dict()
        self.load_py_dict(True)
        
        print("Checking versions of required dependencies...")
        for tool in self.dep_dict:
            if self.dep_dict[tool][3] == "shell":
                print("Checking {0}: Tested version is: {1}; Installed version is: {2}".format(tool, self.dep_dict[tool][0], self.get_version_with_version_string(self.dep_dict[tool][1], self.dep_dict[tool][2])))
            elif self.dep_dict[tool][3] == "python":
                try:
                    print("Checking {0}: Tested version is: {1}; Installed version is: {2}".format(tool, self.dep_dict[tool][0], self.py_dict[tool]))
                except KeyError:
                    print("Error 100: Python module {} is required, but not installed. Please install and try again.".format(tool))
                    print("Installed modules are: {}".format(self.py_dict))
                    sys.exit(100)
            else:
                print("Error 110: unknown problem with dep file")
                sys.exit(110)

        print("Please keep in mind that dependencies of dependencies are neither checked nor listed here!")
        print("Please refer to the manual of the individual tools for their specific requirements!")

    def get_version_with_version_string(self, cmd, grp_str):
        """Return version string from program call on shell"""
        output = ""
        try:
            output = subprocess.Popen(cmd.split(" "), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            (stdoutdata, stderrdata) = output.communicate()
            if stdoutdata and not stderrdata:
                output = stdoutdata
            else:
                output = stderrdata
        except:
            output = "Version NA"
        
        version_string = ""
        try:
            count_out_lines = 0
            for version_line in output.split("\n"):
                count_out_lines += 1
                if grp_str in version_line:
                    version_string = re.search('\d+(\.\d+)*[a-z,+]?', version_line).group(0)
        except AttributeError:
            if count_out_lines == len(output.split("\n")):
                print("No version ID found!")
                version_string = "NA"
        
        return(version_string)

def main():
    """Parse command line arguments and start script"""
    parser = ArgumentParser(description='Version control')
    # required arguments
    parser.add_argument('-i', '--input', dest='input', help='Specify the dependency file.', required=True)
    args = parser.parse_args()
    vc = VersCont(args.input)
    vc.get_and_print_tool_versions()

if __name__ == '__main__':
    main()
