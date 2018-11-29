
import sys,os
#sys.path.insert( 0, os.path.dirname( __file__ ) )
#from common import submit
#from common import display

from argparse import ArgumentParser

class Samples(object):
    def __init__(self,infile):
        self.infile = infile
        if not os.path.exists(self.infile):
            open(self.infile,'a').close()
        self.map = self.read_file()


    def get(self,sample_id):
        state = ''
        if sample_id in self.map:
            state = self.map[sample_id]
        return state

    def add(self,sample_id,state):
        if not sample_id in self.map:
            self.map[sample_id] = state
        self.write_file()

    def update(self,sample_id,state):
        if sample_id in self.map:
            self.map[sample_id] = state
        self.write_file()


    def delete(self,sample_id):
        if sample_id in self.map:
            del self.map[sample_id]
        self.write_file()

    def read_file(self):
        map = {}
        with open(self.infile) as inf:
            for line in inf:
                elements = line.rstrip().split("\t")
                if len(elements) > 0:
                    sample_id = elements[0]
                    state = elements[1]
                    map[sample_id] = state
        return map

    def write_file(self):
        outf = open(self.infile,'w')
        for key in self.map.iterkeys():
            outf.write(key + "\t" + self.map[key]+"\n")
        outf.close()


def main():
    parser = ArgumentParser(description='Handle sample db operations')

    parser.add_argument('-i', '--sample-id', dest='sample_id', help='Specify the sample id to process.',required=True )
    parser.add_argument('-o', '--output-file', dest='output_file', help='Specify the file to save the changes into.',required=True )
    parser.add_argument('-g', '--state', dest='state', help='Specify the new state of your sample id.')
    parser.add_argument('-l', '--action', dest='action', choices=['update','delete','add','get'],help='Select the action to do with the sample id.',required=True)
    args = parser.parse_args()

    
    samples = Samples(args.output_file)
    if args.action == "update":
        samples.update(args.sample_id,args.state)
    elif args.action == "add":
        samples.add(args.sample_id)
    elif args.action == "delete":
        samples.delete(args.sample_id)
    elif args.action == "get":
        samples.get(args.sample_id)
#    api_key = "8299479a72756f05ad267949e299ee65"
#    api_url = "http://192.168.171.195/ngs_lims/api/"
#    proj_id = "959164618"
#    sample_id = "IT_N_230_SL1"
#    samples.get_sample_information(api_key,api_url,proj_id,sample_id)


if __name__ == '__main__':
    main()
