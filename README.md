# EasyFuse 

<img src="img/easyfuse_workflow.png" style="float: right; margin-right: 10px; margin-top: 5px;">

EasyFuse is a pipeline to detect fusion transcripts from RNA-seq data with high accuracy.
EasyFuse uses five fusion gene detection tools, [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki), [InFusion](https://bitbucket.org/kokonech/infusion/src/master/), [MapSplice2](http://www.netlab.uky.edu/p/bioinfo/MapSplice2), [Fusioncatcher](https://github.com/ndaniel/fusioncatcher), and [SoapFuse](https://sourceforge.net/p/soapfuse/wiki/Home/) along with a powerful read filtering strategy, stringent re-quantification of supporting reads, and machine learning for highly accurate predictions.

- Documentation: [EasyFuse Wiki](https://github.com/TRON-Bioinformatics/EasyFuse/wiki)
- Paper: *available soon*

We recommend using EasyFuse with the Docker container.

<br clear="right"/>

## Usage with Docker

### Download the image

The Docker image can be downloaded from [dockerhub](https://hub.docker.com/r/tronbioinformatics/easyfuse) using the following command:

```
docker pull tronbioinformatics/easyfuse:1.3.4
```

### Download references data

Before running EasyFuse the following reference annotation data needs to be downloaded (~92 GB).

```
# Download reference archive
wget ftp://easyfuse.tron-mainz.de/easyfuse_ref_v2.tar.gz
wget ftp://easyfuse.tron-mainz.de/easyfuse_ref_v2.tar.gz.md5

# Check MD5 sums for consistency
md5sum -c easyfuse_ref_v2.tar.gz.md5 easyfuse_ref_v2.tar.gz

# Extract reference archive
tar xvfz easyfuse_ref_v2.tar.gz

```
### Run EasyFuse

The input FASTQ files needs to be provided in a single folder named `input_fastqs`  (here `</path/to/data>/input_fastqs/`). 

Now EasyFuse can be started by mounting the references data folder and the folder with input FASTQ files:

```
docker run \
  --name easyfuse_container \
  -v </path/to/easyfuse_ref>:/ref \
  -v </path/to/data>:/data \
  --rm \
  -it easyfuse:1.3.4 \
  python /code/easyfuse-1.3.4/processing.py -i /data/input_fastqs/ -o /data/results/

```

The output can be found in `</path/to/data>/results/`. The Output format is described in the wiki page [EasyFuse Output](https://github.com/TRON-Bioinformatics/EasyFuse/wiki/EasyFuse-Output)

## Custom Installation

The EasyFuse pipeline depends on multiple external fusion prediction tools and other dependencies. For example:

  - [STAR](https://github.com/alexdobin/STAR) (2.6.1d) 
  - [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki) (2.6.1d) 
  - [Fusioncatcher](https://github.com/ndaniel/fusioncatcher)(1.00)
  - [MapSplice2](https://github.com/davidroberson/MapSplice2) (2.2.1)
  - [InFusion](https://bitbucket.org/kokonech/infusion/src/master/) (0.8)
  - [SOAPfuse](https://sourceforge.net/projects/soapfuse/) (1.2.7) 
  
The custom installation of EasyFuse is described in the [EasyFuse Wiki](https://github.com/TRON-Bioinformatics/EasyFuse/wiki/Installing-EasyFuse)

