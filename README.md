# EasyFuse 

![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/tron-bioinformatics/EasyFuse)
![Docker Image Version (latest semver)](https://img.shields.io/docker/v/tronbioinformatics/easyfuse?label=docker)
[![License](https://img.shields.io/badge/license-GPLv3-green)](https://opensource.org/licenses/GPL-3.0)

EasyFuse is a pipeline to detect fusion transcripts from RNA-seq data with high accuracy.
EasyFuse uses five fusion gene detection tools, [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki), [InFusion](https://bitbucket.org/kokonech/infusion/src/master/), [MapSplice2](http://www.netlab.uky.edu/p/bioinfo/MapSplice2), [Fusioncatcher](https://github.com/ndaniel/fusioncatcher), and [SoapFuse](https://sourceforge.net/p/soapfuse/wiki/Home/) along with a powerful read filtering strategy, stringent re-quantification of supporting reads and machine learning for highly accurate predictions.

<p align="center"><img src="img/easyfuse_workflow.png" width="400"></p>

- Documentation: [EasyFuse Wiki](https://github.com/TRON-Bioinformatics/EasyFuse/wiki)
- Paper: *available soon*

We recommend using EasyFuse with the Docker container.

## Usage with Docker

### Download the image

The Docker image can be downloaded from [dockerhub](https://hub.docker.com/r/tronbioinformatics/easyfuse) using the following command:

```
docker pull tronbioinformatics/easyfuse:1.3.4
```

### Download reference data

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

EasyFuse will require three folders:
* The input data folder containing FASTQ files, in this example `/path/to/input_data`.
* The reference data folder, in this example `/path/to/easyfuse_ref`
* The output folder, in this example `/path/to/output`

Now EasyFuse can be started by mapping the input data, references and output folders.

```
docker run \
  --name easyfuse_container \
  -v </path/to/easyfuse_ref>:/ref \
  -v </path/to/data>:/data \
  -v </path/to/output>:/output \
  --rm \
  -it easyfuse:1.3.4 \
  python /code/easyfuse-1.3.4/processing.py -i /data -o /output
```


The output can be found in `</path/to/output>/results/`. The Output format is described in the wiki page [EasyFuse Output](https://github.com/TRON-Bioinformatics/EasyFuse/wiki/EasyFuse-Output)

## Custom Installation

The EasyFuse pipeline depends on multiple external fusion prediction tools and other dependencies. For example:

  - [STAR](https://github.com/alexdobin/STAR) (2.6.1d) 
  - [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki) (2.6.1d) 
  - [Fusioncatcher](https://github.com/ndaniel/fusioncatcher)(1.00)
  - [MapSplice2](https://github.com/davidroberson/MapSplice2) (2.2.1)
  - [InFusion](https://bitbucket.org/kokonech/infusion/src/master/) (0.8)
  - [SOAPfuse](https://sourceforge.net/projects/soapfuse/) (1.2.7) 
  
The custom installation of EasyFuse is described in the [EasyFuse Wiki](https://github.com/TRON-Bioinformatics/EasyFuse/wiki/Installing-EasyFuse)

