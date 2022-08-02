# EasyFuse 

[![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/tron-bioinformatics/EasyFuse)](https://github.com/TRON-Bioinformatics/EasyFuse/releases)
[![Docker Image Version (latest semver)](https://img.shields.io/docker/v/tronbioinformatics/easyfuse?label=docker)](https://hub.docker.com/r/tronbioinformatics/easyfuse)
[![License](https://img.shields.io/badge/license-GPLv3-green)](https://opensource.org/licenses/GPL-3.0)

EasyFuse is a pipeline to detect fusion transcripts from RNA-seq data with high accuracy.
EasyFuse uses five fusion gene detection tools, [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki), [InFusion](https://bitbucket.org/kokonech/infusion/src/master/), [MapSplice2](http://www.netlab.uky.edu/p/bioinfo/MapSplice2), [Fusioncatcher](https://github.com/ndaniel/fusioncatcher), and [SoapFuse](https://sourceforge.net/p/soapfuse/wiki/Home/) along with a powerful read filtering strategy, stringent re-quantification of supporting reads and machine learning for highly accurate predictions.

<p align="center"><img src="img/easyfuse_workflow.png"></p>

 - Documentation: [EasyFuse Wiki](https://github.com/TRON-Bioinformatics/EasyFuse/wiki)
 - Publication: [Weber D, Ibn-Salem J, Sorn P, et al. Nat Biotechnol. 2022](https://doi.org/10.1038/s41587-022-01247-9)

We recommend using EasyFuse with Docker or Singularity.

## Usage


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

### Run EasyFuse with Docker

The Docker image can be downloaded from [dockerhub](https://hub.docker.com/r/tronbioinformatics/easyfuse) using the following command:

```
docker pull tronbioinformatics/easyfuse:latest
```

EasyFuse will require three folders:

* The input data folder containing FASTQ files, in this example `/path/to/data`.
* The reference data folder, in this example `/path/to/easyfuse_ref`
* The output folder, in this example `/path/to/output`

Now EasyFuse can be started by mapping the input data, references and output folders.

Using Docker:

```
docker run \
  --name easyfuse_container \
  -v </path/to/easyfuse_ref>:/ref \
  -v </path/to/data>:/data \
  -v </path/to/output>:/output \
  --rm \
  -it easyfuse:latest \
  python /code/easyfuse/processing.py -i /data -o /output

```

### Run EasyFuse with Singularity

```
singularity exec 
  --containall \
  --bind </path/to/easyfuse_ref>:/ref \
  --bind </path/to/data>:/data \
  --bind </path/to/output>:/output \  
  docker://tronbioinformatics/easyfuse:latest \
  python /code/easyfuse/processing.py -i /data/ -o /output

```

The output can be found in `</path/to/output>/FusionSummary`. The Output format is described in the wiki page [EasyFuse Output](https://github.com/TRON-Bioinformatics/EasyFuse/wiki/EasyFuse-Output)

### Custom Installation

The EasyFuse pipeline depends on multiple external fusion prediction tools and other dependencies. For example:

  - [STAR](https://github.com/alexdobin/STAR) (2.6.1d) 
  - [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki) (1.5.0) 
  - [Fusioncatcher](https://github.com/ndaniel/fusioncatcher)(1.00)
  - [MapSplice2](https://github.com/davidroberson/MapSplice2) (2.2.1)
  - [InFusion](https://bitbucket.org/kokonech/infusion/src/master/) (0.8)
  - [SOAPfuse](https://sourceforge.net/projects/soapfuse/) (1.2.7) 
  
It is recommended to run EasyFuse with Docker or Singularity.


### EasyFuse output format description

EasyFuse creates three output files per run. The file *Sample_Name*_fusRank_1.csv for each batchsample in the `FusionSummary` folder. Each line describes a candidate fusion transcript. \\


| FGID                                    | context_sequence_id | FTID                                                                    | Fusion_Gene | Breakpoint1   | Breakpoint2   | context_sequence_100_id | type      | exon_nr | exon_starts                  | exon_ends                    | exon_boundary1 | exon_boundary2 | exon_boundary | bp1_frame | bp2_frame | frame     | context_sequence                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          | context_sequence_bp | neo_peptide_sequence                                                                               | neo_peptide_sequence_bp | fusioncatcher_detected | fusioncatcher_junc | fusioncatcher_span | starfusion_detected | starfusion_junc | starfusion_span | infusion_detected | infusion_junc  | infusion_span | mapsplice_detected | mapsplice_junc | mapsplice_span | soapfuse_detected | soapfuse_junc  | soapfuse_span  | tool_count | tool_frac | ft_bp_best | ft_a_best        | ft_b_best     | ft_junc_best   | ft_span_best | ft_anch_best | wt1_bp_best | wt1_a_best       | wt1_b_best       | wt1_junc_best   | wt1_span_best   | wt1_anch_best | wt2_bp_best | wt2_a_best      | wt2_b_best      | wt2_junc_best   | wt2_span_best   | wt2_anch_best | ft_bp_cnt_best | ft_a_cnt_best | ft_b_cnt_best | ft_junc_cnt_best | ft_span_cnt_best | ft_anch_cnt_best | wt1_bp_cnt_best | wt1_a_cnt_best | wt1_b_cnt_best | wt1_junc_cnt_best | wt1_span_cnt_best | wt1_anch_cnt_best | wt2_bp_cnt_best | wt2_a_cnt_best | wt2_b_cnt_best | wt2_junc_cnt_best | wt2_span_cnt_best | wt2_anch_cnt_best |
|-----------------------------------------|---------------------|-------------------------------------------------------------------------|-------------|---------------|---------------|-------------------------|-----------|---------|------------------------------|------------------------------|----------------|----------------|---------------|-----------|-----------|-----------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------|----------------------------------------------------------------------------------------------------|-------------------------|------------------------|--------------------|--------------------|---------------------|-----------------|-----------------|-------------------|----------------|---------------|--------------------|----------------|----------------|-------------------|----------------|----------------|------------|-----------|------------|------------------|---------------|----------------|--------------|--------------|-------------|------------------|------------------|-----------------|-----------------|---------------|-------------|-----------------|-----------------|-----------------|-----------------|---------------|----------------|---------------|---------------|------------------|------------------|------------------|-----------------|----------------|----------------|-------------------|-------------------|-------------------|-----------------|----------------|----------------|-------------------|-------------------|-------------------|
| CYTH1_17:78782202:-_EIF3H_8:116755797:- | f008ea580a54bd40    | CYTH1_17:78782202:-_ENST00000591455_EIF3H_8:116755797:-_ENST00000276682 | CYTH1_EIF3H | 17:78782202:- | 8:116755797:- | 83c8f2d5bcc48ee3        | trans     | 9       | 78782202*116726016\\|116755666 | 78782272*116726172\\|116755797 | left_boundary  | right_boundary | both          | 2         | 0         | out_frame | CGGCGGGACGCGAGCGGCGAGCCGGAGCGCGGAGCCCGGCTCCCGCACCATGGAGGAGGACGACAGCTACGATGGCGTCCCGCAAGGAAGGTACCGGCTCTACTGCCACCTCTTCCAGCTCCACCGCCGGCGCAGCAGGGAAAGGCAAAGGCAAAGGCGGCTCGGGAGATTCAGCCGTGAAGCAAGTGCAGATAGATGGCCTTGTGGTATTAAAGATAATCAAACATTATCAAGAAGAAGGACAAGGAACTGAAGTTGTTCAAGGAGTGCTTTTGGGTCTGGTTGTAGAAGATCGGCTTGAAATTACCAACTGCTTTCCTTTCCCTCAGCACACAGAGGATGATGCTGACTTTGATGAAGTCCAATATCAGATGGAAATGATGCGGAGCCTTCGCCATGTAAACATTGATCATCTTCACGTGGGCTGGTATCAGTCCACATACTATGGCTCATTCGTTACCCGGGCACTCC                                                                                                                                                                                                                                                   | 71                  | MEEDDSYDGVPQGRYRLYCHLFQLHRRRSRERQRQRRLGRFSREASADRWPCGIKDNQTLSRRRTRN                                | 7.3                     | 1                      | 0.0891906638483    | 0.649109831341     | 0                   | 0.0             | 0.0             | 1                 | 0.232886733382 | 0.0           | 1                  | 0.21306658586  | 0.0            | 1                 | 0.148651106414 | 0.163516217055 | 4          | 0.8       | 71         | 0.0              | 0.0           | 0.0            | 0.0          | 0            | 71          | 0.0              | 0.0              | 0.0             | 0.0             | 0             | 400         | 0.0445953319242 | 0.0545054056851 | 0.0346852581632 | 0.0247751844023 | 49            | 71             | 0             | 0             | 0                | 0                | 0                | 71              | 0              | 0              | 0                 | 0                 | 0                 | 400             | 9              | 11             | 7                 | 5                 | 49                |
| UQCC1_20:35409345:-_ERBB2_17:39715286:+ | fdeb8916217eaccc    | UQCC1_20:35409345:-_ENST00000374394_ERBB2_17:39715286:+_ENST00000578502 | UQCC1_ERBB2 | 20:35409345:- | 17:39715286:+ | e099eaa09f63a1b3        | trans_inv | 3       | 35409345*39715286            | 35409467*39715292            | left_boundary  | left_boundary  | both          | -1        | 1         | no_frame  | ATTGAGGAACATGGCGTTGCTGGTGCGAGTCCTTACGGAGTTTTGCTCTTGTTGCGCAAGCTGGAGTGCAATGGCGTGATCTTGGCTCACTGCAACCTCTGCCTCCCGGGTTCAAGCGATTCTCCTGCCTCAGCCTCCTGAGTAGCTGGGATTACAGGGACCCA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | 157                 | 0                                                                                                  | 13.0                    | 1                      | 0.0396402950437    | 0.0396402950437    | 0                   | 0.0             | 0.0             | 0                 | 0.0            | 0.0           | 0                  | 0.0            | 0.0            | 0                 | 0.0            | 0.0            | 1          | 0.2       | 157        | 0.113965848251   | 0.0           | 0.0            | 0.0          | 0            | 157         | 0.307212286589   | 8.73077498338    | 0.0297302212828 | 0.0297302212828 | 43            | 400         | 0.0             | 0.0             | 0.0             | 0.0             | 0             | 157            | 23            | 0             | 0                | 0                | 0                | 157             | 62             | 1762           | 6                 | 6                 | 43                | 400             | 0              | 0              | 0                 | 0                 | 0                 |
| ACTB_7:5540676:-_KRT8_12:52898899:-     | 631ae5dab7164484    | ACTB_7:5540676:-_ENST00000414620_KRT8_12:52898899:-_ENST00000547031     | ACTB_KRT8   | 7:5540676:-   | 12:52898899:- | 3bb70dce2c7af939        | trans     | 1       | 5540676*52898469             | 5540771*52898899             | left_boundary  | right_boundary | both          | -1        | -1        | no_frame  | ATTCTGCTGAGTCACTTCGGAGGCCACAGCTGCAACCTTGGGAACAATACGAGTCATTCAATTCACTAGGCAGCGTCCCAGCGCTCACCTGGGGGACTTGGCCACGGTCTCTTGTTAGAAGCCTCTTCATGGACAACGTGGCCAACCCCACCTCACCCCTATCGTCTCCCCCGGCCACACCCTGGGCCAATCGGGCTCCCACGTCAGGGTCTTTGCATAGAGATGGAGTTTCAACATGTTGTCCAGGCTGGTCTTGAACTCGTGACCACAAGTGATCCACCCACCATGGCCTCCCAAAGTGCTGGGATTACAGAGGGCTTCCCTGGAGGCCGCCATTGCAGATGCCGAGCAGCGTGGAGAGCTGGCCATTAAGGATGCCAACGCCAAGTTGTCCGAGCTGGAGGCCGCCCTGCAGCGGGCCAAGCAGGACATGGCGCGGCAGCTGCGTGAGTACCAGGAGCTGATGAACGTCAAGCTGGCCCTGGACATCGAGATCGCCACCTACAGGAAGCTGCTGGAGGGCGAGGAGAGCCGGTGGGTGTGGGTACCTCTGACCGGACCTGCTTCCCTATCCCTGGGACCTGGGGTGGGGACGGTGGGAGCCCCCTGAAGCCCCTTGGGACTTGGGGTCCTGTTGTTCTGGGCCAAGAAGGGCTAGGAGTTGGTCCTGACACCCCATTTGTCTGGGTACAGGCTGGAGTCTGGGATGCAGA | 313                 | RASLEAAIADAEQRGELAIKDANAKLSELEAALQRAKQDMARQLREYQELMNVKLALDIEIATYRKLLEGEESRWVWVPLTGPASLSLGPGVGTVGAP | 0.0                     | 1                      | 0.0346852581632    | 0.53018894621      | 0                   | 0.0             | 0.0             | 0                 | 0.0            | 0.0           | 0                  | 0.0            | 0.0            | 0                 | 0.0            | 0.0            | 1          | 0.2       | 313        | 0.00495503688046 | 0.10405577449 | 0.178381327697 | 0.0          | 34           | 313         | 0.00495503688046 | 0.00991007376093 | 0.0             | 0.0             | 0             | 400         | 0.0             | 0.0             | 0.0             | 0.0             | 0             | 313            | 1             | 21            | 36               | 0                | 34               | 313             | 1              | 2              | 0                 | 0                 | 0                 | 400             | 0              | 0              | 0                 | 0                 | 0                 |


**FGID:** The FGID is an identifier composed of GeneName1_chr1:position1:strand1_GeneName2_chr2:position2:strand2  
**context_sequence_id:** The context sequence id is a unique identifier calculated from 800 bp context sequence (400 upstream and 400 bp downstream from the breakpoint)  
**FTID:** The FTID is a unique identifier composed of GeneName1_chr1:position1:strand1_transcript1_GeneName2_chr2:position2:strand2_transcript2. All transcript combinations are considered.  
**Fusion_Gene:** Fusion Gene is a combination of GeneName1_GeneName2.  
**Breakpoint1:** Breakpoint1 is a combination of chr1:position1:strand1.  
**Breakpoint2:** Breakpoint2 is a combination of chr2:position2:strand2.  
**context_sequence_100_id:** The context sequence 100 id is a unique identifier calculated from 200 bp context sequence (100 upstream and 100 bp downstream from the breakpoint). This shorter sequence is used for primer design.  
**type:** EasyFuse identifies six different types of fusion genes  

| type      | comment                                                                                                                       |
|-----------|-------------------------------------------------------------------------------------------------------------------------------|
| cis_near  | Fusion Genes from neighbouring genes on the same chromosome < 1MB apart                                                       |
| cis_far   | Fusion Genes from neighbouring genes the same chromosome > 1MB apart                                                          |
| cis_inv   | Fusion genes with both partners on different strands of the same chromosome > 1MB apart                                       |
| cis_trans | Fusion Genes on the same strand but both partners not neighbouring, probably a result of chromosomal rearrangement (Deletion) |
| trans     | Fusion Genes with partners on different chromosomes but on the same strand                                                    |
| trans_inv | Fusion Genes with partners on different chromosomes but on different strands                                                  |

**exon_nr:** Number of exons involved in fusion transcript  
**exon_starts:** Genomic starting positions of involved exons according to Ensembl  
**exon_ends:** Genomic ending positions of involved exons according to Ensembl  
**exon_boundary1:** Exon boundary of breakpoint in Gene1 (left boundary is 5' in strand orientation, right boundary is 3' in strand orientation, within means breakpoint is inside exon)  
**exon_boundary2:** Exon boundary of breakpoint in Gene2 (left boundary is 5' in strand orientation, right boundary is 3' in strand orientation,  within means breakpoint is inside exon)  
**exon_boundary:** Which of the partner genes has their breakpoint on an exon boundary

| exon_boundary1 | exon_boundary2 | exon_boundary |
|----------------|----------------|---------------|
| left_boundary  | right_boundary | both          |
| left_boundary  | within         | 5prime        |
| within         | right_boundary | 3prime        |
| within         | within         | no_match      |

**bp1_frame:** Frameshift of Gene1 translation because of breakpoint (0 = no change, -1 is one base less, 2 is two bases more etc.)  
**bp2_frame:** Frameshift of Gene2 translation because of breakpoint (0 = no change, -1 is one base less, 2 is two bases more etc.)  
**frame:** Type of frame for translation of fusion gene  

| frame     | description                                             |
|-----------|---------------------------------------------------------|
| in_frame  | translation not affected (both frames = 0)              |
| neo_frame | translation leads to new, unknown peptide sequence, since one fusion partner is in untranslated region      |
| no_frame  | no translation or translation stops after first partner |
| out_frame | one partner has a frameshift, the other has not         |

**context_sequence:** Genomic sequence 400 bp downstream and upstream from the breakpoint (max 800 bp, shorter if one or both breakpoints are not on exon boundary)  
**context_sequence_bp:** Location of breakpoint on context sequence (400 for an 800 bp context sequence)  
**neo_peptide_sequence:** Translated peptide sequence of context sequence until next stop codon (no_frame=0)   
**neo_peptide_sequence_bp:** Breakpoint on translated peptide sequence  
***tool*_detected:** 1 if Breakpoint was detected by *tool*, 0 if not (fusioncatcher,starfusion,infusion, mapsplice, or soapfuse)    
***tool*_junc:** Junction Reads (reads covering breakpoint) detected by *tool*  
***tool*_span:** Spanning Reads ( read pairs with each partner on one side of breakpoint) detected by *tool*  
**tool_count:** Number of tools detecting fusion gene (1-5)  
**tool_frac:** Fraction of tools out of 5  
***category*_bp_best:** Location of breakpoint on context sequence (400 for an 800 bp context sequence)   
***category*_a_best:**  Fraction of read counts from 1 million reads that map to either context sequence (ft) or wildtype sequence (wt) 400 left of breakpoint  
***category*_b_best:**  Fraction of read counts from 1 million reads that map to either context sequence (ft) or wildtype sequence (wt) 400 right of breakpoint  
***category*_junc_best:** Fraction of read counts from 1 million reads that map to sequence and overlap breakpoint by at least 10 bp  
***category*_span_best:** Fraction of read pairs from 1 million sequenced read pairs, that map to both sides of breakpoint position  
***category*_anch_best:** Maximal read anchor size across all junction reads, where the anchor size for a given read is defined as the minimum distance between read start and breakpoint or read end and the breakpoint.   
***category*_bp_cnt_best:** Location of breakpoint on context sequence ( 400 for an 800 bp context sequence)  
***category*_a_cnt_best:** Number of reads, that map to either context sequence (ft) or wildtype sequence (wt) 400 left of breakpoint  
***category*_b_cnt_best:** Number of reads, that map to either context sequence (ft) or wildtype sequence (wt) 400 right of breakpoint  
***category*_junc_cnt_best:** Number of reads that map to sequence and overlap breakpoint by at least 10 bp  
***category*_span_cnt_best:** Number of read pairs, that map to both sides of breakpoint position  
***category*_anch_cnt_best:** Maximal read anchor size across all junction reads, where the anchor size for a given read is defined as the minimum distance between read start and breakpoint or read end and the breakpoint.

| category | comment                     |
|----------|-----------------------------|
| ft       | Context_sequence            |
| wt1      | Fusion Partner 1 (Wildtype) |
| wt2      | Fusion Partner 2 (Wildtype) |

**prediction_prob:** The predicted probability according to the machine learning model that the fusion candidate is a true positive. 
**prediction_class:** The predicted class (`negative` or `positive`) according to the machine learning model. This classification relies on a user-defined threshold on the `precition_prob` column. 
 

