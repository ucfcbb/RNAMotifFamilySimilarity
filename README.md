## RNAMotifComp

**A comprehensive method to analyze and identify structurally similar RNA motif families**

* Md Mahfuzur Rahaman<sup>†</sup>, mahfuzur at knights dot ucf dot edu
* Nabila Shahnaz Khan<sup>†</sup>, nabilakhan at knights dot ucf dot edu
* Shaojie Zhang*<sup>†</sup>, shzhang at cs dot ucf dot edu

<sup>†</sup>Department Computer Science, University of Central Florida, Orlando, FL 32816 USA \
*To whom correspondence should be addressed.

---

RNAMotifComp is a method to generate similarity graph from RNA structural motifs using either RNAMotifScanX or RNA-align alignment tool. It is developed and tested in a **64-bit Linux** machine. For the basic features, only python (3.x recommended) is required to be installed on the computer with **64-bit Linux** environment. The preliminary files to run this code is included here.

### 1. Installation

#### 1.1 Install prelimineries

All recent Linux systems normally come with python installed by default. If not, please install `python` and `pip` before proceeding to the next step.

#### 1.2: Install required Python libraries

It is required to install several python libraries to run RNAMotifComp. These libraries are included in the [requirements.txt](requirements.txt) file. To install all required python libraries, please navigate to the RNAMotifComp home directory in the terminal and execute the following command.

```
pip install -r requirements.txt
```

### 2. Input Specifications

Motif Family Similarity Analysis tool takes input from any text file (e.g. ‘*.in’) which needs to be in the [RNAMotifFamilySimilarity/data/](data) directory or any subdirectory inside. Each line in that file represents a motif family. Each line starts with the motif family name, followed by a comma-separated list of motifs (the indices for motifs are PDB index, but FASTA index can also be used by setting a parameter in the configuration file). Please check the sample input file ([sample1.in](data/sample1.in)) provided in the [data](data) directory to look into the formatting of input data in detail.

### 3. Commands for usage

```
usage: python3 run.py [-h] -i I [-r [R]][-p [P]][-t [T]]
I - <input_file_name> [Input file name under the data directory (Required). e.g.: sample1.in]
R - <rmsd_threshold> [RMSD threshold to consider as similar. e.g.: 1.0]
P - <participating_motif_instance_threshold> [Percentage of participating motif instances threshold. e.g.: 20.0]
T - <alignment_tool> [Alignment tool to be used, default - ScanX, options: ScanX / RNAalign]
```

**Examples:**

To generate a similarity graph along with the text files containing similar motif instances from [sample1.in](data/sample1.in) using default alignment tool (RNAMotifScanX), use the following command:

```
python3 run.py -i sample1.in
```

To generate a similarity graph along with the text files containing similar motif instances from [sample1.in](data/sample1.in) using ‘RNAalign’ alignment tool, use the following command:

```
python3 run.py -i sample1.in -t RNAalign
```

To generate a similarity graph along with the text files containing similar motif instances from [sample1.in](data/sample1.in) using ‘ScanX’ alignment tool while setting the RMSD threshold as 2.0 and participating motif instances threshold as 10%, use the following command:

```
python3 run.py -i sample1.in -t ScanX -r 2.0 -p 10.0
```

We provided pre-generated data for 360 internal loop motifs from 11 families. For any new dataset, it will automatically download and/or generate required data files (e.g. *.cif, *.fasta, *.aln, etc.) which might take some time. Please make sure to provide valid (not obsolete) PDB number in the input data.

### 4. Output specification

A subdirectory will be created inside the ‘output’ directory using the name of the alignment tool used. Inside that subdirectory there will be an image file (‘&ast;.png’) that represents the generated graph after analyzing the similarities among the input motif families. There will also be a couple of text files (‘&ast;.csv’), each of which represents the similar motif instances for a motif family pair.

### ACKNOWLEDGEMENTS

Motif Family Similarity Analysis is developed for an NIH funded project (R01GM102515).
  
### CONTACTS

The author of the source code is Md Mahfuzur Rahaman. The benchmarking using ML part is authored by Nabila Shahnaz Khan. For bug reports or comments please contact mahfuz@ucf.edu.
