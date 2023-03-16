## RNAMotifComp

**A comprehensive method to analyze and identify structurally similar RNA motif families**

* Md Mahfuzur Rahaman<sup>†</sup>, mahfuzur at knights dot ucf dot edu
* Nabila Shahnaz Khan<sup>†</sup>, nabilakhan at knights dot ucf dot edu
* Shaojie Zhang*<sup>†</sup>, shzhang at cs dot ucf dot edu

<sup>†</sup>Department Computer Science, University of Central Florida, Orlando, FL 32816 USA \
*To whom correspondence should be addressed.

---

This is a tool to generate similarity graph from RNA structural motifs using either RNAMotifScanX or RNA-align alignment tool. It is developed and tested in a **64-bit Linux** machine. Currently, this is a beta version. We will update the latest version as soon as it is ready. For the basic features, python (3.x recommended) and PyMOL are required to be installed on the computer with the previously mentioned environments. The preliminary files to run this code is included here.

### 2. Input Specifications

Motif Family Similarity Analysis tool takes input from any text file (e.g. ‘*.in’) which needs to be in the [RNAMotifFamilySimilarity/data/](data) directory or any subdirectory inside. Each line in that file represents a motif family. The motif family starts with a name, followed by a comma-separated list of motifs (the indices for motifs are expected to be in the PBD index, but it can be changed to Fasta index by setting a parameter in the configuration file). To see examples of formats, please check the sample input file ([sample1.in](data/sample1.in)) provided in the [data](data) directory.

### 3. Commands for usage

Note: **MacOS users** might get an error message saying `'align_ga.mac' cannot be opened because it is from an unidentified developer`. To get rid of this error, please navigate to: 'System Preferences > Security & Privacy > General' and set 'Allow apps downloaded from' to 'Anywhere'.

```
usage: run.py [-h] -i I [-t [T]]
I - <input_file_name> [Get text outputs for user input file - required]
R - <rmsd_threshold> [RMSD threshold to consider as similar]
P - <participating_motif_instance_threshold> [Percentage of participating motif instances threshold]
T - <alignment_tool> [Specify alignment tool to be used, default - ScanX]
```

**Example:**
To generate a similarity graph along with the text files containing similar motif instances from [sample1.in](data/sample1.in) using default alignment tool:

```
python run.py -i sample1.in
```

To generate a similarity graph along with the text files containing similar motif instances from [sample1.in](data/sample1.in) using ‘ScanX’ alignment tool:

```
python run.py -i sample1.in -t ScanX
```

To generate a similarity graph along with the text files containing similar motif instances from [sample1.in](data/sample1.in) using ‘ScanX’ alignment tool while setting the RMSD threshold as 2.0 and participating motif instances threshold as 10%:

```
python run.py -i sample1.in -t ScanX -r 2.0 -p 10.0
```

We provided pregenerated data for 360 internal loop motifs from 11 families. For any new dataset, it will automatically download and/or generate required data files (e.g. *.cif, *.fasta, *.aln, etc.) which might take some time. Please make sure to provide valid (not obsolete) PDB number in the input data, otherwise it might end up getting 404 error while downloading PDB and FASTA files.

### 4. Output specification

A subdirectory will be created inside the ‘output’ directory using the name of the alignment tool used. Inside that subdirectory there will be an image file (‘&ast;.png’) that represents the generated graph after analyzing the similarities among the input motif families. There will also be a couple of text files (‘&ast;.csv’), each of which represents the similar motif instances for a motif family pair.

### ACKNOWLEDGEMENTS

Motif Family Similarity Analysis is developed for an NIH funded project (R01GM102515).
  
### CONTACTS

The author of the source code is Md Mahfuzur Rahaman. The benchmarking using ML part is authored by Nabila Shahnaz Khan. For bug reports or comments please contact mahfuzur@knights.ucf.edu
