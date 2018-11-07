# PhyloRecomb

Pipeline to independent contrast recombination events by topological tests (ELW) and phylogenetic signal (likelihood mapping analysis). 

### Prerequisites

Install the following dependencies, before you begin:

* [Python 2.7](https://www.python.org/downloads/) 
* [Biopython 1.72](https://biopython.org/wiki/Download) 
* [R (≥ 3.2.0)](https://cran.r-project.org/mirrors.html) 


On Ubuntu, you can use apt to install the required packages:
```
sudo apt install seqtk iqtree mafft
```

Lastly, PhyloRecomb relies on some R packages which can be installed from the terminal or manually. 

```
sudo su - -c "R -e \"install.packages('ape', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('sqldf', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('phytools', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('optparse', repos='http://cran.rstudio.com/')\""
```

Optionally, you can download and install jpHMM (http://jphmm.gobics.de/) to detect recombination in the first place. 


### Installing

Download PhyloRecomb and assign execution permissions. 

```
git clone https://github.com/BBeamud/PhyloRecomb/
cd PhyloRecomb/
chmod +x PhyloRecomb_v.alfa.sh &&  chmod +x scripts/*
```

The downloaded directory also contains things that might be useful for PhyloRecomb users:
* input/ contains different HIV alignments that can be used for recombination and phylogenetic analyses 
* scripts/ contains some scripts that can be used for further analyses and graphic representations 

### Usage

```./PhyloRecomb_v.alfa.sh

usage: ./PhyloRecomb_v.alfa.sh  -f <fasta> -d <y|n> -c <.|file> [-ba balignment] [-ra ralignment ] [-o outdir] [-t threads] [-s] \
  -f|--fasta        FASTA file of your recombinant query sequence \
  -d|--detection    'yes' or 'y' for previous recombination detection by jpHMM \
  -c|--coords       file with recombination coordinates \
  -ba|--balignment  path to the alignment for congruence testing (background alignment) \
  -ra|--ralignment  path to the alignment for relatives analyses (representative alignment) \
  -t|--threads      threads to use \
  -o|--outdir       outdir to save all results \
  -s|--save         save intermediate files \
```

### Input

#### With recombination detection 

#### Without recombination detection

##### FASTA file

##### COORDS file

A BED tab-separated file should be passed with the recombination breakpoints/fragments to be evaluated. In the first field should be the FASTA sequence ID, in the second the position of START of the fragment and in the third field the END of it. 

Example: Suppose our sequence name is 3164 the COORDS file would be: 

```
3164	1	3910
3164	3911	4131
3164	4132	6180
3164	6181	6405
3164	6406	6509
3164	6510	6719
3164	6720	6852
3164	6853	7153
3164	7154	7266
```
