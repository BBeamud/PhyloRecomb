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

### Input files 

#### With recombination detection 

* ##### FASTA file


#### Without recombination detection

* ##### FASTA file

* ##### COORDS file

  A BED tab-separated file should be passed with the recombination breakpoints/fragments to be evaluated. In the first field    should be the FASTA sequence ID, in the second the START position of the fragment and in the third field the END of it. 

  Example: Suppose our sequence name is 2011 the COORDS file would be: 

  ```
  2011	1	1472
  2011	1473	2378
  2011	2379	5234
  2011	5235	5378
  2011	5379	8731
  2011	8732	8852
  ```
  
### Output files 

The main results are shown in the **[seqid].summary** tab-delim file. The file looks like this:
```
# CL: ./PhyloRecomb.sh -f 2011.fasta -c 2011.bed 
# SEQID: 2011
START END	SUBTYPE	LMAP(UQ)	CRs
1	1472	A1.CONSENSUS_A1	1.00%	
1473	2378	D.CONSENSUS_D|b.consensus_b	14.70%	
2379	5234	A1.CONSENSUS_A1	1.30%	
5235	5378	G.CONSENSUS_G|a1.consensus_a1|c.consensus_c	37.10%	
5379	8731	D.CONSENSUS_D	0.00%	
8732	8852	A1.CONSENSUS_A1|A2.CONSENSUS_A2	40.60%	
```

The first and second field correspond to the fragment evaluated (same than above). The third field is the subtype(s) assigned at that fragment based on ELW test. The fourth field highlights the phylogenetic signal of that fragment, concretely the unresolved quartets(UQ) of likelihooh mapping analysis. Lastly, CRs are shown if analysis was set. 

These results, in concrete the coordinates, will be also obtained in a reference-based way (in this case, HXB2) 
* **[seqid]_hxb2.summary**

If the user wants to explore each fragment evaluation separately, several files will be available for this purposes:
* **[seqid]\_[fragment].fragment** -> FASTA sequence of the fragment
* **[seqid]\_[fragment]\_aln.fasta** --> Subalignment of the fragment with background alignment
* **[seqid]\_[fragment]\_raln.fasta** --> Subalignment of the fragment with representative alignment (if provided)
* **[seqid]\_[fragment].treefile** --> Original ML tree inferred from above alignment with IQ-TREE in NWK
* **[seqid]\_[fragment].treefile** --> List of all alternative trees in NWK
* **[seqid]\_[fragment].results** --> ELW values for each alternative tree evaluated (in the order of the sequences in the background alignment) 





