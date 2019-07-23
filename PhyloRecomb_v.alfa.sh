#!/bin/bash
# PhyloRecomb: A pipeline to detect and further evaluated recombination events in HIV-1 (ALFA VERSION)
# Author: Beatriz Beamud Aranguren (beatriz.beamud@uv.es) 
## Don't work with sh execution! & Bash version >=4 is required.  


######################
# USAGE & ARGUMENTS ##
######################
if [ ${#@} == 0 ]; then # By default, usage is show. 
programname=$0
function usage {
    echo "usage: $programname  -f <fasta> -d <y|n> -c <.|file> [-ba balignment] [-ra ralignment ] [-o outdir] [-t threads] [-s]"
    echo "  -f|--fasta        FASTA file of our sequence of interest"
    echo "  -d|--detection    'yes' or 'y' for previous recombination detection by jpHMM"
    echo "  -c|--coords        file with recombination coordinates" 
    echo "  -ba|--balignment    path to the alignment for congruence testing (background alignment)"
    echo "  -ra|--ralignment    path to the alignment for relatives analyses (representative alignment)"
    echo "  -t|--threads      threads to use. DEFAULT: 1/auto?" 
    echo "  -o|--outdir       outdir to save all results" 
    echo "  -s|--save         save intermediate files. DEFAULT: NO" 
    exit 1
}
usage
fi

#if [ ${#@} == 0 ]; then
#	usage
#fi	

# Define arguments
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -d|--detection)
    DETECTION="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output)
    OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    -c|--coordenates)
    COORDS="$2"
    shift
    shift
    ;;	
    -f|--fasta)
    FASTA="$2"
    shift # past argument
    shift # past value
    ;;
    -ba|--B_alignment)
    ALN="$2"
    shift # past argument
    shift # past value
    ;;
    -ra|--R_alignment)
    ALN_R="$2"
    shift # past argument
    shift # past value
    ;;
    -t|--threads)
    T="$2"
    shift
    shift
    ;;
    --default)
    DEFAULT=YES
    shift # past argument
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

echo DETECTION  = "${DETECTION}"
echo COORDS = "${COORDS}"
echo FASTA     = "${FASTA}"
echo ALN    = "${ALN}"
echo DEFAULT         = "${DEFAULT}"

# Obtain BED file from COORDS, if provided (Sample_name Init End)
if [ -n "$COORDS" ]; then 
BED=`cat "$COORDS" | cut -f 1,2,3`
SUBTYPES=`cat "$COORDS" | cut -f 4`
# ASSOCIATIVE ARRAY #
# Save the subtype of each fragment 
declare -A SUB_MAP
while read line; do
echo "-----------------------------------------------------------" 
start=`echo $line | cut -d ' ' -f2 | sed 's/ //g'`
#  seqtk coordenate correction 
start=$((start +1))
end=`echo $line | cut -d ' ' -f3 | sed 's/ //g'`
sub=`echo $line | cut -d ' ' -f4 | sed 's/ //g'`
SUB_MAP[$start"-"$end]=$sub
done < $COORDS 
fi 

echo "$BED" > coordenates.bed

# Default options:
if [ -z "$T" ]
then
 T=1
fi



#############################
# RECOMBINATION DETECTION  ##
#############################
# Save jpHMM results
if [ "$DETECTION" == 'yes' ]||[ "$DETECTION" == 'y' ]; then
	jpHMM -s $FASTA -i ./input/HIV_alignment.fas -P ./priors/ -v HIV -o $OUTPUT
	RESULTS="./$OUTPUT/recombination.txt"
	# We delete comment lines, blank lines and identifiers
	PRECOORDS=$(grep -v "#" "$RESULTS" | sed '/^\s*$/d' | grep -v ">" | cut -f 1,2,3)
	TAXA=`grep ">" "$FASTA" | cut -d " " -f1 | sed 's/>//g'`
	COORDS=`while read -r line; do echo -e $TAXA '\t' $line; done <<< "$PRECOORDS"`
	# 5' or 3' regions that are exclusive of the sample are removed for recombinarion analysis. 
	# We create a bed file with identifier and intervals. 
	echo "$COORDS" | grep -v "Insertion" | sed 's/	 //g'| awk -v myvar="$TAXA" -F ' ' '{ print myvar "\t" $2 "\t" $3 }'  > ./$OUTPUT/coordenates.bed 
	#echo "$COORDS" | grep -v "Insertion" | sed 's/	 //g'| cut -d " " -f1,2,3 > ./$OUTPUT/coordenates.bed 
	SUBTYPES=`echo "$COORDS" | grep -v "Insertion"| grep -v "deletion" | sed 's/.*	[0-9]* [0-9]* [0-9]* //g'| sed '/^\s*$/d'`
	# ASSOCIATIVE ARRAY #	
	# Save the subtype of each fragment 
	declare -A SUB_MAP
	while read line; do
	start=`echo $line | cut -d ' ' -f2 | sed 's/ //g'`
	# Correction for seqtk
	start=$((start +1))
	end=`echo $line | cut -d ' ' -f3 | sed 's/ //g'`
	sub=`echo $line | cut -d ' ' -f4 | sed 's/ //g'`
	SUB_MAP[$start"-"$end]=$sub
	done <<< "$COORDS"
	rm posterior_probabilities_for_seq_1.txt alignment_to_msa.txt 
	fi


 
#################################
# ALIGNMENTS  + CONGRUENCE TEST #
#################################

# Extract fragments from sequence
seqtk subseq $FASTA ./$OUTPUT/coordenates.bed > ./$OUTPUT/subseqs.fasta
SEQID=`echo $FASTA | sed 's/\.[^.]*$//'`
touch ./$OUTPUT/summary.txt
echo "# CL: "$@ >> ./$OUTPUT/summary.txt
echo "# SEQID: "$SEQID >> ./$OUTPUT/summary.txt
echo -e "START\tEND\tSUBTYPE\tLMAP(U)\tC.RELATIVES\t" >> ./$OUTPUT/summary.txt
# One file for each fragment (multifasta_2_fasta)
cd $OUTPUT
sed 's/:/_/g' subseqs.fasta | awk -F '>' '/^>/ {F=sprintf("%s.fragment", $2); print > F;next;} {print >>F; close(F)}'
# Align sample with reference consensus alignment 
for fragment in *.fragment; do
name=`echo $fragment| sed 's/.fragment//g'`
mafft --addfragments $fragment --thread -$T ../$ALN > $name"_aln.fasta"
# Extract the region of the alignment where the recombinant fragment is present 
python ../scripts/parse_msa.py $name"_aln.fasta"
iqtree -s $name"_subaln.fasta" -st DNA -m GTR+F+I+G4 -djc -redo -nt AUTO
tree=$name"_subaln.fasta.treefile"
name=`echo $tree | sed 's/.fasta.treefile//g'` 
cat $tree | sed 's/:0.[0-9]*//g' > $name.tre
tre=$name".tre"
taxa_id=`echo $name.tre | sed 's/_subaln.tre//g'`
ref_sub=`grep -m 1 '>' ../$ALN | tr -d '>'` 
#alt=`echo "$SUBTYPES" | sort | uniq | grep -v $ref_sub` # (Option if we want only to test against sample subtypes)
#int=`echo "$alt"| tr '\n' ','`
# Comparison against all subypes present in the alignment
#alt_def="A1,A2,A3,A4,A6,C,D,F1,F2,G,H,J,K,01_AE,A,F"
alt_def =`grep '>' ../$ALN | tr -d '>' | paste -sd "," - |sed  "s/,$ref_sub//g"`  
Rscript --vanilla ../scripts/alternative_trees2.R -t $tre -s $taxa_id -r $ref_sub -a $alt_def -p $name
iqtree -s $name".fasta" -z $name.trees -st DNA -n 0 -zb 10000 -m GTR+F+I+G4 -redo -wsl -pre topotest'_'$name -nt AUTO
# AU from consel
makermt --puzzle topotest'_'$name".sitelh" $name.rmt > $name'.consel.log'
consel $name.rmt >> $name'.consel.log'
rm $name.rmt
consel=`catpv $name.pv | tail -n +3 | sed 's/#//g' | awk '{ print $2 "\t" $4}' | sort -n | sed -e '1d' | cut -f2`
n_trees=`wc -l $name.trees| cut -d " " -f1`
results=`cat "topotest_"$name".iqtree" | grep -A $((n_trees+1)) "Tree      logL    deltaL  bp-RELL    p-KH     p-SH    c-ELW" | sed -e '2d' | sed 's/ - / /g' | sed 's/ + / /g'`
paste -d ' ' <(echo "$results") <(echo "$consel") | awk '{ print $(NF-1) "\t" $NF}' > $name".results"
interval=`echo $name | awk -F "_" '{ print $(NF-1)}' | sed 's/.fsa//g'`
start=`echo $interval | awk -F "-" '{ print $1}'`
end=`echo $interval | awk -F "-" '{ print $2}'`
# We do likelihood mapping, adding a prefix to the filex to no rewritting 
iqtree -s $name".fasta" -lmap 1000 -n 0 -m GTR+F+I+G4 -pre lmap'_'$name -redo -nt AUTO
unresolved=`grep 'Number of unresolved' "lmap_"$name".iqtree" | awk '{ print $NF }' | tr -d '(|=|)'`
# Representative analyses
assign=`Rscript --vanilla ../scripts/parse_congruence4.r -s $name".results" -r $ref_sub -a $alt_def`
assign_def=`echo $assign | sed 's/\[1\] //' | sed 's/"//g' | tr ' ' ','`
# Subset of alignment
# all=`echo $assign_def | awk -F '|' '{for(i=1;i<=NF;i++) print $i}'`
# touch 'crlist_'$fragment 
# while IFS= read -r line; do
#    grep -i -E '
#done <<< "$the_list"
#grep -i -E '>A|>B|>[0-9]*_A|cpx|[0-9]*_B' HIV2_RIP_2017_URFs2018.msa | tr -d '>' | wc -l
mafft --addfragments $fragment --thread -$T ../$ALN_R > $name"_raln.fasta"
iqtree  -s $name"_raln.fasta" -st DNA -m GTR+F+I+G4 -bb 1000 -pre cr'_'$name -djc -nt AUTO -redo
cr_list=`Rscript --vanilla ../scripts/get_brothers_bs70.r -t 'cr_'$name".treefile" -s $name` 
cr=`echo $cr_list | sed 's/\[1\] //' | sed 's/"//g' | tr ' ' ','`
echo -e $start"\t"$end"\t"$assign_def"\t"$unresolved"\t"$cr >> summary.txt
rm *.sitelh *.svg *.gz *.eps *.nex *.contree;done
cd .. 
