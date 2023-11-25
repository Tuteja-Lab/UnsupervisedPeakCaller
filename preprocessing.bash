#!/bin/bash
# Preprocessing as per Section "Prediction region selection" in Genome Research manuscript [10.1101/gr.277677.123].
# Steps are as in publication. Stages are for debugging.
#
# Input: threshold t or "median", R replicate BAM files, one merged BAM file, size region L
#
## Stage 0: Create merged BAM file if it does not exist
## Stage 1 (-bga.bedGraph): Compute replicate and merged file coverage.
## Stage 2 (_median.txt): Compute threshold if user requests "median".
# -  For each replicate, compute median of nonzero coverage values on each chromsome. Minimum across replicates is threshold t_c for chromosome c.
## Stage 3 (regions1)
# Step 1:
# - Retain genome positions with coverage >t_c in all R individual BAM files.
## Stage 4 (regions2)
# Step 2:
# - Merge retained sites within 90 bp.
# - Retain all regions longer than 100 bp
## Stage 5  (regions3)
# Step 3.1: Extend regions short to L.
# - Extending (L/2) bp up- and down-stream of the midpoint.
# Step 3.2: Split regions larger than L.
# - Find positions with coverage summed across replicates $\geq 0.95$ quantile of the region (obtained from the merged BAM file).
# - Positions within L bp are merged
# - Extended (L/2) bp up- and down-stream of each merged region's midpoint
## Stage 6 (regions4)
# - Remove blacklist regions.
## Stage 7 (covBga.txt): Get coverage of selected regions.
## Stage 8 (bigInputs.txt): Get biginput regions (the original regions)
last_stage=10


### Process command line.
## default values
pdir="."
indir="example"
outdir=""
genome="hg38"
cutoff="median"
threads=1
prefix="out"
inputLn=1000
quiet=0
verbose=0
overwrite=0
zerowrite=0
save=0
ref_prefix="chr"
exclude="MATCHES_NOTHING_ON_PLANET_EARTH"

helpFunction()
{
	echo ""
	echo "Usage: $0 -d INPUT_DIR -b \"BAM_REP1 BAM_REP2[ BAM_REP3...]\" -n NAME"
	echo -e "\t-b STR BAM files of replicate experiments. Provide as space-separated strings surround by double quotes." 
	echo -e "\t       For example: -b \"MCF7_chr10_rep1.bam MCF7_chr10_rep2.bam\""
		#\n\t	Must be sorted and indexed, though if index is missing, we will sort and index for you."
	echo -e "\t-c STR Cutoff for coverage. Either \"median\" or integer (DEFAULT: cutoff)."
	echo -e "\t-d STR Input directory where BAM and BAI input files are stored (DEFAULT: $indir)."
	echo -e "\t-e STR Exclude reference sequence names matching grep extended regexp (DEFAULT: not set)."
	echo -e "\t       For example, set to \"_\" to get rid of unplaced scaffolds."
	echo -e "\t-g STR Genome the data are aligned to. Currently support mm10 (Ensembl), hg38 (Ensembl), or any unique id for another genome (DEFAULT: $genome)."
	echo -e "\t-L INT Length of RCL input segments (DEFAULT: $inputLn)."
	echo -e "\t-m STR Provide the merged BAM file if you have already created it. (DEFAULT: not named)"
	echo -e "\t-n STR File name prefix (DEFAULT: $prefix)."
	echo -e "\t-o STR Output directory where output files will go (DEFAULT: same as input directory)."
	echo -e "\t-p STR Program directory where preprocessing scripts are located (DEFAULT: $pdir)."
	echo -e "\t-t INT Number of threads to use (DEFAULT: $threads)."
	echo -e "\t-q     Quiet, no output of progress (DEFAULT: no)."
	echo -e "\t-r STR Reference sequence name prefix (DEAFULT: $ref_prefix)."
	echo -e "\t-s     Save intermediate files (DEFAULT: no)."
	echo -e "\t-v     Verbose debugging output (DEFAULT: no)."
	echo -e "\t-w INT Overwrite existing files from stage INT (DEFAULT: no)."
	echo -e "\t       Stage 1: Coverage bedgraphs per replicate, merged per reference sequence (*-bga.bedGraph)."
	echo -e "\t       Stage 2: Median coverage per replicate/reference sequence (*_median.txt)."
	echo -e "\t       Stage 3: Sites exceeding threshold (*-regions1.txt: reference, start, end, coverage)."
	echo -e "\t       Stage 4: Candidate regions (*-regions2.txt: ref, reg_start, reg_end, reg_name, ref, cov_start, cov_end, sum_coverage)."
	echo -e "\t       Stage 5: Bed file of initial L-length segments (step 3) (*-regions3.txt)."
	echo -e "\t       Stage 6: Bed file of blacklist filtered L-length segments (*-regions4.txt)."
	echo -e "\t       Stage 7: Bedgraph of candidate segments (*covBga.txt)."
	echo -e "\t       Stage 8: Final candidate regions (*bigInputs.txt)."
	echo -e "\t-z INT Overwrite only zero-length existing files from stage INT (DEFAULT: no)."
	echo -e "\t-?     This help."

	exit 1 # Exit script after printing help
}

stop=0

## read command line
while getopts "?12345678qsvb:c:d:e:g:L:m:n:p:r:o:t:w:z:" opt
do
	case "$opt" in
		b) indivBam="$OPTARG";;
		c) cutoff="$OPTARG";;
		d) indir="$OPTARG";;
		e) exclude="$OPTARG";;
		g) genome="$OPTARG";;
		L) inputLn="$OPTARG";;
		m) mergedBam="$OPTARG";;
		n) prefix="$OPTARG";;
		o) outdir="$OPTARG";;
		p) pdir="$OPTARG";;
		q) quiet=1;;
		r) ref_prefix="$OPTARG";;
		s) save=1;;
		t) threads="$OPTARG";;
		v) verbose=1;;
		w) overwrite=$OPTARG;;
		z) zerowrite=$OPTARG;;
		1) stop=1;;
		2) stop=2;;
		3) stop=3;;
		4) stop=4;;
		5) stop=5;;
		6) stop=6;;
		7) stop=7;;
		8) stop=8;;
		?) helpFunction ;; # Print helpFunction in case parameter is non-existent
	esac
done

if [ $quiet -eq 1 ]; then
	verbose=0
fi

if [ $overwrite -eq 0 ]; then
	overwrite=$((last_stage+1))
fi

## reset defaults based on user input
indir=${indir%/}
if [ $verbose -eq 1 ]; then echo "[STAGE 0] ($0:${LINENO}) Setting input directory to \"$indir\"."; fi

if [ -z "$outdir" ]; then
	outdir=$indir
	if [ $quiet -eq 0 ]; then echo "[STAGE 0] ($0:${LINENO}) Setting output directory to \"$outdir.\""; fi
fi

if [[ $genome == "hg38" ]]; then
	genome="hg"
fi

if [[ $genome == "mm10" ]]; then
	genome="mm"
fi

basename=""

## guess the merged bam file name
if [ -z $mergedBam ]; then
	for rep in $indivBam; do
		mergedBam=${rep/_rep1.bam/_merged.bam}
		if [[ "$mergedBam" == "$rep" ]]; then
			mergedBam=${rep/.bam/_merged.bam}
			basename=${rep/.bam/}
		else
			basename=${rep/_rep1.bam/}
		fi
		if [ $quiet -eq 0 ]; then echo "[STAGE 0] ($0:${LINENO}) Guessing the merged BAM file is \"$indir/$mergedBam\"."; fi
		break
	done
fi

user_prefix=$prefix	# user's prefix
if [ ! -z $basename ]; then
	prefix="$prefix"_"$basename"
fi

## check for missing user input and help user
if [ -z "$pdir" ] || [ -z "$indir" ] || [ -z "$outdir" ] || [ -z "$genome" ]  || [ -z "$cutoff" ] || #|| [ -z "$mergedBam" ]
	[ -z "$indivBam" ] || [ -z "$threads" ] || [ -z "$prefix" ] || [ -z "$inputLn" ] #|| [ -z "$skip" ]
then
   echo "[ERROR] ($0:${LINENO}) Some or all of the parameters are empty";
   helpFunction
fi

### Stage 0 Create any missing files.

## sort and index individual BAM files if they are not yet indexed
# cannot tell if they are not sorted, so will always sort (again) if there is no index
nindivBam=""
for bam in $indivBam; do
	# BAM file missing
	if [ ! -s "$indir/$bam" ]; then
		echo "[ERROR] BAM file $indir/$bam not found or empty."
		cnt=0
		for bam in $indivBam; do
			cnt=$((cnt+1))
		done
		echo -e "[INFO] You named $cnt BAM files \"$indivBam\" on your command line.\n[INFO] There can be no spaces in the BAM file names."

		exit 1
	fi

	# create missing BAM file index
	if [ ! -s "$indir/${bam}.bai" ]; then
		if [[ "$bam" =~ ".*\.bam" ]]; then
			echo "[ERROR] BAM file $indir/$bam must use extension .bam."
			exit 1
		fi
		echo "[STAGE 0] ($0:${LINENO}) Detected no index for $indir/$bam..."
		# no easy way to know if sorted, so sort first
		obam=${bam/.bam/.srt.bam}
		if [ $quiet -eq 0 ]; then echo "[STAGE 0] ($0:${LINENO}) Sorting $indir/$bam and writing to $indir/$obam..."; fi
		if [ $verbose -eq 1 ]; then set -x; fi
		samtools sort -@ $threads -o "$indir"/"$obam" "$indir"/"$bam"
		set +x
		if [ $quiet -eq 0 ]; then echo "[STAGE 0] ($0:${LINENO}) Indexing $indir/$obam..."; fi
		if [ $verbose -eq 1 ]; then set -x; fi
		samtools index -@ $threads "$indir"/"$obam"
		set +x
		if [ -z $nindivBam ]; then
			nindivBam="$obam"
		else
			nindivBam="$nindivBam $obam"
		fi
	else	# this one does have an index!
		if [ -z $nindivBam ]; then
			nindivBam="$bam"
		else
			nindivBam="$nindivBam $bam"
		fi
	fi
done

# overwrite BAM file list
if [ -z "$nindivBam" ]; then
	if [ $verbose -eq 1 ]; then echo "[STAGE 0] ($0:${LINENO}) Resetting individual BAM files to $nindivBam"; fi
	indivBam=$nindivBam
fi

## create and index merged BAM file if it does not exist
if [ ! -s "$indir"/"$mergedBam" ]; then
	if [ $quiet -eq 0 ]; then echo "[STAGE 0] ($0:${LINENO}) Creating merged bam file $indir/$mergedBam..."; fi
	if [ $verbose -eq 1 ]; then set -x; fi
	pushd "$indir"; samtools merge -@ $threads "$mergedBam" $indivBam; popd
	set +x
fi
if [ ! -s "$indir"/"$mergedBam".bai ]; then
	if [ $verbose -eq 1 ]; then set -x; fi
	samtools index -@ $threads "$indir"/"$mergedBam"
	set +x
fi


### STAGE 1: get necessary coverage files (Step 0)

## Step 0.0: read chromosome names

if [ $quiet -eq 0 ]; then echo -n "[STAGE 0] ($0:${LINENO}) References sequences found in mergedBam file: "; fi
if [ $verbose -eq 1 ]; then set -x; fi
chrList=`samtools view -H "$indir"/"$mergedBam" | grep "SN:" | cut -f 2 | grep -v -E $exclude | sed 's/SN://g' | tr "\n" " "`
set +x
chrList=${chrList% }
if [ $quiet -eq 0 ]; then echo $chrList; fi


## Step 0.1: get bga genomecov of the data for each chromosome (-bga.bedGraph)

if [ $quiet -eq 0 ]; then echo "[STAGE 1] ($0:${LINENO}) Getting coverage data..."; fi

# first for merged BAM file
nchrList=""
for chr in $chrList
do
	OUT="$outdir"/$ref_prefix"$chr"/"$prefix"_"$chr"_merged-bga.bedGraph
	if [[ ! -e "$OUT" || $overwrite -le 1 || ( $zerowrite -le 1 && ! -s "$OUT") ]]; then
		# check for reads aligned to this chromosome
		nreads=`samtools view "$indir"/"$mergedBam" $chr | head | wc --lines`
		if [ $nreads -gt 0 ]; then
			if [ $quiet -eq 0 ]; then echo "[STAGE 1] ($0:${LINENO}) Writing merged coverage for reference sequence $ref_prefix$chr to file \"$OUT\"."; fi
			if [ $verbose -eq 1 ]; then set -x; fi
			mkdir --parents "$outdir"/$ref_prefix$chr
			samtools view -hb -@ $threads "$indir"/"$mergedBam" $chr |
				bedtools genomecov -ibam - -pc -bga |
				grep -w "^$chr" \
				> "$OUT"
			set +x
			if [ ! -s "$OUT" ]; then
				echo "[ERROR] Failed to create \"$OUT\" with coverage data."
				exit 1
			fi
			if [ -z "$nchrList" ]; then
				nchrList="$chr"
			else
				nchrList="$nchrList $chr"
			fi
		elif [ $quiet -eq 0 ]; then
			echo "[STAGE 1] ($0:${LINENO}) No reads aligned to reference sequence $chr in \"$indir/$mergedBam\"."
		fi
	else
		if [ $quiet -eq 0 ]; then echo "[STAGE 1] Using existing coverage file $OUT"; fi
		if [ -z "$nchrList" ]; then
			nchrList="$chr"
		else
			nchrList="$nchrList $chr"
		fi
	fi
done
if [ -z "$nchrList" ]; then
	echo "[ERROR] There were no aligned reads in $indir/$mergedBam"
	exit 1
else
	if [ $verbose -eq 1 ]; then echo "[STAGE 1] ($0:${LINENO}) Resetting list of reference names to: $nchrList"; fi
	chrList=$nchrList
fi

# then for individual BAM file
for ifile in $indivBam
do	
	ofile=${ifile/.bam/}
	OUT="$outdir"/"$user_prefix"_"$ofile"-bga.bedGraph
	if [[ ! -e "$OUT" || $overwrite -le 1 || ( $zerowrite -le 1 && ! -s "$OUT") ]]; then
		if [ $quiet -eq 0 ]; then echo "[STAGE 1] ($0:${LINENO}) Writing coverage in replicate $ifile to \"$OUT\"."; fi
		if [ $verbose -eq 1 ]; then set -x; fi
		bedtools genomecov -ibam "$indir"/"$ifile" -pc -bga > "$OUT"
		set +x
		if [ ! -s "$OUT" ]; then
			echo "[ERROR] Failed to create \"$OUT\" with coverage data."
			exit 1
		fi
	else
		if [ $quiet -eq 0 ]; then echo "[STAGE 1] ($0:${LINENO}) Using existing coverage file \"$OUT\""; fi
	fi
done

if [ $stop -eq 1 ]; then 
	if [ $quiet -eq 0 ]; then echo "[STAGE 1] ($0:${LINENO}) Exiting after stage 1."; fi
	exit 1;
fi

## Step 0.2: compute medians for each chromosome, each replicate

OUT="$outdir"/"$prefix"_median.txt
if [[ ! -e "$OUT" || $overwrite -le 2 || ( $zerowrite -le 2 && ! -s "$OUT") ]]; then
	if [ $quiet -eq 0 ]; then echo "[STAGE 2] ($0:${LINENO}) Computing coverage threshold(s)..."; fi
	if [[ $cutoff == "median" ]]
	then
		if [ -s "$OUT" ]; then rm "$OUT"; fi	# since $pdir/getMedian.bash appends, need to truncate first
		for bam in $indivBam
		do
			obam=${bam/.bam/}
			if [ $verbose -eq 1 ]; then set -x; fi
			bash "$pdir"/getMedian.bash "$chrList" "$outdir"/"$user_prefix"_"$obam"-bga.bedGraph "$OUT" $threads
			set +x
		done
		if [ ! -s "$OUT" ]; then
			echo "[ERROR] Failed to create \"$OUT\" with median data."
			exit 1
		fi
	fi
else
	if [ $quiet -eq 0 ]; then echo "[STAGE 2] ($0:${LINENO}) Using existing median file \"$OUT\""; fi
fi

if [ $stop -eq 2 ]; then 
	if [ $quiet -eq 0 ]; then echo "[STAGE 2] ($0:${LINENO}) Exiting after stage 2."; fi
	exit 1;
fi

### Step 1: produce regions

## Step 1.2: get regions surpassing cutoff in every replicate (STAGE 3 regions1)
nchrList=""
for chr in $chrList
do
	# check if file already exists
	todo=0
	for bam in $indivBam; do
		obam=${bam/.bam/}
		OUT="$outdir"/$ref_prefix"$chr"/"$user_prefix"_"$obam"_"$chr"-regions1.txt
		if [[ ! -e "$OUT" || $overwrite -le 3 || ($zerowrite -le 3 && ! -s "$OUT" ) ]]; then
			todo=1
			break
		else
			if [ $quiet -eq 0 ]; then echo "[STAGE 3] ($0:${LINENO}) Using existing regions1 file \"$OUT\""; fi
		fi
	done
	if [ $todo -eq 0 ]; then
		if [ -z "$nchrList" ]; then nchrList="$chr"; else nchrList="$nchrList $chr"; fi
		continue;
	fi

	# obtain (chromosome-specific) cutoff: minimum of medians already computed
	if [[ $cutoff == "median" ]]; then
		if [ $verbose -eq 1 ]; then set -x; fi
		chr_cutoff=`grep -w "^$chr" "$outdir"/"$prefix"_median.txt | sort -k2,2n | cut -f 2 | head -n 1`
		set +x
		if [ $chr_cutoff -eq 0 ]; then
			echo "[WARNING] Minimum median is 0 for reference sequence $chr. Skipping."
			continue
		fi
	else
		chr_cutoff=$cutoff
	fi

	# create bed file with selected regions
	if [ $quiet -eq 0 ]; then echo "[STAGE 3] ($0:${LINENO}) Extracting regions with coverage >$chr_cutoff to file \"$outdir/$ref_prefix$chr/$user_prefix_*_$chr-regions1.txt\""; fi
	if [ $verbose -eq 1 ]; then set -x; fi
	echo $indivBam |
		sed 's/.bam//g' |
		sed 's/ /\n/g' |
		parallel "samtools view -hb \"$indir\"/{}.bam \"$chr\" -@ $threads | bedtools genomecov -ibam - -pc -bga | grep -w \"^$chr\" | awk -v c=$chr_cutoff '{if(\$4>c){print}}' > \"$outdir\"/$ref_prefix\"$chr\"/\"$user_prefix\"_{}_\"$chr\"-regions1.txt"
		# TODO: getAboveThreshold.bash script no longer needed
		#parallel "bash \"$pdir\"/getAboveThreshold.bash \"$indir\"/{} \"$chr\" $chr_cutoff \"$outdir\"/$ref_prefix\"$chr\"/\"$user_prefix\"_{}_\"$chr\"-regions1.txt $threads"
	set +x
	chr_added=0
	for bam in $indivBam; do
		obam=${bam/.bam/}
		OUT="$outdir"/$ref_prefix"$chr"/"$user_prefix"_"$obam"_"$chr"-regions1.txt
		if [ ! -e "$OUT" ]; then
			echo "[ERROR] Failed to create \"$OUT\" with median data."
			exit 1
		elif [ ! -s "$OUT" ]; then
			echo "[WARNING] No regions found for replicate $obam, reference sequence $chr."
		elif [ $chr_added -eq 0 ]; then
			if [ -z "$nchrList" ]; then nchrList="$chr"; else nchrList="$nchrList $chr"; fi
			chr_added=1
		fi
	done
	#echo $indivBam | sed 's/ /\n/g' | parallel "bash \"$pdir\"/getAboveThreshold.bash \"$indir\"/{} \"$i\" \"$c\" \"$outdir\"/chr\"$i\"/\"$i\"_{}-temp.txt \"$threads\""
done

if [ -z "$nchrList" ]; then
	echo "[ERROR] There are no reference sequences with median >0."
	exit 1
else
	chrList=$nchrList
fi

if [ $stop -eq 3 ]; then
	if [ $quiet -eq 0 ]; then echo "[STAGE 3] ($0:${LINENO}) Exiting after stage 3."; fi
	exit 1
fi

## Step 2.1: (STAGE 4 regions2)
 # get positions that pass the threshold in every replicate,
 # merge if within 90 bp, keep if longer than 100 bp,
 # then record merged coverage at every retained position for later analysis.
if [ $quiet -eq 0 ]; then echo "[STAGE 4] ($0:${LINENO}) Finding candidate regions (regions2)..."; fi
nchrList=""
for chr in $chrList
do
	OUT="$outdir"/$ref_prefix"$chr"/"$prefix"_"$chr"-regions2.txt
	if [[ ! -e "$OUT" || $overwrite -le 4 || ($zerowrite -le 4 && ! -s "$OUT") ]]; then
		if [ $quiet -eq 0 ]; then echo "[STAGE 4] ($0:${LINENO}) Writing initial candidate regions to file \"$OUT\""; fi
		if [ $verbose -eq 1 ]; then set -x; fi
		bedops --intersect "$outdir"/$ref_prefix"$chr"/"$user_prefix"_*_"$chr"-regions1.txt |
			bedtools merge -d 90 -i - |
			awk '{if ($3-$2 >= 100) {print}}' |
			awk '{$(NF+1)="segment"NR}1' |
			sed 's/ /\t/g' |
			bedtools intersect -wa -wb -a - -b "$outdir"/$ref_prefix"$chr"/"$prefix"_"$chr"_merged-bga.bedGraph > "$OUT"
		set +x
		if [ ! -e "$OUT" ]; then
			echo "[ERROR] Failed to create \"$OUT\" with region2 data."
			exit 1
		elif [ ! -s "$OUT" ]; then
			echo "[WARNING] No surviving regions in \"$chr\"."
		else
			if [ -z "$nchrList" ]; then nchrList="$chr"; else nchrList="$nchrList $chr"; fi
		fi
	else
		if [ $quiet -eq 0 ]; then echo "[STAGE 4] ($0:${LINENO}) Using existing regions2 file \"$OUT\"."; fi
		if [ -z "$nchrList" ]; then nchrList="$chr"; else nchrList="$nchrList $chr"; fi
	fi
done

if [ -z "$nchrList" ]; then
	echo "[ERROR] There are no reference sequences with median >0."
	exit 1
else
	chrList=$nchrList
fi

if [ $stop -eq 4 ]; then 
	if [ $quiet -eq 0 ]; then echo "[STAGE 4] ($0:${LINENO}) Exiting after stage 4."; fi
	exit 1
fi

## Step 3.1 and 3.2: (STAGE 5 regions3)
 # get $inputLn bp regions
if [ $quiet -eq 0 ]; then echo "[STAGE 5] ($0:${LINENO}) Extracting initial $inputLn bp segments from candidate regions (regions3)...."; fi
arr=($chrList)
m=${#arr[@]}
n=0
nchrList=""
while [ $n -lt $m ]; do
	# gather next up to $threads reference names to process
	subchrList=""
	ncpu=0
	while [ $n -lt $m -a $ncpu -lt $threads ]; do
		chr=${arr[$n]}
		OUT="$outdir"/$ref_prefix"$chr"/"$prefix"_"$chr"-regions3.txt
		if [[ ! -e "$OUT" || $overwrite -le 5 || ($zerowrite -le 5 && ! -s "$OUT") ]]; then 
			if [ -z "$subchrList" ]; then
				subchrList="$chr"
				ncpu=$((ncpu+1))
			else
				subchrList="$subchrList $chr"
			fi
		else
			if [ $quiet -eq 0 ]; then echo "[STAGE 5] ($0:${LINENO}) Using existing regions3 file \"$OUT\""; fi
			if [ -z "$nchrList" ]; then nchrList="$chr"; else nchrList="$nchrList $chr"; fi
		fi
		n=$((n+1))
	done
	
	if [ $ncpu -eq 0 ]; then
		continue
	fi

	# run next up to $threads reference names: add chr to reference name files
	if [ $quiet -eq 0 ]; then echo "[STAGE 5] ($0:${LINENO}) Writing initial $inputLn-length segments to file \"$outdir/${ref_prefix}*/${prefix}_*-regions3.txt\""; fi
	if [ $verbose -eq 1 ]; then set -x; fi
	echo "$subchrList" |
		sed 's/ /\n/g' |
		parallel "Rscript --vanilla \"$pdir\"/getCountFiles.R \"$outdir\"/$ref_prefix{}/\"$prefix\"_{}-regions2.txt \"$inputLn\" $threads \"$outdir\"/$ref_prefix{}/\"$prefix\"_{}-regions3.txt"
	set +x
	for chr in $subchrList; do
		OUT="$outdir"/$ref_prefix"$chr"/"$prefix"_"$chr"-regions3.txt
		if [ ! -e "$OUT" ]; then
			echo "[ERROR] Failed to write \"$OUT\" regions3 file."
			exit 1
		elif [ ! -s "$OUT" ]; then
			echo "[WARNING] Reference sequence \"$ref_prefix$chr\" has no surviving regions in \"$OUT\"."
		else
			if [ -z "$nchrList" ]; then nchrList="$chr"; else nchrList="$nchrList $chr"; fi
		fi
	done
done

if [ -z "$nchrList" ]; then
	echo "[ERROR] There are no reference sequences with median >0."
	exit 1
else
	chrList=$nchrList
fi

if [ $stop -eq 5 ]; then 
	if [ $quiet -eq 0 ]; then echo "[STAGE 5] ($0:${LINENO}) Exiting after stage 5."; fi
	exit 1
fi

## Step 3.3: remove blacklist regions (STAGE 6 regions4)
if [ $quiet -eq 0 ]; then echo "[STAGE 6] ($0:${LINENO}) Removing blacklist regions if possible..."; fi
if [[ $genome == "hg" ]]
then
	echo "[WARNING] Blacklist regions are from Ensembl genome. This will not work if the genome used for alignment is from UCSC."
	bl="$pdir/hg38-blacklist.v2.ensembl.bed"
elif [[ $genome == "mm" ]]
then
	echo "[WARNING] Blacklist regions are from Ensembl genome. This will not work if the genome used for alignment is from UCSC."
	bl="$pdir/mm10-blacklist.v2.ensembl.bed"
else
	echo "[WARNING] Only human and mouse genome blacklists are supported for now. Blacklist regions are Ensembl genome. This will not work if the genome used for alignment is from UCSC. If not mouse or human, will skip removing blacklist regions."
fi

nchrList=""
for chr in $chrList
do
	IN="$outdir"/$ref_prefix"$chr"/"$prefix"_"$chr"-regions3.txt
	if [ ! -s $IN ]; then
		echo "[ERROR] Failed to create \"$IN\" region3 data."
		exit 1
	fi
	OUT="$outdir"/$ref_prefix"$chr"/"$prefix"_"$chr"-regions4.txt
	if [[ ! -e "$OUT" || $overwrite -le 6 || ( $zerowrite -le 6 && ! -s "$OUT") ]]; then
		if [ $quiet -eq 0 ]; then echo "[STAGE 6] ($0:${LINENO}) Writing $inputLn-length without blacklist regions to file \"$OUT\""; fi
		if [[ $genome == "hg" ]] || [[ $genome == "mm" ]]; then
			if [ $verbose -eq 1 ]; then set -x; fi
			if [[ "$chr" == chr* ]]; then
				tail -n +2 "$IN" |
					sed 's/^chr//' |
					bedtools intersect -v -a - -b "$bl" |
					sed 's/^/chr/' > "$OUT"
			else
				tail -n +2 "$IN" |
					bedtools intersect -v -a - -b "$bl" > "$OUT"
			fi
			set +x
		else
			if [ $verbose -eq 1 ]; then set -x; fi
			if [[ "$chr" == chr* ]]; then
	                	tail -n +2 "$IN" > "$OUT"
			else
	                	tail -n +2 "$IN" |
					sed 's/chr//g' > "$OUT"
			fi
			set +x
		fi
	elif [ $quiet -eq 0 ]; then
		echo "[STAGE 6] ($0:${LINENO}) Using existing regions4 file \"$OUT\"."
	fi

	if [ ! -e "$OUT" ]; then
		echo "[ERROR] Failed to create \"$OUT\" with region4 data."
		exit 1
	elif [ ! -s "$OUT" ]; then
		echo "[WARNING] No surviving regions in \"$chr\" after intersection with blacklists."
	else
		if [ -z "$nchrList" ]; then nchrList="$chr"; else nchrList="$nchrList $chr"; fi
	fi
done

if [ -z "$nchrList" ]; then
	echo "[ERROR] There are no reference sequences with peaks outside blacklist regions."
	exit 1
else
	chrList=$nchrList
fi

if [ $stop -eq 6 ]; then
	if [ $quiet -eq 0 ]; then echo "[STAGE 6] ($0:${LINENO}) Exiting after stage 6."; fi
	exit 1
fi

### Step 4: get counts (STAGE 7 *.covBga.txt)
if [ $quiet -eq 0 ]; then echo "[STAGE 7] ($0:${LINENO}) Extracting segment coverage per replicate..."; fi
declare -a rcl_input_files
arr=($chrList)
m=${#arr[@]}
for bam in $indivBam; do
	obam=${bam/.bam/}
	n=0
	while [ $n -lt $m ]; do
		subchrList=""
		ncpu=0
		while [ $n -lt $m -a $ncpu -lt $threads ]; do
			chr=${arr[$n]}
			OUT="$outdir"/$ref_prefix"$chr"/"$user_prefix"_"$obam"_"$chr".covBga.txt
			rcl_input_files+=("$OUT")
			if [[ $overwrite -le 7 || ! -s "$OUT" ]]; then	# there should be no empty files at this stage
				to_run=1
				if [ -z "$subchrList" ]; then
					subchrList="$chr"
				else
					subchrList="$subchrList $chr"
				fi
				ncpu=$((ncpu+1))
			else
				if [ $quiet -eq 0 ]; then echo "[STAGE 7] ($0:${LINENO}) Using existing final RCL input file \"$OUT\""; fi
			fi
			n=$((n+1))
		done
	
		if [ $ncpu -eq 0 ]; then
			continue
		fi

		if [ $quiet -eq 0 ]; then echo -e "[STAGE 7] ($0:${LINENO}) Writing RCL $inputLn-length coverage files to file \"$outdir/${ref_prefix}*/${user_prefix}_${obam}_*.covBga.txt\" for reference sequences\n\t$subchrList"; fi
		if [ $verbose -eq 1 ]; then set -x; fi
		echo "$subchrList" |
			sed 's/ /\n/g' |
			parallel "bedtools intersect -wb -a \"$outdir\"/$ref_prefix{}/\"$prefix\"_{}-regions4.txt -b \"$outdir\"/\"$user_prefix\"_\"$obam\"-bga.bedGraph | cut -f1-4,8 > \"$outdir\"/$ref_prefix{}/\"$user_prefix\"_\"$obam\"_{}.covBga.txt"
		set +x
                #echo ${arr[n]} ${arr[n+1]} ${arr[n+2]} | sed 's/ /\n/g' | parallel "bedtools intersect -wb -a \"$outdir\"/chr\"$i\"/\"$i\"regions4.txt -b \"$outdir\"/\"$prefix\"_{}-bga.bedGraph | cut -f1-4,8 > \"$outdir\"/chr\"$i\"/\"$i\".\"$prefix\"_{}.covBga.txt"
	done
done

if [ $stop -eq 7 ]; then
	if [ $quiet -eq 0 ]; then echo "[STAGE 7] ($0:${LINENO}) Exiting after stage 7."; fi
	exit 1
fi

### Step 5: get big inputs (STAGE 8)
echo $chrList | sed 's/ /\n/g' > "$outdir"/chrList.txt
if [ $quiet -eq 0 ]; then echo "[STAGE 8] ($0:${LINENO}) Getting coordinates of candidate regions..."; fi
candidate_regions_file="$outdir"/"$prefix"_bigInputs.txt
if [ $overwrite -le 8 -o ! -s "$candidate_regions_file" ]; then	# there should be no empty files at this stage
	if [ $quiet -eq 0 ]; then echo "[STAGE 8] ($0:${LINENO}) Writing RCL candidate regions to file \"$candidate_regions_file\""; fi
	if [ $verbose -eq 1 ]; then set -x; fi
	Rscript --vanilla "$pdir"/bigInputs.R "$outdir" "chrList.txt" "$prefix" "$candidate_regions_file" "$ref_prefix"
	set +x
	if [ ! -s "$candidate_regions_file" ]; then
		echo "[ERROR] Failed to create \"$candidate_regions_file\" with candidate regions for RCL."
		exit 1
	fi
else
	if [ $quiet -eq 0 ]; then echo "[STAGE 8] ($0:${LINENO}) Using existing candidate regions file \"$candidate_regions_file\""; fi
fi

if [ $stop -eq 8 ]; then
	if [ $quiet -eq 0 ]; then echo "[STAGE 8] ($0:${LINENO}) Exiting after stage 8."; fi
	exit 1
fi

### STAGE 9: copy files to RCL expected locations
# Copy final files to place consistent with run_rcl.bash. 
# Force overwrite because these files are not linked to this run and can be overwritten any time.
if [ $quiet -eq 0 ]; then echo "[STAGE 9] ($0:${LINENO}) RCL input files being copied to name that run_rcl.bash expects."; fi
for bam in $indivBam; do
	obam=${bam/.bam/}
	for chr in $chrList; do
		IN="$outdir"/$ref_prefix"$chr"/"$user_prefix"_"$obam"_"$chr".covBga.txt
		if [ ! -e $IN ]; then
			echo "[ERROR] Failed to create RCL input data \"$IN\"."
			exit 1
		fi
		cp "$IN" "$outdir"/$ref_prefix"$chr"/"$obam".covBga.txt
	done
done
cp "$candidate_regions_file" "$outdir"/bigInputs.txt

if [ $quiet -eq 0 ]; then
	echo -e "[STAGE 9] ($0:${LINENO}) Data preprocessing has completed succesfully.\n\n"
fi
echo -e "RCL input data are in files:"
for rcl in "${rcl_input_files[@]}"; do
	echo -e "\t$rcl"
done
echo -e "Candidate regions are in file:\n\t$candidate_regions_file"

### STAGE 10: clean up temporary files
if [ $save -eq 0 ]; then
	if [ $quiet -eq 0 ]; then echo "[STAGE 10] ($0:${LINENO}) Cleaning up temporary files..."; fi
	if [ $verbose -eq 1 ]; then set -x; fi
	rm "$outdir"/"$prefix"*-bga.bedGraph
	rm "$outdir"/$ref_prefix*/"$prefix"*-bga.bedGraph
	rm "$outdir"/$ref_prefix*/"$prefix"*-regions[1234].txt
	set +x
fi
