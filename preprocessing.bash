#!/bin/bash
# Preprocessing as per Section "Prediction region selection" in Genome Research manuscript [10.1101/gr.277677.123].
#
# Input: threshold t or "median", R replicate BAM files, one merged BAM file, size region L
#
# Step 0: Compute threshold if user requests "median".
# -  For each replicate, compute median of nonzero coverage values on each chromsome. Minimum across replicates is threshold t_c for chromosome c.
# Step 1: produces regions1
# - Retain genome positions with coverage >t_c in all R individual BAM files.
# Step 2.1: produce regions2
# - Merge retained sites within 90 bp.
# - Retain all regions longer than 100 bp
# Step 3.1: Extend regions short to L.
# - Extending (L/2) bp up- and down-stream of the midpoint.
# Step 3.2: Split regions larger than L.
# - Find positions with coverage summed across replicates $\geq 0.95$ quantile of the region (obtained from the merged BAM file).
# - Positions within L bp are merged
# - Extended (L/2) bp up- and down-stream of each merged region's midpoint
# Step 3.3: Remove blacklist regions.
# Step 4: Get coverage of selected regions.
# Step 5: Get biginput regions (the original regions)


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
save=0
ref_prefix="chr"

helpFunction()
{
	echo ""
	echo "Usage: $0 -d INPUT_DIR -b \"BAM_REP1 BAM_REP2[ BAM_REP3...]\" -n NAME"
	echo -e "\t-p Program directory where preprocessing scripts are located (DEFAULT: $pdir)."
	echo -e "\t-d Input directory where BAM and BAI input files are stored (DEFAULT: $indir)."
	echo -e "\t-o Output directory where output files will go (DEFAULT: same as input directory)."
	echo -e "\t-g Genome the data are aligned to. Currently support mm10 (Ensembl), hg38 (Ensembl), or any unique id for another genome (DEFAULT: $genome)."
	echo -e "\t-c Cutoff for coverage. Either \"median\" or integer (DEFAULT: cutoff)."
	echo -e "\t-m Provide the merged BAM file if you have already created it. (DEFAULT: not named)"
	echo -e "\t-b BAM files of replicate experiments. Provide as space-separated strings surround by double quotes." 
	echo -e "\t   For example: -b \"MCF7_chr10_rep1.bam MCF7_chr10_rep2.bam\""
		#\n\t	Must be sorted and indexed, though if index is missing, we will sort and index for you."
	echo -e "\t-t Number of threads to use (DEFAULT: $threads)."
	echo -e "\t-n File name prefix (DEFAULT: $prefix)."
	echo -e "\t-L Length of RCL input segments (DEFAULT: $inputLn)."
	echo -e "\t-r Reference sequence name prefix (DEAFULT: $ref_prefix)."
	echo -e "\t-s Save intermediate files (DEAFULT: no)."
	echo -e "\t-w Overwrite existing files (DEFAULT: no)."
	echo -e "\t-q Quiet, no output of progress (DEFAULT: no)."
	echo -e "\t-v Verbose debugging output (DEFAULT: no)."
	echo -e "\t-? This help."

	exit 1 # Exit script after printing help
}

stop0=0
stop1=0
stop2=0
stop3=0
stop4=0

## read command line
while getopts "svqw01234?p:d:o:g:c:m:b:t:n:L:r:" opt
do
	case "$opt" in
		p) pdir="$OPTARG";;
		d) indir="$OPTARG";;
		o) outdir="$OPTARG";;
		g) genome="$OPTARG";;
		c) cutoff="$OPTARG";;
		m) mergedBam="$OPTARG";;
		b) indivBam="$OPTARG";;
		t) threads="$OPTARG";;
		n) prefix="$OPTARG";;
		L) inputLn="$OPTARG";;
		q) quiet=1;;
		r) ref_prefix="$OPTARG";;
		s) save=1;;
		v) verbose=1;;
		w) overwrite=1;;
		s) skip="$OPTARG";;
		0) stop0=1;;
		1) stop1=1;;
		2) stop2=1;;
		3) stop3=1;;
		4) stop4=1;;
		?) helpFunction ;; # Print helpFunction in case parameter is non-existent
	esac
done

if [ $quiet -eq 1 ]; then
	verbose=0
fi


## reset defaults based on user input
indir=${indir%/}
if [ $verbose -eq 1 ]; then echo "[STATUS] Setting input directory to \"$indir\"."; fi

if [ -z "$outdir" ]; then
	outdir=$indir
	if [ $quiet -eq 0 ]; then echo "[STATUS] Setting output directory to \"$outdir.\""; fi
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
		if [ $quiet -eq 0 ]; then echo "[STATUS] Guessing the merged BAM file is \"$mergedBam\"."; fi
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
   echo "[STATUS] Some or all of the parameters are empty";
   helpFunction
fi

### Step 0 Create any missing files.

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
		echo -e "[STATUS] Detected $cnt BAM files. If this is not correct, check your command line.\nNOTE: There can be no spaces in the BAM file names."

		exit 1
	fi

	# create missing BAM file index
	if [ ! -s "$indir/${bam}.bai" ]; then
		if [[ "$bam" =~ ".*\.bam" ]]; then
			echo "[ERROR] BAM file $indir/$bam must use extension .bam."
			exit 1
		fi
		echo "[STATUS] Detected no index for $indir/$bam..."
		# no easy way to know if sorted, so sort first
		obam=${bam/.bam/.srt.bam}
		if [ $quiet -eq 0 ]; then echo "[STATUS] Sorting $indir/$bam and writing to $indir/$obam..."; fi
		if [ $verbose -eq 1 ]; then set -x; fi
		samtools sort -@ $threads -o "$indir"/"$obam" "$indir"/"$bam"
		set +x
		if [ $quiet -eq 0 ]; then echo "[STATUS] Indexing $indir/$obam..."; fi
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
	if [ $verbose -eq 1 ]; then echo "[STATUS] Resetting individual BAM files to $nindivBam"; fi
	indivBam=$nindivBam
fi

## create and index merged BAM file if it does not exist
if [ ! -s "$indir"/"$mergedBam" ]; then
	if [ $quiet -eq 0 ]; then echo "[STATUS] Creating merged bam file $indir/$mergedBam..."; fi
	if [ $verbose -eq 1 ]; then set -x; fi
	pushd "$indir"; samtools merge -@ $threads "$mergedBam" $indivBam; popd
	set +x
fi
if [ ! -s "$indir"/"$mergedBam".bai ]; then
	if [ $verbose -eq 1 ]; then set -x; fi
	samtools index -@ $threads "$indir"/"$mergedBam"
	set +x
fi


### Step 0: get necessary coverage files

## Step 0.0: read chromosome names

if [ $quiet -eq 0 ]; then echo -n "[STATUS] Getting chromosomes from mergedBam file: "; fi
if [ $verbose -eq 1 ]; then set -x; fi
chrList=`samtools view -H "$indir"/"$mergedBam" | grep "SN:" | cut -f 2 | sed 's/SN://g' | tr "\n" " "`
set +x
chrList=${chrList% }
if [ $quiet -eq 0 ]; then echo $chrList; fi


## Step 0.1: get bga genomecov of the data for each chromosome (-bga.bedGraph)

if [ $quiet -eq 0 ]; then echo "[STATUS] Getting coverage data..."; fi

# first for merged BAM file
nchrList=""
for chr in $chrList
do
	if [ $verbose -eq 1 ]; then set -x; fi
	mkdir --parents "$outdir"/$ref_prefix$chr
	set +x
	OUT="$outdir"/$ref_prefix"$chr"/"$prefix"_"$chr"_merged-bga.bedGraph
	if [ $overwrite -eq 1 -o ! -s "$OUT" ]; then
		# check for reads aligned to this chromosome
		nreads=`samtools view "$indir"/"$mergedBam" $chr | head | wc --lines`
		if [ $nreads -gt 0 ]; then
			if [ $quiet -eq 0 ]; then echo "[STATUS] Computing merged coverage for reference sequence $chr."; fi
			if [ $verbose -eq 1 ]; then set -x; fi
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
			echo "[STATUS] No reads aligned to reference sequence $chr in \"$indir/$mergedBam\"."
		fi
	else
		if [ $quiet -eq 0 ]; then echo "[STATUS] Using existing coverage file $OUT"; fi
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
	if [ $verbose -eq 1 ]; then echo "[STATUS] Resetting list of reference names to: $nchrList"; fi
	chrList=$nchrList
fi

# then for individual BAM file
for ifile in $indivBam
do	
	ofile=${ifile/.bam/}
	OUT="$outdir"/"$user_prefix"_"$ofile"-bga.bedGraph
	if [ $overwrite -eq 1 -o ! -s "$OUT" ]; then
		if [ $quiet -eq 0 ]; then echo "[STATUS] Computing coverage in replicate $ifile."; fi
		if [ $verbose -eq 1 ]; then set -x; fi
		bedtools genomecov -ibam "$indir"/"$ifile" -pc -bga > "$OUT"
		set +x
		if [ ! -s "$OUT" ]; then
			echo "[ERROR] Failed to create \"$OUT\" with coverage data."
			exit 1
		fi
	else
		if [ $quiet -eq 0 ]; then echo "[STATUS] Using existing coverage file \"$OUT\""; fi
	fi
done

## Step 0.2: compute medians for each chromosome, each replicate

OUT="$outdir"/"$prefix"_median.txt
if [ $overwrite -eq 1 -o ! -s "$OUT" ]; then
	if [ $quiet -eq 0 ]; then echo "[STATUS] Getting threshold -t..."; fi
	if [[ $cutoff == "median" ]]
	then
		if [ ! -s "$OUT" ]; then rm "$OUT"; fi	# since $pdir/getMedian.bash appends, need to truncate first
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
	if [ $quiet -eq 0 ]; then echo "[STATUS] Using existing median file \"$OUT\""; fi
fi

if [ $stop0 -eq 1 ]; then 
	if [ $quiet -eq 0 ]; then echo "[STATUS] Exiting after step 0."; fi
	exit 1;
fi

### Step 1: produce regions

## Step 1.2: get regions surpassing cutoff in every replicate (regions1)
nchrList=""
for chr in $chrList
do
	# check if file already exists
	todo=0
	for bam in $indivBam; do
		obam=${bam/.bam/}
		OUT="$outdir"/$ref_prefix"$chr"/"$user_prefix"_"$obam"_"$chr"-regions1.txt
		if [ ! -s "$OUT" -o $overwrite -eq 1 ]; then
			todo=1
			break
		else
			if [ $quiet -eq 0 ]; then echo "[STATUS] Using existing regions1 file \"$OUT\""; fi
		fi
	done
	if [ $todo -eq 0 ]; then
		if [ -z "$nchrList" ]; then nchrList="$chr"; else nchrList="$nchrList $chr"; fi
		continue;
	fi

	echo "[DEBUGGING] nchrList = $nchrList"

	# obtain (chromosome-specific) cutoff: minimum of medians already computed
	if [[ $cutoff == "median" ]]; then
		if [ $verbose -eq 1 ]; then set -x; fi
		chr_cutoff=`grep -w "^$chr" "$outdir"/"$prefix"_median.txt | sort -k2,2n | cut -f 2 | head -n 1`
		set +x
		if [ $chr_cutoff -eq 0 ]; then
			echo "[WARNING] Minimum median is 0 for reference sequence $chr. Skipping."
			continue
		else
			if [ -z "$nchrList" ]; then nchrList="$chr"; else nchrList="$nchrList $chr"; fi
		fi
	else
		chr_cutoff=$cutoff
	fi

	# create bed file with selected regions
	if [ $verbose -eq 1 ]; then set -x; fi
	echo $indivBam |
		sed 's/.bam//g' |
		sed 's/ /\n/g' |
		parallel "samtools view -hb \"$indir\"/{}.bam \"$chr\" -@ $threads | bedtools genomecov -ibam - -pc -bga | grep -w \"^$chr\" | awk -v c=$chr_cutoff '{if(\$4>c){print}}' > \"$outdir\"/$ref_prefix\"$chr\"/\"$user_prefix\"_{}_\"$chr\"-regions1.txt"
		# TODO: getAboveThreshold.bash script no longer needed
		#parallel "bash \"$pdir\"/getAboveThreshold.bash \"$indir\"/{} \"$chr\" $chr_cutoff \"$outdir\"/$ref_prefix\"$chr\"/\"$user_prefix\"_{}_\"$chr\"-regions1.txt $threads"
	set +x
	for bam in $indivBam; do
		obam=${bam/.bam/}
		OUT="$outdir"/$ref_prefix"$chr"/"$user_prefix"_"$obam"_"$chr"-regions1.txt
		if [ ! -s "$OUT" ]; then
			echo "[ERROR] Failed to create \"$OUT\" with median data."
			exit 1
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

if [ $stop1 -eq 1 ]; then
	if [ $quiet -eq 0 ]; then echo "[STATUS] Exiting after step 1."; fi
	exit 1
fi

## Step 2.1: (regions2)
 # get positions that pass the threshold in every replicate,
 # merge if within 90 bp, keep if longer than 100 bp,
 # then record merged coverage at every retained position for later analysis.
if [ $quiet -eq 0 ]; then echo "[STATUS] Getting bases that pass threshold..."; fi
nchrList=""
for chr in $chrList
do
	OUT="$outdir"/$ref_prefix"$chr"/"$prefix"_"$chr"-regions2.txt
	if [ ! -s "$OUT" -o $overwrite -eq 1 ]; then
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
		if [ $quiet -eq 0 ]; then echo "[STATUS] Using existing regions2 file \"$OUT\"."; fi
		if [ -z "$nchrList" ]; then nchrList="$chr"; else nchrList="$nchrList $chr"; fi
	fi
done

if [ -z "$nchrList" ]; then
	echo "[ERROR] There are no reference sequences with median >0."
	exit 1
else
	chrList=$nchrList
fi

if [ $stop2 -eq 1 ]; then 
	if [ $quiet -eq 0 ]; then echo "[STATUS] Exiting after step 2."; fi
	exit 1
fi

## Step 3.1 and 3.2: (regions3)
 # get $inputLn bp regions
if [ $quiet -eq 0 ]; then echo "[STATUS] Getting $inputLn bp regions..."; fi
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
		if [ $overwrite -eq 1 -o ! -s "$OUT" ]; then
			if [ -z "$subchrList" ]; then
				subchrList="$chr"
				ncpu=$((ncpu+1))
			else
				subchrList="$subchrList $chr"
			fi
		else
			if [ $quiet -eq 0 ]; then echo "[STATUS] Using existing regions3 file \"$OUT\""; fi
			if [ -z "$nchrList" ]; then nchrList="$chr"; else nchrList="$nchrList $chr"; fi
		fi
		n=$((n+1))
	done
	
	if [ $ncpu -eq 0 ]; then
		continue
	fi

	# run next up to $threads reference names: add chr to reference name files
	if [ $verbose -eq 1 ]; then set -x; fi
	echo "$subchrList" |
		sed 's/ /\n/g' |
		parallel "Rscript --vanilla \"$pdir\"/getCountFiles.R \"$outdir\"/$ref_prefix{}/\"$prefix\"_{}-regions2.txt \"$inputLn\" $threads \"$outdir\"/$ref_prefix{}/\"$prefix\"_{}-regions3.txt"
	set +x
	for chr in $subchrList; do
		OUT="$outdir"/$ref_prefix"$chr"/"$prefix"_"$chr"-regions3.txt
		if [ ! -s "$OUT" ]; then
			if [ ! -e "$OUT" ]; then
				echo "[ERROR] Failed to write \"$OUT\" regions3 file."
				exit 1
			else
				echo "[WARNING] Reference sequence \"$chr\" has no surviving regions."
			fi
		else
			if [ -z "$nchrList" ]; then nchrList="$chr"; else nchrList="$nchrList $chr"; fi
		fi
	done
done

## Step 3.3: remove blacklist regions (regions4)
if [ $quiet -eq 0 ]; then echo "[STATUS] Removing blacklist regions if possible..."; fi
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

for chr in $chrList
do
	IN="$outdir"/$ref_prefix"$chr"/"$prefix"_"$chr"-regions3.txt
	if [ ! -s $IN ]; then
		echo "[ERROR] Failed to create \"$IN\" region3 data."
		exit 1
	fi
	OUT="$outdir"/$ref_prefix"$chr"/"$prefix"_"$chr"-regions4.txt
	if [ $overwrite -eq 1 -o ! -s "$OUT" ]; then
		if [[ $genome == "hg" ]] || [[ $genome == "mm" ]]; then
			if [ $verbose -eq 1 ]; then set -x; fi
			tail -n +2 "$IN" |
				sed 's/chr//g' |
				bedtools intersect -v -a - -b "$bl" > "$OUT"
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
		echo "[STATUS] Using existing regions4 file \"$OUT\"."
	fi
done

if [ $stop3 -eq 1 ]; then
	if [ $quiet -eq 0 ]; then echo "[STATUS] Exiting after step 3."; fi
	exit 1
fi

### Step 4: get counts (*.covBga.txt)
if [ $quiet -eq 0 ]; then echo "[STATUS] Getting counts per replicate..."; fi
declare -a rcl_input_files
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
			if [ $overwrite -eq 1 -o ! -s "$OUT" ]; then
				to_run=1
				if [ -z "$subchrList" ]; then
					subchrList="$chr"
				else
					subchrList="$subchrList $chr"
				fi
				ncpu=$((ncpu+1))
			else
				if [ $quiet -eq 0 ]; then echo "[STATUS] Using existing final RCL input file \"$OUT\""; fi
			fi
			n=$((n+1))
		done
	
		if [ $ncpu -eq 0 ]; then
			continue
		fi

		if [ $verbose -eq 1 ]; then set -x; fi
		echo "$subchrList" |
			sed 's/ /\n/g' |
			parallel "bedtools intersect -wb -a \"$outdir\"/$ref_prefix{}/\"$prefix\"_{}-regions4.txt -b \"$outdir\"/\"$user_prefix\"_\"$obam\"-bga.bedGraph | cut -f1-4,8 > \"$outdir\"/$ref_prefix{}/\"$user_prefix\"_\"$obam\"_{}.covBga.txt"
		set +x
                #echo ${arr[n]} ${arr[n+1]} ${arr[n+2]} | sed 's/ /\n/g' | parallel "bedtools intersect -wb -a \"$outdir\"/chr\"$i\"/\"$i\"regions4.txt -b \"$outdir\"/\"$prefix\"_{}-bga.bedGraph | cut -f1-4,8 > \"$outdir\"/chr\"$i\"/\"$i\".\"$prefix\"_{}.covBga.txt"
	done
done

if [ $stop4 -eq 1 ]; then
	if [ $quiet -eq 0 ]; then echo "[STATUS] Exiting after step 4."; fi
	exit 1
fi

### Step 5: get big inputs
if [ $quiet -eq 0 ]; then echo "[STATUS] Getting coordinates of candidate regions..."; fi
candidate_regions_file="$outdir"/"$prefix"_bigInputs.txt
if [ $overwrite -eq 1 -o ! -s "$candidate_regions_file" ]; then
	echo $chrList | sed 's/ /\n/g' > "$outdir"/chrList.txt
	if [ $verbose -eq 1 ]; then set -x; fi
	Rscript --vanilla "$pdir"/bigInputs.R "$outdir" "chrList.txt" "$prefix" "$candidate_regions_file" "$ref_prefix"
	set +x
	if [ ! -s "$candidate_regions_file" ]; then
		echo "[ERROR] Failed to create \"$candidate_regions_file\" with candidate regions for RCL."
		exit 1
	fi
else
	if [ $quiet -eq 0 ]; then echo "[STATUS] Using existing candidate regions file \"$candidate_regions_file\""; fi
fi

# Copy final files to place consistent with run_rcl.bash. 
# Force overwrite because these files are not linked to this run and can be overwritten any time.
if [ $quiet -eq 0 ]; then echo "[STATUS] RCL input files being copied to name that run_rcl.bash expects."; fi
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
	echo -e "[STATUS] Data preprocessing has completed succesfully.\n\n"
fi
echo -e "RCL input data are in files:"
for rcl in "${rcl_input_files[@]}"; do
	echo -e "\t$rcl"
done
echo -e "Candidate regions are in file:\n\t$candidate_regions_file"

### Step 6: clean up temporary files
if [ $save -eq 0 ]; then
	if [ $quiet -eq 0 ]; then echo "[STATUS] Cleaning up temporary files..."; fi
	if [ $verbose -eq 1 ]; then set -x; fi
	rm "$outdir"/"$prefix"*-bga.bedGraph
	rm "$outdir"/$ref_prefix*/"$prefix"*-bga.bedGraph
	rm "$outdir"/$ref_prefix*/"$prefix"*-regions[1234].txt
	set +x
fi
