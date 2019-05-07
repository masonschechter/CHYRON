#module load anaconda/2-2.3.0
## You need regex library in your python distribution
#source activate my_root
module load enthought_python/7.3.2
# module load python/2.7.15

forward=$1
reverse=$2
barcodes=$3
refs=$4
exp=$5
out=$6

umi="False"
umilength=20

outdir="$out/results/"
sep_dir="$outdir""separated_by_barcode"
comb_dir="$outdir""combined_by_pear"
align_dir="$outdir""aligned_to_ref"
merg_dir="$outdir""merged"
temp_dir="$outdir""temp_for_umi"
ins_dir="$outdir""insertions"

mkdir "$outdir"
mkdir "$sep_dir"  
mkdir "$comb_dir" 
mkdir "$align_dir"  
mkdir "$merg_dir" 
mkdir "$ins_dir"
mkdir "$temp_dir"

# # Barcode.txt File Format
# # "Sample name" "forward barcode" "reverse barcode" "20bp around cutsite" "Exp"

# ############################################################
#Separate
if [ "$umi" = "False" ]
then
python separate_by_barcodes_modified.py $forward $reverse fastq $barcodes $sep_dir
else
python separate_by_barcodes_modified_UMIversion.py $forward $reverse fastq $barcodes $sep_dir $umilength
fi
echo "separation by barcodes was done"

# # ############################################################
# # Combine by Pear
while read c1 c2 c3 c4 c5; do
../bin/pear-0.9.10-bin-64 -f $sep_dir/"$c1"_forward_*.fastq -r $sep_dir/"$c1"_reverse_*.fastq -o $comb_dir/"$c1".fastq -j 16 ; echo $c1 "is done"
done < $barcodes
echo "combining by PEAR was done"

# ############################################################
# Map and get stats
while read c1 c2; do
# mkdir $outdir$c1
echo $c2 | tr /a-z/ /A-Z/ > $outdir$exp"_ref.txt"
# grep "$c1" $barcodes >> $outdir$exp"_barcodes.txt"
# done < $refs
if [ "$umi" = "False" ]
then
# while read f1 f2 f3 f4 f5; do
# python multithread_maptogetDXMIseq.py $outdir$exp $barcodes $sep_dir $align_dir ##skips PEAR
python multithread_maptogetDXMIseq.py $outdir$exp $barcodes $comb_dir $align_dir
# echo $f1 $f2 $f3 $f4 $f5
# done < $barcodes
else
module load mafft/7.305
#Initial mapping - Not that efficient
# python map_togetDXMIseq_umi.py $outdir$exp $f1 $comb_dir $temp_dir
python multithread_maptogetDXMI_umi.py $outdir$exp $barcodes $comb_dir $temp_dir
python multithread_merge_umi_barcodes.py $barcodes $temp_dir

while read f1 f2 f3 f4 f5; do
#Extract similar barcodes
# python merge-umi-barcodes.py $f1 $temp_dir $temp_dir
#Multi Align and Merge
python multi_align_by_umibarcode.py $temp_dir/$f1".new.aligned.txt" $temp_dir/$f1".newi.barcodes.txt" $temp_dir/$f1".txt" $temp_dir/$f1".seq" "$exp" "$f1" $outdir
python add_mapped_info.py $temp_dir/$f1".seq" $temp_dir/$f1".aligned.txt" $align_dir/$f1".aligned.txt"
# mv temp_for_umi/* $outdir"temp_for_umi/"
done < $barcodes
fi

##Final stats
module unload enthought_python/7.3.2
#module load anaconda/2.7-4.3.1
module load python/2.7.15
#source activate my_root
## You need regex library in your python distribution
python statisticsfrommapped.py $barcodes $align_dir $outdir$exp"_ref.txt" $outdir$exp $ins_dir
done < $refs
###########################################################
mv insertions $outdir"insertions"