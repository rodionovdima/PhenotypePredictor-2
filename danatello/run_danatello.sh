#!/usr/bin/bash

# mkdir -p DIAMOND
# Format of the DIAMON table
fmt="qseqid sseqid qlen qstart qend evalue pident bitscore"
# Path to the mcSEED-DIAMOND database directory
# Directory should contain:
#    - all_no_func_representative_d.dmnd - DIAMOND sequence database
#    - all_subsystems.txt - table of subsystems genes and their annotations
db_dir="database/"
# Number of parallel processes to use. Set in a range of 1-16
ncpu=16

echo "Danatello pipeline"

# Run DIAMOND
for genome in ./faa/*.fa*
do
   echo $genome
   name=$(basename $genome)
   echo Processing $name
   if [ ! -f DIAMOND/"$name" ]
   then
       diamond blastp \
           -q $genome \
           -d "$db_dir"/all_no_func_representative_d.dmnd \
           -o DIAMOND/"$name" \
           -f 6 $fmt \
           -p $ncpu \
           -k 100
   else
      echo "File DIAMOND/"$name" exists"
   fi
done

# Annotating DIAMOND output using Danatello
mkdir -p annotation
ident=0.8
echo "Annotating subsystems set with KDE threshold selection"

python annotate.py \
    -d DIAMOND \
    -a "$db_dir"/all_subsystems.txt \
    -f $fmt \
    -p $ncpu \
    -i $ident \
    -m 1 \
    -o annotation \
    --dir \
    --kde
