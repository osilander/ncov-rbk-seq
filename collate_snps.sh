#!/bin/bash

# this is the script to collate and calculate fractions of correctly called SNPs
echo "fraction" > all_collated.txt
for P in 40 80 120 160 200 300 400 500 600 800 1000 1200 1400 1600 1800 2000 2400 2800 3200 3600 4000 5000 6000 7000 8000 9000 10000
do
    # first we add to the column indicating the number of "contaminating" reads
    echo -n $P >> all_collated.txt
    echo -n -e "\t" >> all_collated.txt
    # now we loop over all BC files having contaminating reads
    for F in 1 2 4 5
    do
    # and loop again
      for S in 1 2 4 5
        do
            if [ $F != $S ]
            then
                # we cut out *only* the 2nd column from the *true* pass.txt and get rid of the header bit
                cut -f2 2020-05-05_2246_barcode0${S}.pass.txt | tail -n +2 | uniq > tempf.cut.txt
                # this lets us know the total number of SNPs called
                D=($(wc -l tempf.cut.txt))
                # we do the same for the contaminating file
                cut -f2 2020-05-05_2246_barcode0${F}_0${S}_${P}.pass.txt | tail -n +2 | uniq > tempfs.cut.txt
                ## IGNORE grep -f tempf.cut.txt tempfs.cut.txt | wc -l | tr '\n' '/' >> all_collated.txt

                # here we use grep to check whether the locations of the SNPs found in the 
                # contaminated file are the same as the true positive SNPs
                T=($(grep -f tempf.cut.txt tempfs.cut.txt | wc -l))
                ## IGNORE echo -n -e $D >> all_collated.txt
                # now we just get the fraction of TP SNPs called and output to the file.
                R=($(echo "scale=2 ; $T / $D" | bc))
                echo -n -e $R >> all_collated.txt
                # and add a tab for reading into R
                echo -n -e "\t" >> all_collated.txt
            fi
        done
    done
    echo "" >> all_collated.txt
done