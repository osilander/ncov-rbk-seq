#!/bin/bash

# this is the read number of the "contaminating" barcode
for P in 40 80 120 160 200 300 400 500 600 800 1000 1200 1400 1600 1800 2000 2400 2800 3200 3600 4000 5000 6000 7000 8000 9000 10000
do
    # we use a *total* of 20,000 reads, so that's a combination of
    # the "true" reads and the contaminant
    U=$((20000-$P))

    # this loops over the different barcodes (in this case, BC01, BC02, BC04, BC05)
    for F in 1 2 4 5
    do
        # and again we loop, and do nothing if we're "contaminating" from the same barcode
        for S in 1 2 4 5
        do
            if [ $F != $S ]
            then
                # and here we make the new fastq. It is named according the number of contaminating reads (P)
                seqtk sample 2020-05-05_2246_barcode0${F}.fastq $P > 2020-05-05_2246_barcode0${F}_0${S}_${P}.fastq
                seqtk sample 2020-05-05_2246_barcode0${S}.fastq $U >> 2020-05-05_2246_barcode0${F}_0${S}_${P}.fastq
            fi
        done
    done
done