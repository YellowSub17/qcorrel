#!/bin/bash


prots=("1al1" "1cos" "1mft")
rmaxs=("120" "110" "100" "90" "80" "70 ""60" "50" "40" "30")
plot_types=("r1r2" "sumax0" "sumax2")



for prot in ${prots[@]}; do
    #mkdir ${prot}

    for rmax in ${rmaxs[@]}; do

        #mkdir ${prot}/${rmax}_rmax
        for plot_type in ${plot_types[@]}; do
#            mkdir ${prot}/${rmax}_rmax/${plot_type} 

            mv ${prot}*rmax${rmax}*${plot_type}.png ${prot}/${rmax}_rmax/${plot_type}
        done        

    done
done
