#!/bin/bash

# -n: Make afresh
# -p: Production
# -g: GPU

n_flag=0
p_flag=0
g_flag=0

while getopts 'npgh' flag
do
    case "${flag}" in
        n) n_flag=1 ;;
        p) p_flag=1 ;;
        g) g_flag=1 ;;
        h) echo "-n: build anew; -p: production; -g: GPU acceration"
            exit 1 ;;
    esac
done

if [ $p_flag -ne 1 ]
then
    cmakeopt='-DCMAKE_BUILD_TYPE=Debug'
else
    cmakeopt=''
fi

if [ $g_flag -eq 1 ]
then
    export FC=nvfortran
fi

rm src/*.old

fppfiles="src/linalg.fpp"

for fppfile in ${fppfiles}
do
    fypp -m itertools ${fppfile} "${fppfile%.fpp}.f90"
done

if [ $n_flag -eq 1 ]
then
    rm -r build
    mkdir build; cd $_
    cmake .. ${cmakeopt}
    make
    cp els.x ..
    cd ..
else
    cd build
    make
    cp els.x ..
    cd ..
fi

for fppfile in ${fppfiles}
do
    mv "${fppfile%.fpp}.f90" "${fppfile%.fpp}.old"
done
