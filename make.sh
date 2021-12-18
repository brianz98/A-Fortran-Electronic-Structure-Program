#!/bin/bash

# -n: Make afresh
# -p: Production

n_flag=0
p_flag=0

while getopts 'np' flag
do
    case "${flag}" in
        n) n_flag=1 ;;
        p) p_flag=1 ;;
        *) n_flag=0; p_flag=0
            exit 1 ;;
    esac
done

if [ $p_flag -ne 1 ]
then
    cmakeopt='-DCMAKE_BUILD_TYPE=Debug'
else
    cmakeopt=''
fi

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
