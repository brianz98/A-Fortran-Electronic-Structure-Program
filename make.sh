#!/bin/bash

n_flag=0

while getopts 'n' flag
do
    case "${flag}" in
        n) n_flag=1 ;;
        *) n_flag=0
            exit 1 ;;
    esac
done

if [ $n_flag -eq 1 ]
then
    rm -r build
    mkdir build; cd $_
    cmake .. -DCMAKE_BUILD_TYPE=Debug
    make
    cp els.x ..
    cd ..
else
    cd build
    make
    cp els.x ..
    cd ..
fi
