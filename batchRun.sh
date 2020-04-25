#!/bin/bash
# Execute Fit_COVID_models.m for all possible combinations of i_country and i_model

declare -i i=0
declare -i j=0

for((i=1;i<=9;i=i+1))
do
    for((j=1;j<=3;j=j+1))
    do
        echo "i_country=$i,i_model=$j"
        nohup matlab -r "i_country=$i;i_model=$j;Fit_COVID_models" > out.$i.$j.txt &
    done
done
