#!/bin/bash

declare -i i=0

for((i=1;i<=51;i=i+1))
do
    echo "i_state=$i"
    nohup matlab -r "i_state=$i;State_wise_models" > out.$i.txt &
done
