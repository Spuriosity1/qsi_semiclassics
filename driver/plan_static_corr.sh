#!/bin/bash


PROG=../bin/qsi_static
DIR="$HOME/out/static_20"

mkdir -p $DIR

INFILE=../input/SSF/medium.toml

PLANFILE=static_corr.plan

# Temperature points from file
read -d '' -r -a temps < "temps"

echo '' > $PLANFILE

for i in `seq 0 18`;
do
    mkdir "$DIR/$i"
    for j in `seq 17 32`;
    do
        echo "$PROG $INFILE \"$DIR/$i\" $((64*$j+$i)) --temperature=${temps[$i]} >\"$DIR/$i/$j.list\" 2>\"$DIR/$i/$j.track\"" >> $PLANFILE
    done
done
