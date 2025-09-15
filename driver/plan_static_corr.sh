#!/bin/bash

PROG=../bin/qsi_static
DIR="$HOME/out/static_20"

mkdir -p $DIR


POSITIONAL_ARGS=()
TEMPFILE="temps"


while [[ $# -gt 0 ]]; do
  case $1 in
    -f|--force)
      FORCEMKDIR="-p"
      shift # past argument
      ;;
    -t|--temps)
      TEMPFILE="$2"
      shift
      shift
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters


INFILE=../input/SSF/medium.toml

PLANFILE=static_corr.plan

# Temperature points from file
echo "Importing temperature points from file '${TEMPFILE}'"
read -d '' -r -a temps < "${TEMPFILE}"


for i in `seq 0 18`;
do
    mkdir $FORCEMKDIR "$DIR/$i"
    for j in `seq 17 32`;
    do
        echo "$PROG $INFILE \"$DIR/$i\" $((64*$j+$i)) --temperature=${temps[$i]} >\"$DIR/$i/$j.list\" 2>\"$DIR/$i/$j.track\"" >> $PLANFILE
    done
done
