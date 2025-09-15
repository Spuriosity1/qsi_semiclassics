#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Check if the user provided an input file
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <command_file> <cpus_per_task>"
    exit 1
elif [ "$#" -eq 1 ]; then
    CPUS_PER_TASK=1
else
    CPUS_PER_TASK="$2"
fi

COMMAND_FILE="$1"

# Check if the file exists
if [ ! -f "$COMMAND_FILE" ]; then
    echo "Error: File '$COMMAND_FILE' not found."
    exit 1
fi


# Determine the template file based on the first line of the command file
FIRST_LINE=$(head -n 1 "$COMMAND_FILE")
if [[ "$FIRST_LINE" == python* ]]; then
    TEMPLATE_FILE="julia.sh"
elif [[ "$FIRST_LINE" == python* ]]; then
    TEMPLATE_FILE="python.sh"
else
    TEMPLATE_FILE="default.sh"
fi

HEADER_FILE="$SCRIPT_DIR/slurm_templates/header.slurm"
TEMPLATE_FILE="$SCRIPT_DIR/slurm_templates/$TEMPLATE_FILE"

if [ ! -f "$TEMPLATE_FILE" ]; then
    echo "Error: File '$TEMPLATE_FILE' not found."
    exit 1
fi

STEM="${COMMAND_FILE%.plan}"


# Get the number of lines (commands) in the file
NUM_LINES=$(wc -l < "$COMMAND_FILE" | xargs)

# Create a temporary job script
# JOB_SCRIPT=$(mktemp)
mkdir -p tmp
mkdir -p track
JOB_SCRIPT="tmp/`date +%s`.slurm"

# Copy template content to job script
cp "$HEADER_FILE" "$JOB_SCRIPT"

echo "#!
#SBATCH -J $STEM
#SBATCH --output=track/$STEM.%A_%a.out
#SBATCH --error=track/$STEM.%A_%a.error
#SBATCH --cpus-per-task=$CPUS_PER_TASK
#!
#SBATCH --array=1-$NUM_LINES 
" >> $JOB_SCRIPT

cat "$TEMPLATE_FILE" >> $JOB_SCRIPT

echo "
workdir="$SCRIPT_DIR"

cd \$workdir
echo -e \"Changed directory to \`pwd\`.\n\"


# Load the command from the file
COMMAND=\$(sed \"\${SLURM_ARRAY_TASK_ID}q;d\" \"$COMMAND_FILE\")

echo \"Running: \$COMMAND\"
eval \"\$COMMAND\"
" >> "$JOB_SCRIPT"

# Submit the array job

echo Wrote job script $JOB_SCRIPT
echo "Run with"
# echo sbatch --array=1-$NUM_LINES  "$JOB_SCRIPT"
echo sbatch "$JOB_SCRIPT"


