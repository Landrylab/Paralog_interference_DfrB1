#!/bin/bash

# This script will run Rosetta DDG prediction on SLURM servers. 
# Flexddg compute ddG of binding
# Your input directory should have :
#     1) PDB file
#     2) Mutations list file (or positions list + residue types for saturation mutagenesis)
#     3) Params files of ligand in your structure (Usually generate with mol_to_params.py script)
#     4) Configuration run file (default flexddg_talaris2014)
#     5) Configuration settings file (default rosettampi)
# Output will be in the flexddg folder created on in your running directory

# Define paths to rosetta intallation directory and the required virtual environment
ROSETTADDGPRED_VENV=/path/to/virtualenv
ROSETTA_INSTALL_DIR=/path/Rosetta/install/

# Function to display usage instructions
usage() {
  echo "Usage: $(basename $0) [-h] [-i] [-o] -p PDBFILE [--config-run CONFIG_RUN_FILE] [--config-set CONFIGFILE_SETTINGS] [--config-agg CONFIG_AGG_FILE] [-l LIST_FILE] [--saturation] [--reslistfile RESLIST_FILE] [--extra-res-path EXTRA_RES_PATH] [--nstruct NSTRUCT] [--auto_setup_metals True] [--extract-pdb True] [-n NUM_PROCESSES] [--partition ] [--mem 1G] [--time 1-00:00:00]"
  echo "Options:"
  echo "  -i, --in                 Input directory (default: current directory)"
  echo "  -o, --out                Output directory (default: current directory)"
  echo "  -p, --pdbfile            PDB file"
  echo "  -l, --listfile           File with list of mutations (or position for saturation mutagenesis) (default: mutation_list.txt)"
  echo "      --saturation         Perform saturation mutagenesis on selected positions (no argument required)"
  echo "      --reslistfile        File with list of residue types to perform saturation mutagenesis"
  echo "      --config-run         Configuration run file (default: flexddg_talaris2014.yaml)"
  echo "      --config-set         Configuration settings file (default: rosettampi.yaml)"  
  echo "      --config-agg         Configuration file for aggregation (default: aggregate.yaml)"
  echo "      --extra-res-path     Path to extra residues params files (default: input directory)"
  echo "      --nstruct            Number of structures to be generated for each mutation (default: 35)"
  echo "      --auto_setup_metals  Whether rosetta should setup metals ions automatically (True or False) (default: True)"
  echo "      --extract-pdb        Whether to extract the PDB structures from the database file (True or False) (default: True)"
  echo "  -n, --np                 SLURM - Number of processes (default: 1)"
  echo "      --partition          SLURM - Disk partition to use (default:)"
  echo "      --mem                SLURM - Memory to use (default: 1G per processes)"
  echo "      --time               SLURM - Time limit for computation (default: 1-00:00:00 (1 Day))"
  echo "  -s, --batch_folder       A subdirectory to keep all the slurm batch files"
  echo "  -h, --help               Display this help message"
  exit 1
}

# Parse command-line arguments using getopt
VALID_ARGS=$(getopt -o i:o:p:l:n:s:h --long in:,out:,pdbfile:,listfile:,saturation,reslistfile:,config-run:,config-set:,config-agg:,extra-res-path:,nstruct:,auto_setup_metals:,extract-pdb:,np:,partition:,mem:,time:,batch_folder:,help -- "$@")
if [[ $? != 0 ]]; then
    usage
fi

eval set -- "$VALID_ARGS"

# Default values for variables
INPUT_DIR=$PWD
OUTPUT_DIR=$PWD
LIST_FILE=mutations_list.txt
SATURATION=false
RESLIST_FILE=reslist.txt
CONFIG_RUN_FILE=flexddg_talaris2014.yaml
CONFIG_SETTINGS_FILE=rosettampi.yaml
CONFIG_AGG_FILE=aggregate.yaml
NUM_PROCESSES=1
EXTRA_RES_PATH=$INPUT_DIR/params
AUTO_SETUP_METALS=True
EXTRACT_PDB=True
NSTRUCT=35
TIME=1-04:00:00

# Process command-line options
while true; do
    case "$1" in
        -i | --in ) INPUT_DIR="$2"; shift 2;;
        -o | --out ) OUTPUT_DIR="$2"; shift 2;;
        -p | --pdbfile ) PDB_FILE="$2"; shift 2;;
        -l | --listfile ) LIST_FILE="$2"; shift 2;;
	-s | --batch_folder ) BATCH_FOLDER="$2"; shift 2;;
        --saturation ) SATURATION=true; shift;;
        --reslistfile ) RESLIST_FILE="$2"; shift 2;;
        --config-run ) CONFIG_RUN_FILE="$2"; shift 2;;
        --config-set ) CONFIG_SETTINGS_FILE="$2"; shift 2;;
	--config-agg ) CONFIG_AGG_FILE="$2"; shift 2;;
        --extra-res-path ) EXTRA_RES_PATH="$2"; shift 2;;
        --nstruct ) NSTRUCT="$2"; shift 2;;
        --auto_setup_metals ) AUTO_SETUP_METALS="$2"; shift 2;;
        --extract-pdb ) EXTRACT_PDB="$2"; shift 2;;
        -n | --np ) NUM_PROCESSES="$2"; shift 2;;
        --partition ) PARTITION="$2"; shift 2;;
        --mem ) MEM="$2"; shift 2;;
        --time ) TIME="$2"; shift 2;;
        -h | --help ) usage; exit 0;;
        -- ) shift; break;;
        * ) usage; exit 1;;
    esac
done

MEM=$((1 * $NUM_PROCESSES))G

# Extract the protein name from the PDB file name and date for job name
PROT=$(basename $PDB_FILE .pdb)
DATE=$(date +'%d%m%Y_%H-%M')

echo Checking output subdirectories

mkdir -p ${BATCH_FOLDER}

for subdirectory in "$OUTPUT_DIR"/*/; do
    if [ -d "$subdirectory" ]; then
    
	# Create a Slurm job script for Rosetta DDG prediction
        cat > ${BATCH_FOLDER}/flexddg_${PROT}_$(basename "${INPUT_DIR}").sbatch << EOF
#!/bin/bash

#SBATCH -D $subdirectory
#SBATCH -J flexddg_$(basename "$subdirectory")_${PROT}
#SBATCH -o flexddg_${PROT}_${DATE}.out
#SBATCH --account=def-<account>
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=$NUM_PROCESSES
#SBATCH --mem=$MEM
#SBATCH --time=$TIME
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<email>

## <load virtualenv>

cd $PWD

# Use yq to set the option in your YAML file
#yq e -i '.steps.flexddg.options."-in:file:extra_res_path" = "'"$EXTRA_RES_PATH"'"' $INPUT_DIR/$CONFIG_RUN_FILE
#yq e -i '.steps.flexddg.extract_structures.options."-in:file:extra_res_path" = "'"$EXTRA_RES_PATH"'"' $INPUT_DIR/$CONFIG_RUN_FILE
#yq e -i '.steps.flexddg.options."-in:auto_setup_metals" = $AUTO_SETUP_METALS' $INPUT_DIR/$CONFIG_RUN_FILE
#yq e -i '.steps.flexddg.extract_structures.options."-in:auto_setup_metals" = $AUTO_SETUP_METALS' $INPUT_DIR/$CONFIG_RUN_FILE
#yq e -i '.steps.flexddg.extract_structures.options."-corrections:restore_talaris_behavior" = True' $INPUT_DIR/$CONFIG_RUN_FILE
#yq e -i '.steps.flexddg.extract_structures.extract = $EXTRACT_PDB' $INPUT_DIR/$CONFIG_RUN_FILE
#yq e -i ".mutations.nstruct = $NSTRUCT" $INPUT_DIR/$CONFIG_RUN_FILE

# Run rosetta DDG prediction
if [ $SATURATION == false ]; then
    rosetta_ddg_run -d $subdirectory \
        -p $INPUT_DIR/$PDB_FILE \
        -cr $INPUT_DIR/$CONFIG_RUN_FILE \
        -cs $INPUT_DIR/$CONFIG_SETTINGS_FILE \
        -r $ROSETTA_INSTALL_DIR \
        -l $INPUT_DIR/$LIST_FILE \
        -n $NUM_PROCESSES
elif [ $SATURATION == true ]; then
    rosetta_ddg_run -d $subdirectory \
        -p $INPUT_DIR/$PDB_FILE \
        -cr $INPUT_DIR/$CONFIG_RUN_FILE \
        -cs $INPUT_DIR/$CONFIG_SETTINGS_FILE \
        -r $ROSETTA_INSTALL_DIR \
        -l $subdirectory/$LIST_FILE \
        --saturation \
        --reslistfile $INPUT_DIR/$RESLIST_FILE \
        -n $NUM_PROCESSES
fi

wait

rosetta_ddg_aggregate -cr $INPUT_DIR/$CONFIG_RUN_FILE \
    -cs $INPUT_DIR/$CONFIG_SETTINGS_FILE \
    -ca $INPUT_DIR/$CONFIG_AGG_FILE \
    -mf $subdirectory/flexddg/mutinfo.txt \
    -n $NUM_PROCESSES \
    -d $subdirectory \
    -od $subdirectory

EOF
    fi
done


