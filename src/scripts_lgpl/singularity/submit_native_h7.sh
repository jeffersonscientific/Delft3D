#! /bin/bash
  
# Usage:
#   - Modify this script where needed (e.g. number of nodes, number of tasks per node).
#   - Execute this script from the command line of H7 using:
#     sbatch ./submit_native_h7.sh
#
# This is an h7 specific script for single or multi-node simulations
  
# Specify Slurm SBATCH directives 
#SBATCH --nodes=1              # Number of nodes.
#SBATCH --ntasks-per-node=4     # The number of tasks to be invoked on each node.
                                # For sequential runs, the number of tasks should be '1'.
                                # Note: SLURM_NTASKS is equal to "--nodes" multiplied by "--ntasks-per-node".
#SBATCH --job-name=test_model   # Specify a name for the job allocation.
#SBATCH --time=00:15:00         # Set a limit on the total run time of the job allocation.
#SBATCH --partition=4vcpu       # Request a specific partition for the resource allocation.
                                # See: https://publicwiki.deltares.nl/display/Deltareken/Compute+nodes.
##SBATCH --exclusive            # The job allocation can not share nodes with other running jobs.
                                # In many cases this option can be omitted.
##SBATCH --contiguous           # The allocated nodes must form a contiguous set, i.e. next to each other.
                                # In many cases this option can be omitted.

# Set bash options. Exit on failures (and propagate errors in pipes).
set -eo pipefail

# These variables should be modified.
DIMR_FOLDER=/p/d-hydro/development/alma8/dimrsets/2.26.04.78838
DIMR_CONFIG_FOLDER=$PWD
DIMR_FILE=dimr_config.xml

# Set MPI/OpenMP options. Uncomment to override default settings.
# Reference on intel MPI environment variables: 
# https://www.intel.com/content/www/us/en/docs/mpi-library/developer-reference-linux/2021-8/environment-variable-reference.html

# export I_MPI_DEBUG=5
# export I_MPI_FABRICS=ofi
# export FI_PROVIDER=tcp
# export I_MPI_OFI_PROVIDER=tcp
# export OMP_NUM_THREADS=1

# You shouldn't need to modify the script below this line.
CONFIG_PATH="${DIMR_CONFIG_FOLDER}/${DIMR_FILE}"

# Load modules.
module purge
module load intelmpi/2021.11.0   # Load the message-passing library for parallel simulations.

# Set the maximum stacksize to 'unlimited'. 
# In some cases the process simulating the model will run out of stack memory and crash.
ulimit -s unlimited
  
# Modify the value of the process tag in dimr_config.xml.
PROCESS_STR="$(seq -s " " 0 $((SLURM_NTASKS-1)))"
sed -i "s/\(<process.*>\)[^<>]*\(<\/process.*\)/\1$PROCESS_STR\2/" $CONFIG_PATH

# The name of the MDU file is read from the DIMR configuration file.
MDU_FILENAME=$(sed -n 's/\r//; s/<inputFile>\(.*[.]mdu\)<\/inputFile>/\1/p' $CONFIG_PATH | sed -n 's/^\s*\(.*\)\s*$/\1/p')

# Set the shared library path to the `lib` folder in the dimrset
export LD_LIBRARY_PATH=$DIMR_FOLDER/lnx64/lib:$LD_LIBRARY_PATH

# Partition data using dflowfm.
if [[ "$SLURM_NTASKS" -gt 1 ]]; then
    echo "Partitioning parallel model..."
    echo "Run dflowfm on $MDU_FILENAME with $SLURM_NTASKS partitions."
    srun -n 1 -N 1 ${DIMR_FOLDER}/lnx64/bin/dflowfm --nodisplay --autostartstop --partition:ndomains=${SLURM_NTASKS}:icgsolver=6 $MDU_FILENAME
else
    echo "Sequential model..."
fi

# Run simulation using dimr.
echo "Run Simulation..."
srun ${DIMR_FOLDER}/lnx64/bin/dimr $CONFIG_PATH
