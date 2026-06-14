#!/bin/bash
#SBATCH --job-name=amrex_lite_scale
#SBATCH --output=scaling_lite_%j.out
#SBATCH --error=scaling_lite_%j.err
#SBATCH --time=00:15:00               # 15 minutes! Scheduler loves this.
#SBATCH --nodes=1                     # Only 1 node required
#SBATCH --ntasks=32                   # Max cores we will test
#SBATCH --partition=normal            

# Load required modules
module load gcc/11.3.0
module load openmpi/4.1.4

EXE="./weno_3d_amrex"

# Scale across CORES on a single node, rather than multiple nodes
CORE_COUNTS=(1 2 4 8 16 32)

# Keep the grid moderate so it fits in memory but has enough work
# N=256 in your code (256^3) is perfect for a 64-core test
AMREX_ARGS="amrex.tiny_profiler=1 amrex.v=1"

echo "=========================================================="
echo "Starting Low-Intensity Scaling Test (Single Node)"
echo "Job ID: $SLURM_JOB_ID"
echo "=========================================================="

for cores in "${CORE_COUNTS[@]}"; do
    
    echo ""
    echo "----------------------------------------------------------"
    echo "Running on $cores MPI tasks (cores)..."
    echo "----------------------------------------------------------"
    
    # Restrict the run to $cores tasks
    srun --exact --nodes=1 --ntasks=$cores --mpi=pmix_v3 \
        $EXE $AMREX_ARGS > "scaling_out_${cores}_cores.log" 2>&1
    
    echo "Done. Profiling snippet:"
    # Print just the important lines to the main output file
    grep "Evolution_Loop" "scaling_out_${cores}_cores.log"
    grep "compute_wb_flux_3D" "scaling_out_${cores}_cores.log" 
    
done

echo ""
echo "=========================================================="
echo "Finished!"
echo "=========================================================="
