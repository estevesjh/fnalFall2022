#!/bin/bash
# working directory, use high performance shared area on perlmutter
# Now set up CosmoSIS
export TOP_DIR=/global/common/software/des/annis
export COSMOSIS_REPO_DIR=${TOP_DIR}/cosmosis
export CSL_DIR=${TOP_DIR}/cosmosis-standard-library

export INTEGRATION_TOOLS_DIR=${TOP_DIR}/y3_pipe_under
export Y3PIPE_DIR=${TOP_DIR}/y3_cluster_cpp
export Y3_CLUSTER_WORK_DIR=${Y3PIPE_DIR}/release-build
export Y3_CLUSTER_CPP_DIR=${Y3PIPE_DIR}
export COSMOSIS_STANDARD_LIBRARY=${CSL_DIR}

export CUBA_DIR=${INTEGRATION_TOOLS_DIR}/cuba
export CUBA_CPP_DIR=$INTEGRATION_TOOLS_DIR/cubacpp
export GPU_INT_DIR=${INTEGRATION_TOOLS_DIR}/gpuintegration 

# We set OMP_NUM_THREADS to avoid oversubscribing the CPU cores.
# This value has not been carefully tuned.
export OMP_NUM_THREADS=4

# Now set up CosmoSIS
source ${COSMOSIS_REPO_DIR}/setup-cosmosis-nersc /global/common/software/des/common/Conda_Envs/cosmosis-global
#module load python/3.9-anaconda-2021.11  # to get condaconda activate ./env
#conda activate /global/common/software/des/common/Conda_Envs/cosmosis-global
exec python -m ipykernel "$@"

