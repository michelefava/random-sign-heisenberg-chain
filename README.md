# random-sign-heisenberg-chain

Samples of Julia scripts used to generate the data of the paper "Heisenberg spin chain with random-sign couplings".

sample_spin_wave_script.jl can be run directly.

The DMRG scripts are meant to be used as follows. The scripts read the simulation parameters from the file input.json in the current directory (a sample input file is included).
First run write_random_cj.jl to save a file containing a random c_j configuration. Then run one of the DMRG scripts to compute the ground state property of a spin chain with that disorder realization.
