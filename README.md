# random-sign-heisenberg-chain

Julia scripts used to generate the data of the paper ["Heisenberg spin chain with random-sign couplings"](https://arxiv.org/abs/2312.13452), together with DMRG data repported in the same paper.

`sample_spin_wave_script.jl` can be run directly.

The DMRG scripts are meant to be used as follows. The scripts read the simulation parameters from the `file input.json` in the current directory (a sample input file is included).
First run `write_random_cj.jl` to save a file containing a random c_j configuration. Then run one of the DMRG scripts to compute the ground state property of a spin chain with that disorder realization.

DMRG data is stored in HDF5 files. A sample script showing how to load the data into Pandas Dataframe is provided in the `DMRG_data` folder.
