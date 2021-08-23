# Turning PathFinder Upside-Down: Exploring FPGA Switch-Blocks by Negotiating Switch Presence


This repository holds the source code of the FPGA switch-block pattern exploration tool presented in the paper entitled "Turning PathFinder Upside-Down: Exploring FPGA
Switch-Blocks by Negotiating Switch Presence".


## Setting up VPR

Extensions to VTR version 8.0 (the latest stable release at the time of writing the paper) are distributed as a patch in [vpr/vtr8_avalanche.patch](https://github.com/EPFL-LAP/fpl21-avalanche/edit/master/vpr/vtr8_avalanche.patch).
Running [vpr/get_avalanche_vpr.sh](https://github.com/EPFL-LAP/fpl21-avalanche/edit/master/vpr/get_avalanche_vpr.sh) will download VTR 8.0 and apply the patch, after which VTR can be built as usual.
In order for the scripts to be able to actually use it, please update [setenv.py](https://github.com/EPFL-LAP/fpl21-avalanche/edit/master/setenv.py) appropriately.

## SPICE Models

Please follow the instructions from [EPFL-LAP/fpga21-scaled-tech](https://github.com/EPFL-LAP/fpga21-scaled-tech) to set up the SPICE models.

## Running the Flow

Scripts that are part of the exploration flow can be found in the [avalanche/](https://github.com/EPFL-LAP/fpl21-avalanche/edit/master/avalanche/) directory.
The results presented in the paper were obtained by the running [avalanche/explore_avalanche.py](https://github.com/EPFL-LAP/fpl21-avalanche/edit/master/avalanche/explore_avalanche.py) with the following parameters:

`python -u explore_avalanche.py --base_cost 1.0 --scaling_factor 9 --avalanche_iter 25 --adoption_threshold 1.1 --wd res_dir`.

The meaning of the parameters is defined in the script itself as well as Sections III-A and V-C of the paper.
To run the simple greedy algorithm presented in Section III, please specify `--greedy 1`. The set of circuits used for exploration can be changed by modifying lines 60 and 61 of the script and adding the .blif files of the circuits to [avalanche/benchmarks](https://github.com/EPFL-LAP/fpl21-avalanche/edit/master/avalanche/) directory.


## Contact

If you find any bugs please open an issue. For all other questions, including getting access to the development branch, please contact Stefan Nikolic (firstname dot lastname at epfl dot ch).
