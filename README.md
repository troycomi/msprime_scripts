[![Build Status](https://travis-ci.com/troycomi/msprime_scripts.svg?branch=master)](https://travis-ci.com/troycomi/msprime_scripts)

# msprime_scripts
> Untangling human migration through simulated admixture

Collection of python and shell scripts for performing msprime simulations 
relevant for evaluating admixture and back migration events.

## Installation
All required packages are listed in environment.txt.  
To create a conda environment for the project on Della, first activate anaconda3:
```
module load anaconda3
```
The conda environment is then created with:
```
conda create -n msprime_scripts --file environment.txt
```

## Usage
### Environment Setup
Prior to running on Della, need to activate anaconda3:
```
module load anaconda3
```
The conda environment is then activated with:
```
conda activate msprime_scripts
```

### Normal Operation
Coming soon

### Tests
The unit tests for msprime_scripts are run with `pytest` 
from anywhere in the directory.

## License
TBD
