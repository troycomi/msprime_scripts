[![Build Status](https://travis-ci.com/troycomi/msprime_scripts.svg?branch=master)](https://travis-ci.com/troycomi/msprime_scripts)
[![codecov](https://codecov.io/gh/troycomi/msprime_scripts/branch/master/graph/badge.svg)](https://codecov.io/gh/troycomi/msprime_scripts)

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
### Command Line
#### Environment Setup
Prior to running on Della, you must activate anaconda3:
```
module load anaconda3
```
The conda environment is then activated with:
```
conda activate msprime_scripts
```

#### Normal Operation
With the correct environment, an admixture simulation is started simply with
```
python Admixture_Simulation.py
```
to use all default values.  Several options are available and can be displayed
with the -h flag.  Most affect parameters of the msprime simulation.  Output
options include the debug output, options listing, haplotype call file, ILS
call file, vcf and f4dstat information.  By default, the debug information
will be printed to stdout and the script will exit.  Flags can activate other
file types including the output file names.  For debug, options, haplo and ils,
excluding a file name will cause them to output to stdout.  Finally, specifying
only an output directory will generate all file types with preset names.

#### Tests
The unit tests for msprime_scripts are run with `pytest` 
from anywhere in the directory.

### Slurm Jobs
Slurm jobs are submitted using the corresponding submit script within the
slurmJobs directory.  The behavior of the submission is customized within
the submit script which performs some basic operations before submitting
the slurm array job.

### Seff Examination
The reportSeff and formatSeff scripts can be used to automatically
view seff information from the slurmout files.  From the snakemake directory
```
./reportSeff.sh slurm_out/ | python formatSeff.py
```
Depending on the size of the output, pipe through less or write to an output
file.

### Snakemake
A snakemake workflow is available for performing numerous simulations and 
analyzing the results with [S star](https://github.com/bvernot/freezing-archer)
and [match p value](https://github.com/lparsons/archaic_match).
Two configurations are provided for the Princeton clusters
della and gencomp1.  Execution is controlled through config files in the 
snakefiles directory.  The run scripts will begin a workflow execution.  The
loop\_run script gives an example of how to run several instances with different
parameters in a sequential mode.

#### Bugs in other code
Dependencies of the snakemake workflow contain bugs which should be corrected
for best performance.  

The archaic match \_\_main\_\_.py has a bug where Nonetype callsets are not
handled properly.  After the conda environments are created, find the offending
file and change line 390 to:
```
if callset is None or 'calldata/GT' not in callset:
```

Snakemake, as of version 5.4.0, has errors when attempting to restart group
jobs which are partially completed.  A key error is generated in the dag.py
file which can be corrected by changing line 818 to:
```
stop = lambda j: j.group != job.group or j not in self.needrun_jobs
```

## License
TBD
