# Installation
The docker of MetaGP is available in, `https://hub.docker.com/r/aritramhp/metagp`
Use `docker pull aritramhp/metagp:<tag>` to download and install the image

# Execution of the pipeline:
All the functions are implemented in the `metagp_core.py` file. 
Each function performs a specific task present in the docker. 
For easy to execute, `wrapper_metagp.py` file has been provided. 
> `python wrapper_metagp.py -i <input_directory> [--pre|--qc|--taxo|--func]`

Use `python wrapper_metagp.py -h` for the list of all arguments.
