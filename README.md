# myBrain-Seq [![license](https://img.shields.io/badge/license-MIT-brightgreen)](https://github.com/sing-group/my-brain-seq) [![dockerhub](https://img.shields.io/badge/hub-docker-blue)](https://hub.docker.com/r/sing-group/my-brain-seq) [![compihub](https://img.shields.io/badge/hub-compi-blue)](https://www.sing-group.org/compihub/explore/625e719acc1507001943ab7f)
> **myBrain-Seq** is a [Compi](https://www.sing-group.org/compi/) pipeline for miRNA-Seq data analysis in neuropsychiatry . A Docker image is available for this pipeline in [this Docker Hub repository](https://hub.docker.com/r/sing-group/my-brain-seq).

## IPSSA repositories

- [GitHub](https://github.com/sing-group/my-brain-seq)
- [Docker Hub](https://hub.docker.com/r/sing-group/my-brain-seq)
- [CompiHub](https://www.sing-group.org/compihub/explore/625e719acc1507001943ab7f)

# What does IPSSA do?

**IPSSA** (Integrated Positively Selected Sites Analyses) is a [Compi](https://www.sing-group.org/compi/) pipeline to automatically identify positively selected amino acid sites using three different methods, namely CodeML, omegaMap, and FUBAR. Moreover, it looks for evidence of recombination in the data.
 
IPSSA applies the same steps to each input FASTA file separately. This process comprises:
1. Checking if the input FASTA file contains ambiguous nucleotide positions or non-multiple of three sequences. If so, the pipeline stops at this point and the files must be fixed.
2. Extract a random subset of sequences according to the sequence limit specified to create the master set of sequences.
3. Translate and align the master set of sequences.
4. The master protein alignment is then backtranslated to produce a master DNA alignment, used to:
    1. Run phipack.
    2. Create the PSS subsets for CodeML, omegaMap, and FUBAR, according to the number of sequences and replicas specified for each method.
5. The master protein alignment is also filtered to remove low confidence positions, according to the value specified. These filtered files are then converted into DNA files, which are split into the same subsets used by CodeML, omegaMap, and FUBAR. These are the files used by MrBayes to produce a Bayesian phylogenetic tree that is used by FUBAR and CodeML.
6. Run phipack for each one of the PSS subsets.
7. Run MrBayes for each one of the PSS subsets (using the filtered DNA files produced in step 5).
8. Run CodeML, omegaMap, and FUBAR, using their corresponding PSS subsets.
9. Finally, gather the results of all PSS methods into a tabular format.

# Using the IPSSA image in Linux
In order to use the IPSSA image, create first a directory in your local file system (`ipssa_project` in the example) with the following structure: 

```bash
ipssa_project/
├── input
│   ├── 1.fasta
│   ├── 2.fasta
│   ├── .
│   ├── .
│   ├── .
│   └── n.fasta
└── ipssa-project.params
```

Where:
- The input FASTA files to be analized must be placed in the `ipssa_project/input` directory.
- If neccessary, the Compi parameters file is located at `ipssa_project/ipssa-project.params`.

## Running IPSSA without a parameters file

Once this structure and files are ready, you should run and adapt the following commands to run the entire pipeline using the default parameters (i.e. without a Compi parameters file). Here, you only need to set `PROJECT_DIR` to the right path in your local file system and `COMPI_NUM_TASKS` to the maximum number of parallel tasks that can be run. Pipeline parameters can be added at the end of the pipeline (e.g. `--sequence_limit 10`). Note that the `--host_working_dir` is mandatory and must point to the pipeline working directory in the host machine.

```bash
PROJECT_DIR=/path/to/ipssa_project
COMPI_NUM_TASKS=6

PIPELINE_WORKING_DIR=${PROJECT_DIR}/pipeline_working_dir
INPUT_DIR=${PROJECT_DIR}/input

# Run with default parameter values
docker run -v /tmp:/tmp -v /var/run/docker.sock:/var/run/docker.sock -v ${PIPELINE_WORKING_DIR}:/working_dir -v ${INPUT_DIR}:/input --rm pegi3s/ipssa -o --logs /working_dir/logs --num-tasks ${COMPI_NUM_TASKS} -- --host_working_dir ${PIPELINE_WORKING_DIR}
```

## Running IPSSA with a parameters file

If you want to specify the pipeline parameters using a Compi parameters file, you should run and adapt the following commands. These are the same commands as above but with the addition of the `PARAMS_DIR` variable.

An example of a Compi parameters file can be obtained running the following command: `docker run --rm --entrypoint cat pegi3s/ipssa /resources/ipssa-project.params`.

This parameters file contains the default values recommended for running IPSSA. Please, note that you must update the value of the `host_working_dir` parameter in this file before using it.

```bash
PROJECT_DIR=/path/to/ipssa_project
COMPI_NUM_TASKS=6

PIPELINE_WORKING_DIR=${PROJECT_DIR}/pipeline_working_dir
INPUT_DIR=${PROJECT_DIR}/input
PARAMS_DIR=${PROJECT_DIR}

docker run -v /tmp:/tmp -v /var/run/docker.sock:/var/run/docker.sock -v ${PIPELINE_WORKING_DIR}:/working_dir -v ${INPUT_DIR}:/input -v ${PARAMS_DIR}:/params --rm pegi3s/ipssa -o --logs /working_dir/logs --num-tasks ${COMPI_NUM_TASKS} -pa /params/ipssa-project.params
```

## Extra

### Find out tasks with errors

Some tasks may produce errors that do not cause the pipeline to fail, but they can be important. Such errors are reported in the log files produced in the `logs` directory of the pipeline working directory. The `find-error-tasks.sh` script of the `pegi3s/ipssa` Docker image displays the errored tasks (i.e. those containing the word error in their log files) along with the names of the corresponding input files. Run the following command to find them (assuming the environment variables with the project and working directory paths have been declared):

```bash
docker run --entrypoint /opt/scripts/find-error-tasks.sh -v ${PIPELINE_WORKING_DIR}:/working_dir -v ${INPUT_DIR}:/input --rm pegi3s/ipssa /working_dir/logs /input /working_dir/run_lists
```

### Re-running the pipeline

To re-run the pipeline in the same project directory, run the following command first in order to clean the pipeline working directory:

```bash
sudo rm -rf ${PIPELINE_WORKING_DIR}
``` 

# Pipeline parameters

These are the pipeline parameters:
		
- `sequence_limit`: The maximum number of sequences to use for the master file. The default value is `90`.
- `random_seed`: The random seed.
- `align_method`: The alignment method to use, one of: `clustalw`, `muscle`, `kalign`, `t_coffee`, or `amap`. The default value is `muscle`.
- `tcoffee_min_score`: The minimum support value for alignment positions. The default value is `3`.
- `mrbayes_generations`: The number of iterations in MrBayes. The default value is `1000000`.
- `mrbayes_burnin`: The MrBayes burnin. The default value is `2500`.
- `fubar_sequence_limit`: The maximum number of sequences to be used by FUBAR. The default value is `90`.
- `fubar_runs`: The number of independent replicas for FUBAR. The default value is `1`.
- `codeml_sequence_limit`: The maximum number of sequences to be used by CodeML. The default value is `30`.
- `codeml_runs`: The number of independent replicas for CodeML. The default value is `1`.
- `codeml_models`: The CodeML models to be run, one or more of: '1', '2', '7', and/or '8'. To declare more than one model use a blank space between models. The default value is `1 2 7 8`.
- `omegamap_sequence_limit`: The maximum number of sequences to be used by omegaMap. The default value is `90`.
- `omegamap_iterations`: The number of omegaMap iterations. the default value is `1`.
- `omegamap_runs`: The number of independent replicas for omegaMap. The default value is `2500`.
- `omegamap_recomb`: A flag indicating if omegaMap must be executed only if recombination is detected in the master file. By default, the flag is not present and thus omegaMap is executed (if `omegamap_iterations` > 0).

# Test data

The sample data is available [here](https://github.com/pegi3s/ipssa/raw/master/resources/test-data/ipssa-m-leprae.zip). Download and uncompress it, and move to the directory named `ipssa-m-leprae`, where you will find:

- A directory called `ipssa-m-leprae-project`, that contains the structure described previously.
- A file called `run.sh`, that contains the following commands (where you should adapt the `PROJECT_DIR` path) to test the pipeline:

```bash
PROJECT_DIR=/path/to/ipssa-m-leprae-project
COMPI_NUM_TASKS=8

PIPELINE_WORKING_DIR=${PROJECT_DIR}/pipeline_working_dir
INPUT_DIR=${PROJECT_DIR}/input
PARAMS_DIR=${PROJECT_DIR}

docker run -v /tmp:/tmp -v /var/run/docker.sock:/var/run/docker.sock -v ${PIPELINE_WORKING_DIR}:/working_dir -v ${INPUT_DIR}:/input -v ${PARAMS_DIR}:/params --rm pegi3s/ipssa -o --logs /working_dir/logs --num-tasks ${COMPI_NUM_TASKS} -pa /params/ipssa-project.params
```

## Running times

- ≈ 207 minutes - 50 parallel tasks - Ubuntu 18.04.2 LTS, 96 CPUs (AMD EPYC™ 7401 @ 2GHz), 1TB of RAM and SSD disk.
- ≈ 345 minutes - 16 parallel tasks - Ubuntu 18.04.3 LTS, 12 CPUs (AMD Ryzen 5 2600 @ 3.40GHz), 16GB of RAM and SSD disk.

# For Developers

## Building the Docker image

To build the Docker image, [`compi-dk`](https://www.sing-group.org/compi/#downloads) is required. Once you have it installed, simply run `compi-dk build` from the project directory to build the Docker image. The image will be created with the name specified in the `compi.project` file (i.e. `pegi3s/ipssa:latest`). This file also specifies the version of compi that goes into the Docker image.
