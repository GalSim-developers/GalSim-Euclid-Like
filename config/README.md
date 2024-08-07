## Config imsim

`was.yaml`: _euclidlike_imsim exemple_ config file.
You will want to update the following entries:

- `input.obseq_data.file_name`: path to the observing sequence. Link to `euclidlike/data/euclid_obseq.fits`
- `input.sky_catalog.file_name`: path to the skyCatalog to use
- `output.dir`: path to the output directory for the simulated images
- `output.truth.dir`: path to the ouput directory for the true catalogs

To run the code:
```bash
galsim was.yaml
```

You might want to specify some config entries in the command line, like:
```bash
galsim was.yaml input.obseq_data.visit=33690 image.CCD=1
```

## Config SLURM

`slurm_runner.sh` contain the SLURM configuration to run "large scale" simulations.
You will want to update the following line:

- `#SBATCH --output=/path/to/slurm-%A-%a.out`: slurm stdout file
- `#SBATCH --error=/path/to/slurm-%A-%a.err`: slurm stderr file
- `source activate [env_name]`: conda environment to use
- `file_list='/path/to/run_list.txt'`: file containing the pointings to simulate (see note below)

the `run_list.txt` is a 2 columns file with the pointing and the CCD_ID to simulate. It should lookq like:
```
33688 0
33688 1
33688 2
[...]
33688 35
33689 0
33689 1
[...]
```
