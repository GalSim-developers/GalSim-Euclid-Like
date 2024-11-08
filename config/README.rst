Config imsim
============

``was.yaml``: *euclidlike_imsim example* config file.  
You will want to update the following entries:

- ``input.obseq_data.file_name``: path to the observing sequence. Link to ``euclidlike/data/euclid_obseq.fits``
- ``input.sky_catalog.file_name``: path to the skyCatalog to use
- ``output.dir``: path to the output directory for the simulated images
- ``output.truth.dir``: path to the output directory for the true catalogs

To run the code:

.. code-block:: bash

   galsim was.yaml

You might want to specify some config entries on the command line, like:

.. code-block:: bash

   galsim was.yaml input.obseq_data.visit=33690 image.CCD=1


Config SLURM
============

``slurm_runner.sh`` contains the SLURM configuration to run "large scale" simulations.  
You will want to update the following lines:

- ``#!/bin/zsh``: depending on the shell you are using, you might want to change it to: ``#!/bin/bash``
- ``#SBATCH --output=/path/to/slurm-%A-%a.out``: SLURM stdout file
- ``#SBATCH --error=/path/to/slurm-%A-%a.err``: SLURM stderr file
- ``source activate [env_name]``: conda environment to use
- ``file_list='/path/to/run_list.txt'``: file containing the pointings to simulate (see note below)

The ``run_list.txt`` is a 2-column file with the pointing and the ``CCD_ID`` to simulate. It should look like:

.. code-block:: text

   33688 0
   33688 1
   33688 2
   [...]
   33688 35
   33689 0
   33689 1
   [...]

