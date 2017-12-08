# GROMACS Parser
[NOMAD Laboratory CoE](http://nomad-coe.eu) parser for [GROMACS](http://www.gromacs.org/)
## Version 0.0.2

This is the parser for GROMACS gmx Molecular Dynamics (mdrun) in [GROMACS](http://www.gromacs.org).
The official version lives at:

    git@gitlab.mpcdf.mpg.de:nomad-lab/parser-gromacs.git

You can browse it at:

    https://gitlab.rzg.mpg.de/nomad-lab/parser-gromacs

It relies on having the nomad-meta-info and the python-common repositories one level higher.
The simplest way to have this is to check out nomad-lab-base recursively:

    git clone --recursive git@gitlab.mpcdf.mpg.de:nomad-lab/nomad-lab-base.git

This parser will be in the directory parsers/gromacs of this repository.

## Running and Testing the Parser
### Requirements
The required python packages can be installed with (see [python-common](https://gitlab.rzg.mpg.de/nomad-lab/python-common)):

    pip install -r nomad-lab-base/python-common/requirements.txt

### Usage
GROMACS gmx mdrun log output files can be parsed with:

    python GromacsParser.py [path/toFile]

### Test Files
Example log output files of GROMACS gmx can be found in the directory test/examples.
More details about the calculations and files will be explained in this README file.


