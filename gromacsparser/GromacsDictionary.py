# Copyright 2016-2018 Berk Onat, Fawzi Mohamed
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import numpy as np
from .GromacsCommon import PARSERNAME, PROGRAMNAME, PARSERTAG
from nomadcore.smart_parser.SmartParserDictionary import metaNameConverter, metaNameConverter_UnderscoreSpaceDash
from nomadcore.smart_parser.SmartParserDictionary import MetaInfoMap, FileInfoMap, MapDictionary
from nomadcore.smart_parser.SmartParserDictionary import getDict_MetaStrInDict, getList_MetaStrInDict, get_unitDict
import nomadcore.md_data_access.MDDataAccess as MDData

def get_fileListDict():
    """Loads dictionary for file namelists.

    Returns:
        the list of defaults file namelists
    """

    startpage = {
        'nameTranslate'   : metaNameConverter,
        'metaHeader'      : PARSERTAG,
        'metaNameTag'     : 'inout_file',
        'metaInfoType'    : 'C',
        'activeMetaNames' : [],
        'activeSections'  : [PARSERTAG+'_section_input_output_files']
        }
    namelist = {
        'topoltpr'     : FileInfoMap(startpage, fileFormat=['.tpr', '.tpb'],
                                     fileInterface=['mdanalysis','pymolfile','mdtraj'],
                                     infoPurpose=['topology', 'unitcell'], matchNames="-s"),
        'trajtrr'      : FileInfoMap(startpage, fileFormat=['.trr', '.trj', '.xdr', '.cpt'],
                                     #fileInterface=['mdanalysis','pymolfile','mdtraj'],
                                     infoPurpose=['trajectory', 'unitcell'], matchNames="-o"),
        'traj_compxtc' : FileInfoMap(startpage, fileFormat=['.xtc', '.xdr'],
                                     #fileInterface=['mdanalysis','pymolfile','mdtraj'],
                                     infoPurpose=['trajectory', 'unitcell'], matchNames="-x"),
        #'statecpt'     : FileInfoMap(startpage, activeInfo=True, fileFormat=['.cpt'],
        #                             fileInterface=['mdanalysis','mdtraj'],
        #                             infoPurpose=['inputcoordinates'],
        #                             matchNames="-cpi"),
        #'confoutgro'   : FileInfoMap(startpage, activeInfo=True, fileFormat=['.gro', '.g96', '.pdb'],
        #                             fileInterface=['mdanalysis','mdtraj'],
        #                             infoPurpose=['inputcoordinates'],
        #                             matchNames="-c"),
        'eneredr'      : FileInfoMap(startpage, activeInfo=True, fileFormat=['.edr'],
                                     infoPurpose=['thermostats'], matchNames="-e"),
        'mdlog'        : FileInfoMap(startpage),
        'rerunxtc'     : FileInfoMap(startpage),
        'rotationlog'  : FileInfoMap(startpage),
        'rotangleslog' : FileInfoMap(startpage),
        'rotslabslog'  : FileInfoMap(startpage),
        'rottorquelog' : FileInfoMap(startpage),
        'nmmtx'        : FileInfoMap(startpage),
        'dipolendx'    : FileInfoMap(startpage),
        'membeddat'    : FileInfoMap(startpage),
        'membedtop'    : FileInfoMap(startpage),
        'membedndx'    : FileInfoMap(startpage)
        }
    return MapDictionary(namelist)

def get_nameListDict(self, deflist):
    """Loads control in data of GROMACS.

    Args:
        deflist: name list definition (cntrl/arch/parm/qmopts/grpopts/anneal).

    matchWith parameters:
        EOL = End of line
        PL  = Match with previous line
        PW  = Previous word
        NW  = Next word (until space, comma) (DEFAULT)
        UD  = until delimeter (can be any string/character)

        The text matched inside apostrophe or qoutation is extracted

    Returns:
        the list of namelists
    """
    startpage = {
        'nameTranslate'   : metaNameConverter_UnderscoreSpaceDash,
        'metaHeader'      : PARSERTAG,
        'metaNameTag'     : 'inout_control',
        'metaInfoType'    : 'C',
        'activeMetaNames' : [],
        'activeSections'  : [PARSERTAG+'_section_control_parameters'],
        'matchWith'       : 'EOL'
        }

    archlist = {
        'GROMACS version'    : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'Precision'          : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'Memory model'       : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'MPI library'        : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'OpenMP support'     : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'GPU support'        : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'OpenCL support'     : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'invsqrt routine'    : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'SIMD instructions'  : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'FFT library'        : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'RDTSCP usage'       : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'C[+][+]11 compilation'  : MetaInfoMap(startpage, defaultValue='',
                                               replaceTag='Cxx11 compilation', matchWith='EOL'),
        'TNG support'        : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'Tracing support'    : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'Built on'           : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'Built by'           : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'Build OS/arch'      : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'Build CPU vendor'   : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'Build CPU brand'    : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'Build CPU family'   : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'Build CPU features' : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'C compiler'         : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'C compiler flags'   : MetaInfoMap(startpage, defaultValue='', matchWith='EOL'),
        'C[+][+] compiler'       : MetaInfoMap(startpage, defaultValue='',
                                               replaceTag='Cxx compiler', matchWith='EOL'),
        'C[+][+] compiler flags' : MetaInfoMap(startpage, defaultValue='',
                                               replaceTag='Cxx compiler flags', matchWith='EOL'),
        'Boost version'      : MetaInfoMap(startpage, defaultValue='', matchWith='EOL')
        }

    startpage = {
        'nameTranslate'   : metaNameConverter,
        'metaHeader'      : PARSERTAG,
        'metaNameTag'     : 'inout_control',
        'metaInfoType'    : 'C',
        'activeMetaNames' : [],
        'activeSections'  : [PARSERTAG+'_section_control_parameters'],
        'matchWith'       : 'EOL'
        }

    cntrllist = {
        'integrator'         : MetaInfoMap(startpage, defaultValue=''),
        'tinit'              : MetaInfoMap(startpage, defaultValue=0),
        'dt'                 : MetaInfoMap(startpage, defaultValue=0.0001),
        'nsteps'             : MetaInfoMap(startpage, defaultValue=0),
        'init-step'          : MetaInfoMap(startpage, defaultValue=0),
        'simulation-part'    : MetaInfoMap(startpage, defaultValue=0),
        'comm-mode'          : MetaInfoMap(startpage, defaultValue=0),
        'nstcomm'            : MetaInfoMap(startpage, defaultValue=0),
        'bd-fric'            : MetaInfoMap(startpage, defaultValue=0),
        'ld-seed'            : MetaInfoMap(startpage, defaultValue=0),
        'emtol'              : MetaInfoMap(startpage, defaultValue=0),
        'emstep'             : MetaInfoMap(startpage, defaultValue=0),
        'niter'              : MetaInfoMap(startpage, defaultValue=0),
        'fcstep'             : MetaInfoMap(startpage, defaultValue=0),
        'nstcgsteep'         : MetaInfoMap(startpage, defaultValue=0),
        'nbfgscorr'          : MetaInfoMap(startpage, defaultValue=0),
        'rtpi'               : MetaInfoMap(startpage, defaultValue=0),
        'nstxout'            : MetaInfoMap(startpage, defaultValue=0),
        'nstvout'            : MetaInfoMap(startpage, defaultValue=0),
        'nstfout'            : MetaInfoMap(startpage, defaultValue=0),
        'nstlog'             : MetaInfoMap(startpage, defaultValue=0),
        'nstcalcenergy'      : MetaInfoMap(startpage, defaultValue=0),
        'nstenergy'          : MetaInfoMap(startpage, defaultValue=0),
        'nstxout-compressed' : MetaInfoMap(startpage, defaultValue=0),
        'compressed-x-precision' : MetaInfoMap(startpage, defaultValue=0),
        'cutoff-scheme'      : MetaInfoMap(startpage, defaultValue=0),
        'nstlist'            : MetaInfoMap(startpage, defaultValue=0),
        'ns-type'            : MetaInfoMap(startpage, defaultValue=0),
        'pbc'                : MetaInfoMap(startpage, defaultValue=0),
        'periodic-molecules' : MetaInfoMap(startpage, defaultValue=0),
        'verlet-buffer-tolerance': MetaInfoMap(startpage, defaultValue=0),
        'rlist'              : MetaInfoMap(startpage, defaultValue=0),
        'rlistlong'          : MetaInfoMap(startpage, defaultValue=0),
        'nstcalclr'          : MetaInfoMap(startpage, defaultValue=0),
        'coulombtype'        : MetaInfoMap(startpage, defaultValue=0),
        'coulomb-modifier'   : MetaInfoMap(startpage, defaultValue=0),
        'rcoulomb-switch'    : MetaInfoMap(startpage, defaultValue=0),
        'rcoulomb'           : MetaInfoMap(startpage, defaultValue=0),
        'epsilon-r'          : MetaInfoMap(startpage, defaultValue=0),
        'epsilon-rf'         : MetaInfoMap(startpage, defaultValue=0),
        'vdw-type'           : MetaInfoMap(startpage, defaultValue=0),
        'vdw-modifier'       : MetaInfoMap(startpage, defaultValue=0),
        'rvdw-switch'        : MetaInfoMap(startpage, defaultValue=0),
        'rvdw'               : MetaInfoMap(startpage, defaultValue=0),
        'DispCorr'           : MetaInfoMap(startpage, defaultValue=0),
        'table-extension'    : MetaInfoMap(startpage, defaultValue=0),
        'fourierspacing'     : MetaInfoMap(startpage, defaultValue=0),
        'fourier-nx'         : MetaInfoMap(startpage, defaultValue=0),
        'fourier-ny'         : MetaInfoMap(startpage, defaultValue=0),
        'fourier-nz'         : MetaInfoMap(startpage, defaultValue=0),
        'pme-order'          : MetaInfoMap(startpage, defaultValue=0),
        'ewald-rtol'         : MetaInfoMap(startpage, defaultValue=0),
        'ewald-rtol-lj'      : MetaInfoMap(startpage, defaultValue=0),
        'lj-pme-comb-rule'   : MetaInfoMap(startpage, defaultValue=0),
        'ewald-geometry'     : MetaInfoMap(startpage, defaultValue=0),
        'epsilon-surface'    : MetaInfoMap(startpage, defaultValue=0),
        'implicit-solvent'   : MetaInfoMap(startpage, defaultValue=0),
        'gb-algorithm'       : MetaInfoMap(startpage, defaultValue=0),
        'nstgbradii'         : MetaInfoMap(startpage, defaultValue=0),
        'rgbradii'           : MetaInfoMap(startpage, defaultValue=0),
        'gb-epsilon-solvent' : MetaInfoMap(startpage, defaultValue=0),
        'gb-saltconc'        : MetaInfoMap(startpage, defaultValue=0),
        'gb-obc-alpha'       : MetaInfoMap(startpage, defaultValue=0),
        'gb-obc-beta'        : MetaInfoMap(startpage, defaultValue=0),
        'gb-obc-gamma'       : MetaInfoMap(startpage, defaultValue=0),
        'gb-dielectric-offset' : MetaInfoMap(startpage, defaultValue=0),
        'sa-algorithm'       : MetaInfoMap(startpage, defaultValue=0),
        'sa-surface-tension' : MetaInfoMap(startpage, defaultValue=0),
        'tcoupl'             : MetaInfoMap(startpage, defaultValue=0),
        'nsttcouple'         : MetaInfoMap(startpage, defaultValue=0),
        'nh-chain-length'    : MetaInfoMap(startpage, defaultValue=0),
        'print-nose-hoover-chain-variables' : MetaInfoMap(startpage, defaultValue=0),
        'pcoupl'             : MetaInfoMap(startpage, defaultValue=0),
        'pcoupltype'         : MetaInfoMap(startpage, defaultValue=0),
        'nstpcouple'         : MetaInfoMap(startpage, defaultValue=0),
        'tau-p'              : MetaInfoMap(startpage, defaultValue=0),
        'compressibility\[    0\]' : MetaInfoMap(startpage, defaultValue=0,
                replaceDict = {'{' : '[', '}' : ']'},
                subFunc=lambda x: np.asarray(
                    self.literal_eval(x)) if x is not None else None),
        'compressibility\[    1\]' : MetaInfoMap(startpage, defaultValue=0,
                replaceDict = {'{' : '[', '}' : ']'},
                subFunc=lambda x: np.asarray(
                    self.literal_eval(x)) if x is not None else None),
        'compressibility\[    2\]' : MetaInfoMap(startpage, defaultValue=0,
                replaceDict = {'{' : '[', '}' : ']'},
                subFunc=lambda x: np.asarray(
                    self.literal_eval(x)) if x is not None else None),
        'ref-p\[    0\]'       : MetaInfoMap(startpage, defaultValue=0,
                replaceDict = {'{' : '[', '}' : ']'},
                subFunc=lambda x: np.asarray(
                    self.literal_eval(x)) if x is not None else None),
        'ref-p\[    1\]'       : MetaInfoMap(startpage, defaultValue=0,
                replaceDict = {'{' : '[', '}' : ']'},
                subFunc=lambda x: np.asarray(
                    self.literal_eval(x)) if x is not None else None),
        'ref-p\[    2\]'       : MetaInfoMap(startpage, defaultValue=0,
                replaceDict = {'{' : '[', '}' : ']'},
                subFunc=lambda x: np.asarray(
                    self.literal_eval(x)) if x is not None else None),
        'refcoord-scaling'   : MetaInfoMap(startpage, defaultValue=0),
        'posres-com\[0\]'      : MetaInfoMap(startpage, defaultValue=0),
        'posres-com\[1\]'      : MetaInfoMap(startpage, defaultValue=0),
        'posres-com\[2\]'      : MetaInfoMap(startpage, defaultValue=0),
        'posres-comB\[0\]'     : MetaInfoMap(startpage, defaultValue=0),
        'posres-comB\[1\]'     : MetaInfoMap(startpage, defaultValue=0),
        'posres-comB\[2\]'     : MetaInfoMap(startpage, defaultValue=0),
        'QMMM'               : MetaInfoMap(startpage, defaultValue=0),
        'QMconstraints'      : MetaInfoMap(startpage, defaultValue=0),
        'QMMMscheme'         : MetaInfoMap(startpage, defaultValue=0),
        'MMChargeScaleFactor' : MetaInfoMap(startpage, defaultValue=0),
        }

    qmlist = {
        'ngQM' : MetaInfoMap(startpage, defaultValue=0),
        'constraint-algorithm' : MetaInfoMap(startpage, defaultValue=0),
        'continuation' : MetaInfoMap(startpage, defaultValue=0),
        'Shake-SOR' : MetaInfoMap(startpage, defaultValue=0),
        'shake-tol' : MetaInfoMap(startpage, defaultValue=0),
        'lincs-order' : MetaInfoMap(startpage, defaultValue=0),
        'lincs-iter' : MetaInfoMap(startpage, defaultValue=0),
        'lincs-warnangle' : MetaInfoMap(startpage, defaultValue=0),
        'nwall' : MetaInfoMap(startpage, defaultValue=0),
        'wall-type' : MetaInfoMap(startpage, defaultValue=0),
        'wall-r-linpot' : MetaInfoMap(startpage, defaultValue=0),
        'wall-atomtype\[0\]' : MetaInfoMap(startpage, defaultValue=0),
        'wall-atomtype\[1\]' : MetaInfoMap(startpage, defaultValue=0),
        'wall-density\[0\]' : MetaInfoMap(startpage, defaultValue=0),
        'wall-density\[1\]' : MetaInfoMap(startpage, defaultValue=0),
        'wall-ewald-zfac' : MetaInfoMap(startpage, defaultValue=0),
        'pull' : MetaInfoMap(startpage, defaultValue=0),
        'rotation' : MetaInfoMap(startpage, defaultValue=0),
        'interactiveMD' : MetaInfoMap(startpage, defaultValue=0),
        'disre' : MetaInfoMap(startpage, defaultValue=0),
        'disre-weighting' : MetaInfoMap(startpage, defaultValue=0),
        'disre-mixed' : MetaInfoMap(startpage, defaultValue=0),
        'dr-fc' : MetaInfoMap(startpage, defaultValue=0),
        'dr-tau' : MetaInfoMap(startpage, defaultValue=0),
        'nstdisreout' : MetaInfoMap(startpage, defaultValue=0),
        'orire-fc' : MetaInfoMap(startpage, defaultValue=0),
        'orire-tau' : MetaInfoMap(startpage, defaultValue=0),
        'nstorireout' : MetaInfoMap(startpage, defaultValue=0),
        'free-energy' : MetaInfoMap(startpage, defaultValue=0),
        'cos-acceleration' : MetaInfoMap(startpage, defaultValue=0),
        'deform\[    0\]' : MetaInfoMap(startpage, defaultValue=0,
                replaceDict = {'{' : '[', '}' : ']'},
                subFunc=lambda x: np.asarray(
                    self.literal_eval(x)) if x is not None else None),
        'deform\[    1\]' : MetaInfoMap(startpage, defaultValue=0,
                replaceDict = {'{' : '[', '}' : ']'},
                subFunc=lambda x: np.asarray(
                    self.literal_eval(x)) if x is not None else None),
        'deform\[    2\]' : MetaInfoMap(startpage, defaultValue=0,
                replaceDict = {'{' : '[', '}' : ']'},
                subFunc=lambda x: np.asarray(
                    self.literal_eval(x)) if x is not None else None),
        'simulated-tempering' : MetaInfoMap(startpage, defaultValue=0),
        'E-x' : MetaInfoMap(startpage, defaultValue=0, nextMatch='n', concatMatch=True),
        'E-xt' : MetaInfoMap(startpage, defaultValue=0, nextMatch='n', concatMatch=True),
        'E-y' : MetaInfoMap(startpage, defaultValue=0, nextMatch='n', concatMatch=True),
        'E-yt' : MetaInfoMap(startpage, defaultValue=0, nextMatch='n', concatMatch=True),
        'E-z' : MetaInfoMap(startpage, defaultValue=0, nextMatch='n', concatMatch=True),
        'E-zt' : MetaInfoMap(startpage, defaultValue=0, nextMatch='n', concatMatch=True),
        'swapcoords' : MetaInfoMap(startpage, defaultValue=0),
        'adress' : MetaInfoMap(startpage, defaultValue=0),
        'userint1' : MetaInfoMap(startpage, defaultValue=0),
        'userint2' : MetaInfoMap(startpage, defaultValue=0),
        'userint3' : MetaInfoMap(startpage, defaultValue=0),
        'userint4' : MetaInfoMap(startpage, defaultValue=0),
        'userreal1' : MetaInfoMap(startpage, defaultValue=0),
        'userreal2' : MetaInfoMap(startpage, defaultValue=0),
        'userreal3' : MetaInfoMap(startpage, defaultValue=0),
        'userreal4' : MetaInfoMap(startpage, defaultValue=0),
        }

    grplist = {
        'nrdf' : MetaInfoMap(startpage, defaultValue=0),
        'ref-t' : MetaInfoMap(startpage, defaultValue=0),
        'tau-t' : MetaInfoMap(startpage, defaultValue=0),
        }

    anneallist = {
        'annealing' : MetaInfoMap(startpage, defaultValue=0),
        'annealing-npoints' : MetaInfoMap(startpage, defaultValue=0),
        'acc' : MetaInfoMap(startpage, defaultValue=0,
            subFunc=lambda x: [float(i) for i in x.split()] if(x is not None and
                isinstance(x,str) is True) else None),
        'nfreeze' : MetaInfoMap(startpage, defaultValue=0,
            subFunc=lambda x: x.split() if(x is not None and
                isinstance(x,str) is True) else None),
        'energygrp-flags\[  0\]' : MetaInfoMap(startpage, defaultValue=0),
        'energygrp-flags\[  1\]' : MetaInfoMap(startpage, defaultValue=0),
        'energygrp-flags\[  2\]' : MetaInfoMap(startpage, defaultValue=0),
        }

    startpage.update({
        'nameTranslate'   : metaNameConverter_UnderscoreSpaceDash,
        'metaNameTag'     : 'parm'
        })
    parmlist = {
        'Domain Decomposition on' : MetaInfoMap(startpage, defaultValue=0, replaceTag='Domain Decomposition Ranks', matchWith='EOL'),
        'Dynamic load balancing' : MetaInfoMap(startpage, defaultValue=0),
        '\(option -dds\)' : MetaInfoMap(startpage, defaultValue=0, matchWith='NW'),
        'Optimizing the DD grid for' : MetaInfoMap(startpage, defaultValue=0, replaceTag='Number of cells'),
        'with a minimum initial size of' : MetaInfoMap(startpage, defaultValue=0, replaceTag='Min init size cells', matchWith='EOL'),
        'Domain decomposition grid' : MetaInfoMap(startpage, defaultValue=0),
        'separate PME ranks' : MetaInfoMap(startpage, defaultValue=0),
        'Domain decomposition rank' : MetaInfoMap(startpage, defaultValue=0),
        'coordinates' : MetaInfoMap(startpage, defaultValue=0, replaceTag='Domain Decomposition coordinates'),
        'Table routines are used for coulomb' : MetaInfoMap(startpage, defaultValue=0, replaceTag='use coulomb tables'),
        'Table routines are used for vdw' : MetaInfoMap(startpage, defaultValue=0, replaceTag='use vdw tables'),
        'NS' : MetaInfoMap(startpage, defaultValue=0, alsoMatch=["Cut-off's:"], replaceTag='cutoff ns'),
        'Coulomb' : MetaInfoMap(startpage, defaultValue=0, alsoMatch=["Cut-off's:"], replaceTag='cutoff coulomb'),
        'LJ' : MetaInfoMap(startpage, defaultValue=0, alsoMatch=["Cut-off's:"], replaceTag='cutoff lj'),
        'System total charge' : MetaInfoMap(startpage, defaultValue=0),
        'data points for 1-4 COUL' : MetaInfoMap(startpage, defaultValue=0, matchWith='PW'),
        'data points for 1-4 LJ6' : MetaInfoMap(startpage, defaultValue=0, matchWith='PW'),
        'data points for 1-4 LJ12' : MetaInfoMap(startpage, defaultValue=0, matchWith='PW'),
        'Potential shift: LJ r\^-12:' : MetaInfoMap(startpage, defaultValue=0, replaceTag='Potential shift LJ r12'),
        'r\^-6:' : MetaInfoMap(startpage, defaultValue=0, replaceTag='Potential shift LJ r6'),
        ', Coulomb' : MetaInfoMap(startpage, defaultValue=0, replaceTag='Potential shift Coulomb'),
        'Initializing Parallel LINear Constraint' : MetaInfoMap(startpage, defaultValue='yes', changeTags={'PLINCS' : 'Yes'}, mute=True),
        'PLINCS' : MetaInfoMap(startpage, ignoreMatch=True),
        'number of constraints is' : MetaInfoMap(startpage, defaultValue=0, removeText=' is'),
        'domain decomposition cell size is' : MetaInfoMap(startpage, defaultValue=0, replaceTag='Min size cells'),
        'shrink of DD cells \(option -dds\) is' : MetaInfoMap(startpage, defaultValue=0, removeText=' \(option -dds\) is'),
        'allowed shrink of domain decomposition cells is' : MetaInfoMap(startpage, defaultValue=0, removeText=' is'),
        'non-bonded interactions' : MetaInfoMap(startpage, defaultValue=0),
        'two-body bonded interactions  \(-rdd\)' : MetaInfoMap(startpage, defaultValue=0, replaceTag='two-body bonded interactions'),
        'multi-body bonded interactions  \(-rdd\)' : MetaInfoMap(startpage, defaultValue=0, replaceTag='multi-body bonded interactions'),
        'atoms separated by up to 4 constraints  \(-rcon\)' : MetaInfoMap(startpage, defaultValue=0, replaceTag='atoms separated 4 constraints'),
        'domain decomposition grid' : MetaInfoMap(startpage, defaultValue=0),
        'Center of mass motion removal mode is' : MetaInfoMap(startpage, defaultValue=0, removeText=' is'),
        'Atom distribution over' : MetaInfoMap(startpage, defaultValue=0, replaceTag='Atom distribution domains'),
        'Grid' : MetaInfoMap(startpage, defaultValue=0),
        'Initial temperature' : MetaInfoMap(startpage, defaultValue=0),
        }

    startpage.update({
        'metaNameTag'     : 'mdout',
        'activeSections'  : ['section_single_configuration_calculation']
        })
    mddatalist = {
        'Step' : MetaInfoMap(startpage),
        'Time' : MetaInfoMap(startpage),
        'Lambda' : MetaInfoMap(startpage),
        'Bond' : MetaInfoMap(startpage),
        'Angle' : MetaInfoMap(startpage),
        'G96Bond' : MetaInfoMap(startpage),
        'G96Angle' : MetaInfoMap(startpage),
        'U-B' : MetaInfoMap(startpage),
        'Proper Dih.' : MetaInfoMap(startpage),
        'Improper Dih.' : MetaInfoMap(startpage),
        'CMAP Dih.' : MetaInfoMap(startpage),
        'LJ (SR)' : MetaInfoMap(startpage),
        'LJ (LR)' : MetaInfoMap(startpage),
        'Disper. corr.' : MetaInfoMap(startpage),
        'Coulomb (SR)' : MetaInfoMap(startpage),
        'Coul. recip.' : MetaInfoMap(startpage),
        'LJ-14' : MetaInfoMap(startpage),
        'Coulomb-14' : MetaInfoMap(startpage),
        'Potential' : MetaInfoMap(startpage),
        'Kinetic En.' : MetaInfoMap(startpage),
        'Total Energy' : MetaInfoMap(startpage),
        'Conserved En.' : MetaInfoMap(startpage),
        'Temperature' : MetaInfoMap(startpage),
        'Pres. DC (bar)' : MetaInfoMap(startpage),
        'Pressure (bar)' : MetaInfoMap(startpage),
        'Constr. rmsd' : MetaInfoMap(startpage),
        'dVremain/dl' : MetaInfoMap(startpage),
        'Virial-Tensor' : MetaInfoMap(startpage),
        'Pressure-Tensor' : MetaInfoMap(startpage),
        }

    startpage.update({
        'nameTranslate'   : metaNameConverter_UnderscoreSpaceDash,
        'metaNameTag'     : 'mdout',
        'activeSections'  : ['section_single_configuration_calculation']
        })
    edrdatalist = {
        'Time' : MetaInfoMap(startpage),
        'Bond' : MetaInfoMap(startpage),
        'Angle' : MetaInfoMap(startpage),
        'G96Bond' : MetaInfoMap(startpage),
        'G96Angle' : MetaInfoMap(startpage),
        'U-B' : MetaInfoMap(startpage),
        'Proper Dih.' : MetaInfoMap(startpage),
        'Improper Dih.' : MetaInfoMap(startpage),
        'CMAP Dih.' : MetaInfoMap(startpage),
        'LJ (SR)' : MetaInfoMap(startpage),
        'LJ (LR)' : MetaInfoMap(startpage),
        'Disper. corr.' : MetaInfoMap(startpage),
        'Coulomb (SR)' : MetaInfoMap(startpage),
        'Coul. recip.' : MetaInfoMap(startpage),
        'LJ-14' : MetaInfoMap(startpage),
        'Coulomb-14' : MetaInfoMap(startpage),
        'Potential' : MetaInfoMap(startpage),
        'Kinetic En.' : MetaInfoMap(startpage),
        'Total Energy' : MetaInfoMap(startpage),
        'Conserved En.' : MetaInfoMap(startpage),
        'Temperature' : MetaInfoMap(startpage),
        'Pres. DC (bar)' : MetaInfoMap(startpage),
        'Pressure (bar)' : MetaInfoMap(startpage),
        'Constr. rmsd' : MetaInfoMap(startpage),
        'dVremain/dl' : MetaInfoMap(startpage),
        'Virial-Tensor' : MetaInfoMap(startpage),
        'Pressure-Tensor' : MetaInfoMap(startpage),
        }

    startpage.update({
        'metaNameTag'     : 'parm',
        'activeSections'  : [PARSERTAG+'_mdin_method']
        })
    extralist = {
        'flags' : MetaInfoMap(startpage),
        'number_of_atoms' : MetaInfoMap(startpage),
        'box_info' : MetaInfoMap(startpage),
        'unitcell_radius' : MetaInfoMap(startpage),
        'total_memory' : MetaInfoMap(startpage),
        'file_format' : MetaInfoMap(startpage),
        'file_version' : MetaInfoMap(startpage),
        'file_date' : MetaInfoMap(startpage),
        'file_time' : MetaInfoMap(startpage)
        }

    startpage.update({
        'metaNameTag'     : 'mdout',
        'activeSections'  : ['section_single_configuration_calculation']
        })
    # Additional dictionary to store stepping parameters
    # such as output log frequency, trajectory writing
    # freuquency and the actually needed number of
    # section_single_configuration_calculation
    #stepcontrol = {
    #    'logsteps'       : MetaInfoMap(startpage),
    #    'trajsteps'      : MetaInfoMap(startpage),
    #    'enersteps'      : MetaInfoMap(startpage),
    #    'steps'          : MetaInfoMap(startpage),
    #    'follow'         : MetaInfoMap(startpage),
    #    'sectioncontrol' : MetaInfoMap(startpage)
    #    }

    if deflist == 'mddata':
        namelist = mddatalist
    elif deflist == 'edrdata':
        namelist = edrdatalist
    elif deflist == 'arch':
        namelist = archlist
    elif deflist == 'qm':
        namelist = qmlist
    elif deflist == 'grp':
        namelist = grplist
    elif deflist == 'anneal':
        namelist = anneallist
    elif deflist == 'parm':
        namelist = parmlist
    elif deflist == 'extra':
        namelist = extralist
    elif deflist == 'stepcontrol':
        namelist = stepcontrol
    else:
        namelist = cntrllist
    return MapDictionary(namelist)

def set_Dictionaries(self):
        self.unitDict = get_unitDict('si')
        self.fileDict = get_fileListDict()
        self.cntrlDict = get_nameListDict(self, 'cntrl')
        self.archDict = get_nameListDict(self, 'arch')
        self.qmDict = get_nameListDict(self, 'qm')
        self.grpDict = get_nameListDict(self, 'grp')
        self.edrDict = get_nameListDict(self, 'edrdata')
        self.annealDict = get_nameListDict(self, 'anneal')
        self.parmDict = get_nameListDict(self, 'parm')
        self.mddataDict = get_nameListDict(self, 'mddata')
        self.extraDict = get_nameListDict(self, 'extra')
        #self.stepcontrolDict = get_nameListDict('stepcontrol')
        self.stepcontrolDict = {
            'logsteps'       : None,
            'nextlogsteps'   : None,
            'trajsteps'      : None,
            'enersteps'      : None,
            'steps'          : None,
            'follow'         : None,
            'sectioncontrol' : None,
        }
        self.metaDicts = {
            'file'   : self.fileDict,
            'cntrl'  : self.cntrlDict,
            'arch'   : self.archDict,
            'qm'     : self.qmDict,
            'grp'    : self.grpDict,
            'anneal' : self.annealDict,
            'parm'   : self.parmDict,
            'mddata' : self.mddataDict,
            'edr'    : self.edrDict,
            'extra'  : self.extraDict,
            }
        self.sectionDict = {
                'section' : {
                    "metaNameTag"      : ['input_output_files',
                                          'control_parameters'],
                    "listTypStr"       : 'type_section',
                    "repeatingSection" : False,
                    "supraNames"       : ["section_run"]
                    }
                }


def get_updateDictionary(self, defname):

    def value_convert_array(itemdict):
        val = itemdict["value"]
        if(val is not None and isinstance(val, str)):
            if("[" in val or "(" in val):
                return False, np.asarray(val), itemdict
        return False, val, itemdict

    def value_convert_array_norm(itemdict):
        val = []
        depdict = itemdict["depends"]
        for item in depdict["value"]:
            st=self.metaStorage.fetchAttrValue(item)
            if(st is not None and isinstance(st, str)):
                if("[" in st or "(" in st):
                    val.append(list(storedtext))
        if val:
            return False, np.linalg.norm(np.asarray(val)), itemdict
        else:
            return False, itemdict["value"], itemdict

    startpage = {
        'metaHeader'      : '',
        'metaNameTag'     : '',
        'activeMetaNames' : [],
        'activeSections'  : ['section_sampling_method','section_single_configuration_calculation',PARSERTAG+'_mdout']
        }

    # ---------------------------------------------------------------
    #   Definitions of meta data values for section_sampling_method
    # ---------------------------------------------------------------
    sampling = {
        'ensemble_type' : MetaInfoMap(startpage, activeInfo=True,
            depends=[
                {'test' : [['integrator', ' in [\"steep\",'
                                          '\"cg\",'
                                          '\"l-bfgs\"]']],
                 'assign' : 'minimization'},
                {'test' : [['integrator', ' in [\"md\",'
                                          '\"md-vv\",'
                                          '\"sd\",'
                                          '\"sd2\",'
                                          '\"md-vv-avek\"]'],
                           ['tcoupl', ' in [\"no\",\"No\",\"NO\"]']],
                 'assign' : 'NVE'},
                {'test' : [['integrator', ' in [\"md\",'
                                          '\"sd\",'
                                          '\"sd2\",'
                                          '\"md-vv\",'
                                          '\"md-vv-avek\"]'],
                           ['tcouple', ' in [\"berendsen\",'
                                       '\"nose-hoover\",'
                                       '\"andersen\",'
                                       '\"anderson-massive\",'
                                       '\"v-rescale\"]'],
                           ['pcoupl', ' in [\"no\",\"No\",\"NO\"]']],
                 'assign' : 'NVT'},
                {'test' : [['integrator', ' in [\"md\",'
                                          '\"sd\",'
                                          '\"sd2\",'
                                          '\"md-vv\",'
                                          '\"md-vv-avek\"]'],
                           ['tcouple', ' in [\"berendsen\",'
                                       '\"nose-hoover\",'
                                       '\"andersen\",'
                                       '\"anderson-massive\",'
                                       '\"v-rescale\"]'],
                           ['pcoupl', ' in [\"berendsen\",'
                                      '\"Parinello-Rahman\"]']],
                 'assign' : 'NPT'},
                {'test' : [['integrator', '== \"bd\"']],
                 'assign' : 'Langevin'}
                ],
            lookupdict=self.cntrlDict
            ),
        'sampling_method' : MetaInfoMap(startpage, activeInfo=True,
            depends=[
                {'test' : [['integrator', ' in [\"steep\",'
                                          '\"cg\",'
                                          '\"l-bfgs\"]']],
                 'assign' : 'geometry_optimization'},
                {'test' : [['integrator', ' in [\"md\",'
                                          '\"sd\",'
                                          '\"sd2\",'
                                          '\"bd\",'
                                          '\"md-vv\",'
                                          '\"md-vv-avek\"]']],
                 'assign' : 'molecular_dynamics'}
                ],
            lookupdict=self.cntrlDict
            ),
#            'settings_geometry_optimization' : MetaInfoMap(startpage),
#            'settings_metadynamics' : MetaInfoMap(startpage),
#            'settings_molecular_dynamics' : MetaInfoMap(startpage),
#            'settings_Monte_Carlo' : MetaInfoMap(startpage),
#        'geometry_optimization_energy_change' : MetaInfoMap(startpage,
#            depends={
#                '' : {'imin' : '1'},
#                },
#            lookupdict=self.cntrlDict
#            ),
#       'geometry_optimization_geometry_change' : MetaInfoMap(startpage),
        'geometry_optimization_method' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['integrator', '== \"steep\"']],
                 'assign' : 'SD'},
                {'test' : [['integrator', '== \"cg\"']],
                 'assign' : 'CG'},
                {'test' : [['integrator', '== \"l-bfgs\"']],
                 'assign' : 'L-BFGS'},
                ],
            lookupdict=self.cntrlDict,
            activeInfo=True,
            ),
#       'geometry_optimization_threshold_force' : MetaInfoMap(startpage),
        'x_gromacs_barostat_type' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['integrator', ' in [\"md\",'
                                          '\"sd\",'
                                          '\"sd2\",'
                                          '\"bd\",'
                                          '\"md-vv\",'
                                          '\"md-vv-avek\"]'],
                           ['pcoupl', '== \"berendsen\"']],
                 'assign' : 'Berendsen'},
                {'test' : [['integrator', ' in [\"md\",'
                                          '\"sd\",'
                                          '\"sd2\",'
                                          '\"bd\",'
                                          '\"md-vv\",'
                                          '\"md-vv-avek\"]'],
                           ['pcoupl', '== \"Parinello-Rahman\"']],
                 'assign' : 'Parinello-Rahman'}
                ],
            lookupdict=self.cntrlDict,
            #autoSections=True,
            activeInfo=True,
            activeSections=['settings_barostat']
            ),
        'x_gromacs_barostat_target_pressure' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['integrator', ' in [\"md\",'
                                          '\"sd\",'
                                          '\"sd2\",'
                                          '\"bd\",'
                                          '\"md-vv\",'
                                          '\"md-vv-avek\"]'],
                           ['pcoupl', 'not in [\"no\",\"No\",\"NO\"]']],
                 'value' : 'ref-p\[    0\]'}
                ],
            value=(lambda x: np.array(x['val']).flatten().shape[0] if(
                x is not None and x['val'] is not None) else None)(
                self.metaStorage.fetchAttrValue(
                    'x_gromacs_inout_control_refp0'
                    )
                ),
            lookupdict=self.cntrlDict,
            valtype='float',
            unitdict=self.unitDict,
            unit='bar',
            #autoSections=True,
            activeInfo=True,
            activeSections=['settings_barostat']
            ),
        'x_gromacs_barostat_tau' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['integrator', ' in [\"md\",'
                                          '\"sd\",'
                                          '\"sd2\",'
                                          '\"bd\",'
                                          '\"md-vv\",'
                                          '\"md-vv-avek\"]'],
                           ['pcoupl', ' not in [\"no\",\"No\",\"NO\"]']],
                 'value' : 'tau-p'}
                ],
            lookupdict=self.cntrlDict,
            valtype='float',
            unitdict=self.unitDict,
            unit='pico-second',
            #autoSections=True,
            activeInfo=True,
            activeSections=['settings_barostat']
            ),
        'x_gromacs_integrator_type' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['integrator', ' in [\"md\",'
                                          '\"sd\",'
                                          '\"sd2\",'
                                          '\"md-vv\",'
                                          '\"md-vv-avek\"]']],
                 'assign' : 'verlet'},
                {'test' : [['integrator', ' in [\"steep\",'
                                          '\"cg\",'
                                          '\"l-bfgs\"]']],
                 'assign' : 'minimization'},
                {'test' : [['integrator', '== \"bd\"']],
                 'assign' : 'Langevin'}
                ],
            lookupdict=self.cntrlDict,
            #autoSections=True,
            activeInfo=True,
            activeSections=['settings_integrator']
            ),
        'x_gromacs_integrator_dt' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['integrator', ' in [\"md\",'
                                          '\"sd\",'
                                          '\"sd2\",'
                                          '\"bd\",'
                                          '\"md-vv\",'
                                          '\"md-vv-avek\"]']],
                 'value' : 'dt'}
                ],
            lookupdict=self.cntrlDict,
            valtype='float',
            unitdict=self.unitDict,
            unit='pico-second',
            #autoSections=True,
            activeInfo=True,
            activeSections=['settings_integrator']
            ),
        'x_gromacs_number_of_steps_requested' : MetaInfoMap(startpage,
            depends=[{'value' : 'nsteps'}],
            lookupdict=self.cntrlDict,
            valtype='int',
            #autoSections=True,
            activeInfo=True,
            activeSections=['settings_integrator']
            ),
        'x_gromacs_thermostat_type' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['integrator', ' in [\"md\",'
                                          '\"sd\",'
                                          '\"sd2\",'
                                          '\"bd\",'
                                          '\"md-vv\",'
                                          '\"md-vv-avek\"]'],
                           ['tcouple', '== \"berendsen\"']],
                 'assign' : 'Berendsen'},
                {'test' : [['integrator', ' in [\"md\",'
                                          '\"sd\",'
                                          '\"sd2\",'
                                          '\"bd\",'
                                          '\"md-vv\",'
                                          '\"md-vv-avek\"]'],
                           ['tcouple', '== \"nose-hoover\"']],
                 'assign' : 'Nose-Hoover'},
                {'test' : [['integrator', ' in [\"md\",'
                                          '\"sd\",'
                                          '\"sd2\",'
                                          '\"bd\",'
                                          '\"md-vv\",'
                                          '\"md-vv-avek\"]'],
                           ['tcouple', '== \"andersen\"']],
                 'assign' : 'Andersen'},
                {'test' : [['integrator', ' in [\"md\",'
                                          '\"sd\",'
                                          '\"sd2\",'
                                          '\"bd\",'
                                          '\"md-vv\",'
                                          '\"md-vv-avek\"]'],
                           ['tcouple', '== \"anderson-massive\"']],
                 'assign' : 'Andersen-Massive'},
                {'test' : [['integrator', ' in [\"md\",'
                                          '\"sd\",'
                                          '\"sd2\",'
                                          '\"bd\",'
                                          '\"md-vv\",'
                                          '\"md-vv-avek\"]'],
                           ['tcouple', '== \"v-rescale\"']],
                 'assign' : 'Velocity Rescaling'},
                ],
            lookupdict=self.cntrlDict,
            #autoSections=True,
            activeInfo=True,
            activeSections=['settings_thermostat']
            ),
        'x_gromacs_thermostat_target_temperature' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['integrator', ' in [\"md\",'
                                          '\"sd\",'
                                          '\"sd2\",'
                                          '\"bd\",'
                                          '\"md-vv\",'
                                          '\"md-vv-avek\"]'],
                           ['tcoupl', ' not in [\"no\",\"No\",\"NO\"]']],
                 'value' : 'ref-t'}
                ],
            lookupdict=self.cntrlDict,
            valtype='float',
            unitdict=self.unitDict,
            unit='kelvin',
            #autoSections=True,
            activeInfo=True,
            activeSections=['settings_thermostat']
            ),
        'x_gromacs_thermostat_tau' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['integrator', ' in [\"md\",'
                                          '\"sd\",'
                                          '\"sd2\",'
                                          '\"bd\",'
                                          '\"md-vv\",'
                                          '\"md-vv-avek\"]'],
                           ['tcoupl', ' not in [\"no\",\"No\",\"NO\"]']],
                 'value' : 'tau-t'},
                ],
            lookupdict=self.cntrlDict,
            valtype='float',
            unitdict=self.unitDict,
            unit='pico-second',
            #autoSections=True,
            activeInfo=True,
            activeSections=['settings_thermostat']
            ),
        'x_gromacs_langevin_bd_fric' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['integrator', '== \"bd\"']],
                 'value' : 'bd-fric'}
                ],
            lookupdict=self.cntrlDict,
            valtype='float',
            #autoSections=True,
            activeInfo=True,
            activeSections=['settings_thermostat']
            ),
        'x_gromacs_periodicity_type' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['pbc', ' in [\"no\",\"No\",\"NO\"]']],
                 'assign' : 'no periodic boundaries'},
                {'test' : [['pbc', ' not in [\"no\",\"No\",\"NO\"]']],
                 'assign' : 'periodic boundaries'}
                ],
            lookupdict=self.cntrlDict,
            #autoSections=True,
            activeInfo=True,
            activeSections=['settings_thermostat']
            )
        }

    # ------------------------------------------------------------
    #   Definitions for section_single_configuration_calculation
    # ------------------------------------------------------------
    singleconfcalc = {
        #'atom_forces_type' : MetaInfoMap(startpage,
        #    depends=[{'assign' : 'Force Field'}],
        #    lookupdict=self.mddataDict
        #    ),
        'energy_correction_entropy' : MetaInfoMap(startpage),
        #'energy_current' : MetaInfoMap(startpage,
        #    depends=[{'value' : 'EPtot'}],
        #    lookupdict=self.mddataDict
        #    ),
        'energy_current' : MetaInfoMap(startpage),
        'energy_electrostatic' : MetaInfoMap(startpage,
            depends=[{'value' : 'Potential'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='kilo-joule/mol',
            activeInfo=True,
            lookupdict=self.mddataDict
            ),
        'energy_free_per_atom' : MetaInfoMap(startpage),
        'energy_free' : MetaInfoMap(startpage),
        #'energy_method_current' : MetaInfoMap(startpage,
        #    depends=[{'assign' : 'Force Field'}],
        #    lookupdict=self.mddataDict
        #    ),
        'energy_T0_per_atom' : MetaInfoMap(startpage),
        'energy_total_T0_per_atom' : MetaInfoMap(startpage),
        'energy_total_T0' : MetaInfoMap(startpage),
        'energy_total' : MetaInfoMap(startpage,
            depends=[{'value' : 'Total Energy'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='kilo-joule/mol',
            activeInfo=True,
            lookupdict=self.mddataDict
            ),
        'hessian_matrix' : MetaInfoMap(startpage),
        'single_configuration_calculation_converged' : MetaInfoMap(startpage),
        'single_configuration_calculation_to_system_ref' : MetaInfoMap(startpage),
        'single_configuration_to_calculation_method_ref' : MetaInfoMap(startpage),
        'time_calculation' : MetaInfoMap(startpage),
        'time_single_configuration_calculation_cpu1_end' : MetaInfoMap(startpage),
        'time_single_configuration_calculation_cpu1_start' : MetaInfoMap(startpage),
        'time_single_configuration_calculation_date_end' : MetaInfoMap(startpage),
        'time_single_configuration_calculation_date_start' : MetaInfoMap(startpage),
        'time_single_configuration_calculation_wall_end' : MetaInfoMap(startpage),
        'time_single_configuration_calculation_wall_start' : MetaInfoMap(startpage),
#        'stress_tensor_kind' : MetaInfoMap(startpage),
#        'stress_tensor_kind' : MetaInfoMap(startpage,
#            depends=[
#                {'test' : [['integrator', '== \"steep\"']],
#                 'assign' : 'geometry_optimization'},
#                {'test' : [['integrator', '== \"cg\"']],
#                 'assign' : 'geometry_optimization'},
#                {'test' : [['integrator', '== \"l-bfgs\"']],
#                 'assign' : 'geometry_optimization'},
#                {'test' : [['integrator', '== \"md\"']],
#                 'assign' : 'molecular_dynamics'},
#                {'test' : [['integrator', '== \"md-vv\"']],
#                 'assign' : 'molecular_dynamics'},
#                {'test' : [['integrator', '== \"md-vv-avek\"']],
#                 'assign' : 'molecular_dynamics'}
#                ],
#            activeInfo=True,
#            lookupdict=self.cntrlDict
#            ),
        'stress_tensor_value' : MetaInfoMap(startpage)
        }

    # section_single_energy_van_der_Waals of section_single_configuration_calculation
    singlevdw = {
        'energy_van_der_Waals_value' : MetaInfoMap(startpage,
            depends=[{'value' : 'VDWAALS'}],
            lookupdict=self.mddataDict,
            valtype='float',
            unitdict=self.unitDict,
            unit='kilo-joule/mol',
            #autoSections=True,
            activeSections=['section_energy_van_der_Waals']
            ),
        'x_gromacs_energy_lj_14_value' : MetaInfoMap(startpage,
            depends=[{'value' : 'LJ-14'}],
            lookupdict=self.mddataDict,
            valtype='float',
            unitdict=self.unitDict,
            unit='kilo-joule/mol',
            #autoSections=True,
            activeSections=['section_energy_van_der_Waals']
            ),
        }

    # ------------------------------------------
    #   Definitions for section_frame_sequence
    # ------------------------------------------
    frameseq = {
        'frame_sequence_conserved_quantity_frames' : MetaInfoMap(startpage,
            depends=[{'store' : 'Step'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'frame_sequence_conserved_quantity_stats' : MetaInfoMap(startpage),
        'frame_sequence_conserved_quantity' : MetaInfoMap(startpage,
            depends=[{'store' : 'Conserved En.'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='kilo-joule/mol',
            lookupdict=self.mddataDict
            ),
        'frame_sequence_continuation_kind' : MetaInfoMap(startpage),
        'frame_sequence_external_url' : MetaInfoMap(startpage),
        'frame_sequence_kinetic_energy_frames' : MetaInfoMap(startpage,
            depends=[{'store' : 'Step'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'frame_sequence_kinetic_energy_stats' : MetaInfoMap(startpage),
        'frame_sequence_kinetic_energy' : MetaInfoMap(startpage,
            depends=[{'store' : 'Kinetic En.'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='kilo-joule/mol',
            lookupdict=self.mddataDict
            ),
        'frame_sequence_local_frames_ref' : MetaInfoMap(startpage),
        'frame_sequence_potential_energy_frames' : MetaInfoMap(startpage,
            depends=[{'store' : 'Step'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'frame_sequence_potential_energy_stats' : MetaInfoMap(startpage),
        'frame_sequence_potential_energy' : MetaInfoMap(startpage,
            depends=[{'store' : 'Potential'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='kilo-joule/mol',
            lookupdict=self.mddataDict
            ),
        'frame_sequence_pressure_frames' : MetaInfoMap(startpage,
            depends=[{'store' : 'Step'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'frame_sequence_pressure_stats' : MetaInfoMap(startpage),
        'frame_sequence_pressure' : MetaInfoMap(startpage,
            depends=[{'store' : 'Pressure (bar)'}],
            valtype='float',
            lookupdict=self.mddataDict
            ),
        'frame_sequence_temperature_frames' : MetaInfoMap(startpage,
            depends=[{'test'  : [['Volume', ' is not None']],
                      'store' : 'Step'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'frame_sequence_temperature_stats' : MetaInfoMap(startpage),
        'frame_sequence_temperature' : MetaInfoMap(startpage,
            depends=[{'store' : 'Temperature'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='Kelvin',
            lookupdict=self.mddataDict
            ),
        'frame_sequence_time' : MetaInfoMap(startpage,
            depends=[{'store' : 'Time'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='pico-second',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_volume_frames' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'test'  : [['Volume', ' is not None']],
                      'store' : 'Step'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_volume' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'Volume'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='Angstrom**3',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_density_frames' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'Step'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_density' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'Density'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='1.0/Angstrom**3',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_bond_energy_frames' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'Step'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_bond_energy' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'Bond'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='kilo-joule/mol',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_ubond_energy_frames' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'Step'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_ubond_energy' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'U-B'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='kilo-joule/mol',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_proper_dihedral_energy_frames' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'Step'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_proper_dihedral_energy' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'Proper Dih.'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='kilo-joule/mol',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_improper_dihedral_energy_frames' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'Step'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_improper_dihedral_energy' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'Improper Dih.'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='kilo-joule/mol',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_cmap_dihedral_energy_frames' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'Step'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_cmap_dihedral_energy' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'CMAP Dih.'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='kilo-joule/mol',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_lj_14_energy_frames' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'Step'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_lj_14_energy' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'LJ-14'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='kilo-joule/mol',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_coulomb_14_energy_frames' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'Step'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_coulomb_14_energy' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'Coulomb-14'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='kilo-joule/mol',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_lj_sr_energy_frames' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'Step'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_lj_sr_energy' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'LJ (SR)'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='kilo-joule/mol',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_coulomb_sr_energy_frames' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'Step'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_coulomb_sr_energy' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'Coulomb (SR)'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='kilo-joule/mol',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_constrain_rmsd_frames' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'Step'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'x_gromacs_frame_sequence_constrain_rmsd' : MetaInfoMap(startpage,
            activeSections=['section_frame_sequence'],
            depends=[{'store' : 'Constr. rmsd'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='kilo-joule/mol',
            lookupdict=self.mddataDict
            ),
        #'frame_sequence_to_sampling_ref' : MetaInfoMap(startpage),
        'geometry_optimization_converged' : MetaInfoMap(startpage,
            value=self.minConverged
            ),
        #'previous_sequence_ref' : MetaInfoMap(startpage)
        }

    frameseqend = {
        'number_of_conserved_quantity_evaluations_in_sequence' : MetaInfoMap(startpage,
            activeInfo=True,
            value=(lambda x: np.array(x['val']).flatten().shape[0] if(
                x is not None and x['val'] is not None) else None)(
                self.metaStorage.fetchAttrValue(
                    'frame_sequence_conserved_quantity_frames'
                    )
                )
            ),
        'number_of_frames_in_sequence' : MetaInfoMap(startpage, activeInfo=True,
            value=(lambda x: np.array(x['val']).flatten().shape[0] if(
                x is not None and x['val'] is not None) else None)(
                self.metaStorage.fetchAttrValue(
                    'frame_sequence_potential_energy_frames'
                    )
                )
            ),
        'number_of_kinetic_energies_in_sequence' : MetaInfoMap(startpage, activeInfo=True,
            value=(lambda x: np.array(x['val']).flatten().shape[0] if(
                x is not None and x['val'] is not None) else None)(
                self.metaStorage.fetchAttrValue(
                    'frame_sequence_kinetic_energy_frames'
                    )
                )
            ),
        'number_of_potential_energies_in_sequence' : MetaInfoMap(startpage,
            activeInfo=True,
            value=(lambda x: np.array(x['val']).flatten().shape[0] if(
                x is not None and x['val'] is not None) else None)(
                self.metaStorage.fetchAttrValue(
                    'frame_sequence_potential_energy_frames'
                    )
                )
            ),
        'number_of_pressure_evaluations_in_sequence' : MetaInfoMap(startpage,
            activeInfo=True,
            value=(lambda x: np.array(x['val']).flatten().shape[0] if(
                x is not None and x['val'] is not None) else None)(
                self.metaStorage.fetchAttrValue(
                    'frame_sequence_pressure_frames'
                    )
                )
            ),
        'number_of_temperatures_in_sequence' : MetaInfoMap(startpage,
            activeInfo=True,
            value=(lambda x: np.array(x['val']).flatten().shape[0] if(
                x is not None and x['val'] is not None) else None)(
                self.metaStorage.fetchAttrValue(
                    'frame_sequence_temperature_frames'
                    )
                )
            ),
        'x_gromacs_number_of_volumes_in_sequence' : MetaInfoMap(startpage,
            activeInfo=True,
            activeSections=['section_frame_sequence'],
            value=(lambda x: np.array(x['val']).flatten().shape[0] if(
                x is not None and x['val'] is not None) else None)(
                self.metaStorage.fetchAttrValue(
                    'x_gromacs_frame_sequence_volume_frames'
                    )
                )
            ),
        'x_gromacs_number_of_densities_in_sequence' : MetaInfoMap(startpage,
            activeInfo=True,
            activeSections=['section_frame_sequence'],
            value=(lambda x: np.array(x['val']).flatten().shape[0] if(
                x is not None and x['val'] is not None) else None)(
                self.metaStorage.fetchAttrValue(
                    'x_gromacs_frame_sequence_density_frames'
                    )
                )
            ),
        'x_gromacs_number_of_bond_energies_in_sequence' : MetaInfoMap(startpage,
            activeInfo=True,
            activeSections=['section_frame_sequence'],
            value=(lambda x: np.array(x['val']).flatten().shape[0] if(
                x is not None and x['val'] is not None) else None)(
                self.metaStorage.fetchAttrValue(
                    'x_gromacs_frame_sequence_bond_energy_frames'
                    )
                )
            ),
        'x_gromacs_number_of_ubond_energies_in_sequence' : MetaInfoMap(startpage,
            activeInfo=True,
            activeSections=['section_frame_sequence'],
            value=(lambda x: np.array(x['val']).flatten().shape[0] if(
                x is not None and x['val'] is not None) else None)(
                self.metaStorage.fetchAttrValue(
                    'x_gromacs_frame_sequence_ubond_energy_frames'
                    )
                )
            ),
        'x_gromacs_number_of_coulomb_sr_energies_in_sequence' : MetaInfoMap(startpage,
            activeInfo=True,
            activeSections=['section_frame_sequence'],
            value=(lambda x: np.array(x['val']).flatten().shape[0] if(
                x is not None and x['val'] is not None) else None)(
                self.metaStorage.fetchAttrValue(
                    'x_gromacs_frame_sequence_coulomb_sr_energy_frames'
                    )
                )
            ),
        'x_gromacs_number_of_coulomb_14_energies_in_sequence' : MetaInfoMap(startpage,
            activeInfo=True,
            activeSections=['section_frame_sequence'],
            value=(lambda x: np.array(x['val']).flatten().shape[0] if(
                x is not None and x['val'] is not None) else None)(
                self.metaStorage.fetchAttrValue(
                    'x_gromacs_frame_sequence_coulomb_14_energy_frames'
                    )
                )
            ),
        'x_gromacs_number_of_lj_14_energies_in_sequence' : MetaInfoMap(startpage,
            activeInfo=True,
            activeSections=['section_frame_sequence'],
            value=(lambda x: np.array(x['val']).flatten().shape[0] if(
                x is not None and x['val'] is not None) else None)(
                self.metaStorage.fetchAttrValue(
                    'x_gromacs_frame_sequence_lj_14_energy_frames'
                    )
                )
            ),
        'x_gromacs_number_of_lj_sr_energies_in_sequence' : MetaInfoMap(startpage,
            activeInfo=True,
            activeSections=['section_frame_sequence'],
            value=(lambda x: np.array(x['val']).flatten().shape[0] if(
                x is not None and x['val'] is not None) else None)(
                self.metaStorage.fetchAttrValue(
                    'x_gromacs_frame_sequence_lj_sr_energy_frames'
                    )
                )
            ),
        'x_gromacs_number_of_proper_dihedral_energies_in_sequence' : MetaInfoMap(startpage,
            activeInfo=True,
            activeSections=['section_frame_sequence'],
            value=(lambda x: np.array(x['val']).flatten().shape[0] if(
                x is not None and x['val'] is not None) else None)(
                self.metaStorage.fetchAttrValue(
                    'x_gromacs_frame_sequence_proper_dihedral_energy_frames'
                    )
                )
            ),
        'x_gromacs_number_of_improper_dihedral_energies_in_sequence' : MetaInfoMap(startpage,
            activeInfo=True,
            activeSections=['section_frame_sequence'],
            value=(lambda x: np.array(x['val']).flatten().shape[0] if(
                x is not None and x['val'] is not None) else None)(
                self.metaStorage.fetchAttrValue(
                    'x_gromacs_frame_sequence_improper_dihedral_energy_frames'
                    )
                )
            ),
        'x_gromacs_number_of_cmap_dihedral_energies_in_sequence' : MetaInfoMap(startpage,
            activeInfo=True,
            activeSections=['section_frame_sequence'],
            value=(lambda x: np.array(x['val']).flatten().shape[0] if(
                x is not None and x['val'] is not None) else None)(
                self.metaStorage.fetchAttrValue(
                    'x_gromacs_frame_sequence_cmap_dihedral_energy_frames'
                    )
                )
            ),
        'x_gromacs_number_of_constrain_rmsd_in_sequence' : MetaInfoMap(startpage,
            activeInfo=True,
            activeSections=['section_frame_sequence'],
            value=(lambda x: np.array(x['val']).flatten().shape[0] if(
                x is not None and x['val'] is not None) else None)(
                self.metaStorage.fetchAttrValue(
                    'x_gromacs_frame_sequence_constrain_rmsd_frames'
                    )
                )
            ),
        #'previous_sequence_ref' : MetaInfoMap(startpage)
        }

    # ----------------------------------
    #   Definitions for section_system
    # ----------------------------------
    # section_system
    sec_system = {
        #'topology_ref' : MetaInfoMap(startpage),
        'atom_velocities' : MetaInfoMap(startpage),
        'configuration_raw_gid' : MetaInfoMap(startpage),
        'local_rotations' : MetaInfoMap(startpage),
        'number_of_atoms' : MetaInfoMap(startpage,
            depends=[{'value' : 'number_of_atoms'}],
            valtype='int',
            lookupdict=self.extraDict
            ),
        'number_of_sites' : MetaInfoMap(startpage),
        'number_of_symmetry_operations' : MetaInfoMap(startpage),
        'reduced_symmetry_matrices' : MetaInfoMap(startpage),
        'reduced_symmetry_translations' : MetaInfoMap(startpage),
        'SC_matrix' : MetaInfoMap(startpage),
        'spacegroup_3D_choice' : MetaInfoMap(startpage),
        'spacegroup_3D_hall' : MetaInfoMap(startpage),
        'spacegroup_3D_international' : MetaInfoMap(startpage),
        'spacegroup_3D_number' : MetaInfoMap(startpage),
        'spacegroup_3D_origin_shift' : MetaInfoMap(startpage),
        'spacegroup_3D_pointgroup' : MetaInfoMap(startpage),
        'spacegroup_3D_std_lattice' : MetaInfoMap(startpage),
        'spacegroup_3D_std_positions' : MetaInfoMap(startpage),
        'spacegroup_3D_std_types' : MetaInfoMap(startpage),
        'spacegroup_3D_trasformation_matrix' : MetaInfoMap(startpage),
        'spacegroup_3D_wyckoff' : MetaInfoMap(startpage),
        'symmorphic' : MetaInfoMap(startpage),
        'system_name' : MetaInfoMap(startpage,
            subfunction={
                'function' : MDData.MDDataConverter.topology_system_name,
                'supportDict' : self.topology.topoDict,
                },
            #functionbase=self
            ),
        'time_reversal_symmetry' : MetaInfoMap(startpage)
        }

    # section_configuration_core of section_system
    configuration_core = {
        'number_of_electrons' : MetaInfoMap(startpage,
            value=0,
            valtype='float',
            ),
        'atom_labels' : MetaInfoMap(startpage,
            #subfunction=self.system_atom_labels()
            ),
        'atom_positions' : MetaInfoMap(startpage,
            #subfunction=self.system_atom_positions()
            ),
        'configuration_periodic_dimensions' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['pbc', '== "no"']],
                 'assign' : np.asarray([False, False, False])},
                {'test' : [['pbc', '== "xy"']],
                 'assign' : np.asarray([True, True, False])},
                {'test' : [['pbc', '== "xyz"']],
                 'assign' : np.asarray([True, True, True])}
                ],
            lookupdict=self.cntrlDict,
            ),
        'embedded_system' : MetaInfoMap(startpage),
        'lattice_vectors' : MetaInfoMap(startpage,
            #subfunction=self.system_lattice_vectors()
            ),
        'simulation_cell' : MetaInfoMap(startpage,
            #subfunction=self.system_simulation_cell()
            )
        }

    # section_spacegroup_3D_operation of section_system
    spacegroup_op = {
            'spacegroup_3D_rotation' : MetaInfoMap(startpage),
            'spacegroup_3D_translation' : MetaInfoMap(startpage)
        }

    # section_system_to_system_refs of section_system
    sys_to_sys = {
            'system_to_system_kind' : MetaInfoMap(startpage),
            'system_to_system_ref' : MetaInfoMap(startpage)
        }

    # --------------------------------------------------------
    #   Definitions of meta data values for section_topology
    # --------------------------------------------------------
    # section_topology of section_run
    topology = {
        'atom_to_molecule' : MetaInfoMap(startpage,
            subfunction={
                'function' : MDData.MDDataConverter.topology_atom_to_mol,
                'supportDict' : self.topology.topoDict,
                },
            #functionbase=self
            ),
        'molecule_to_molecule_type_map' : MetaInfoMap(startpage),
        'number_of_topology_atoms' : MetaInfoMap(startpage,
            depends=[{'value' : 'NATOM'}],
            valtype='int',
            lookupdict=self.parmDict
            ),
        'number_of_topology_molecules' : MetaInfoMap(startpage,
            subfunction={
                'function' : MDData.MDDataConverter.topology_num_topo_mol,
                'supportDict' : self.topology.topoDict,
                },
            valtype='int',
            ),
        'topology_force_field_name' : MetaInfoMap(startpage,
            value='Force Field',
            )
        }

    # section_atom_type of section_topology
    atom_type = {
        'atom_type_charge' : MetaInfoMap(startpage),
        'atom_type_mass' : MetaInfoMap(startpage),
        'atom_type_name' : MetaInfoMap(startpage)
        }

    # section_constraint of section_topology
    constraint = {
        'constraint_atoms' : MetaInfoMap(startpage),
        'constraint_kind' : MetaInfoMap(startpage),
        'constraint_parameters' : MetaInfoMap(startpage),
        'number_of_atoms_per_constraint' : MetaInfoMap(startpage),
        'number_of_constraints' : MetaInfoMap(startpage)
        }

    # section_interaction of section_topology
    interaction = {
        'interaction_atoms' : MetaInfoMap(startpage),
        'interaction_kind' : MetaInfoMap(startpage),
        'interaction_parameters' : MetaInfoMap(startpage),
        'number_of_atoms_per_interaction' : MetaInfoMap(startpage),
        'number_of_interactions' : MetaInfoMap(startpage)
        }

    # -------------------------------------------------------------
    #   Definitions of meta data values for section_molecule_type
    # -------------------------------------------------------------
    # section_molecule_type of section_topology
    mol_type = {
        'molecule_type_name' : MetaInfoMap(startpage),
        'number_of_atoms_in_molecule' : MetaInfoMap(startpage),
        'settings_atom_in_molecule' : MetaInfoMap(startpage)
        }

    # section_molecule_constraint of section_molecule_type
    mol_constraint = {
        'molecule_constraint_atoms' : MetaInfoMap(startpage),
        'molecule_constraint_kind' : MetaInfoMap(startpage),
        'molecule_constraint_parameters' : MetaInfoMap(startpage),
        'number_of_atoms_per_molecule_constraint' : MetaInfoMap(startpage),
        'number_of_molecule_constraints' : MetaInfoMap(startpage)
        }

    # section_molecule_interaction of section_molecule_type
    mol_interaction = {
        'molecule_interaction_atoms' : MetaInfoMap(startpage),
        'molecule_interaction_kind' : MetaInfoMap(startpage),
        'molecule_interaction_parameters' : MetaInfoMap(startpage),
        'number_of_atoms_per_molecule_interaction' : MetaInfoMap(startpage),
        'number_of_molecule_interactions' : MetaInfoMap(startpage)
        }

    # section_atom_in_molecule of section_molecule_type
    atom_in_mol = {
        'atom_in_molecule_charge' : MetaInfoMap(startpage),
        'atom_in_molecule_name' : MetaInfoMap(startpage),
        'atom_in_molecule_to_atom_type_ref' : MetaInfoMap(startpage)
        }


    if defname == 'system':
        dictionary = systemDict
    elif defname == 'topology':
        dictionary = topology
    elif defname == 'singleconfcalc':
        dictionary = singleconfcalc
    elif defname == 'frameseq':
        dictionary = frameseq
    elif defname == 'frameseqend':
        dictionary = frameseqend
    elif defname == 'atom_type':
        dictionary = atom_type
    elif defname == 'molecule_type':
        dictionary = moltypeDict
    elif defname == 'interaction':
        dictionary = interaction
    elif defname == 'sampling':
        dictionary = sampling
    elif defname == 'singlevdw':
        dictionary = singlevdw
    else:
        dictionary = singleconfcalclist
    return MapDictionary(dictionary)

