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

import argparse

def get_commandLineArguments(cmdline):
    """Parse command line arguments of GROMACS
       using argparse

    Returns:
        the dictionary of arguments
    """
    cmdlineparser = argparse.ArgumentParser()
    cmdlineparser.add_argument('-s', nargs='?', const='topol.tpr', default='topol.tpr',
        help='[<.tpr/.tpb/...>](Input)  Run input file: tpr tpb tpa')
    cmdlineparser.add_argument('-o', nargs='?', const='traj.trr', default='traj.trr',
        help='[<.trr/.cpt/...>](Output)  Full precision trajectory: trr cpt trj tng')
    cmdlineparser.add_argument('-x', nargs='?', const='traj_comp.xtc', default='traj_comp.xtc',
        help='[<.xtc/.tng>](Output, Optional) Compressed trajectory (tng format or portable xdr format)')
    cmdlineparser.add_argument('-cpi', nargs='?', const='state.cpt', default='state.cpt',
        help='[<.cpt>](Input, Optional) Checkpoint file')
    cmdlineparser.add_argument('-cpo', nargs='?', const='state.cpt', default='state.cpt',
        help='[<.cpt>](Output, Optional) Checkpoint file')
    cmdlineparser.add_argument('-c', nargs='?', const='confout.gro', default='confout.gro',
        help='[<.gro/.g96/...>](Output)  Structure file: gro g96 pdb brk ent esp')
    cmdlineparser.add_argument('-e', nargs='?', const='ener.edr', default='ener.edr',
        help='[<.edr>](Output)  Energy file')
    cmdlineparser.add_argument('-g', nargs='?', const='md.log', default='md.log',
        help='[<.log>](Output)  Log file')
    cmdlineparser.add_argument('-dhdl', nargs='?', const='dhdl.xvg', default='dhdl.xvg',
        help='[<.xvg>](Output, Optional) xvgr/xmgr file')
    cmdlineparser.add_argument('-field', nargs='?', const='field.xvg', default='field.xvg',
        help='[<.xvg>](Output, Optional) xvgr/xmgr file')
    cmdlineparser.add_argument('-table', nargs='?', const='table.xvg', default='table.xvg',
        help='[<.xvg>](Input, Optional) xvgr/xmgr file')
    cmdlineparser.add_argument('-tabletf', nargs='?', const='tabletf.xvg', default='tabletf.xvg',
        help='[<.xvg>](Input, Optional) xvgr/xmgr file')
    cmdlineparser.add_argument('-tablep', nargs='?', const='tablep.xvg', default='tablep.xvg',
        help='[<.xvg>](Input, Optional) xvgr/xmgr file')
    cmdlineparser.add_argument('-tableb', nargs='?', const='table.xvg', default='table.xvg',
        help='[<.xvg>](Input, Optional) xvgr/xmgr file')
    cmdlineparser.add_argument('-rerun', nargs='?', const='rerun.xtc', default='rerun.xtc',
        help='[<.xtc/.trr/...>](Input, Optional) Trajectory: xtc trr cpt trj gro g96 pdb tng')
    cmdlineparser.add_argument('-tpi', nargs='?', const='tpi.xvg', default='tpi.xvg',
        help='[<.xvg>](Output, Optional) xvgr/xmgr file')
    cmdlineparser.add_argument('-tpid', nargs='?', const='tpidist.xvg', default='tpidist.xvg',
        help='[<.xvg>](Output, Optional) xvgr/xmgr file')
    cmdlineparser.add_argument('-ei', nargs='?', const='sam.edi', default='sam.edi',
        help='[<.edi>](Input, Optional) ED sampling input')
    cmdlineparser.add_argument('-eo', nargs='?', const='edsam.xvg', default='edsam.xvg',
        help='[<.xvg>](Output, Optional) xvgr/xmgr file')
    cmdlineparser.add_argument('-devout', nargs='?', const='deviatie.xvg', default='deviatie.xvg',
        help='[<.xvg>](Output, Optional) xvgr/xmgr file')
    cmdlineparser.add_argument('-runav', nargs='?', const='runaver.xvg', default='runaver.xvg',
        help='[<.xvg>](Output, Optional) xvgr/xmgr file')
    cmdlineparser.add_argument('-px', nargs='?', const='pullx.xvg', default='pullx.xvg',
        help='[<.xvg>](Output, Optional) xvgr/xmgr file')
    cmdlineparser.add_argument('-pf', nargs='?', const='pullf.xvg', default='pullf.xvg',
        help='[<.xvg>](Output, Optional) xvgr/xmgr file')
    cmdlineparser.add_argument('-ro', nargs='?', const='rotation.xvg', default='rotation.xvg',
        help='[<.xvg>](Output, Optional) xvgr/xmgr file')
    cmdlineparser.add_argument('-ra', nargs='?', const='rotangles.log', default='rotangles.log',
        help='[<.log>](Output, Optional) Log file')
    cmdlineparser.add_argument('-rs', nargs='?', const='rotslabs.log', default='rotslabs.log',
        help='[<.log>](Output, Optional) Log file')
    cmdlineparser.add_argument('-rt', nargs='?', const='rottorque.log', default='rottorque.log',
        help='[<.log>](Output, Optional) Log file')
    cmdlineparser.add_argument('-mtx', nargs='?', const='nm.mtx', default='nm.mtx',
        help='[<.mtx>](Output, Optional) Hessian matrix')
    cmdlineparser.add_argument('-dn', nargs='?', const='dipole.ndx', default='dipole.ndx',
        help='[<.ndx>](Output, Optional) Index file')
    cmdlineparser.add_argument('-multidir', nargs='?', const='rundir', default='rundir',
        help='[<dir> [...]](Input, Optional) Run directory')
    cmdlineparser.add_argument('-membed', nargs='?', const='membed.dat', default='membed.dat',
        help='[<.dat>](Input, Optional) Generic data file')
    cmdlineparser.add_argument('-mp', nargs='?', const='membed.top', default='membed.top',
        help='[<.top>](Input, Optional) Topology file')
    cmdlineparser.add_argument('-mn', nargs='?', const='membed.ndx', default='membed.ndx',
        help='[<.ndx>](Input, Optional) Index file')
    cmdlineparser.add_argument('-if', nargs='?', const='imdforces.xvg', default='imdforces.xvg',
        help='[<.xvg>](Output, Optional) xvgr/xmgr file')
    cmdlineparser.add_argument('-swap', nargs='?', const='swapions.xvg', default='swapions.xvg',
        help='[<.xvg>](Output, Optional) xvgr/xmgr file')
    cmdlineparser.add_argument('-nice', nargs='?', const='0', default='0',
        help='<int>  Set the nicelevel')
    cmdlineparser.add_argument('-deffnm', nargs='?', const='', default='',
        help='<string>  Set the default filename for all file options')
    cmdlineparser.add_argument('-xvg', nargs='?', const='xmgrace', default='xmgrace',
        help='<enum>  xvg plot formatting: xmgrace, xmgr, none')
    cmdlineparser.add_argument('-dd', nargs='?', const='0 0 0', default='0 0 0',
        help='<vector> Domain decomposition grid, 0 is optimize')
    cmdlineparser.add_argument('-ddorder', nargs='?', const='interleave', default='interleave',
        help='<enum>  DD rank order: interleave, pp_pme, cartesian')
    cmdlineparser.add_argument('-npme', nargs='?', const='-1', default='-1',
        help='<int>  Number of separate ranks to be used for PME, -1 is guess')
    cmdlineparser.add_argument('-nt', nargs='?', const='0', default='0',
        help='<int>  Total number of threads to start (0 is guess)')
    cmdlineparser.add_argument('-ntmpi', nargs='?', const='0', default='0',
        help='<int>  Number of thread-MPI threads to start (0 is guess)')
    cmdlineparser.add_argument('-ntomp', nargs='?', const='0', default='0',
        help='<int>  Number of OpenMP threads per MPI rank to start (0 is guess)')
    cmdlineparser.add_argument('-ntomp_pme', nargs='?', const='0', default='0',
        help='<int>  Number of OpenMP threads per MPI rank to start (0 is -ntomp)')
    cmdlineparser.add_argument('-pin', nargs='?', const='auto', default='auto',
        help='<enum>  Set thread affinities: auto, on, off')
    cmdlineparser.add_argument('-pinoffset', nargs='?', const='0', default='0',
        help='<int>  The starting logical core number for pinning to cores; used to avoid pinning threads from different mdrun instances to the same core')
    cmdlineparser.add_argument('-pinstride', nargs='?', const='0', default='0',
        help='<int>  Pinning distance in logical cores for threads, use 0 to minimize the number of threads per physical core')
    cmdlineparser.add_argument('-gpu_id', nargs='?', const='', default='',
        help='<string>  List of GPU device id-s to use, specifies the per-node PP rank to GPU mapping')
    cmdlineparser.add_argument('-ddcheck', nargs='?', const=True, default='yes',
        help='  Check for all bonded interactions with DD')
    cmdlineparser.add_argument('-noddcheck', nargs='?', const=False, default='yes',
        help='  Check for all bonded interactions with DD')
    cmdlineparser.add_argument('-rdd', nargs='?', const='0', default='0',
        help='<real>  The maximum distance for bonded interactions with DD (nm), 0 is determine from initial coordinates')
    cmdlineparser.add_argument('-rcon', nargs='?', const='0', default='0',
        help='<real>  Maximum distance for P-LINCS (nm), 0 is estimate')
    cmdlineparser.add_argument('-dlb', nargs='?', const='auto', default='auto',
        help='<enum>  Dynamic load balancing (with DD): auto, no, yes')
    cmdlineparser.add_argument('-dds', nargs='?', const='0.8', default='0.8',
        help='<real>  Fraction in (0,1) by whose reciprocal the initial DD cell size will be increased in order to provide a margin in which dynamic load balancing can act while preserving the minimum cell size.')
    cmdlineparser.add_argument('-gcom', nargs='?', const='-1', default='-1',
        help='<int>  Global communication frequency')
    cmdlineparser.add_argument('-nb', nargs='?', const='auto', default='auto',
        help='<enum>  Calculate non-bonded interactions on: auto, cpu, gpu, gpu_cpu')
    cmdlineparser.add_argument('-nstlist', nargs='?', const='0', default='0',
        help='<int>  Set nstlist when using a Verlet buffer tolerance (0 is guess)')
    cmdlineparser.add_argument('-tunepme', nargs='?', const=True, default='yes',
        help='  Optimize PME load between PP/PME ranks or GPU/CPU')
    cmdlineparser.add_argument('-notunepme', nargs='?', const=False, default='yes',
        help='  Optimize PME load between PP/PME ranks or GPU/CPU')
    cmdlineparser.add_argument('-testverlet', nargs='?', const=True, default='no',
        help='  Test the Verlet non-bonded scheme')
    cmdlineparser.add_argument('-notestverlet', nargs='?', const=False, default='no',
        help='  Test the Verlet non-bonded scheme')
    cmdlineparser.add_argument('-v', nargs='?', const=True, default='no',
        help='  Be loud and noisy')
    cmdlineparser.add_argument('-nov', nargs='?', const=False, default='no',
        help='  Be loud and noisy')
    cmdlineparser.add_argument('-compact', nargs='?', const=True, default='yes',
        help='  Write a compact log file')
    cmdlineparser.add_argument('-nocompact', nargs='?', const=False, default='yes',
        help='  Write a compact log file')
    cmdlineparser.add_argument('-seppot', nargs='?', const=True, default='no',
        help='  Write separate V and dVdl terms for each interaction type and rank to the log file(s)')
    cmdlineparser.add_argument('-noseppot', nargs='?', const=False, default='no',
        help='  Write separate V and dVdl terms for each interaction type and rank to the log file(s)')
    cmdlineparser.add_argument('-pforce', nargs='?', const='-1', default='-1',
        help='<real>  Print all forces larger than this (kJ/mol nm)')
    cmdlineparser.add_argument('-reprod', nargs='?', const=True, default='no',
        help='  Try to avoid optimizations that affect binary reproducibility')
    cmdlineparser.add_argument('-noreprod', nargs='?', const=False, default='no',
        help='  Try to avoid optimizations that affect binary reproducibility')
    cmdlineparser.add_argument('-cpt', nargs='?', const='15', default='15',
        help='<real>  Checkpoint interval (minutes)')
    cmdlineparser.add_argument('-cpnum', nargs='?', const=True, default='no',
        help='  Keep and number checkpoint files')
    cmdlineparser.add_argument('-nocpnum', nargs='?', const=False, default='no',
        help='  Keep and number checkpoint files')
    cmdlineparser.add_argument('-append', nargs='?', const=True, default='yes',
        help='  Append to previous output files when continuing from checkpoint instead of adding the simulation part number to all file names')
    cmdlineparser.add_argument('-noappend', nargs='?', const=False, default='yes',
        help='  Append to previous output files when continuing from checkpoint instead of adding the simulation part number to all file names')
    cmdlineparser.add_argument('-nsteps', nargs='?', const='-2', default='-2',
        help='<int>  Run this number of steps, overrides .mdp file option')
    cmdlineparser.add_argument('-maxh', nargs='?', const='-1', default='-1',
        help='<real>  Terminate after 0.99 times this time (hours)')
    cmdlineparser.add_argument('-multi', nargs='?', const='0', default='0',
        help='<int>  Do multiple simulations in parallel')
    cmdlineparser.add_argument('-replex', nargs='?', const='0', default='0',
        help='<int>  Attempt replica exchange periodically with this period (steps)')
    cmdlineparser.add_argument('-nex', nargs='?', const='0', default='0',
        help='<int>  Number of random exchanges to carry out each exchange interval (N^3 is one suggestion). -nex zero or not specified gives neighbor replica exchange.')
    cmdlineparser.add_argument('-reseed', nargs='?', const='-1', default='-1',
        help='<int>  Seed for replica exchange, -1 is generate a seed')
    args = cmdlineparser.parse_args(str(cmdline).split())
    return vars(args)

