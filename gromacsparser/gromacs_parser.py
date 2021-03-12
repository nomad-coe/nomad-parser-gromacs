#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD.
# See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import os
import numpy as np
import pint
import logging

import panedr
try:
    import MDAnalysis
except Exception:
    logging.warn('Required module MDAnalysis not found.')
    MDAnalysis = False

from .metainfo import m_env
from nomad.parsing.parser import FairdiParser

from nomad.parsing.file_parser import TextParser, Quantity, FileParser
from nomad.datamodel.metainfo.common_dft import Run, SamplingMethod, System,\
    EnergyContribution, SingleConfigurationCalculation, Topology, Interaction
from .metainfo.gromacs import x_gromacs_section_control_parameters, x_gromacs_section_input_output_files

MOL = 6.022140857e+23


class GromacsLogParser(TextParser):
    def __init__(self):
        run_control = [
            'integrator', 'tinit', 'dt', 'nsteps', 'init-step', 'simulation-part',
            'comm-mode', 'nstcomm', 'comm-grps']

        langevin_dynamics = ['bd-fric', 'ld-seed']

        energy_minimization = ['emtol', 'emstep', 'nstcgsteep', 'nbfgscorr']

        shell_molecular_dynamics = ['niter', 'fcstep']

        test_particle_insertion = ['rtpi']

        output_control = [
            'nstxout', 'nstvout', 'nstfout', 'nstlog', 'nstcalcenergy', 'nstenergy',
            'nstxout-compressed', 'compressed-x-precision', 'compressed-x-grps',
            'energygrps']

        neigbor_searching = [
            'cutoff-scheme', 'nstlist', 'pbc', 'periodic-molecules', 'verlet-buffer-tolerance',
            'rlist']

        electrostatics = [
            'coulombtype', 'coulomb-modifier', 'rcoulomb-switch', 'rcoulomb', 'epsilon-r',
            'epsilon-rf']

        van_der_waals = ['vdwtype', 'vdw-modifier', 'rvdw-switch', 'rvdw', 'DispCorr']

        tables = ['table-extension', 'energygrp-table']

        ewald = [
            'fourierspacing', 'fourier-nx', 'fourier-ny', 'fourier-nz', 'pme-order',
            'ewald-rtol', 'ewald-rtol-lj', 'lj-pme-comb-rule', 'ewald-geometry',
            'epsilon-surface']

        temperature_coupling = [
            'tcoupl', 'nsttcouple', 'nh-chain-length', 'print-nose-hoover-chain-variables',
            'tc-grps', 'tau-t', 'ref-t']

        pressure_coupling = [
            'pcoupl', 'pcoupltype', 'nstpcouple', 'tau-p', 'compressibility', 'ref-p',
            'refcoord-scaling']

        simulated_annealing = [
            'annealing', 'annealing-npoints', 'annealing-time', 'annealing-temp']

        velocity_generation = ['gen-vel', 'gen-temp', 'gen-seed']

        bonds = [
            'constraints', 'constraint-algorithm', 'continuation', 'shake-tol', 'lincs-order',
            'lincs-iter', 'lincs-warnangle', 'morse']

        energy_group_exclusions = ['energygrp-excl']

        walls = [
            'nwall', 'wall-atomtype', 'wall-type', 'table', 'wall-r-linpot', 'wall-density',
            'wall-ewald-zfac']

        com_pulling = [
            'pull', 'pull-cylinder-r', 'pull-constr-tol', 'pull-print-com', 'pull-print-ref-value',
            'pull-print-components', 'pull-nstxout', 'pull-nstfout', 'pull-pbc-ref-prev-step-com',
            'pull-xout-average', 'pull-ngroups', 'pull-ncoords', 'pull-group1-name',
            'pull-group1-weights', 'pull-group1-pbcatom', 'pull-coord1-type',
            'pull-coord1-potential-provider', 'pull-coord1-geometry', 'pull-coord1-groups',
            'pull-coord1-dim', 'pull-coord1-origin', ' pull-coord1-vec', 'pull-coord1-start',
            'pull-coord1-init', 'pull-coord1-rate', 'pull-coord1-k', 'pull-coord1-kB']

        awh_adaptive_biasing = [
            'awh', 'awh-potential', 'awh-share-multisim', 'awh-seed', 'awh-nstout',
            'awh-nstsample', 'awh-nsamples-update', 'awh-nbias', 'awh1-error-init',
            'awh1-growth', 'awh1-equilibrate-histogram', 'awh1-target', 'awh1-target-beta-scaling',
            'awh1-target-cutoff', 'awh1-user-data', 'awh1-share-group', 'awh1-ndim',
            'awh1-dim1-coord-provider', 'awh1-dim1-coord-index', 'awh1-dim1-force-constant',
            'awh1-dim1-start', 'awh1-dim1-end', 'awh1-dim1-diffusion', 'awh1-dim1-cover-diameter']

        enforced_rotation = [
            'rotation', 'rot-ngroups', 'rot-group0', 'rot-type0', 'rot-massw0', 'rot-vec0',
            'rot-pivot0', 'rot-rate0', 'rot-k0', 'rot-slab-dist0', 'rot-min-gauss0',
            'rot-eps0', 'rot-fit-method0', 'rot-potfit-nsteps0', 'rot-potfit-step0',
            'rot-nstrout', 'rot-nstsout']

        nmr_refinement = [
            'disre', 'disre-weighting', 'disre-mixed', 'disre-fc', 'disre-tau', 'nstdisreout',
            'orire', 'orire-fc', 'orire-tau', 'orire-fitgrp', 'nstorireout']

        free_energy_calculations = [
            'free-energy', 'expanded', 'init-lambda', 'delta-lambda', 'init-lambda-state',
            'fep-lambdas', 'coul-lambdas', 'vdw-lambdas', 'bonded-lambdas', 'restraint-lambdas',
            'mass-lambdas', 'temperature-lambdas', 'calc-lambda-neighbors', 'sc-alpha',
            'sc-r-power', 'sc-coul', 'sc-power', 'sc-sigma', 'couple-moltype', 'couple-lambda0',
            'couple-lambda1', 'couple-intramol', 'nstdhdl', 'dhdl-derivatives',
            'dhdl-print-energy', 'separate-dhdl-file', 'dh-hist-size', 'dh-hist-spacing']

        expanded_ensemble_calculations = [
            'nstexpanded', 'lmc-stats', 'lmc-mc-move', 'lmc-seed', 'mc-temperature', 'wl-ratio',
            'wl-scale', 'init-wl-delta', 'wl-oneovert', 'lmc-repeats', 'lmc-gibbsdelta',
            'lmc-forced-nstart', 'nst-transition-matrix', 'symmetrized-transition-matrix',
            'mininum-var-min', 'init-lambda-weights', 'lmc-weights-equil', 'simulated-tempering',
            'sim-temp-low', 'sim-temp-high', 'simulated-tempering-scaling']

        non_equilibrium_md = [
            'acc-grps', 'accelerate', 'freezegrps', 'freezedim', 'cos-acceleration',
            'deform']

        electric_fields = ['electric-field-x', 'electric-field-y', 'electric-field-z']

        mixed_quantum_classical_molecular_dynamics = [
            'QMMM', 'QMMM-grps', 'QMMMscheme', 'QMmethod', 'QMbasis', 'QMcharge', 'QMmult',
            'CASorbitals', 'CASelectrons', 'SH']

        computational_electrophysiology = [
            'swapcoords', 'swap-frequency', 'split-group0', 'split-group1', 'massw-split0',
            'massw-split1', 'solvent-group', 'coupl-steps', 'iontypes', 'iontype0-name',
            'iontype0-in-A', 'iontype0-in-B', 'bulk-offsetA', 'bulk-offsetB', 'threshold',
            'cyl0-r', 'cyl0-up', 'cyl0-down', 'cyl1-r', 'cyl1-up', 'cyl1-down']

        density_guided_simulations = [
            'density-guided-simulation-active', 'density-guided-simulation-group',
            'density-guided-simulation-similarity-measure', 'density-guided-simulation-atom-spreading-weight',
            'density-guided-simulation-force-constant', 'density-guided-simulation-gaussian-transform-spreading-width',
            'density-guided-simulation-gaussian-transform-spreading-range-in-multiples-of-width',
            'density-guided-simulation-reference-density-filename', 'density-guided-simulation-nst',
            'density-guided-simulation-normalize-densities', 'density-guided-simulation-adaptive-force-scaling',
            'density-guided-simulation-adaptive-force-scaling-time-constant']

        self._commands = ['version']

        self._commands +=\
            run_control + langevin_dynamics + energy_minimization + shell_molecular_dynamics +\
            test_particle_insertion + output_control + neigbor_searching + electrostatics +\
            van_der_waals + tables + ewald + temperature_coupling + pressure_coupling +\
            simulated_annealing + velocity_generation + bonds + energy_group_exclusions +\
            walls + com_pulling + awh_adaptive_biasing + enforced_rotation + nmr_refinement +\
            free_energy_calculations + expanded_ensemble_calculations + non_equilibrium_md +\
            electric_fields + mixed_quantum_classical_molecular_dynamics + computational_electrophysiology +\
            density_guided_simulations

        super().__init__(None)

    def init_quantities(self):
        def str_op(val):
            val = val.strip().replace(',', '').replace('{', '').replace('}', '')
            val = val.split()
            val = val[0] if len(val) == 1 else val
            return val

        self._quantities = []
        for command in self._commands:
            self._quantities.append(
                Quantity(
                    command, r'\s*%s\s*\[*\s*\d*\s*\]*\s*[=:]+([\s\S]*?\n)' % command,
                    str_operation=str_op, repeats=False)
            )

    def get_pbc(self):
        pbc = self.get('pbc', 'xyz')
        return ['x' in pbc, 'y' in pbc, 'z' in pbc]

    def get_sampling_settings(self):
        integrator = self.get('integrator', 'md').lower()
        if integrator in ['l-bfgs', 'cg', 'steep']:
            sampling_method = 'geometry_optimization'
        elif integrator in ['bd']:
            sampling_method = 'langevin_dynamics'
        else:
            sampling_method = 'molecular_dynamics'

        ensemble_type = 'NVE' if sampling_method == 'molecular_dynamics' else None
        tcoupl = self.get('tcoupl', 'no').lower()
        if tcoupl != 'no':
            ensemble_type = 'NVT'
            pcoupl = self.get('pcoupl', 'no').lower()
            if pcoupl != 'no':
                ensemble_type = 'NPT'

        return dict(
            sampling_method=sampling_method, integrator_type=integrator,
            ensemble_type=ensemble_type)

    def get_tpstat_settings(self):
        target_T = self.get('ref-t', 0, unit='K')

        thermostat_type = None
        tcoupl = self.get('tcoupl', 'no').lower()
        if tcoupl != 'no':
            thermostat_type = 'Velocity Rescaling' if tcoupl == 'v-rescale' else tcoupl.title()

        thermostat_tau = self.get('tau-t', 0, unit='ps')

        # TODO infer langevin_gamma [s] from bd_fric
        # bd_fric = self.get('bd-fric', 0, unit='amu/ps')
        langevin_gamma = None

        target_P = self.get('ref-p', 0, unit='bar')
        # if P is array e.g. for non-isotropic pressures, get average since metainfo is float
        if hasattr(target_P, 'shape'):
            target_P = np.average(target_P)

        barostat_type = None
        pcoupl = self.get('pcoupl', 'no').lower()
        if pcoupl != 'no':
            barostat_type = pcoupl.title()

        barostat_tau = self.get('tau-p', 0, unit='ps')

        return dict(
            target_T=target_T, thermostat_type=thermostat_type, thermostat_tau=thermostat_tau,
            target_P=target_P, barostat_type=barostat_type, barostat_tau=barostat_tau,
            langevin_gamma=langevin_gamma)


class GromacsEDRParser(FileParser):
    def __init__(self):
        super().__init__(None)
        self._energy_keys = [
            'LJ (SR)', 'Coulomb (SR)', 'Potential', 'Kinetic En.', 'Total Energy',
            'Vir-XX', 'Vir-XY', 'Vir-XZ', 'Vir-YX', 'Vir-YY', 'Vir-YZ', 'Vir-ZX', 'Vir-ZY',
            'Vir-ZZ']
        self._pressure_keys = [
            'Pressure', 'Pres-XX', 'Pres-XY', 'Pres-XZ', 'Pres-YX', 'Pres-YY', 'Pres-YZ',
            'Pres-ZX', 'Pres-ZY', 'Pres-ZZ']
        self._temperature_keys = ['Temperature']
        self._time_keys = ['Time']

    @property
    def fileedr(self):
        if self._file_handler is None:
            self._file_handler = panedr.edr_to_df(self.mainfile)

        return self._file_handler

    def parse(self, key):
        val = self.fileedr.get(key, None)
        if self._results is None:
            self._results = dict()

        if val is not None:
            val = np.asarray(val)
        if key in self._energy_keys:
            val = pint.Quantity(val / MOL, 'kJ')
        elif key in self._temperature_keys:
            val = pint.Quantity(val, 'K')
        elif key in self._pressure_keys:
            val = pint.Quantity(val, 'bar')
        elif key in self._time_keys:
            val = pint.Quantity(val, 'ps')

        self._results[key] = val

    def get_keys(self):
        return list(self.fileedr.keys())

    def get_length(self):
        return self.fileedr.shape[0]


class MDAnalysisParser(FileParser):
    def __init__(self):
        super().__init__(None)

    @property
    def trajectory_file(self):
        return self._trajectory_file

    @trajectory_file.setter
    def trajectory_file(self, val):
        self._file_handler = None
        self._trajectory_file = val

    def get_interactions(self):
        interactions = self.get('interactions', None)

        if interactions is not None:
            return interactions

        interaction_types = ['angles', 'bonds', 'dihedrals', 'impropers']
        interactions = []
        for interaction_type in interaction_types:
            try:
                interaction = getattr(self.universe, interaction_type)
            except Exception:
                continue

            for i in range(len(interaction)):
                interactions.append(
                    (str(interaction[i].type), [interaction[i].value()]))

        self._results['interactions'] = interactions

        return interactions

    def get_n_atoms(self, frame_index):
        return self.get('n_atoms', [0] * frame_index)[frame_index]

    def get_cell(self, frame_index):
        return self.get('cell', [np.zeros((3, 3))] * frame_index)[frame_index]

    def get_atom_labels(self, frame_index):
        return self.get('atom_labels', None)

    def get_positions(self, frame_index):
        return self.get('positions', [None] * frame_index)[frame_index]

    def get_velocities(self, frame_index):
        return self.get('velocities', [None] * frame_index)[frame_index]

    def get_forces(self, frame_index):
        return self.get('forces', [None] * frame_index)[frame_index]

    @property
    def universe(self):
        if self._file_handler is None:
            args = [f for f in [self.trajectory_file] if f is not None]
            self._file_handler = MDAnalysis.Universe(self.mainfile, *args)
        return self._file_handler

    def parse(self, key):
        if not MDAnalysis:
            return

        if self._results is None:
            self._results = dict()

        atoms = list(self.universe.atoms)
        try:
            trajectory = self.universe.trajectory
        except Exception:
            trajectory = []

        unit = None
        val = None
        if key == 'timestep':
            val = trajectory.dt
            unit = 'ps'
        elif key == 'atom_labels':
            val = [
                MDAnalysis.topology.guessers.guess_atom_element(atom.name)
                for atom in atoms]
        elif key == 'n_atoms':
            val = [traj.n_atoms for traj in trajectory] if trajectory else [len(atoms)]
        elif key == 'n_frames':
            val = len(trajectory)
        elif key == 'positions':
            val = [traj.positions if traj.has_positions else None for traj in trajectory]
            unit = 'angstrom'
        elif key == 'velocities':
            val = [traj.velocities if traj.has_velocities else None for traj in trajectory]
            unit = 'angstrom/ps'
        elif key == 'forces':
            val = [traj.forces / MOL if traj.has_forces else None for traj in trajectory]
            unit = 'kJ/angstrom'
        elif key == 'cell':
            val = [traj.triclinic_dimensions for traj in trajectory]
            unit = 'angstrom'

        if unit is not None:
            if isinstance(val, list):
                val = [pint.Quantity(v, unit) if v is not None else v for v in val]
            else:
                val = pint.Quantity(val, unit)

        self._results[key] = val


class GromacsParser(FairdiParser):
    def __init__(self):
        super().__init__(
            name='parsers/gromacs', code_name='Gromacs', code_homepage='http://www.gromacs.org/',
            domain='dft', mainfile_contents_re=r'gmx mdrun, VERSION')
        self._metainfo_env = m_env
        self.log_parser = GromacsLogParser()
        self.traj_parser = MDAnalysisParser()
        self.energy_parser = GromacsEDRParser()

    def get_gromacs_file(self, ext):
        files = [d for d in self._gromacs_files if d.endswith(ext)]

        if len(files) == 1:
            return os.path.join(self._maindir, files[0])

        # we assume that the file has the same basename as the log file e.g.
        # out.log would correspond to out.tpr and out.trr and out.edr
        for f in files:
            if f.split('.')[0] == self._basename:
                return os.path.join(self._maindir, f)

        # if the files are all named differently, we guess that the one that does not
        # share the same basename would be file we are interested in
        # e.g. in a list of files out.log someout.log out.tpr out.trr another.tpr file.trr
        # we guess that the out.* files belong together and the rest that does not share
        # a basename would be grouped together
        counts = []
        for f in files:
            count = 0
            for reff in self._gromacs_files:
                if f.split('.')[0] == reff.split('.')[0]:
                    count += 1
            if count == 1:
                return os.path.join(self._maindir, f)
            counts.append(count)

        return os.path.join(self._maindir, files[counts.index(min(counts))])

    def parse_thermodynamic_data(self):
        sec_run = self.archive.section_run[-1]
        sec_sccs = sec_run.section_single_configuration_calculation

        energy_keys_mapping = {
            'LJ (SR)': 'Leonard-Jones', 'Coulomb (SR)': 'coulomb',
            'Potential': 'potential', 'Kinetic En.': 'kinetic'}

        n_evaluations = self.energy_parser.get_length()
        create_scc = True
        if sec_sccs:
            if len(sec_sccs) != n_evaluations:
                self.logger.warn(
                    '''Mismatch in number of calculations and number of thermodynamic
                    evaluations, will create new sections''',
                    data=dict(n_calculations=len(sec_sccs), n_evaluations=n_evaluations))

            else:
                create_scc = False

        timestep = self.traj_parser.get('timestep', 1.0, unit='ps')

        keys = self.energy_parser.get_keys()
        for n in range(n_evaluations):
            if create_scc:
                sec_scc = sec_run.m_create(SingleConfigurationCalculation)
            else:
                sec_scc = sec_sccs[n]

            for key in keys:
                val = self.energy_parser.get(key)[n]
                if key in energy_keys_mapping:
                    sec_energy = sec_scc.m_create(EnergyContribution)
                    sec_energy.energy_contibution_kind = energy_keys_mapping[key]
                    sec_energy.energy_contribution_value = val

                elif key == 'Total Energy':
                    sec_scc.energy_method_current = val

                elif key == 'Pressure':
                    sec_scc.pressure = val

                elif key == 'Temperature':
                    sec_scc.temperature = val

                elif key == 'Time':
                    sec_scc.time_step = int((val / timestep).magnitude)

    def parse_topology(self):
        sec_run = self.archive.section_run[-1]

        sec_topology = sec_run.m_create(Topology)
        try:
            n_atoms = self.traj_parser.get('n_atoms', [0])[0]
        except Exception:
            gro_file = self.get_gromacs_file('gro')
            self.traj_parser.mainfile = gro_file
            n_atoms = self.traj_parser.get('n_atoms', [0])[0]

        sec_topology.number_of_topology_atoms = n_atoms

        interactions = self.traj_parser.get_interactions()
        for interaction in interactions:
            if not interaction[0] or not interaction[1]:
                continue
            sec_interaction = sec_topology.m_create(Interaction)
            sec_interaction.interaction_kind = interaction[0]
            sec_interaction.interaction_parameters = interaction[1]

    def parse_system(self):
        sec_run = self.archive.section_run[-1]

        n_frames = self.traj_parser.get('n_frames', 0)

        create_scc = True
        sec_sccs = sec_run.section_single_configuration_calculation
        if sec_sccs:
            if len(sec_sccs) != n_frames:
                self.logger.warn(
                    '''Mismatch in number of calculations and number of structures,
                    will create new sections''',
                    data=dict(n_calculations=len(sec_sccs), n_structures=n_frames))

            else:
                create_scc = False

        pbc = self.log_parser.get_pbc()
        for n in range(n_frames):
            positions = self.traj_parser.get_positions(n)
            if positions is None:
                continue

            sec_system = sec_run.m_create(System)
            sec_system.number_of_atoms = self.traj_parser.get_n_atoms(n)
            sec_system.configuration_periodic_dimensions = pbc
            sec_system.simulation_cell = self.traj_parser.get_cell(n)
            sec_system.atom_labels = self.traj_parser.get_atom_labels(n)
            sec_system.atom_positions = positions

            velocities = self.traj_parser.get_velocities(n)
            if velocities is not None:
                sec_system.atom_velocities = velocities

            forces = self.traj_parser.get_forces(n)
            if forces is not None:
                if create_scc:
                    sec_scc = sec_run.m_create(SingleConfigurationCalculation)
                else:
                    sec_scc = sec_sccs[n]
                sec_scc.atom_forces = forces

    def parse_sampling_method(self):
        sec_run = self.archive.section_run[-1]
        sec_sampling_method = sec_run.m_create(SamplingMethod)

        sampling_settings = self.log_parser.get_sampling_settings()

        sec_sampling_method.sampling_method = sampling_settings['sampling_method']
        sec_sampling_method.ensemble_type = sampling_settings['ensemble_type']
        sec_sampling_method.x_gromacs_integrator_type = sampling_settings['integrator_type']

        timestep = self.log_parser.get('dt', 0)
        sec_sampling_method.x_gromacs_integrator_dt = timestep

        nsteps = self.log_parser.get('nsteps', 0)
        sec_sampling_method.x_gromacs_number_of_steps_requested = nsteps

        tp_settings = self.log_parser.get_tpstat_settings()

        target_T = tp_settings.get('target_T', None)
        if target_T is not None:
            sec_sampling_method.x_gromacs_thermostat_target_temperature = target_T
        thermostat_tau = tp_settings.get('thermostat_tau', None)
        if thermostat_tau is not None:
            sec_sampling_method.x_gromacs_thermostat_tau = thermostat_tau
        target_P = tp_settings.get('target_P', None)
        if target_P is not None:
            sec_sampling_method.x_gromacs_barostat_target_pressure = target_P
        barostat_tau = tp_settings.get('barostat_P', None)
        if barostat_tau is not None:
            sec_sampling_method.x_gromacs_barostat_tau = barostat_tau
        langevin_gamma = tp_settings.get('langevin_gamma', None)
        if langevin_gamma is not None:
            sec_sampling_method.x_gromacs_langevin_gamma = langevin_gamma

    def parse_input(self):
        sec_run = self.archive.section_run[-1]
        sec_input_output_files = sec_run.m_create(x_gromacs_section_input_output_files)

        topology_file = os.path.basename(self.traj_parser.mainfile)
        if topology_file.endswith('tpr'):
            sec_input_output_files.x_gromacs_inout_file_topoltpr = topology_file
        elif topology_file.endswith('gro'):
            sec_input_output_files.x_gromacs_inout_file_confoutgro = topology_file

        trajectory_file = os.path.basename(self.traj_parser.trajectory_file)
        sec_input_output_files.x_gromacs_inout_file_trajtrr = trajectory_file

        edr_file = os.path.basename(self.energy_parser.mainfile)
        sec_input_output_files.x_gromacs_inout_file_eneredr = edr_file

        sec_control_parameters = sec_run.m_create(x_gromacs_section_control_parameters)
        keys = self.log_parser.keys()
        for key in keys:
            val = self.log_parser.get(key)
            if val is None:
                continue
            key = 'x_gromacs_inout_control_%s' % key.replace('-', '').lower()
            if hasattr(sec_control_parameters, key):
                setattr(sec_control_parameters, key, str(val))

    def init_parser(self):
        self.log_parser.mainfile = self.filepath
        self.log_parser.logger = self.logger
        self.traj_parser.logger = self.logger
        self.energy_parser.logger = self.logger

    def reuse_parser(self, parser):
        self.log_parser.quantities = parser.log_parser.quantities

    def parse(self, filepath, archive, logger):
        self.filepath = os.path.abspath(filepath)
        self.archive = archive
        self.logger = logging if logger is None else logger
        self._maindir = os.path.dirname(self.filepath)
        self._gromacs_files = os.listdir(self._maindir)
        self._basename = os.path.basename(filepath).split('.')[0]

        self.init_parser()

        sec_run = self.archive.m_create(Run)

        version = self.log_parser.get('version', ['VERSION', 'unknown'])
        sec_run.program_name = 'GROMACS'
        sec_run.program_version = version[1]

        self.parse_sampling_method()

        topology_file = self.get_gromacs_file('tpr')
        trajectory_file = self.get_gromacs_file('trr')
        self.traj_parser.mainfile = topology_file
        self.traj_parser.trajectory_file = trajectory_file

        self.parse_topology()

        self.parse_system()

        edr_file = self.get_gromacs_file('edr')
        self.energy_parser.mainfile = edr_file

        self.parse_thermodynamic_data()

        self.parse_input()
