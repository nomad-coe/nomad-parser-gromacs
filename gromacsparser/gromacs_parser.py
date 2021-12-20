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
import logging
import re
import datetime

import panedr
try:
    import MDAnalysis
except Exception:
    logging.warn('Required module MDAnalysis not found.')
    MDAnalysis = False

from .metainfo import m_env
from nomad.units import ureg
from nomad.parsing.parser import FairdiParser

from nomad.parsing.file_parser import TextParser, Quantity, FileParser
from nomad.datamodel.metainfo.common_dft import Run, SamplingMethod, System,\
    EnergyContribution, SingleConfigurationCalculation, Topology, Interaction, Method
from .metainfo.gromacs import x_gromacs_section_control_parameters, x_gromacs_section_input_output_files

MOL = 6.022140857e+23


class GromacsLogParser(TextParser):
    def __init__(self):
        super().__init__(None)

    def init_quantities(self):
        def str_to_header(val_in):
            val = [v.split(':', 1) for v in val_in.strip().split('\n')]
            return {v[0].strip(): v[1].strip() for v in val if len(v) == 2}

        def str_to_input_parameters(val_in):
            re_array = re.compile(r'\s*([\w\-]+)\[[\d ]+\]\s*=\s*\{*(.+)')
            re_scalar = re.compile(r'\s*([\w\-]+)\s*=\s*(.+)')
            parameters = dict()
            val = val_in.strip().split('\n')
            for val_n in val:
                val_scalar = re_scalar.match(val_n)
                if val_scalar:
                    parameters[val_scalar.group(1)] = val_scalar.group(2)
                    continue
                val_array = re_array.match(val_n)
                if val_array:
                    parameters.setdefault(val_array.group(1), [])
                    value = [float(v) for v in val_array.group(2).rstrip('}').split(',')]
                    parameters[val_array.group(1)].append(value[0] if len(value) == 1 else value)
            return parameters

        def str_to_energies(val_in):
            energy_keys_re = re.compile(r'(.+?)(?:  |\Z| P)')
            keys = []
            values = []
            energies = dict()
            for val in val_in.strip().split('\n'):
                val = val.strip()
                if val[0].isalpha():
                    keys = [k.strip() for k in energy_keys_re.findall(val)]
                    keys = ['P%s' % k if k.startswith('res') else k for k in keys if k]
                else:
                    values = val.split()
                    for n, key in enumerate(keys):

                        if key == 'Temperature':
                            energies[key] = float(values[n]) * ureg.kelvin
                        elif key.startswith('Pres'):
                            key = key.rstrip(' (bar)')
                            energies[key] = float(values[n]) * ureg.bar
                        else:
                            energies[key] = float(values[n]) / MOL * ureg.kJ
            return energies

        def str_to_step_info(val_in):
            val = val_in.strip().split('\n')
            keys = val[0].split()
            values = [float(v) for v in val[1].split()]
            return {key: values[n] for n, key in enumerate(keys)}

        thermo_quantities = [
            Quantity(
                'energies',
                r'Energies \(kJ/mol\)\s*([\s\S]+?)\n\n',
                str_operation=str_to_energies, convert=False),
            Quantity(
                'step_info',
                r'(Step.+\n[\d\.\- ]+)',
                str_operation=str_to_step_info, convert=False)]

        self._quantities = [
            Quantity('time_start', r'Log file opened on (.+)', flatten=False),
            Quantity(
                'host_info',
                r'Host:\s*(\S+)\s*pid:\s*(\d+)\s*rank ID:\s*(\d+)\s*number of ranks:\s*(\d*)'),
            Quantity('module_version', r'GROMACS:\s*(.+?),\s*VERSION\s*(\S+)', flatten=False),
            Quantity('execution_path', r'Executable:\s*(.+)'),
            Quantity('working_path', r'Data prefix:\s*(.+)'),
            # TODO cannot understand treatment of the command line in the old parser
            Quantity(
                'header',
                r'(GROMACS version:[\s\S]+?)\n\n', str_operation=str_to_header),
            Quantity(
                'input_parameters',
                r'Input Parameters:\s*([\s\S]+?)\n\n', str_operation=str_to_input_parameters),
            Quantity(
                'step',
                r'(Step\s*Time[\s\S]+?Energies[\s\S]+?\n\n)',
                repeats=True, sub_parser=TextParser(quantities=thermo_quantities)),
            Quantity(
                'averages',
                r'A V E R A G E S  ====>([\s\S]+?\n\n\n)',
                sub_parser=TextParser(quantities=thermo_quantities)),
            Quantity('time_end', r'Finished \S+ on rank \d+ (.+)', flatten=False)]

    def get_pbc(self):
        pbc = self.get('input_parameters', {}).get('pbc', 'xyz')
        return ['x' in pbc, 'y' in pbc, 'z' in pbc]

    def get_sampling_settings(self):
        input_parameters = self.get('input_parameters', {})
        integrator = input_parameters.get('integrator', 'md').lower()
        if integrator in ['l-bfgs', 'cg', 'steep']:
            sampling_method = 'geometry_optimization'
        elif integrator in ['bd']:
            sampling_method = 'langevin_dynamics'
        else:
            sampling_method = 'molecular_dynamics'

        ensemble_type = 'NVE' if sampling_method == 'molecular_dynamics' else None
        tcoupl = input_parameters.get('tcoupl', 'no').lower()
        if tcoupl != 'no':
            ensemble_type = 'NVT'
            pcoupl = input_parameters.get('pcoupl', 'no').lower()
            if pcoupl != 'no':
                ensemble_type = 'NPT'

        return dict(
            sampling_method=sampling_method, integrator_type=integrator,
            ensemble_type=ensemble_type)

    def get_tpstat_settings(self):
        input_parameters = self.get('input_parameters', {})
        target_T = input_parameters.get('ref-t', 0) * ureg.kelvin

        thermostat_type = None
        tcoupl = input_parameters.get('tcoupl', 'no').lower()
        if tcoupl != 'no':
            thermostat_type = 'Velocity Rescaling' if tcoupl == 'v-rescale' else tcoupl.title()

        thermostat_tau = input_parameters.get('tau-t', 0) * ureg.ps

        # TODO infer langevin_gamma [s] from bd_fric
        # bd_fric = self.get('bd-fric', 0, unit='amu/ps')
        langevin_gamma = None

        target_P = input_parameters.get('ref-p', 0) * ureg.bar
        # if P is array e.g. for non-isotropic pressures, get average since metainfo is float
        if hasattr(target_P, 'shape'):
            target_P = np.average(target_P)

        barostat_type = None
        pcoupl = input_parameters.get('pcoupl', 'no').lower()
        if pcoupl != 'no':
            barostat_type = pcoupl.title()

        barostat_tau = input_parameters.get('tau-p', 0) * ureg.ps

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
            try:
                self._file_handler = panedr.edr_to_df(self.mainfile)
            except Exception:
                self.logger.error('Error reading edr file.')

        return self._file_handler

    def parse(self, key):
        if self.fileedr is None:
            return

        val = self.fileedr.get(key, None)
        if self._results is None:
            self._results = dict()

        if val is not None:
            val = np.asarray(val)
        if key in self._energy_keys:
            val = val / MOL * ureg.kJ
        elif key in self._temperature_keys:
            val = val * ureg.kelvin
        elif key in self._pressure_keys:
            val = val * ureg.bar
        elif key in self._time_keys:
            val = val * ureg.ps

        self._results[key] = val

    def keys(self):
        return list(self.fileedr.keys())

    @property
    def length(self):
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
            try:
                args = [f for f in [self.trajectory_file] if f is not None]
                self._file_handler = MDAnalysis.Universe(self.mainfile, *args)
            except Exception:
                self.logger.error('Error setting up MDAnalysis.')
        return self._file_handler

    def parse(self, key):
        if self._results is None:
            self._results = dict()

        if self.universe is None:
            return

        if not MDAnalysis:
            return

        atoms = list(self.universe.atoms)
        try:
            trajectory = self.universe.trajectory
        except Exception:
            trajectory = []

        unit = None
        val = None
        if key == 'timestep':
            val = trajectory.dt
            unit = ureg.ps
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
            unit = ureg.angstrom
        elif key == 'velocities':
            val = [traj.velocities if traj.has_velocities else None for traj in trajectory]
            unit = ureg.angstrom / ureg.ps
        elif key == 'forces':
            val = [traj.forces / MOL if traj.has_forces else None for traj in trajectory]
            unit = ureg.kJ / ureg.angstrom
        elif key == 'cell':
            val = [traj.triclinic_dimensions for traj in trajectory]
            unit = ureg.angstrom

        if unit is not None:
            if isinstance(val, list):
                val = [v * unit if v is not None else v for v in val]
            else:
                val = val * unit

        self._results[key] = val


class GromacsParser(FairdiParser):
    def __init__(self):
        super().__init__(
            name='parsers/gromacs', code_name='Gromacs', code_homepage='http://www.gromacs.org/',
            domain='dft', mainfile_contents_re=r'gmx mdrun, (VERSION|version)')
        self._metainfo_env = m_env
        self.log_parser = GromacsLogParser()
        self.traj_parser = MDAnalysisParser()
        self.energy_parser = GromacsEDRParser()
        self._metainfo_mapping = {
            'LJ (SR)': 'Leonard-Jones', 'Coulomb (SR)': 'coulomb',
            'Potential': 'potential', 'Kinetic En.': 'kinetic'}

    def get_gromacs_file(self, ext):
        files = [d for d in self._gromacs_files if d.endswith(ext)]

        if len(files) == 1:
            return os.path.join(self._maindir, files[0])

        # we assume that the file has the same basename as the log file e.g.
        # out.log would correspond to out.tpr and out.trr and out.edr
        for f in files:
            if f.rsplit('.', 1)[0] == self._basename:
                return os.path.join(self._maindir, f)

        for f in files:
            if f.rsplit('.', 1)[0].startswith(self._basename):
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
                if f.rsplit('.', 1)[0] == reff.rsplit('.', 1)[0]:
                    count += 1
            if count == 1:
                return os.path.join(self._maindir, f)
            counts.append(count)

        return os.path.join(self._maindir, files[counts.index(min(counts))])

    def parse_thermodynamic_data(self):
        sec_run = self.archive.section_run[-1]

        forces = self.traj_parser.get('forces')
        for n, forces_n in enumerate(forces):
            sec_scc = sec_run.m_create(SingleConfigurationCalculation)
            sec_scc.atom_forces = forces_n
            sec_scc.single_configuration_calculation_to_system_ref = sec_run.section_system[n]
            sec_scc.single_configuration_to_calculation_method_ref = sec_run.section_method[-1]

        # get it from edr file
        if len(self.energy_parser.keys()) > 0:
            thermo_data = self.energy_parser
            n_evaluations = self.energy_parser.length
        else:
            # try to get it from log file
            steps = self.log_parser.get('step', [])
            thermo_data = dict()
            for n, step in enumerate(steps):
                n = int(step.get('step_info', {}).get('Step', n))
                if step.energies is None:
                    continue
                keys = step.energies.keys()
                for key in keys:
                    thermo_data.setdefault(key, [None] * len(forces))
                    thermo_data[key][n] = step.energies.get(key)
                info = step.get('step_info', {})
                thermo_data.setdefault('Time', [None] * len(forces))
                thermo_data['Time'][n] = info.get('Time', None)
            n_evaluations = n + 1

        create_scc = False
        if len(forces) != n_evaluations:
            self.logger.warn(
                'Mismatch in number of calculations and number of thermodynamic'
                'evaluations, will create new sections')
            create_scc = True

        timestep = self.traj_parser.get('timestep')
        if timestep is None:
            timestep = self.log_parser.get('input_parameters', {}).get('dt', 1.0) * ureg.ps

        # TODO add other energy contributions, properties
        energy_keys = ['LJ (SR)', 'Coulomb (SR)', 'Potential', 'Kinetic En.']

        for n in range(n_evaluations):
            if create_scc:
                sec_scc = sec_run.m_create(SingleConfigurationCalculation)
            else:
                sec_scc = sec_run.section_single_configuration_calculation[n]

            for key in thermo_data.keys():
                val = thermo_data.get(key)[n]
                if val is None:
                    continue

                if key == 'Total Energy':
                    sec_scc.energy_total = val

                elif key == 'Pressure':
                    sec_scc.pressure = val

                elif key == 'Temperature':
                    sec_scc.temperature = val

                elif key == 'Time':
                    sec_scc.time_step = int((val / timestep).magnitude)
                if key in energy_keys:
                    sec_energy = sec_scc.m_create(EnergyContribution)
                    sec_energy.energy_contribution_kind = self._metainfo_mapping[key]
                    sec_energy.energy_contribution_value = val

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

        sec_run.section_method[-1].method_to_topology_ref = sec_topology

    def parse_system(self):
        sec_run = self.archive.section_run[-1]

        n_frames = self.traj_parser.get('n_frames', 0)

        pbc = self.log_parser.get_pbc()
        for n in range(n_frames):
            positions = self.traj_parser.get_positions(n)
            sec_system = sec_run.m_create(System)
            if positions is None:
                continue

            sec_system.number_of_atoms = self.traj_parser.get_n_atoms(n)
            sec_system.configuration_periodic_dimensions = pbc
            sec_system.simulation_cell = self.traj_parser.get_cell(n)
            sec_system.lattice_vectors = self.traj_parser.get_cell(n)
            sec_system.atom_labels = self.traj_parser.get_atom_labels(n)
            sec_system.atom_positions = positions

            velocities = self.traj_parser.get_velocities(n)
            if velocities is not None:
                sec_system.atom_velocities = velocities

    def parse_method(self):
        self.archive.section_run[-1].m_create(Method)
        # TODO what to put here

    def parse_sampling_method(self):
        sec_run = self.archive.section_run[-1]
        sec_sampling_method = sec_run.m_create(SamplingMethod)

        sampling_settings = self.log_parser.get_sampling_settings()

        sec_sampling_method.sampling_method = sampling_settings['sampling_method']
        sec_sampling_method.ensemble_type = sampling_settings['ensemble_type']
        sec_sampling_method.x_gromacs_integrator_type = sampling_settings['integrator_type']

        input_parameters = self.log_parser.get('input_parameters', {})
        timestep = input_parameters.get('dt', 0)
        sec_sampling_method.x_gromacs_integrator_dt = timestep

        nsteps = input_parameters.get('nsteps', 0)
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
        input_parameters = self.log_parser.get('input_parameters', {})
        input_parameters.update(self.log_parser.get('header', {}))
        for key, val in input_parameters.items():
            key = 'x_gromacs_inout_control_%s' % key.replace('-', '').replace(' ', '_').lower()
            if hasattr(sec_control_parameters, key):
                val = str(val) if not isinstance(val, np.ndarray) else val
                setattr(sec_control_parameters, key, val)

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
        self.logger = logging.getLogger(__name__) if logger is None else logger
        self._maindir = os.path.dirname(self.filepath)
        self._gromacs_files = os.listdir(self._maindir)
        self._basename = os.path.basename(filepath).rsplit('.', 1)[0]

        self.init_parser()

        sec_run = self.archive.m_create(Run)

        sec_run.program_name = 'GROMACS'

        header = self.log_parser.get('header', {})
        sec_run.program_version = header.get('GROMACS version', 'unknown').lstrip('VERSION ')

        for key in ['start', 'end']:
            time = self.log_parser.get('time_%s')
            if time is None:
                continue
            setattr(sec_run, 'time_run_date_%s' % key, datetime.datetime.strptime(
                time, '%a %b %d %H:%M:%S %Y').timestamp())

        host_info = self.log_parser.get('host_info')
        if host_info is not None:
            sec_run.x_gromacs_program_execution_host = host_info[0]
            sec_run.x_gromacs_parallel_task_nr = host_info[1]
            sec_run.x_gromacs_number_of_tasks = host_info[2]

        self.parse_method()

        self.parse_sampling_method()

        topology_file = self.get_gromacs_file('tpr')
        # I have no idea if output trajectory file can be specified in input
        trajectory_file = self.get_gromacs_file('trr')
        self.traj_parser.mainfile = topology_file
        self.traj_parser.trajectory_file = trajectory_file

        self.parse_topology()

        self.parse_system()

        # TODO read also from ene
        edr_file = self.get_gromacs_file('edr')
        self.energy_parser.mainfile = edr_file

        self.parse_thermodynamic_data()

        self.parse_input()
