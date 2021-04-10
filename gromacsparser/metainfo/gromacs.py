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
import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference
)
from nomad.metainfo.legacy import LegacyDefinition

from nomad.datamodel.metainfo import common
from nomad.datamodel.metainfo import public

m_package = Package(
    name='gromacs_nomadmetainfo_json',
    description='None',
    a_legacy=LegacyDefinition(name='gromacs.nomadmetainfo.json'))


class x_gromacs_mdin_input_output_files(MCategory):
    '''
    Parameters of mdin belonging to x_gromacs_section_control_parameters.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='x_gromacs_mdin_input_output_files'))


class x_gromacs_mdin_control_parameters(MCategory):
    '''
    Parameters of mdin belonging to x_gromacs_section_control_parameters.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='x_gromacs_mdin_control_parameters'))


class x_gromacs_mdin_method(MCategory):
    '''
    Parameters of mdin belonging to section method.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='x_gromacs_mdin_method'))


class x_gromacs_mdout_single_configuration_calculation(MCategory):
    '''
    Parameters of mdout belonging to section_single_configuration_calculation.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='x_gromacs_mdout_single_configuration_calculation'))


class x_gromacs_mdout_method(MCategory):
    '''
    Parameters of mdin belonging to section method.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='x_gromacs_mdout_method'))


class x_gromacs_mdout_run(MCategory):
    '''
    Parameters of mdin belonging to settings run.
    '''

    m_def = Category(
        categories=[public.settings_run],
        a_legacy=LegacyDefinition(name='x_gromacs_mdout_run'))


class x_gromacs_mdin_run(MCategory):
    '''
    Parameters of mdin belonging to settings run.
    '''

    m_def = Category(
        categories=[public.settings_run],
        a_legacy=LegacyDefinition(name='x_gromacs_mdin_run'))


class x_gromacs_section_input_output_files(MSection):
    '''
    Section to store input and output file names
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_gromacs_section_input_output_files'))

    x_gromacs_inout_file_topoltpr = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs input topology file.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_file_topoltpr'))

    x_gromacs_inout_file_trajtrr = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs input trajectory file.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_file_trajtrr'))

    x_gromacs_inout_file_trajcompxtc = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs input compressed trajectory file.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_file_trajcompxtc'))

    x_gromacs_inout_file_statecpt = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs input coordinates and state file.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_file_statecpt'))

    x_gromacs_inout_file_confoutgro = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs output configuration file.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_file_confoutgro'))

    x_gromacs_inout_file_eneredr = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs output energies file.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_file_eneredr'))


class x_gromacs_section_control_parameters(MSection):
    '''
    Section to store the input and output control parameters
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_gromacs_section_control_parameters'))

    x_gromacs_inout_control_gromacs_version = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_gromacs_version'))

    x_gromacs_inout_control_precision = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_precision'))

    x_gromacs_inout_control_memory_model = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_memory_model'))

    x_gromacs_inout_control_mpi_library = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_mpi_library'))

    x_gromacs_inout_control_openmp_support = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_openmp_support'))

    x_gromacs_inout_control_gpu_support = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_gpu_support'))

    x_gromacs_inout_control_opencl_support = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_opencl_support'))

    x_gromacs_inout_control_invsqrt_routine = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_invsqrt_routine'))

    x_gromacs_inout_control_simd_instructions = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_simd_instructions'))

    x_gromacs_inout_control_fft_library = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_fft_library'))

    x_gromacs_inout_control_rdtscp_usage = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_rdtscp_usage'))

    x_gromacs_inout_control_cxx11_compilation = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_cxx11_compilation'))

    x_gromacs_inout_control_tng_support = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_tng_support'))

    x_gromacs_inout_control_tracing_support = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_tracing_support'))

    x_gromacs_inout_control_built_on = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_built_on'))

    x_gromacs_inout_control_built_by = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_built_by'))

    x_gromacs_inout_control_build_osarch = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_build_osarch'))

    x_gromacs_inout_control_build_cpu_vendor = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_build_cpu_vendor'))

    x_gromacs_inout_control_build_cpu_brand = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_build_cpu_brand'))

    x_gromacs_inout_control_build_cpu_family = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_build_cpu_family'))

    x_gromacs_inout_control_build_cpu_features = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_build_cpu_features'))

    x_gromacs_inout_control_c_compiler = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_c_compiler'))

    x_gromacs_inout_control_c_compiler_flags = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_c_compiler_flags'))

    x_gromacs_inout_control_cxx_compiler = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_cxx_compiler'))

    x_gromacs_inout_control_cxx_compiler_flags = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_cxx_compiler_flags'))

    x_gromacs_inout_control_boost_version = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_boost_version'))

    x_gromacs_inout_control_integrator = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_integrator'))

    x_gromacs_inout_control_tinit = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_tinit'))

    x_gromacs_inout_control_dt = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_dt'))

    x_gromacs_inout_control_nsteps = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nsteps'))

    x_gromacs_inout_control_initstep = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_initstep'))

    x_gromacs_inout_control_simulationpart = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_simulationpart'))

    x_gromacs_inout_control_commmode = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_commmode'))

    x_gromacs_inout_control_nstcomm = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nstcomm'))

    x_gromacs_inout_control_bdfric = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_bdfric'))

    x_gromacs_inout_control_ldseed = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_ldseed'))

    x_gromacs_inout_control_emtol = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_emtol'))

    x_gromacs_inout_control_emstep = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_emstep'))

    x_gromacs_inout_control_niter = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_niter'))

    x_gromacs_inout_control_fcstep = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_fcstep'))

    x_gromacs_inout_control_nstcgsteep = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nstcgsteep'))

    x_gromacs_inout_control_nbfgscorr = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nbfgscorr'))

    x_gromacs_inout_control_rtpi = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_rtpi'))

    x_gromacs_inout_control_nstxout = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nstxout'))

    x_gromacs_inout_control_nstvout = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nstvout'))

    x_gromacs_inout_control_nstfout = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nstfout'))

    x_gromacs_inout_control_nstlog = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nstlog'))

    x_gromacs_inout_control_nstcalcenergy = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nstcalcenergy'))

    x_gromacs_inout_control_nstenergy = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nstenergy'))

    x_gromacs_inout_control_nstxoutcompressed = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nstxoutcompressed'))

    x_gromacs_inout_control_compressedxprecision = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_compressedxprecision'))

    x_gromacs_inout_control_cutoffscheme = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_cutoffscheme'))

    x_gromacs_inout_control_nstlist = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nstlist'))

    x_gromacs_inout_control_nstype = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nstype'))

    x_gromacs_inout_control_pbc = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_pbc'))

    x_gromacs_inout_control_periodicmolecules = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_periodicmolecules'))

    x_gromacs_inout_control_verletbuffertolerance = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_verletbuffertolerance'))

    x_gromacs_inout_control_rlist = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_rlist'))

    x_gromacs_inout_control_rlistlong = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_rlistlong'))

    x_gromacs_inout_control_nstcalclr = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nstcalclr'))

    x_gromacs_inout_control_coulombtype = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_coulombtype'))

    x_gromacs_inout_control_coulombmodifier = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_coulombmodifier'))

    x_gromacs_inout_control_rcoulombswitch = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_rcoulombswitch'))

    x_gromacs_inout_control_rcoulomb = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_rcoulomb'))

    x_gromacs_inout_control_epsilonr = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_epsilonr'))

    x_gromacs_inout_control_epsilonrf = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_epsilonrf'))

    x_gromacs_inout_control_vdwtype = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_vdwtype'))

    x_gromacs_inout_control_vdwmodifier = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_vdwmodifier'))

    x_gromacs_inout_control_rvdwswitch = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_rvdwswitch'))

    x_gromacs_inout_control_rvdw = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_rvdw'))

    x_gromacs_inout_control_dispcorr = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_dispcorr'))

    x_gromacs_inout_control_tableextension = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_tableextension'))

    x_gromacs_inout_control_fourierspacing = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_fourierspacing'))

    x_gromacs_inout_control_fouriernx = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_fouriernx'))

    x_gromacs_inout_control_fourierny = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_fourierny'))

    x_gromacs_inout_control_fouriernz = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_fouriernz'))

    x_gromacs_inout_control_pmeorder = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_pmeorder'))

    x_gromacs_inout_control_ewaldrtol = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_ewaldrtol'))

    x_gromacs_inout_control_ewaldrtollj = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_ewaldrtollj'))

    x_gromacs_inout_control_ljpmecombrule = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_ljpmecombrule'))

    x_gromacs_inout_control_ewaldgeometry = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_ewaldgeometry'))

    x_gromacs_inout_control_epsilonsurface = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_epsilonsurface'))

    x_gromacs_inout_control_implicitsolvent = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_implicitsolvent'))

    x_gromacs_inout_control_gbalgorithm = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_gbalgorithm'))

    x_gromacs_inout_control_nstgbradii = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nstgbradii'))

    x_gromacs_inout_control_rgbradii = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_rgbradii'))

    x_gromacs_inout_control_gbepsilonsolvent = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_gbepsilonsolvent'))

    x_gromacs_inout_control_gbsaltconc = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_gbsaltconc'))

    x_gromacs_inout_control_gbobcalpha = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_gbobcalpha'))

    x_gromacs_inout_control_gbobcbeta = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_gbobcbeta'))

    x_gromacs_inout_control_gbobcgamma = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_gbobcgamma'))

    x_gromacs_inout_control_gbdielectricoffset = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_gbdielectricoffset'))

    x_gromacs_inout_control_saalgorithm = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_saalgorithm'))

    x_gromacs_inout_control_sasurfacetension = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_sasurfacetension'))

    x_gromacs_inout_control_tcoupl = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_tcoupl'))

    x_gromacs_inout_control_nsttcouple = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nsttcouple'))

    x_gromacs_inout_control_nhchainlength = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nhchainlength'))

    x_gromacs_inout_control_printnosehooverchainvariables = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_printnosehooverchainvariables'))

    x_gromacs_inout_control_pcoupl = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_pcoupl'))

    x_gromacs_inout_control_pcoupltype = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_pcoupltype'))

    x_gromacs_inout_control_nstpcouple = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nstpcouple'))

    x_gromacs_inout_control_taup = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_taup'))

    x_gromacs_inout_control_compressibility = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_compressibility'))

    x_gromacs_inout_control_compressibility0 = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_compressibility0'))

    x_gromacs_inout_control_compressibility1 = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_compressibility1'))

    x_gromacs_inout_control_compressibility2 = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_compressibility2'))

    x_gromacs_inout_control_refp = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_refp'))

    x_gromacs_inout_control_refp0 = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_refp0'))

    x_gromacs_inout_control_refp1 = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_refp1'))

    x_gromacs_inout_control_refp2 = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_refp2'))

    x_gromacs_inout_control_refcoordscaling = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_refcoordscaling'))

    x_gromacs_inout_control_posrescom = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_posrescom'))

    x_gromacs_inout_control_posrescom0 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_posrescom0'))

    x_gromacs_inout_control_posrescom1 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_posrescom1'))

    x_gromacs_inout_control_posrescom2 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_posrescom2'))

    x_gromacs_inout_control_posrescomb = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_posrescomb'))

    x_gromacs_inout_control_posrescomb0 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_posrescomb0'))

    x_gromacs_inout_control_posrescomb1 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_posrescomb1'))

    x_gromacs_inout_control_posrescomb2 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_posrescomb2'))

    x_gromacs_inout_control_qmmm = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_qmmm'))

    x_gromacs_inout_control_qmconstraints = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_qmconstraints'))

    x_gromacs_inout_control_qmmmscheme = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_qmmmscheme'))

    x_gromacs_inout_control_mmchargescalefactor = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_mmchargescalefactor'))

    x_gromacs_inout_control_ngqm = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_ngqm'))

    x_gromacs_inout_control_constraintalgorithm = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_constraintalgorithm'))

    x_gromacs_inout_control_continuation = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_continuation'))

    x_gromacs_inout_control_shakesor = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_shakesor'))

    x_gromacs_inout_control_shaketol = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_shaketol'))

    x_gromacs_inout_control_lincsorder = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_lincsorder'))

    x_gromacs_inout_control_lincsiter = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_lincsiter'))

    x_gromacs_inout_control_lincswarnangle = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_lincswarnangle'))

    x_gromacs_inout_control_nwall = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nwall'))

    x_gromacs_inout_control_walltype = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_walltype'))

    x_gromacs_inout_control_wallrlinpot = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_wallrlinpot'))

    x_gromacs_inout_control_wallatomtype = Quantity(
        type=np.dtype(np.float64),
        shape=[2],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_wallatomtype'))

    x_gromacs_inout_control_wallatomtype0 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_wallatomtype0'))

    x_gromacs_inout_control_wallatomtype1 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_wallatomtype1'))

    x_gromacs_inout_control_walldensity = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_walldensity'))

    x_gromacs_inout_control_walldensity0 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_walldensity0'))

    x_gromacs_inout_control_walldensity1 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_walldensity1'))

    x_gromacs_inout_control_wallewaldzfac = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_wallewaldzfac'))

    x_gromacs_inout_control_pull = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_pull'))

    x_gromacs_inout_control_rotation = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_rotation'))

    x_gromacs_inout_control_interactivemd = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_interactivemd'))

    x_gromacs_inout_control_disre = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_disre'))

    x_gromacs_inout_control_disreweighting = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_disreweighting'))

    x_gromacs_inout_control_disremixed = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_disremixed'))

    x_gromacs_inout_control_drfc = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_drfc'))

    x_gromacs_inout_control_drtau = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_drtau'))

    x_gromacs_inout_control_nstdisreout = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nstdisreout'))

    x_gromacs_inout_control_orirefc = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_orirefc'))

    x_gromacs_inout_control_oriretau = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_oriretau'))

    x_gromacs_inout_control_nstorireout = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nstorireout'))

    x_gromacs_inout_control_freeenergy = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_freeenergy'))

    x_gromacs_inout_control_cosacceleration = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_cosacceleration'))

    x_gromacs_inout_control_deform = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_deform'))

    x_gromacs_inout_control_deform0 = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_deform0'))

    x_gromacs_inout_control_deform1 = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_deform1'))

    x_gromacs_inout_control_deform2 = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_deform2'))

    x_gromacs_inout_control_simulatedtempering = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_simulatedtempering'))

    x_gromacs_inout_control_ex = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_ex'))

    x_gromacs_inout_control_ext = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_ext'))

    x_gromacs_inout_control_ey = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_ey'))

    x_gromacs_inout_control_eyt = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_eyt'))

    x_gromacs_inout_control_ez = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_ez'))

    x_gromacs_inout_control_ezt = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_ezt'))

    x_gromacs_inout_control_swapcoords = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_swapcoords'))

    x_gromacs_inout_control_adress = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_adress'))

    x_gromacs_inout_control_userint1 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_userint1'))

    x_gromacs_inout_control_userint2 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_userint2'))

    x_gromacs_inout_control_userint3 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_userint3'))

    x_gromacs_inout_control_userint4 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_userint4'))

    x_gromacs_inout_control_userreal1 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_userreal1'))

    x_gromacs_inout_control_userreal2 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_userreal2'))

    x_gromacs_inout_control_userreal3 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_userreal3'))

    x_gromacs_inout_control_userreal4 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_userreal4'))

    x_gromacs_inout_control_nrdf = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nrdf'))

    x_gromacs_inout_control_reft = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_reft'))

    x_gromacs_inout_control_taut = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_taut'))

    x_gromacs_inout_control_annealing = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_annealing'))

    x_gromacs_inout_control_annealingnpoints = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_annealingnpoints'))

    x_gromacs_inout_control_acc = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_acc'))

    x_gromacs_inout_control_nfreeze = Quantity(
        type=str,
        shape=[3],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_nfreeze'))

    x_gromacs_inout_control_energygrpflags = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 2],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_energygrpflags'))

    x_gromacs_inout_control_energygrpflags0 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_energygrpflags0'))

    x_gromacs_inout_control_energygrpflags1 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_energygrpflags1'))

    x_gromacs_inout_control_energygrpflags2 = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs running environment and control parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_inout_control_energygrpflags2'))


class x_gromacs_section_atom_to_atom_type_ref(MSection):
    '''
    Section to store atom label to atom type definition list
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_gromacs_section_atom_to_atom_type_ref'))

    x_gromacs_atom_to_atom_type_ref = Quantity(
        type=np.dtype(np.int64),
        shape=['number_of_atoms_per_type'],
        description='''
        Reference to the atoms of each atom type.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_atom_to_atom_type_ref'))


class x_gromacs_section_single_configuration_calculation(MSection):
    '''
    section for gathering values for MD steps
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_gromacs_section_single_configuration_calculation'))


class section_system(public.section_system):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_system'))

    x_gromacs_atom_positions_image_index = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_atoms', 3],
        unit='dimensionless',
        description='''
        PBC image flag index.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_atom_positions_image_index'))

    x_gromacs_atom_positions_scaled = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        unit='dimensionless',
        description='''
        Position of the atoms in a scaled format [0, 1].
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_atom_positions_scaled'))

    x_gromacs_atom_positions_wrapped = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        unit='meter',
        description='''
        Position of the atoms wrapped back to the periodic box.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_atom_positions_wrapped'))

    x_gromacs_lattice_lengths = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        description='''
        Lattice dimensions in a vector. Vector includes [a, b, c] lengths.
        ''',
        categories=[public.configuration_core],
        a_legacy=LegacyDefinition(name='x_gromacs_lattice_lengths'))

    x_gromacs_lattice_angles = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        description='''
        Angles of lattice vectors. Vector includes [alpha, beta, gamma] in degrees.
        ''',
        categories=[public.configuration_core],
        a_legacy=LegacyDefinition(name='x_gromacs_lattice_angles'))

    x_gromacs_dummy = Quantity(
        type=str,
        shape=[],
        description='''
        dummy
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_dummy'))

    x_gromacs_mdin_finline = Quantity(
        type=str,
        shape=[],
        description='''
        finline in mdin
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_mdin_finline'))

    x_gromacs_traj_timestep_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_traj_timestep_store'))

    x_gromacs_traj_number_of_atoms_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_traj_number_of_atoms_store'))

    x_gromacs_traj_box_bound_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_traj_box_bound_store'))

    x_gromacs_traj_box_bounds_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_traj_box_bounds_store'))

    x_gromacs_traj_variables_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_traj_variables_store'))

    x_gromacs_traj_atoms_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_traj_atoms_store'))


class section_sampling_method(public.section_sampling_method):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_sampling_method'))

    x_gromacs_barostat_target_pressure = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='pascal',
        description='''
        MD barostat target pressure.
        ''',
        categories=[public.settings_sampling, public.settings_molecular_dynamics, public.settings_barostat],
        a_legacy=LegacyDefinition(name='x_gromacs_barostat_target_pressure'))

    x_gromacs_barostat_tau = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        MD barostat relaxation time.
        ''',
        categories=[public.settings_sampling, public.settings_molecular_dynamics, public.settings_barostat],
        a_legacy=LegacyDefinition(name='x_gromacs_barostat_tau'))

    x_gromacs_barostat_type = Quantity(
        type=str,
        shape=[],
        description='''
        MD barostat type, valid values are defined in the barostat_type wiki page.
        ''',
        categories=[public.settings_sampling, public.settings_molecular_dynamics, public.settings_barostat],
        a_legacy=LegacyDefinition(name='x_gromacs_barostat_type'))

    x_gromacs_integrator_dt = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        MD integration time step.
        ''',
        categories=[public.settings_sampling, public.settings_molecular_dynamics, public.settings_integrator],
        a_legacy=LegacyDefinition(name='x_gromacs_integrator_dt'))

    x_gromacs_integrator_type = Quantity(
        type=str,
        shape=[],
        description='''
        MD integrator type, valid values are defined in the integrator_type wiki page.
        ''',
        categories=[public.settings_sampling, public.settings_molecular_dynamics, public.settings_integrator],
        a_legacy=LegacyDefinition(name='x_gromacs_integrator_type'))

    x_gromacs_periodicity_type = Quantity(
        type=str,
        shape=[],
        description='''
        Periodic boundary condition type in the sampling (non-PBC or PBC).
        ''',
        categories=[public.settings_sampling, public.settings_molecular_dynamics, public.settings_integrator],
        a_legacy=LegacyDefinition(name='x_gromacs_periodicity_type'))

    x_gromacs_langevin_gamma = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Langevin thermostat damping factor.
        ''',
        categories=[public.settings_thermostat, public.settings_sampling, public.settings_molecular_dynamics],
        a_legacy=LegacyDefinition(name='x_gromacs_langevin_gamma'))

    x_gromacs_number_of_steps_requested = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Number of requested MD integration time steps.
        ''',
        categories=[public.settings_sampling, public.settings_molecular_dynamics, public.settings_integrator],
        a_legacy=LegacyDefinition(name='x_gromacs_number_of_steps_requested'))

    x_gromacs_thermostat_level = Quantity(
        type=str,
        shape=[],
        description='''
        MD thermostat level (see wiki: single, multiple, regional).
        ''',
        categories=[public.settings_thermostat, public.settings_sampling, public.settings_molecular_dynamics],
        a_legacy=LegacyDefinition(name='x_gromacs_thermostat_level'))

    x_gromacs_thermostat_target_temperature = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='kelvin',
        description='''
        MD thermostat target temperature.
        ''',
        categories=[public.settings_thermostat, public.settings_sampling, public.settings_molecular_dynamics],
        a_legacy=LegacyDefinition(name='x_gromacs_thermostat_target_temperature'))

    x_gromacs_thermostat_tau = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        MD thermostat relaxation time.
        ''',
        categories=[public.settings_thermostat, public.settings_sampling, public.settings_molecular_dynamics],
        a_legacy=LegacyDefinition(name='x_gromacs_thermostat_tau'))

    x_gromacs_thermostat_type = Quantity(
        type=str,
        shape=[],
        description='''
        MD thermostat type, valid values are defined in the thermostat_type wiki page.
        ''',
        categories=[public.settings_thermostat, public.settings_sampling, public.settings_molecular_dynamics],
        a_legacy=LegacyDefinition(name='x_gromacs_thermostat_type'))


class section_atom_type(common.section_atom_type):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_atom_type'))

    x_gromacs_atom_name = Quantity(
        type=str,
        shape=[],
        description='''
        Atom name of an atom in topology definition.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_atom_name'))

    x_gromacs_atom_type = Quantity(
        type=str,
        shape=[],
        description='''
        Atom type of an atom in topology definition.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_atom_type'))

    x_gromacs_atom_element = Quantity(
        type=str,
        shape=[],
        description='''
        Atom type of an atom in topology definition.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_atom_element'))

    x_gromacs_atom_type_element = Quantity(
        type=str,
        shape=[],
        description='''
        Element symbol of an atom type.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_atom_type_element'))

    x_gromacs_atom_type_radius = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        van der Waals radius of an atom type.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_atom_type_radius'))

    number_of_atoms_per_type = Quantity(
        type=int,
        shape=[],
        description='''
        Number of atoms involved in this type.
        ''',
        a_legacy=LegacyDefinition(name='number_of_atoms_per_type'))


class section_interaction(common.section_interaction):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_interaction'))

    x_gromacs_interaction_atom_to_atom_type_ref = Quantity(
        type=common.section_atom_type,
        shape=['number_of_atoms_per_interaction'],
        description='''
        Reference to the atom type of each interaction atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_interaction_atom_to_atom_type_ref'))

    x_gromacs_number_of_defined_pair_interactions = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of defined pair interactions (L-J pairs).
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_number_of_defined_pair_interactions'))

    x_gromacs_pair_interaction_atom_type_ref = Quantity(
        type=common.section_atom_type,
        shape=['x_gromacs_number_of_defined_pair_interactions', 'number_of_atoms_per_interaction'],
        description='''
        Reference to the atom type for pair interactions.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_pair_interaction_atom_type_ref'))

    x_gromacs_pair_interaction_parameters = Quantity(
        type=np.dtype(np.float64),
        shape=['x_gromacs_number_of_defined_pair_interactions', 2],
        description='''
        Pair interactions parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_pair_interaction_parameters'))


class section_molecule_interaction(common.section_molecule_interaction):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_molecule_interaction'))

    x_gromacs_molecule_interaction_atom_to_atom_type_ref = Quantity(
        type=common.section_atom_type,
        shape=['number_of_atoms_per_interaction'],
        description='''
        Reference to the atom type of each molecule interaction atoms.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_molecule_interaction_atom_to_atom_type_ref'))

    x_gromacs_number_of_defined_molecule_pair_interactions = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of defined pair interactions within a molecule (L-J pairs).
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_number_of_defined_molecule_pair_interactions'))

    x_gromacs_pair_molecule_interaction_parameters = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_defined_molecule_pair_interactions', 2],
        description='''
        Molecule pair interactions parameters.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_pair_molecule_interaction_parameters'))

    x_gromacs_pair_molecule_interaction_to_atom_type_ref = Quantity(
        type=common.section_atom_type,
        shape=['x_gromacs_number_of_defined_pair_interactions', 'number_of_atoms_per_interaction'],
        description='''
        Reference to the atom type for pair interactions within a molecule.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_pair_molecule_interaction_to_atom_type_ref'))


class section_run(public.section_run):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_run'))

    x_gromacs_program_version_date = Quantity(
        type=str,
        shape=[],
        description='''
        Program version date.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_program_version_date'))

    x_gromacs_parallel_task_nr = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Program task no.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_parallel_task_nr'))

    x_gromacs_number_of_tasks = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Number of tasks in parallel program (MPI).
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_number_of_tasks'))

    x_gromacs_program_module_version = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs program module (gmx) version.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_program_module_version'))

    x_gromacs_program_license = Quantity(
        type=str,
        shape=[],
        description='''
        Gromacs program license.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_program_license'))

    x_gromacs_xlo_xhi = Quantity(
        type=str,
        shape=[],
        description='''
        test
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_xlo_xhi'))

    x_gromacs_data_file_store = Quantity(
        type=str,
        shape=[],
        description='''
        Filename of data file
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_data_file_store'))

    x_gromacs_program_working_path = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_program_working_path'))

    x_gromacs_program_execution_host = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_program_execution_host'))

    x_gromacs_program_execution_path = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_program_execution_path'))

    x_gromacs_program_module = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_program_module'))

    x_gromacs_program_execution_date = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_program_execution_date'))

    x_gromacs_program_execution_time = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_program_execution_time'))

    x_gromacs_mdin_header = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_mdin_header'))

    x_gromacs_mdin_wt = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_mdin_wt'))

    x_gromacs_section_input_output_files = SubSection(
        sub_section=SectionProxy('x_gromacs_section_input_output_files'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_gromacs_section_input_output_files'))

    x_gromacs_section_control_parameters = SubSection(
        sub_section=SectionProxy('x_gromacs_section_control_parameters'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_gromacs_section_control_parameters'))


class section_topology(common.section_topology):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_topology'))

    x_gromacs_input_units_store = Quantity(
        type=str,
        shape=[],
        description='''
        It determines the units of all quantities specified in the input script and data
        file, as well as quantities output to the screen, log file, and dump files.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_input_units_store'))

    x_gromacs_data_bond_types_store = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        store temporarly
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_data_bond_types_store'))

    x_gromacs_data_bond_count_store = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        store temporarly
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_data_bond_count_store'))

    x_gromacs_data_angle_count_store = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        store temporarly
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_data_angle_count_store'))

    x_gromacs_data_atom_types_store = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        store temporarly
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_data_atom_types_store'))

    x_gromacs_data_dihedral_count_store = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        store temporarly
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_data_dihedral_count_store'))

    x_gromacs_data_angles_store = Quantity(
        type=str,
        shape=[],
        description='''
        store temporarly
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_data_angles_store'))

    x_gromacs_data_angle_list_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_data_angle_list_store'))

    x_gromacs_data_bond_list_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_data_bond_list_store'))

    x_gromacs_data_dihedral_list_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_data_dihedral_list_store'))

    x_gromacs_data_dihedral_coeff_list_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_data_dihedral_coeff_list_store'))

    x_gromacs_masses_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_masses_store'))

    x_gromacs_data_topo_list_store = Quantity(
        type=str,
        shape=[],
        description='''
        tmp
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_data_topo_list_store'))

    x_gromacs_section_atom_to_atom_type_ref = SubSection(
        sub_section=SectionProxy('x_gromacs_section_atom_to_atom_type_ref'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_gromacs_section_atom_to_atom_type_ref'))


class section_frame_sequence(public.section_frame_sequence):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_frame_sequence'))

    x_gromacs_number_of_volumes_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of volumes in this sequence of frames, see
        x_gromacs_frame_sequence_volume.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_number_of_volumes_in_sequence'))

    x_gromacs_number_of_densities_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of densities in this sequence of frames, see
        x_gromacs_frame_sequence_density.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_number_of_densities_in_sequence'))

    x_gromacs_number_of_ubond_energies_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of ubond_energies in this sequence of frames, see
        x_gromacs_frame_sequence_ubond_energy.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_number_of_ubond_energies_in_sequence'))

    x_gromacs_number_of_bond_energies_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of bond_energies in this sequence of frames, see
        x_gromacs_frame_sequence_bond_energy.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_number_of_bond_energies_in_sequence'))

    x_gromacs_number_of_coulomb_sr_energies_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of coulomb_sr_energies in this sequence of frames, see
        x_gromacs_frame_sequence_coulomb_sr_energy.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_number_of_coulomb_sr_energies_in_sequence'))

    x_gromacs_number_of_coulomb_14_energies_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of coulomb_14_energies in this sequence of frames, see
        x_gromacs_frame_sequence_coulomb_14_energy.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_number_of_coulomb_14_energies_in_sequence'))

    x_gromacs_number_of_lj_sr_energies_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of lj_sr_energies in this sequence of frames, see
        x_gromacs_frame_sequence_lj_sr_energy.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_number_of_lj_sr_energies_in_sequence'))

    x_gromacs_number_of_lj_14_energies_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of lj_14_energies in this sequence of frames, see
        x_gromacs_frame_sequence_lj_14_energy.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_number_of_lj_14_energies_in_sequence'))

    x_gromacs_number_of_proper_dihedral_energies_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of proper_dihedral_energies in this sequence of frames, see
        x_gromacs_frame_sequence_proper_dihedral_energy.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_number_of_proper_dihedral_energies_in_sequence'))

    x_gromacs_number_of_improper_dihedral_energies_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of improper_dihedral_energies in this sequence of frames, see
        x_gromacs_frame_sequence_improper_dihedral_energy.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_number_of_improper_dihedral_energies_in_sequence'))

    x_gromacs_number_of_cmap_dihedral_energies_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of cmap_dihedral_energies in this sequence of frames, see
        x_gromacs_frame_sequence_cmap_dihedral_energy.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_number_of_cmap_dihedral_energies_in_sequence'))

    x_gromacs_number_of_constrain_rmsd_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of constrain_rmsd_energies in this sequence of frames, see
        x_gromacs_frame_sequence_constrain_rmsd_energy.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_number_of_constrain_rmsd_in_sequence'))

    x_gromacs_frame_sequence_density_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['x_gromacs_number_of_densities_in_sequence'],
        description='''
        Array containing the strictly increasing indices of the frames the
        x_gromacs_frame_sequence_densities values refers to. If not given it defaults to
        the trivial mapping 0,1,...
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_density_frames'))

    x_gromacs_frame_sequence_density = Quantity(
        type=np.dtype(np.float64),
        shape=['x_gromacs_number_of_densities_in_sequence'],
        description='''
        Array containing the values of the density along this sequence of frames (i.e., a
        trajectory, a frame is one section_single_configuration_calculation). If not all
        frames have a value the indices of the frames that have a value are stored in
        frame_sequence_density_frames.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_density'))

    x_gromacs_frame_sequence_ubond_energy_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['x_gromacs_number_of_ubond_energies_in_sequence'],
        description='''
        Array containing the strictly increasing indices of the frames the
        x_gromacs_frame_sequence_ubond_energy values refers to. If not given it defaults
        to the trivial mapping 0,1,...
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_ubond_energy_frames'))

    x_gromacs_frame_sequence_ubond_energy = Quantity(
        type=np.dtype(np.float64),
        shape=['x_gromacs_number_of_ubond_energies_in_sequence'],
        description='''
        Array containing the values of the ubond_energy along this sequence of frames
        (i.e., a trajectory, a frame is one section_single_configuration_calculation). If
        not all frames have a value the indices of the frames that have a value are stored
        in frame_sequence_ubond_energy_frames.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_ubond_energy'))

    x_gromacs_frame_sequence_coulomb_sr_energy_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['x_gromacs_number_of_coulomb_sr_energy_in_sequence'],
        description='''
        Array containing the strictly increasing indices of the frames the
        x_gromacs_frame_sequence_coulomb_sr_energy values refers to. If not given it
        defaults to the trivial mapping 0,1,...
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_coulomb_sr_energy_frames'))

    x_gromacs_frame_sequence_coulomb_sr_energy = Quantity(
        type=np.dtype(np.float64),
        shape=['x_gromacs_number_of_coulomb_sr_energy_in_sequence'],
        description='''
        Array containing the values of the coulomb_sr_energy along this sequence of frames
        (i.e., a trajectory, a frame is one section_single_configuration_calculation). If
        not all frames have a value the indices of the frames that have a value are stored
        in frame_sequence_coulomb_sr_energy_frames.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_coulomb_sr_energy'))

    x_gromacs_frame_sequence_coulomb_14_energy_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['x_gromacs_number_of_coulomb_14_energy_in_sequence'],
        description='''
        Array containing the strictly increasing indices of the frames the
        x_gromacs_frame_sequence_coulomb_14_energy values refers to. If not given it
        defaults to the trivial mapping 0,1,...
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_coulomb_14_energy_frames'))

    x_gromacs_frame_sequence_coulomb_14_energy = Quantity(
        type=np.dtype(np.float64),
        shape=['x_gromacs_number_of_coulomb_14_energy_in_sequence'],
        description='''
        Array containing the values of the coulomb_14_energy along this sequence of frames
        (i.e., a trajectory, a frame is one section_single_configuration_calculation). If
        not all frames have a value the indices of the frames that have a value are stored
        in frame_sequence_coulomb_14_energy_frames.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_coulomb_14_energy'))

    x_gromacs_frame_sequence_lj_sr_energy_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['x_gromacs_number_of_lj_sr_energy_in_sequence'],
        description='''
        Array containing the strictly increasing indices of the frames the
        x_gromacs_frame_sequence_lj_sr_energy values refers to. If not given it defaults
        to the trivial mapping 0,1,...
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_lj_sr_energy_frames'))

    x_gromacs_frame_sequence_lj_sr_energy = Quantity(
        type=np.dtype(np.float64),
        shape=['x_gromacs_number_of_lj_sr_energy_in_sequence'],
        description='''
        Array containing the values of the lj_sr_energy along this sequence of frames
        (i.e., a trajectory, a frame is one section_single_configuration_calculation). If
        not all frames have a value the indices of the frames that have a value are stored
        in frame_sequence_lj_sr_energy_frames.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_lj_sr_energy'))

    x_gromacs_frame_sequence_lj_14_energy_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['x_gromacs_number_of_lj_14_energy_in_sequence'],
        description='''
        Array containing the strictly increasing indices of the frames the
        x_gromacs_frame_sequence_lj_14_energy values refers to. If not given it defaults
        to the trivial mapping 0,1,...
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_lj_14_energy_frames'))

    x_gromacs_frame_sequence_lj_14_energy = Quantity(
        type=np.dtype(np.float64),
        shape=['x_gromacs_number_of_lj_14_energy_in_sequence'],
        description='''
        Array containing the values of the lj_14_energy along this sequence of frames
        (i.e., a trajectory, a frame is one section_single_configuration_calculation). If
        not all frames have a value the indices of the frames that have a value are stored
        in frame_sequence_lj_14_energy_frames.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_lj_14_energy'))

    x_gromacs_frame_sequence_constrain_rmsd_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['x_gromacs_number_of_constrain_rmsd_energy_in_sequence'],
        description='''
        Array containing the strictly increasing indices of the frames the
        x_gromacs_frame_sequence_constrain_rmsd_energy values refers to. If not given it
        defaults to the trivial mapping 0,1,...
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_constrain_rmsd_frames'))

    x_gromacs_frame_sequence_constrain_rmsd = Quantity(
        type=np.dtype(np.float64),
        shape=['x_gromacs_number_of_constrain_rmsd_in_sequence'],
        description='''
        Array containing the values of the constrain_rmsd_energy along this sequence of
        frames (i.e., a trajectory, a frame is one
        section_single_configuration_calculation). If not all frames have a value the
        indices of the frames that have a value are stored in
        frame_sequence_constrain_rmsd_energy_frames.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_constrain_rmsd'))

    x_gromacs_frame_sequence_cmap_dihedral_energy_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['x_gromacs_number_of_cmap_dihedral_energy_in_sequence'],
        description='''
        Array containing the strictly increasing indices of the frames the
        x_gromacs_frame_sequence_cmap_dihedral_energy values refers to. If not given it
        defaults to the trivial mapping 0,1,...
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_cmap_dihedral_energy_frames'))

    x_gromacs_frame_sequence_cmap_dihedral_energy = Quantity(
        type=np.dtype(np.float64),
        shape=['x_gromacs_number_of_cmap_dihedral_energy_in_sequence'],
        description='''
        Array containing the values of the cmap_dihedral_energy along this sequence of
        frames (i.e., a trajectory, a frame is one
        section_single_configuration_calculation). If not all frames have a value the
        indices of the frames that have a value are stored in
        frame_sequence_cmap_dihedral_energy_frames.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_cmap_dihedral_energy'))

    x_gromacs_frame_sequence_improper_dihedral_energy_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['x_gromacs_number_of_improper_dihedral_energy_in_sequence'],
        description='''
        Array containing the strictly increasing indices of the frames the
        x_gromacs_frame_sequence_improper_dihedral_energy values refers to. If not given
        it defaults to the trivial mapping 0,1,...
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_improper_dihedral_energy_frames'))

    x_gromacs_frame_sequence_improper_dihedral_energy = Quantity(
        type=np.dtype(np.float64),
        shape=['x_gromacs_number_of_improper_dihedral_energy_in_sequence'],
        description='''
        Array containing the values of the improper_dihedral_energy along this sequence of
        frames (i.e., a trajectory, a frame is one
        section_single_configuration_calculation). If not all frames have a value the
        indices of the frames that have a value are stored in
        frame_sequence_improper_dihedral_energy_frames.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_improper_dihedral_energy'))

    x_gromacs_frame_sequence_proper_dihedral_energy_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['x_gromacs_number_of_proper_dihedral_energy_in_sequence'],
        description='''
        Array containing the strictly increasing indices of the frames the
        x_gromacs_frame_sequence_proper_dihedral_energy values refers to. If not given it
        defaults to the trivial mapping 0,1,...
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_proper_dihedral_energy_frames'))

    x_gromacs_frame_sequence_proper_dihedral_energy = Quantity(
        type=np.dtype(np.float64),
        shape=['x_gromacs_number_of_proper_dihedral_energy_in_sequence'],
        description='''
        Array containing the values of the proper_dihedral_energy along this sequence of
        frames (i.e., a trajectory, a frame is one
        section_single_configuration_calculation). If not all frames have a value the
        indices of the frames that have a value are stored in
        frame_sequence_proper_dihedral_energy_frames.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_proper_dihedral_energy'))

    x_gromacs_frame_sequence_bond_energy_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['x_gromacs_number_of_bond_energies_in_sequence'],
        description='''
        Array containing the strictly increasing indices of the frames the
        x_gromacs_frame_sequence_bond_energy values refers to. If not given it defaults to
        the trivial mapping 0,1,...
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_bond_energy_frames'))

    x_gromacs_frame_sequence_bond_energy = Quantity(
        type=np.dtype(np.float64),
        shape=['x_gromacs_number_of_bond_energies_in_sequence'],
        description='''
        Array containing the values of the bond_energy along this sequence of frames
        (i.e., a trajectory, a frame is one section_single_configuration_calculation). If
        not all frames have a value the indices of the frames that have a value are stored
        in frame_sequence_bond_energy_frames.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_bond_energy'))

    x_gromacs_frame_sequence_volume_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['x_gromacs_number_of_volumes_in_sequence'],
        description='''
        Array containing the strictly increasing indices of the frames the
        x_gromacs_frame_sequence_volume values refers to. If not given it defaults to the
        trivial mapping 0,1,...
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_volume_frames'))

    x_gromacs_frame_sequence_volume = Quantity(
        type=np.dtype(np.float64),
        shape=['x_gromacs_number_of_volumes_in_sequence'],
        description='''
        Array containing the values of the volume along this sequence of frames (i.e., a
        trajectory, a frame is one section_single_configuration_calculation). If not all
        frames have a value the indices of the frames that have a value are stored in
        frame_sequence_volume_frames.
        ''',
        a_legacy=LegacyDefinition(name='x_gromacs_frame_sequence_volume'))


class section_single_configuration_calculation(public.section_single_configuration_calculation):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_single_configuration_calculation'))

    x_gromacs_section_single_configuration_calculation = SubSection(
        sub_section=SectionProxy('x_gromacs_section_single_configuration_calculation'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_gromacs_section_single_configuration_calculation'))


m_package.__init_metainfo__()
