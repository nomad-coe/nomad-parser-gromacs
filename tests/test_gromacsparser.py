#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
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

import pytest
import numpy as np

from nomad.datamodel import EntryArchive
from gromacsparser.gromacs_parser import GromacsParser


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


@pytest.fixture(scope='module')
def parser():
    return GromacsParser()


def test_md_verbose(parser):
    archive = EntryArchive()
    parser.parse('tests/data/fe_test/md.log', archive, None)

    sec_run = archive.run[0]
    assert sec_run.program.version == '5.1.4'
    sec_control = sec_run.x_gromacs_section_control_parameters[0]
    assert sec_control.x_gromacs_inout_control_coulombtype == 'PME'
    assert np.shape(sec_control.x_gromacs_inout_control_deform) == (3, 3)

    sec_md = archive.workflow[0].molecular_dynamics
    assert sec_md.ensemble_type == 'NPT'
    assert sec_md.x_gromacs_integrator_dt.magnitude == 0.0005
    assert sec_md.x_gromacs_barostat_target_pressure.magnitude == approx(33333.33)

    sec_sccs = sec_run.calculation
    assert len(sec_sccs) == 9
    assert sec_sccs[2].energy.total.value.magnitude == approx(-3.2711290665182795e-17)
    assert sec_sccs[5].thermodynamics.pressure.magnitude == approx(1.21842e+08)
    assert sec_sccs[7].energy.contributions[0].value.magnitude == approx(7.377359821857783e-18)
    assert sec_sccs[0].forces.total.value[5][2].magnitude == approx(-7.932968909721231e-10)

    sec_systems = sec_run.system
    assert len(sec_systems) == 2
    assert np.shape(sec_systems[0].atoms.positions) == (1516, 3)
    assert sec_systems[1].atoms.positions[800][1].magnitude == approx(2.4609454e-09)
    assert sec_systems[0].atoms.velocities[500][0].magnitude == approx(869.4773)
    assert sec_systems[1].atoms.lattice_vectors[2][2].magnitude == approx(2.469158e-09)


def test_md_edr(parser):
    archive = EntryArchive()
    parser.parse('tests/data/fe_test/mdrun.out', archive, None)

    assert len(archive.run[0].calculation) == 7
