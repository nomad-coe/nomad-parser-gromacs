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

from builtins import map
from builtins import range
from builtins import object
import numpy as np
import nomadcore.ActivateLogging
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import SimpleMatcher as SM
#import nomadcore.smart_parser.SmartParserCommon as SmartParser
from nomadcore.smart_parser import SmartParserCommon as SmartParser
from nomadcore.smart_parser.SmartParserCommon import get_metaInfo, conv_str, conv_int, conv_float, open_section
from nomadcore.smart_parser.SmartParserDictionary import getList_MetaStrInDict, getDict_MetaStrInDict
from nomadcore.smart_parser.SmartParserDictionary import isMetaStrInDict
from .GromacsDictionary import get_updateDictionary, set_Dictionaries
from .GromacsCmdLineArgs import get_commandLineArguments
from .GromacsCommon import PARSERNAME, PROGRAMNAME, PARSERVERSION, PARSERTAG, LOGGER
from .GromacsCommon import PARSER_INFO_DEFAULT, META_INFO_PATH, set_excludeList, set_includeList
from .GromacsCmdLineArgs import get_commandLineArguments
#import nomadcore.md_data_access.MDDataAccess as MDDA
from nomadcore.md_data_access import MDDataAccess as MDDA
import argparse
import logging
import os
import re
import sys
import datetime
import io
from nomadcore.simple_parser import mainFunction

############################################################
# This is the parser for the main file of Gromacs.
############################################################

#PRINTABLE = re.compile(r"\W+")

class GromacsParser(SmartParser.ParserBase):
    """Context for parsing Gromacs main file.

    This class keeps tracks of several Gromacs settings to adjust the parsing to them.
    The onClose_ functions allow processing and writing of cached values after a section is closed.
    They take the following arguments:
        backend: Class that takes care of writing and caching of metadata.
        gIndex: Index of the section that is closed.
        section: The cached values and sections that were found in the section that is closed.
    """
    def __init__(self):
        # dictionary of energy values, which are tracked between SCF iterations and written after convergence
        self.totalEnergyList = {
                                'energy_electrostatic': None,
                                'energy_total_T0_per_atom': None,
                                'energy_free_per_atom': None,
                               }
        SmartParser.ParserBase.__init__(
            self, re_program_name=re.compile(r"\s*"+PROGRAMNAME+"$"),
            parsertag=PARSERTAG, metainfopath=META_INFO_PATH,
            parserinfodef=PARSER_INFO_DEFAULT)

        set_Dictionaries(self)
        self.metaInfoEnv = get_metaInfo(self)
        self.secGIndexDict = {}

        self.cachingLevelForMetaName = {
                               PARSERTAG + '_trajectory_file_detect': CachingLevel.Cache,
                               PARSERTAG + '_geometry_optimization_cdetect': CachingLevel.Cache,
                               PARSERTAG + '_mdin_finline': CachingLevel.Ignore,
                               #PARSERTAG + '_section_input_output_files': CachingLevel.Ignore,
                               PARSERTAG + '_single_configuration_calculation_detect': CachingLevel.Cache,
                              }
        for name in self.metaInfoEnv.infoKinds:
            metaInfo = self.metaInfoEnv.infoKinds[name]
            if (name.startswith(PARSERTAG + '_mdin_') and
                metaInfo.kindStr == "type_document_content" and
                (PARSERTAG + "_mdin_method" in metaInfo.superNames or
                 PARSERTAG + "_mdin_run" in metaInfo.superNames or
                 PARSERTAG + "_mdin_system" in metaInfo.superNames) or
                name.startswith(PARSERTAG + '_parm_') and
                metaInfo.kindStr == "type_document_content" and
                (PARSERTAG + "_mdin_method" in metaInfo.superNames or
                 PARSERTAG + "_mdin_run" in metaInfo.superNames or
                 PARSERTAG + "_mdin_system" in metaInfo.superNames) or
                #name.startswith(PARSERTAG + '_mdin_file_') and
                name.startswith(PARSERTAG + '_inout_file_') or
                #metaInfo.kindStr == "type_document_content" and
                #(PARSERTAG + "_section_input_output_files" in metaInfo.superNames or
                # "section_run" in metaInfo.superNames) or
                name.startswith(PARSERTAG + '_inout_control_') or
                #(PARSERTAG + "_section_control_parameters" in metaInfo.superNames) or
                #name.startswith(PARSERTAG + '_mdin_') and
                #(PARSERTAG + "_section_control_parameters" in metaInfo.superNames) or
                name.startswith(PARSERTAG + '_mdout_') or
                name.startswith(PARSERTAG + '_mdout_') and
                #metaInfo.kindStr == "type_document_content" and
                (PARSERTAG + "_mdout_method" in metaInfo.superNames or
                 PARSERTAG + "_mdout_system" in metaInfo.superNames or
                 "section_run" in metaInfo.superNames or
                 PARSERTAG + "_mdout_single_configuration_calculation" in metaInfo.superNames)
                or name.startswith('section_single_configuration_calculation')
                or PARSERTAG + '_mdout' in metaInfo.superNames
                or 'section_sampling_method' in metaInfo.superNames
                or 'section_single_configuration_calculation' in metaInfo.superNames
               ):
                self.cachingLevelForMetaName[name] = CachingLevel.Cache
            if name in self.extraDict.keys():
                self.cachingLevelForMetaName[name] = CachingLevel.Ignore


    def initialize_values(self):
        """Initializes the values of certain variables.

        This allows a consistent setting and resetting of the variables,
        when the parsing starts and when a section_run closes.
        """
        set_Dictionaries(self)
        self.secGIndexDict.clear()
        self.secMethodGIndex = None
        self.secSystemGIndex = None
        self.secTopologyGIndex = None
        self.secSamplingGIndex = None
        self.secSingleGIndex = None
        self.secVDWGIndex = None
        self.secAtomType = None
        self.inputMethodIndex = None
        self.inputControlIndex = None
        self.mainMethodIndex = None
        self.mainCalcIndex = None
        self.MD = True
        self.topologyDict = None
        self.topologyTable = None
        self.topologyBonds = None
        self.topology = None
        self.topologyFormat = None
        self.topologyFile = None
        self.trajectory = None
        self.trajectoryFormat = None
        self.trajectoryFile = None
        self.readChunk = 300
        self.boxlengths = None
        self.latticevectors = None
        self.atompositions = None
        self.inputcoords = None
        self.inputpositions = None
        self.outputcoords = None
        self.outputpositions = None
        # start with -1 since zeroth iteration is the initialization
        self.MDiter = -1
        self.MDlogsteps = None
        self.MDlogstep = None
        self.MDstep = 1
        self.MDcurrentstep = -1
        self.MDnextstep = 0
        self.MDnextlogstep = 0
        self.singleConfCalcs = []
        self.minConverged = None
        self.parsedLogFile = False
        self.LogSuperContext = None
        self.forces_raw = []
        self.lastCalculationGIndex = None
        self.logFileName = None
        self.lastfInLine = None
        self.lastfInMatcher = None
        self.secOpen = open_section
        self.superP = self.parser.backend.superBackend
        self.MDData = MDDA.MDDataAccess()
        if self.recordList:
            self.recordList.close()
        self.recordList = io.StringIO()

    def startedParsing(self, fInName, parser):
        """Function is called when the parsing starts.

        Get compiled parser, filename and metadata.

        Args:
            fInName: The file name on which the current parser is running.
            parser: The compiled parser. Is an object of the class SimpleParser in nomadcore.simple_parser.py.
        """
        self.parser = parser
        self.fName = fInName
        # save metadata
        self.metaInfoEnv = self.parser.parserBuilder.metaInfoEnv
        # allows to reset values if the same superContext is used to parse different files
        self.initialize_values()

    def peekline(self, parser):
        pos = parser.fIn.fIn.tell()
        line = parser.fIn.fIn.readline()
        parser.fIn.fIn.seek(pos)
        return line

    def onClose_section_run(self, backend, gIndex, section):
        """Trigger called when section_run is closed.

        Write the keywords from control parametres and
        the Gromacs output from the parsed log output,
        which belong to settings_run.
        Variables are reset to ensure clean start for new run.
        """
        frameSequenceGIndex = backend.openSection("section_frame_sequence")
        section_frameseq_Dict = get_updateDictionary(self, 'frameseqend')
        updateFrameDict = {
                'startSection' : [
                    ['section_frame_sequence']],
                #'muteSections' : [['section_sampling_method']],
                'dictionary' : section_frameseq_Dict
                }
        self.metaStorage.update(updateFrameDict)
        self.metaStorage.updateBackend(backend.superBackend,
                startsection=['section_frame_sequence'],
                autoopenclose=False)
        backend.addValue("frame_sequence_to_sampling_ref", self.secSamplingGIndex)
        backend.addArrayValues("frame_sequence_local_frames_ref", np.asarray(self.singleConfCalcs))
        backend.closeSection("section_frame_sequence", frameSequenceGIndex)

        # reset all variables
        self.initialize_values()

    def onClose_x_gromacs_section_input_output_files(self, backend, gIndex, section):
        """Trigger called when x_gromacs_section_input_output_files is closed.

        Determine whether topology, trajectory and input coordinate files are
        supplied to the parser

        Initiates topology and trajectory file handles.

        Captures topology, atomic positions, atom labels, lattice vectors and
        stores them before section_system and
        section_single_configuration_calculation are encountered.
        """
        # Checking whether topology, input
        # coordinates and trajectory files exist
        atLeastOneFileExist = False
        working_dir_name = os.path.dirname(os.path.abspath(self.fName))
        for k,v in self.fileDict.items():
            if k.startswith(PARSERTAG + '_inout_file'):
                if v.value:
                    file_name = os.path.normpath(os.path.join(working_dir_name, v.value))
                    self.fileDict[k].fileSupplied = os.path.isfile(file_name)
                    self.fileDict[k].activeInfo = False
                    if self.fileDict[k].fileSupplied:
                        self.fileDict[k].fileName = file_name
                        self.fileDict[k].activeInfo = True
                        if self.fileDict[k].activeInfo:
                            atLeastOneFileExist = True
                    else:
                        self.fileDict[k].value = None
        updateDict = {
            'startSection'   : [[PARSERTAG+'_section_input_output_files']],
            'dictionary'     : self.fileDict
            }
        self.metaStorage.update(updateDict)
        self.metaStorage.updateBackend(backend.superBackend,
                startsection=[PARSERTAG+'_section_input_output_files'],
                autoopenclose=False)
        if atLeastOneFileExist:
            self.MDData.initializeFileHandlers(self)
            if self.topology:
                self.topoDict = self.topology.topoDict
            if self.trajectory:
                self.atompositions = self.trajectory.positions()
            if self.inputcoords:
                self.inputpositions = self.inputcoords.positions()
            if self.outputcoords:
                self.outputpositions = self.outputcoords.positions()
        #self.MDiter += self.MDstep
        self.MDiter += 1
        self.MDcurrentstep =  int(self.MDiter)
        self.stepcontrolDict.update({"MDcurrentstep" : int(self.MDcurrentstep)})
            #pass

    def onOpen_x_gromacs_section_control_parameters(self, backend, gIndex, section):
        # keep track of the latest section
        if self.inputControlIndex is None:
            self.inputControlIndex = gIndex

    def onClose_x_gromacs_section_control_parameters(self, backend, gIndex, section):
        section_control_Dict = {}
        section_control_Dict.update(self.archDict)
        section_control_Dict.update(self.cntrlDict)
        section_control_Dict.update(self.qmDict)
        section_control_Dict.update(self.grpDict)
        section_control_Dict.update(self.annealDict)
        updateDict = {
            'startSection' : [[PARSERTAG+'_section_control_parameters']],
            'dictionary'   : section_control_Dict
            }
        self.metaStorage.update(updateDict)
        self.metaStorage.updateBackend(backend.superBackend,
                startsection=[PARSERTAG+'_section_control_parameters'],
                autoopenclose=False)
        # Gromacs prints the initial and final energies to the log file.
        # The total number of MD steps in Gromacs is nsteps irrelevant
        # to the number of steps in log file of energy file (.edr)
        nstep = 0
        nstlog = 0
        nstxout = 0
        nstenergy = 0
        nstepKey = isMetaStrInDict("nsteps", self.cntrlDict)
        nstlogKey = isMetaStrInDict("nstlog", self.cntrlDict)
        nstxoutKey = isMetaStrInDict("nstxout", self.cntrlDict)
        nstenergyKey = isMetaStrInDict("nstenergy", self.cntrlDict)
        if nstepKey is not None:
            if self.cntrlDict[nstepKey].activeInfo:
                nsteps = conv_int(self.cntrlDict[nstepKey].value, default=0)
            else:
                nsteps = conv_int(self.cntrlDict[nstepKey].defaultValue, default=0)
        if nstlogKey is not None:
            if self.cntrlDict[nstlogKey].activeInfo:
                nlogsteps = conv_int(self.cntrlDict[nstlogKey].value, default=0)
            else:
                nlogsteps = conv_int(self.cntrlDict[nstlogKey].defaultValue, default=0)
        if nstxoutKey is not None:
            if self.cntrlDict[nstxoutKey].activeInfo:
                ntrajsteps = conv_int(self.cntrlDict[nstxoutKey].value, default=0)
            else:
                ntrajsteps = conv_int(self.cntrlDict[nstxoutKey].defaultValue, default=0)
        if nstenergyKey is not None:
            if self.cntrlDict[nstenergyKey].activeInfo:
                nenersteps = conv_int(self.cntrlDict[nstenergyKey].value, default=0)
            else:
                nenersteps = conv_int(self.cntrlDict[nstenergyKey].defaultValue, default=0)

        if nlogsteps>0:
            logsteps = [i for i in range(0, nsteps, nlogsteps)]
        else:
            logsteps = [0]
        if ntrajsteps>0:
            trajsteps = [i for i in range(0, nsteps, ntrajsteps)]
        else:
            trajsteps = [0]
        if nenersteps>0:
            enersteps = [i for i in range(0, nsteps, nenersteps)]
        else:
            enersteps = [0]
        logsteps.append(nsteps)
        self.MDlogsteps=logsteps
        self.MDlogstep = 0
        trajsteps.append(nsteps)
        enersteps.append(nsteps)
        self.stepcontrolDict.update({"logsteps"     : logsteps})
        self.stepcontrolDict.update({"nextlogsteps" : logsteps})
        self.stepcontrolDict.update({"trajsteps"    : trajsteps})
        self.stepcontrolDict.update({"enersteps"    : enersteps})
        #steps = logsteps if len(logsteps)>len(trajsteps) else trajsteps
        #followsteps = "log" if len(logsteps)>len(trajsteps) else "traj"
        #steps = logsteps
        steps = []
        for step in range(0,nsteps+1):
            if(step in logsteps or
               step in trajsteps or
               step in enersteps):
               steps.append(step)
        followsteps = "log"
        self.stepcontrolDict.update({"steps"  : steps if len(steps)>=len(enersteps) else enersteps})
        self.stepcontrolDict.update({"follow" : followsteps if len(steps)>=len(enersteps) else "ener"})

    def onOpen_section_method(self, backend, gIndex, section):
        # keep track of the latest method section
        self.secMethodGIndex = gIndex
        if self.inputMethodIndex is None:
            self.inputMethodIndex = gIndex
        else:
            backend.openNonOverlappingSection("section_method_to_method_refs")
            backend.addValue("method_to_method_kind", "core_settings")
            backend.addValue("method_to_method_ref", self.inputMethodIndex)
            backend.closeNonOverlappingSection("section_method_to_method_refs")
        if self.mainMethodIndex is None:
            self.mainMethodIndex = gIndex

    def onClose_section_method(self, backend, gIndex, section):
        """Trigger called when section_method is closed.
        """
        # input method
        pass

    def onOpen_section_sampling_method(self, backend, gIndex, section):
        # keep track of the latest sampling description section
        self.secSamplingGIndex = gIndex
        #self.inputControlIndex = backend.superBackend.openSection(PARSERTAG + '_section_control_parameters')

    def onClose_section_sampling_method(self, backend, gIndex, section):
        """Trigger called when section_sampling_method is closed.

        Writes sampling method details for minimization and molecular dynamics.
        """
        # check control keywords were found throguh dictionary support
        #backend.superBackend.closeSection(PARSERTAG + '_section_control_parameters', self.inputControlIndex)
        section_sampling_Dict = get_updateDictionary(self, 'sampling')
        updateDict = {
            #'startSection' : [['section_sampling_method']],
            'startSection' : [['section_run']],
            'dictionary' : section_sampling_Dict
            }
        self.metaStorage.update(updateDict)
        self.metaStorage.updateBackend(backend.superBackend,
                startsection=['section_run'],
                #startsection=['section_sampling_method'],
                autoopenclose=False)

    def onOpen_section_topology(self, backend, gIndex, section):
        # keep track of the latest topology description section
        if (gIndex is None or gIndex == -1 or gIndex == "-1"):
            self.secTopologyGIndex = backend.superBackend.openSection("section_topology")
        else:
            self.secTopologyGIndex = gIndex

    def onClose_section_topology(self, backend, gIndex, section):
        """Trigger called when section_topology is closed.
        """
        section_topology_Dict = get_updateDictionary(self, 'topology')
        updateDict = {
            'startSection' : [
                ['section_topology']],
            #'muteSections' : [['section_sampling_method']],
            'dictionary' : section_topology_Dict
            }
        self.metaStorage.update(updateDict)
        self.metaStorage.updateBackend(backend.superBackend,
                startsection=['section_topology'],
                autoopenclose=False)
        self.topology_atom_type_and_interactions(backend, gIndex)

        if (gIndex is None or gIndex == -1 or gIndex == "-1"):
            backend.superBackend.closeSection("section_topology", self.secTopologyGIndex)

    def onOpen_section_system(self, backend, gIndex, section):
        # keep track of the latest system description section
        if (gIndex is None or gIndex == -1 or gIndex == "-1"):
            self.secSystemGIndex = backend.superBackend.openSection("section_system")
        else:
            self.secSystemGIndex = gIndex

    def onClose_section_system(self, backend, gIndex, section):
        """Trigger called when section_system is closed.

        Writes atomic positions, atom labels and lattice vectors.
        """
        # Our special recipe for the confused backend and more
        # See our other recipes at the back of this parser: MetaInfoStorage! :)
        if (gIndex is None or gIndex == -1 or gIndex == "-1"):
            SloppyBackend = backend.superBackend
        else:
            SloppyBackend = backend

        numatoms = None

        if self.topology:
            if (self.secTopologyGIndex is None or
                (self.secTopologyGIndex == -1 or
                self.secTopologyGIndex == "-1")):
                self.onOpen_section_topology(backend, None, None)
                self.onClose_section_topology(backend, None, None)
            SloppyBackend.addValue("topology_ref", self.secTopologyGIndex)

        coordinates=None
        if self.trajectory is not None:
            coordinates=self.trajectory
            positions=self.atompositions
        elif(self.inputcoords is not None and
             self.MDcurrentstep == steps[0]):
            coordinates=self.inputcoords
            positions=self.inputpositions
        elif(self.outputcoords is not None and
             self.MDcurrentstep == steps[-1]):
            coordinates=self.outputcoords
            positions=self.outputpositions

        if(coordinates is not None and positions is not None):
            self.trajRefSingleConfigurationCalculation = gIndex
            if coordinates.unitcell_vectors is not None:
                unit_cell = np.asarray(self.metaStorage.convertUnits(
                    self.trajectory.unitcell_vectors, "Angstrom", self.unitDict))
                SloppyBackend.addArrayValues('simulation_cell', unit_cell)
                SloppyBackend.addArrayValues('lattice_vectors', unit_cell)
            if coordinates.unitcell_lengths is not None:
                SloppyBackend.addArrayValues(PARSERTAG + '_lattice_lengths', self.trajectory.unitcell_lengths)
            if coordinates.unitcell_angles is not None:
                SloppyBackend.addArrayValues(PARSERTAG + '_lattice_angles', self.trajectory.unitcell_angles)
            if(isinstance(positions, np.ndarray) or isinstance(positions, (tuple, list))):
                SloppyBackend.addArrayValues('atom_positions', np.transpose(np.asarray(
                    self.metaStorage.convertUnits(positions, "Angstrom", self.unitDict))))
            if coordinates.velocities is not None:
                SloppyBackend.addArrayValues('atom_velocities', np.transpose(np.asarray(
                    self.metaStorage.convertUnits(
                        coordinates.velocities, "Angstrom/(pico-second)", self.unitDict))))
            if self.topology is not None:
                atom_labels = self.topology.topoDict["atom_element_list"]
                numatoms = len(atom_labels)
                SloppyBackend.addArrayValues('atom_labels', np.asarray(atom_labels))
                pass

            # Read the next step at trajectory in advance
            # If iread returns None, it will be the last step
            try:
                self.atompositions = self.trajectory.positions()
                if self.atompositions is not None:
                    numatoms = len(self.atompositions)
                #self.MDiter += self.MDstep
                self.MDiter += 1
            except AttributeError:
                self.atompositions = None
                pass

            #if atom_vel:
            # need to transpose array since its shape is [number_of_atoms,3] in the metadata
            #backend.addArrayValues('atom_velocities', np.transpose(np.asarray(atom_vel)))

        if numatoms:
            SloppyBackend.addValue("number_of_atoms", numatoms)

        if (gIndex is None or gIndex == -1 or gIndex == "-1"):
            SloppyBackend.closeSection("section_system", self.secSystemGIndex)

    def onOpen_section_single_configuration_calculation(self, backend, gIndex, section):
        # write the references to section_method and section_system
        self.singleConfCalcs.append(gIndex)
        self.secSingleGIndex = backend.superBackend.openSection("section_single_configuration_calculation")

    def onClose_section_single_configuration_calculation(self, backend, gIndex, section):
        """Trigger called when section_single_configuration_calculation is closed.

        Write number of steps in MD or Minimization.
        Check for convergence of geometry optimization.
        Write energy values at MD and with error in Minimization.
        Write reference to section_method and section_system
        """

        self.lastCalculationGIndex = gIndex
        steps = self.stepcontrolDict["steps"]
        logsteps = self.stepcontrolDict["logsteps"]
        enersteps = self.stepcontrolDict["enersteps"]
        trajsteps = self.stepcontrolDict["trajsteps"]
        #if len(steps)>1:
        #    self.MDstep = steps[1]-steps[0]
        timestep = 0.0001
        istimestep = isMetaStrInDict("dt",self.cntrlDict)
        if istimestep is not None:
            if self.cntrlDict[istimestep].value is not None:
                timestep = float(self.cntrlDict[istimestep].value)
        if self.MDiter<len(steps):
            self.MDcurrentstep = steps[self.MDiter]
        else:
            self.MDcurrentstep += 1
        if self.MDiter+1<len(steps):
            self.MDnextstep = self.MDcurrentstep + steps[self.MDiter+1] - steps[self.MDiter]
        else:
            self.MDnextstep = steps[-1] + 1
        self.stepcontrolDict.update({"MDcurrentstep" : self.MDcurrentstep})
        self.stepcontrolDict.update({"MDcurrenttime" : self.MDcurrentstep*timestep})

        if(self.MDcurrentstep in logsteps or
           self.MDcurrentstep in enersteps):
            section_frameseq_Dict = get_updateDictionary(self, 'frameseq')
            updateFrameDict = {
                'startSection' : [
                    ['section_frame_sequence']],
                #'muteSections' : [['section_sampling_method']],
                'dictionary' : section_frameseq_Dict
                }
            self.metaStorage.update(updateFrameDict)
            section_singlevdw_Dict = get_updateDictionary(self, 'singlevdw')
            updateDictVDW = {
                'startSection' : [
                    ['section_energy_van_der_Waals']],
                #'muteSections' : [['section_sampling_method']],
                'dictionary' : section_singlevdw_Dict
                }
            self.secVDWGIndex = backend.superBackend.openSection("section_energy_van_der_Waals")
            self.metaStorage.update(updateDictVDW)
            self.metaStorage.updateBackend(backend.superBackend,
                    startsection=['section_energy_van_der_Waals'],
                    autoopenclose=False)
            backend.superBackend.closeSection("section_energy_van_der_Waals", self.secVDWGIndex)
            section_singlecalc_Dict = get_updateDictionary(self, 'singleconfcalc')
            updateDict = {
                'startSection' : [
                    ['section_single_configuration_calculation']],
                #'muteSections' : [['section_sampling_method']],
                'dictionary' : section_singlecalc_Dict
                }
            self.metaStorage.update(updateDict)
            self.metaStorage.updateBackend(backend.superBackend,
                    startsection=['section_single_configuration_calculation'],
                    autoopenclose=False)
        if self.MDcurrentstep in trajsteps:
            self.onOpen_section_system(backend, None, None)
            self.onClose_section_system(backend, None, None)
            #backend.addValue('single_configuration_to_calculation_method_ref', self.secMethodGIndex)
            backend.addValue('single_configuration_calculation_to_system_ref', self.secSystemGIndex)
        else:
            if(self.MDcurrentstep in logsteps or
               self.MDcurrentstep in enersteps):
                self.MDiter += 1
        if self.MDiter<len(steps):
            self.MDcurrentstep = steps[self.MDiter]
        else:
            self.MDcurrentstep += 1
        if self.MDiter+1<len(steps):
            self.MDnextstep = self.MDcurrentstep + steps[self.MDiter+1] - steps[self.MDiter]
        else:
            self.MDnextstep = steps[-1] + 1
        if self.MDcurrentstep in logsteps:
            self.MDlogstep = self.MDcurrentstep
            if len(self.MDlogsteps)>1:
                self.MDlogsteps.pop(0)
        if len(self.MDlogsteps)>0:
            self.MDnextlogstep=self.MDlogsteps[0]
        else:
            if len(logsteps)>0:
                self.MDnextlogstep=logsteps[-1]
        self.stepcontrolDict.update({"nextlogsteps" : logsteps})
        self.stepcontrolDict.update({"MDcurrentstep" : self.MDcurrentstep})
        self.stepcontrolDict.update({"MDcurrenttime" : self.MDcurrentstep*timestep})
        backend.superBackend.closeSection("section_single_configuration_calculation", self.secSingleGIndex)

    def setStartingPointCalculation(self, parser):
        backend = parser.backend
        backend.openSection('section_calculation_to_calculation_refs')
        if self.lastCalculationGIndex:
            backend.addValue('calculation_to_calculation_ref', self.lastCalculationGIndex)
            backend.closeSection('section_calculation_to_calculation_refs')
        return None

    def convertInt(self, fname, value):
        keyMapper = {
                "Step" : "MDcurrentstep",
                }
        if value is not None:
            for k,v in keyMapper.items():
                if fname in k:
                    return int(value)
        return value

    def convertFloat(self, fname, value):
        keyMapper = {
                "Step" : "MDcurrentstep",
                }
        if value is not None:
            for k,v in keyMapper.items():
                if fname in k:
                    return float(value)
        return value

    def updateTimeStep(self, fname, value):
        keyMapper = {
                "Time" : "MDcurrentstep",
                }
        timestep = 0.0001
        istimestep = isMetaStrInDict("dt",self.cntrlDict)
        if istimestep is not None:
            if self.cntrlDict[istimestep].value is not None:
                timestep = float(self.cntrlDict[istimestep].value)
        if value is not None:
            for k,v in keyMapper.items():
                if fname in k:
                    if timestep>0:
                        return int(value/timestep)
        return value

    def build_subMatchers(self):
        """Builds the sub matchers to parse the main output file.
        """
        mddataNameList=getList_MetaStrInDict(self.mddataDict)

        energyFilter = {
                "Proper Dih."    : "Proper-Dih.",
                "Improper Dih."  : "Improper-Dih.",
                "CMAP Dih."      : "CMAP-Dih.",
                "LJ (SR)"        : "LJ-(SR)",
                "LJ (LR)"        : "LJ-(LR)",
                "Disper. corr."  : "Disper.-corr.",
                "Coulomb (SR)"   : "Coulomb-(SR)",
                "Coul. recip."   : "Coul.-recip.",
                "Kinetic En."    : "Kinetic-En.",
                "Total Energy"   : "Total-Energy",
                "Conserved En."  : "Conserved-En.",
                "Pres. DC (bar)" : "Pres.-DC-(bar)",
                "Pressure (bar)" : "Pressure-(bar)",
                "Constr. rmsd"   : "Constr.-rmsd",
                }

        mdoutSubParsers = [
            # step Log Parser
            { "startReStr"        : r"\s*Step\s*Time\s*",
              "parser"            : "table_parser",
              "parsername"        : "log_step_parser",
              "waitlist"          : None,
              "stopOnMatchStr"    : r"(?:^\s*$|\s*Step\s*Time\s*|\s*Energies\s*)",
              "quitOnMatchStr"    : r"\s*(?:A\s*V\s*E\s*R\s*A\s*G\s*E\s*S|#####)\s*",
              "metaNameStart"     : PARSERTAG + "_mdout_",
              "matchNameList"     : mddataNameList,
              "matchNameDict"     : "mddataDict",
              "updateMatchDict"   : True,
              "onlyCaseSensitive" : True,
              "stopOnFirstLine"   : False,
              "parserOptions"     : {
                  "header"           : True,
                  "wrap"             : True,
                  "tablelines"       : 2,
                  "tablestartsat"    : r"\s*Step\s*Time\s*",
                  "tableendsat"      : r"^\s*$",
                  "lineFilter"       : None,
                  "movetostopline"   : False,
                  "parsercntrlattr"  : "MDcurrentstep",
                  "parsercntrlin"    : "MDlogsteps",
                  "lookupdict"       : "stepcontrolDict"
                  }
              },
            # energy EDR Parser
            { "startReStr"        : "AlwaysRun",
              "parser"            : "dictionary_parser",
              "parsername"        : "edr_energy_parser",
              "waitlist"          : None,
              "stopOnMatchStr"    : "AlwaysStop",
              "quitOnMatchStr"    : "AlwaysStop",
              "metaNameStart"     : PARSERTAG + "_mdout_",
              "matchNameList"     : mddataNameList,
              "matchNameDict"     : "mddataDict",
              "updateMatchDict"   : True,
              "onlyCaseSensitive" : True,
              "stopOnFirstLine"   : True,
              "parserOptions"     : {
                  "dictionary"       : "thermoDict",
                  "dicttype"         : "standard", # (standard or smartparser)
                  "readwritedict"    : "read",
                  "keyMapper"        : None,
                  "controlsections"  : ["section_single_configuration_calculation"],
                  "controlsave"      : "sectioncontrol",
                  "controldict"      : "stepcontrolDict",
                  "parsercntrlattr"  : "MDcurrentstep",
                  "parsercntrlin"    : "enersteps",
                  "lookupdict"       : "stepcontrolDict"
                  }
              },
            # time to step converter Parser
            { "startReStr"        : "AlwaysRun",
              "parser"            : "dictionary_parser",
              "parsername"        : "frameseq_step_parser",
              "waitlist"          : None,
              "stopOnMatchStr"    : "AlwaysStop",
              "quitOnMatchStr"    : "AlwaysStop",
              "metaNameStart"     : PARSERTAG + "_mdout_",
              "matchNameList"     : mddataNameList,
              "matchNameDict"     : "mddataDict",
              "updateMatchDict"   : True,
              "onlyCaseSensitive" : True,
              "stopOnFirstLine"   : True,
              "parserOptions"     : {
                  "dictionary"       : "mddataDict",
                  "dicttype"         : "smartparser", # (standard or smartparser)
                  "readwritedict"    : "read",
                  "keyMapper"        : {"Step" : "Time"},
                  "updatefunc"       : "max",
                  "updateattrs"      : ["MDcurrentstep"],
                  "preprocess"       : self.convertFloat,
                  "postprocess"      : self.updateTimeStep,
                  "parsercntrlattr"  : "MDcurrentstep",
                  "parsercntrlin"    : "enersteps",
                  "lookupdict"       : "stepcontrolDict"
                  }
              },
            # energy Log Parser
            { "startReStr"        : r"\s*Energies\s*",
              "parser"            : "table_parser",
              "parsername"        : "log_energy_parser",
              "waitlist"          : [["log_step_parser"]],
              "stopOnMatchStr"    : r"\s*Step\s*Time\s*",
              "quitOnMatchStr"    : r"\s*(?:A\s*V\s*E\s*R\s*A\s*G\s*E\s*S|#####)\s*",
              "metaNameStart"     : PARSERTAG + "_mdout_",
              "matchNameList"     : mddataNameList,
              "matchNameDict"     : "mddataDict",
              "updateMatchDict"   : True,
              "onlyCaseSensitive" : True,
              "stopOnFirstLine"   : False,
              "parserOptions"     : {
                  "header"           : True,
                  "wrap"             : True,
                  "tablelines"       : 2,
                  "tablestartsat"    : r"\s*Energies\s*",
                  "tableendsat"      : r"^\s*$",
                  "skiplines"        : 1,
                  "movetostopline"   : False,
                  "lineFilter"       : energyFilter,
                  #"parsercntrlattr"  : "MDcurrentstep",
                  #"parsercntrlin"    : "MDlogsteps",
                  "controlsections"  : ["section_single_configuration_calculation"],
                  "controlsave"      : "sectioncontrol",
                  "controldict"      : "stepcontrolDict",
                  #"lookupdict"       : "stepcontrolDict"
                  }
              },
            # thermostat save Parser
            { "startReStr"        : "AlwaysRun",
              "parser"            : "dictionary_parser",
              "parsername"        : "frameseq_step_parser",
              "waitlist"          : None,
              "stopOnMatchStr"    : "AlwaysStop",
              "quitOnMatchStr"    : "AlwaysStop",
              "metaNameStart"     : PARSERTAG + "_mdout_",
              "matchNameList"     : mddataNameList,
              "matchNameDict"     : "mddataDict",
              "updateMatchDict"   : True,
              "onlyCaseSensitive" : True,
              "stopOnFirstLine"   : True,
              "parserOptions"     : {
                  "dictionary"       : "stepcontrolDict",
                  "dicttype"         : "standard", # (standard or smartparser)
                  "readwritedict"    : "write",
                  "keyMapper"        : {"Step" : "MDcurrentstep"},
                  "updatefunc"       : "max",
                  "updateattrs"      : ["MDcurrentstep"],
                  #"controlsections"  : None,
                  "controlsections"  : ["section_single_configuration_calculation"],
                  "controlsave"      : "sectioncontrol",
                  "controldict"      : "stepcontrolDict",
                  "controlattrs"     : ["MDcurrentstep"],
                  "preprocess"       : self.convertInt,
                  "postprocess"      : self.convertInt,
                  #"parsercntrlattr"  : "MDcurrentstep",
                  #"parsercntrlin"    : "steps",
                  #"lookupdict"       : "stepcontrolDict"
                  }
              },
            # Section Control Parser
            { "startReStr"        : "AlwaysRun",
              "parser"            : "section_control_parser",
              "parsername"        : "section_parser",
              "waitlist"          : [["log_step_parser", "log_energy_parser"],
                                     ["edr_energy_parser"]],
              "stopOnMatchStr"    : r"\s*(?:Step\s*Time|Energies)\s*",
              "quitOnMatchStr"    : r"\s*(?:Step\s*Time|Energies)\s*",
              "metaNameStart"     : None,
              "matchNameList"     : None,
              "matchNameDict"     : None,
              "updateMatchDict"   : False,
              "onlyCaseSensitive" : True,
              "stopOnFirstLine"   : True,
              "parserOptions"     : {
                  "sectionname"      : "section_single_configuration_calculation",
                  "sectionopen"      : True,
                  "sectionopenattr"  : "MDcurrentstep",
                  "sectionopenin"    : "steps",
                  "sectionclose"     : True,
                  "sectioncloseattr" : "MDcurrentstep",
                  "sectionclosein"   : "steps",
                  "activatesection"  : "sectioncontrol",
                  "lookupdict"       : "stepcontrolDict"
                  }
              },
            # Step Control Parser
            { "startReStr"        : "AlwaysRun",
              "parser"            : "readline_control_parser",
              "parsername"        : "readline_parser",
              "waitlist"          : None,
              "stopOnMatchStr"    : "AlwaysStop",
              "quitOnMatchStr"    : "AlwaysStop",
              "metaNameStart"     : None,
              "matchNameList"     : None,
              "matchNameDict"     : None,
              "updateMatchDict"   : False,
              "onlyCaseSensitive" : True,
              "stopOnFirstLine"   : True,
              "parserOptions"     : {
                  "waitatlineStr"    : r"\s*(?:Step\s*Time|Energies)\s*",
                  "controlwait"      : "nextlogsteps",
                  "controlattr"      : "MDlogstep",
                  "controlnextattr"  : "MDnextlogstep",
                  "controllast"      : -1,
                  "controlskip"      : [],
                  "controlin"        : "logsteps",
                  "controlcounter"   : "targetstep",
                  "controldict"      : "stepcontrolDict",
                  "lookupdict"       : "stepcontrolDict"
                  }
              }
            ]

        archNameList=getList_MetaStrInDict(self.metaDicts['arch'])
        cntrlNameList=getList_MetaStrInDict(self.metaDicts['cntrl'])
        qmNameList=getList_MetaStrInDict(self.metaDicts['qm'])
        grpNameList=getList_MetaStrInDict(self.metaDicts['grp'])
        annealNameList=getList_MetaStrInDict(self.metaDicts['anneal'])

        mdinoutSubParsers = [
            # arch Parser
            { "startReStr"        : r"\s*GROMACS\s*version:",
              "parser"            : "namelist_parser",
              "parsername"        : "arch_parser",
              "stopOnMatchStr"    : r"(?:[+]*\s*PLEASE\s*READ\s*AND\s*"
                                     "CITE\s*THE\s*FOLLOWING\s*REFERENCE\s*"
                                     "[+]*|Running\s*on)\s*",
              "quitOnMatchStr"    : None,
              "metaNameStart"     : PARSERTAG + "_inout_control_",
              "matchNameList"     : archNameList,
              "matchNameDict"     : "archDict",
              "updateMatchDict"   : True,
              "onlyCaseSensitive" : True,
              "stopOnFirstLine"   : False,
              "parserOptions"     : { "lineFilter" : None }
              },
            # cntrl Parser
            { "startReStr"        : r"\s*Input\s*Parameters:",
              "parser"            : "namelist_parser",
              "parsername"        : "cntrl_parser",
              "stopOnMatchStr"    : r"(?:qm-opts:|grpopts:|annealing:|Initializing)\s*",
              "quitOnMatchStr"    : None,
              "metaNameStart"     : PARSERTAG + "_inout_control_",
              "matchNameList"     : cntrlNameList,
              "matchNameDict"     : "cntrlDict",
              "updateMatchDict"   : True,
              "onlyCaseSensitive" : True,
              "stopOnFirstLine"   : False,
              "parserOptions"     : { "lineFilter" : None }
              },
            # qm Parser
            { "startReStr"        : r"\s*qm-opts:\s*",
              "parser"            : "namelist_parser",
              "parsername"        : "qm_parser",
              "stopOnMatchStr"    : r"(?:grpopts:|annealing:|Initializing)\s*",
              "quitOnMatchStr"    : None,
              "metaNameStart"     : PARSERTAG + "_inout_control_",
              "matchNameList"     : qmNameList,
              "matchNameDict"     : "qmDict",
              "updateMatchDict"   : True,
              "onlyCaseSensitive" : True,
              "stopOnFirstLine"   : False,
              "parserOptions"     : { "lineFilter" : None }
              },
            # grp Parser
            { "startReStr"        : r"\s*grpopts:\s*",
              "parser"            : "namelist_parser",
              "parsername"        : "grp_parser",
              "stopOnMatchStr"    : r"(?:annealing:|Initializing)\s*",
              "quitOnMatchStr"    : None,
              "metaNameStart"     : PARSERTAG + "_inout_control_",
              "matchNameList"     : grpNameList,
              "matchNameDict"     : "grpDict",
              "updateMatchDict"   : True,
              "onlyCaseSensitive" : True,
              "stopOnFirstLine"   : False,
              "parserOptions"     : { "lineFilter" : None }
              },
            # anneal Parser
            { "startReStr"        : r"\s*annealing:\s*",
              "parser"            : "namelist_parser",
              "parsername"        : "anneal_parser",
              "stopOnMatchStr"    : r"(?:qm-opts:|grpopts:|Initializing)\s*",
              "quitOnMatchStr"    : None,
              "metaNameStart"     : PARSERTAG + "_inout_control_",
              "matchNameList"     : annealNameList,
              "matchNameDict"     : "annealDict",
              "updateMatchDict"   : True,
              "onlyCaseSensitive" : True,
              "stopOnFirstLine"   : False,
              "parserOptions"     : { "lineFilter" : None }
              }
            ]

        ########################################
        # return main Parser
        return [
            SM(name='NewRun',
                startReStr=r"(?:Log\sfile\sopened\son\s*|"
                            "Host:\s*[a-zA-Z0-9:.]+\s*|"
                            "pid:\s*[0-9]+\s*rank\s*ID:\s*[0-9].\s*|"
                            "number\s*of\s*ranks:\s*[0-9].\s*|"
                            "\s*GROMACS\s*-\s*gmx)",
                endReStr=r"\s*Finished\s*mdrun\s*on\s*rank\s*",
                repeats=True,
                required=True,
                forwardMatch=True,
                sections=['section_run'],
                fixedStartValues={'program_name': PROGRAMNAME},
                subMatchers=[
                    SM(name='logruntime',
                       startReStr=r"\s*Log\sfile\sopened\son\s"
                                   "\s*(?P<"+PARSERTAG+"_mdin_finline>[a-zA-Z0-9:.\s]+)\s*$",
                       adHoc=lambda p: p.backend.addValue(
                           "time_run_date_start", datetime.datetime.strptime(
                               p.lastMatch[PARSERTAG+"_mdin_finline"].replace("\n",""),
                           '%a %b %d %H:%M:%S %Y').timestamp())),
                    SM(name='loghostinfo',
                       startReStr=r"\s*Host:\s*(?P<"+PARSERTAG+"_program_execution_host>[a-zA-Z0-9:.]+)\s*"
                                   "pid:\s*[0-9]+\s*rank\s*ID:\s*(?P<"+PARSERTAG+"_parallel_task_nr>[0-9]+)\s*"
                                   "number\s*of\s*ranks:\s*(?P<"+PARSERTAG+"_number_of_tasks>[0-9]+)\s*$"),
                    SM(name='programheader',
                       startReStr=r"\s*(:-\))?\s*GROMACS\s*-\s*(?P<"+PARSERTAG+"_program_module>[0-9a-zA-Z\s]+)\s*"
                                   ",\s*VERSION\s*(?P<program_version>[0-9.]+)\s*"),
                    # header specifing version, compilation info, task assignment
                    SM(name='license',
                       startReStr=r"\s*GROMACS\s*is\s*written\s*by:\s*",
                       coverageIgnore=True,
                       adHoc=lambda p:
                       self.adHoc_read_store_text_stop_parsing(p,
                           stopOnMatchStr=r"\s*License,\s*or\s*[(]at\s*"
                                           "your\s*option[)]\s*any\s*"
                                           "later\s*version[.]\s*",
                           quitOnMatchStr=None,
                           metaNameStart=PARSERTAG+"_",
                           metaNameStore=PARSERTAG+"_program_license",
                           matchNameList=None,
                           matchNameDict=None,
                           onlyCaseSensitive=True,
                           stopOnFirstLine=False,
                           storeFirstLine=True,
                           storeStopQuitLine=True,
                           onQuitRunFunction=lambda p: p.backend.addValue(
                               PARSERTAG+"_program_license",
                               ' '.join(','.join(re.split(r'\s{2,}',
                                   p.lastMatch[PARSERTAG+"_program_license"])
                                   ).replace('\n', ' ').split()))
                           )
                       ),
                    SM(name='module_version',
                       startReStr=r"\s*GROMACS:\s*(?P<"+PARSERTAG+"_mdin_finline>"
                                   "[a-zA-Z0-9,.]+[\s][a-zA-Z0-9,.]+[\s][a-zA-Z0-9,.]+[\s]"
                                   "[a-zA-Z0-9,.]+[\s]?[a-zA-Z0-9,.]+[\s])\s*",
                       adHoc=lambda p: p.backend.addValue(
                           PARSERTAG+"_program_module_version", p.lastMatch[
                               PARSERTAG+"_mdin_finline"].replace("\n",""))),
                    SM(name='executable_path',
                       startReStr=r"\s*Executable:\s*"
                                   "(?P<"+PARSERTAG+"_mdin_finline>"
                                   "[a-zA-Z0-9/._\s]+)\s*",
                       adHoc=lambda p: p.backend.addValue(
                           PARSERTAG+"_program_execution_path", p.lastMatch[
                               PARSERTAG+"_mdin_finline"].replace("\n",""))),
                    SM(name='working_path',
                       startReStr=r"\s*Data\s*prefix:\s*"
                                   "(?P<"+PARSERTAG+"_mdin_finline>"
                                   "[a-zA-Z0-9/._\s]+)\s*",
                       adHoc=lambda p: p.backend.addValue(
                           PARSERTAG+"_program_working_path", p.lastMatch[
                               PARSERTAG+"_mdin_finline"].replace("\n",""))),
                    SM(name='FileNameMatch',
                       startReStr=r"\s*Command\s*line:",
                       forwardMatch=True,
                       sections=[PARSERTAG+'_section_input_output_files'],
                       coverageIgnore=True,
                       adHoc=lambda p:
                       self.adHoc_read_commandline_stop_parsing(p,
                           stopOnMatchStr=r"\s*gmx\s*",
                           quitOnMatchStr=None,
                           metaNameStart=PARSERTAG + "_",
                           matchNameList="matchNames",
                           matchNameDict=self.fileDict,
                           onlyCaseSensitive=True,
                           stopOnFirstLine=False,
                           commandLineMatchStr=r"\s*(?:gmx|gmx_mpi)\s*mdrun\s*(?P<cmdargs>.*)",
                           commandLineFunction=get_commandLineArguments)
                       ),
                    SM(name='SectionControlParm',
                       startReStr=r"(?:\s*GROMACS\s*version:|\s*Input\s*Parameters:)",
                       endReStr=r"\s*Started\s*mdrun\s*on",
                       forwardMatch=True,
                       sections=['section_sampling_method', PARSERTAG + '_section_control_parameters'],
                       adHoc=lambda p:
                       self.adHoc_takingover_parsing(p,
                           stopOnMatchStr=r"\s*Started\s*mdrun\s*on",
                           quitOnMatchStr=None,
                           stopControl="stopControl", # if None or True, stop with quitMatch, else wait
                           record=False, # if False or None, no record, no replay
                           replay=0, # if 0 or None= no replay, if <0 infinite replay
                           parseOnlyRecorded=False, # if True, parsers only work on record
                           ordered=False,
                           subParsers=mdinoutSubParsers)
                       ), # END SectionControlParm
                    SM(name='SingleConfigurationCalculationWithSystemDescription',
                       startReStr=r"\s*Step\s*Time\s*",
                       endReStr=r"\s*(?:A\s*V\s*E\s*R\s*A\s*G\s*E\s*S\s*|#####)\s*",
                       forwardMatch=True,
                       subMatchers=[
                           # the actual section for a single configuration calculation starts here
                           SM(name='SingleConfigurationCalculation',
                              startReStr=r"\s*Step\s*Time\s*",
                              endReStr=r"\s*(?:A\s*V\s*E\s*R\s*A\s*G\s*E\s*S\s*|#####)\s*",
                              repeats=True,
                              forwardMatch=True,
                              #sections=['section_single_configuration_calculation'],
                              subMatchers=[
                                  SM(name='MDStep',
                                     startReStr=r"\s*(?:Step\s*Time|Energies)\s"
                                     "(?P<"+PARSERTAG+"_mdin_finline>.*)(?:'|\")?\s*,?",
                                     forwardMatch=True,
                                     adHoc=lambda p:
                                     self.adHoc_takingover_parsing(p,
                                         stopOnMatchStr=r"\s*(?:A\s*V\s*E\s*R\s*A\s*G\s*E\s*S|#####)\s*",
                                         quitOnMatchStr=r"\s*(?:A\s*V\s*E\s*R\s*A\s*G\s*E\s*S|#####)\s*",
                                         stopControl="stopControl", # if None or True, stop with quitMatch, else wait
                                         record=False, # if False or None, no record, no replay
                                         replay=0, # if 0 or None= no replay, if <0 infinite replay
                                         parseOnlyRecorded=False, # if True, parsers only work on record
                                         ordered=False,
                                         onlySubParsersReadLine=True,
                                         subParsers=mdoutSubParsers)
                                     )]
                              ) # END SingleConfigurationCalculation
                       ]), # END SingleConfigurationCalculationWithSystemDescription
                    SM(name='Dashes', startReStr=r"\s*#####\s*"),
                    SM(name='Mflops', startReStr=r"\s*(?:M\s*E\s*G\s*A\s*-\s*F\s*L\s*O\s*P\s*S\s*"
                                                  "|A\s*C\s*C\s*O\s*U\s*N\s*T\s*I\s*N\s*G\s*)"),
                    SM(name='TimeCount', startReStr=r"\s*R\sE\sA\sL\s*C\sY\sC\sL\sE\s*A\sN\sD\s*"
                                                     "T\sI\sM\sE\s*A\sC\sC\sO\sU\sN\sT\sI\sN\sG\s*"),
                    # summary of computation
                    SM(name='ComputationTimings',
                       startReStr=r"\s*Computing:\s*Num",
                       subMatchers=[
                           SM(r"\s*Time:\s*[0-9:.eEdD]+\s*[0-9.:eEdD]+\s*[0-9:.]+"),
                       ]), # END Computation
                    SM(name='end_run',
                       startReStr=r"\s*Finished\s*mdrun\s*on\s*rank\s*[0-9]+\s*"
                                   "(?P<"+PARSERTAG+"_mdin_finline>[a-zA-Z0-9:.\s]+)\s*",
                       adHoc=lambda p: p.backend.addValue(
                           "time_run_date_end", datetime.datetime.strptime(
                               p.lastMatch[PARSERTAG+"_mdin_finline"].replace("\n",""),
                           '%a %b %d %H:%M:%S %Y').timestamp())),
                ]) # END NewRun
            ]


class GromacsParserInterface():
   """ A proper class envolop for running this parser from within python. """
   def __init__(self, backend, **kwargs):
       self.backend_factory = backend

   def parse(self, mainfile):
        from unittest.mock import patch
        logging.info('gromacs parser started')
        logging.getLogger('nomadcore').setLevel(logging.WARNING)
        backend = self.backend_factory("gromacs.nomadmetainfo.json")
        parserInfo = {'name': 'gromacs-parser', 'version': '1.0'}
        context = GromacsParser()
        context.coverageIgnore = re.compile(r"^(?:" + r"|".join(context.coverageIgnoreList) + r")$")
        with patch.object(sys, 'argv', ['<exe>', '--uri', 'nmd://uri', mainfile]):
            mainFunction(
                mainFileDescription=context.mainFileDescription(),
                metaInfoEnv=None,
                parserInfo=parserInfo,
                cachingLevelForMetaName=context.cachingLevelForMetaName,
                superContext=context,
                superBackend=backend)
        return backend


if __name__ == "__main__":
    parser = GromacsParser()
    parser.parse()
