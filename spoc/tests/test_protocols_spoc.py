# **************************************************************************
# *
# * Authors:    David Herreros Calero (dherreros@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnologia, CSIC, Spain
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************


from pwem.protocols import ProtImportVolumes

from pyworkflow.tests import BaseTest, setupTestProject, DataSet

from spoc.protocols import ProtResolutionAnalysisFSCFDR, ProtConfidenceMap


class TestFscFdrControl(BaseTest):
    """This class check if the FSC-FDR control protocol works properly"""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('relion_tutorial')
        cls.volume = cls.dataset.getFile('volumes/reference_masked.vol')
        cls.halfOne = cls.runImportVolumes(cls, cls.volume, 1, 'Halfmap One')
        cls.halfTwo = cls.runImportVolumes(cls, cls.volume, 1, 'Halfmap One')

    def runImportVolumes(cls, pattern, samplingRate, label):
        """ Run an Import particles protocol. """
        protImport = cls.newProtocol(ProtImportVolumes,
                                     filesPath=pattern, samplingRate=samplingRate, objLabel=label)
        cls.launchProtocol(protImport)
        return protImport.outputVolume

    def runFscFdrControl(self, halfOne, halfTwo, localRes, label):
        prot = self.newProtocol(ProtResolutionAnalysisFSCFDR, halfOne=halfOne, halfTwo=halfTwo,
                                localRes=localRes, objLabel=label)
        self.launchProtocol(prot)
        if localRes:
            self.assertIsNotNone(prot.outputLocalResMap,
                                 "There was a problem with FSC-FdR protocol output (Local resolution)")
        else:
            self.assertIsNotNone(prot.outputFSC,
                                 "There was a problem with FSC-FdR protocol output (Global FSC)")
        return prot

    def runConfidenceMap(self, map, label):
        prot = self.newProtocol(ProtConfidenceMap, inputMap=map, locResFilter=False, objLabel=label)
        self.launchProtocol(prot)
        self.assertIsNotNone(prot.confidenceMap,
                             "There was a problem with Confidence Map protocol output")
        self.assertIsNotNone(prot.confidenceMap_log10FDR,
                             "There was a problem with Confidence Map protocol output")
        return prot

    def test_fsc_fdr_control_fsc(self):
        prot = self.runFscFdrControl(self.halfOne, self.halfTwo, False, 'Global FSC')
        stdout = prot._getLogsPath('run.stdout')
        with open(stdout) as file:
            for num, line in enumerate(file, 1):
                if 'Resolution at 1 % FDR-FSC' in line:
                    res = [float(s) for s in line.split() if s.replace(".", "", 1).isdigit()][1]
        self.assertEqual(res, 2.0, msg="Unexpected resolution value")
        return prot

    def test_fsc_fdr_control_locres(self):
        # TODO: Add extra checks (probably comparing to an already saved result?)
        prot = self.runFscFdrControl(self.halfOne, self.halfTwo, True, 'Local resolution')
        return prot

    def test_confidence_map(self):
        # TODO: Add extra checks (probably comparing to an already saved result?)
        prot = self.runConfidenceMap(self.halfOne, 'Confidence map')
        return prot