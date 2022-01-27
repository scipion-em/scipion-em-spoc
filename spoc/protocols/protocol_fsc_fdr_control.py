# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

from os.path import abspath
import numpy as np

from pwem.emlib.image import ImageHandler
from pwem.objects import FSC
from pwem.protocols import ProtAnalysis3D

from pyworkflow.protocol import PointerParam, BooleanParam, FloatParam, IntParam, StringParam
from pyworkflow import BETA
import pyworkflow.utils as pwutils

import spoc


class ProtFscFdrControl(ProtAnalysis3D):
    """
    Significance analysis of FSC curves
    """
    _label = 'FSC-FDR control'
    _devStatus = BETA

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('halfOne', PointerParam, pointerClass="Volume",
                      label='First half map', important=True)
        form.addParam('halfTwo', PointerParam, pointerClass="Volume",
                      label='Second half map', important=True)
        form.addParam('localRes', BooleanParam, default=False, label='Compute local resolution?')
        form.addParam('lowRes', FloatParam, default=-1,
                      label='Set lowest resolution', help='If set to -1, this parameter will not be used')
        form.addParam('stepSize', IntParam, default=-1,
                      label='Step size',
                      help='Voxels to skip for local resolution estimation. '
                           'If set to -1, this parameter will not be used')
        form.addParam('sym', StringParam, default='c1', label='Volume symmetry',
                      help='Symmetry for correction of symmetry effects')
        form.addParam('numAsymUnits', IntParam, default=-1,
                      label='Number of asymmetric units',
                      help='Number of asymmetric units for correction of symmetry effects. '
                           'If set to -1, this parameter will not be used')
        form.addParam('bfactor', FloatParam, default=-1,
                      label='B-Factor',
                      help='B-Factor for sharpening of the map. '
                           'If set to -1, this parameter will not be used')
        form.addParam('mask', PointerParam, pointerClass="Volume", allowsNull=True,
                      label='Mask for local map-model FSC calculation')

    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.computeControlStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -------------------------------
    def convertInputStep(self):
        ih = ImageHandler()
        file_halfOne = self.halfOne.get().getFileName()
        file_halfTwo = self.halfTwo.get().getFileName()
        if pwutils.getExt(file_halfOne) == '.mrc':
            file_halfOne += ':mrc'
        if pwutils.getExt(file_halfTwo) == '.mrc':
            file_halfTwo += ':mrc'
        ih.convert(file_halfOne, self._getExtraPath('halfone.mrc'))
        ih.convert(file_halfTwo, self._getExtraPath('halftwo.mrc'))

    def computeControlStep(self):
        args = '--halfmap1 halfone.mrc --halfmap2 halftwo.mrc' \
               ' --apix %f --symmetry %s' % (self.halfOne.get().getSamplingRate(),
                                             self.sym.get().upper())

        if self.localRes.get():
            args += ' -localResolutions'

        if self.lowRes.get() >= 0:
            args += ' -lowRes %f' % self.lowRes.get()

        if self.stepSize.get() >= 0:
            args += ' --window_size %d' % self.stepSize.get()

        if self.numAsymUnits.get() >= 0:
            args += ' --numAsymUnits %d' % self.numAsymUnits.get()

        if self.bfactor.get() >= 0:
            args += ' --bFactor %f' % self.bfactor.get()

        if self.mask.get():
            args += ' --mask %s' % abspath(self.mask.get().getFileName())

        program = spoc.Plugin.getProgram("FSC_FDRcontrol.py")
        self.runJob(program, args, cwd=self._getExtraPath())

    def createOutputStep(self):
        _fsc = FSC(objLabel='FSC')
        data = np.loadtxt(self._getExtraPath('FSC.txt'))
        _fsc.setData(data[0].tolist(), data[1].tolist())
        self._defineOutputs(outputFSC=_fsc)
        self._defineSourceRelation(self.halfOne, _fsc)
        self._defineSourceRelation(self.halfTwo, _fsc)

    # --------------------------- INFO functions ------------------------------
    def _methods(self):
        methods = []
        methods.append('Significance analysis of FSC curves')
        return methods

    def _summary(self):
        summary = []
        if not self.isFinished():
            summary.append("FDR-FSC information not ready yet.")

        if self.getOutputsSize() >= 1:
            stdout = self._getLogsPath('run.stdout')
            with open(stdout) as file:
                for num, line in enumerate(file, 1):
                    if 'Resolution at 1 % FDR-FSC' in line:
                        res = [float(s) for s in line.split() if s.replace(".", "", 1).isdigit()][1]
            summary.append('Resolution at 1 %% FDR-FSC: %.2f Angstrom' % res)
        return summary
