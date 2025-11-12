# **************************************************************************
# *
# * Authors:     J.L.Vilas (jlvilas@cnb.csic.es)
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
import os.path
from os.path import abspath
import numpy as np

from pwem.emlib.image import ImageHandler
from pwem.objects import FSC, Volume
from pwem.protocols import ProtAnalysis3D

from pyworkflow.protocol import PointerParam, BooleanParam, FloatParam
from pyworkflow import BETA
import pyworkflow.utils as pwutils

import spoc

INPUT_MAP = 'inputMap.mrc'
OUTPUT_MAP = '_confidenceMap.mrc'
OUTPUT_LOG10 = '_-log10FDR.mrc'
FILTERED_MAP = '_locFilt.mrc'

class ProtConfidenceMap(ProtAnalysis3D):
    """
    Confidence maps of cryoEM reconstructions
    """
    _label = 'confidence maps'
    _devStatus = BETA

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMap', PointerParam, pointerClass='Volume',
                      label="Input Volume", important=True,
                      help='Select a volume for determining its '
                           'confidence map.')
        line = form.addLine('Center of the noise box',
                             help='Size of the noise estimation region can be adjusted '
                                  'with -w sizeOfRegion. If the default regions for '
                                  'noise estimation fall on top of the molecule, '
                                  'you can specify the center of region of your '
                                  'choice with -noiseBox x y z')
        line.addParam('x_center', FloatParam, allowsNull=True, label='x')
        line.addParam('y_center', FloatParam, allowsNull=True, label='y')
        line.addParam('z_center', FloatParam, allowsNull=True, label='z')
        form.addParam('box', FloatParam, allowsNull=True,
                      label="Size of the noise box",
                      help='If the default regions for '
                          'noise estimation fall on top of the molecule, '
                          'you can specify the center of region of your '
                          'choice with -noiseBox x y z')
        form.addParam('locResFilter', BooleanParam,
                      default=True,
                      label="Filter the map with local resolution?",
                      help='A corresponding local resolution map can be supplied '
                           'to the program. The ouput will be a locally filtered '
                           'and locally normalized map together with the diagnostic'
                           ' image and the Confidence Map. Adjustment of the noise'
                           ' estimation region can be done accordingly and is '
                           'important for accurate background noise estimation.')
        form.addParam('resMap', PointerParam, pointerClass='Volume',
                      condition='locResFilter',
                      label="Local Resolution Map", important=True,
                      help='A corresponding local resolution map can be supplied '
                           'to the program. The ouput will be a locally filtered '
                           'and locally normalized map together with the diagnostic'
                           ' image and the Confidence Map. Adjustment of the noise'
                           ' estimation region can be done accordingly and is '
                           'important for accurate background noise estimation.')



    # --------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.computeConfidenceMapStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -------------------------------
    def convertInputStep(self):
        ih = ImageHandler()
        file_map = self.inputMap.get().getFileName()
        if pwutils.getExt(file_map) == '.mrc':
            file_map += ':mrc'
        ih.convert(file_map, self._getTmpPath(INPUT_MAP))


    def computeConfidenceMapStep(self):
        path_inputMap = os.path.abspath(self._getTmpPath(INPUT_MAP))
        args = ' -em %s ' % path_inputMap
        args += ' -p %f ' % self.inputMap.get().getSamplingRate()
        if self.box.hasValue():
            args +=' -w %i' % self.box.get()
            if self.x_center.hasValue() and self.y_center.hasValue() \
                    and self.z_center.hasValue():
                args += ' -noiseBox %i %i %i' % (self.x_center.get(), self.y_center.get(),
                                             self.z_center.get())
        if self.locResFilter:
            args += ' -locResMap %s ' % os.path.abspath(self.resMap.get().getFileName())

        program = spoc.Plugin.getProgram("FDRcontrol.py")
        self.runJob(program, args, cwd=self._getExtraPath())

    def createOutputStep(self):

        confMap = Volume()
        fnbase = pwutils.removeBaseExt(INPUT_MAP)
        confMap.setFileName(self._getExtraPath(fnbase+OUTPUT_MAP))
        confMap.setSamplingRate(self.inputMap.get().getSamplingRate())
        self._defineOutputs(confidenceMap=confMap)
        self._defineSourceRelation(self.inputMap, confMap)

        confMaplog = Volume()
        fnbaseout = pwutils.removeBaseExt(OUTPUT_MAP) + OUTPUT_LOG10
        confMaplog.setFileName(self._getExtraPath(fnbase+fnbaseout))
        confMaplog.setSamplingRate(self.inputMap.get().getSamplingRate())
        self._defineOutputs(confidenceMap_log10FDR=confMaplog)
        self._defineSourceRelation(self.inputMap, confMaplog)

        if self.locResFilter:
            filtMap = Volume()
            fnbaseout = pwutils.removeBaseExt(OUTPUT_MAP)

            filtMap.setFileName(self._getExtraPath(fnbase + FILTERED_MAP))
            filtMap.setSamplingRate(self.inputMap.get().getSamplingRate())
            self._defineOutputs(localfiltMap=filtMap)
            self._defineSourceRelation(self.inputMap, filtMap)

    # --------------------------- INFO functions ------------------------------
    def _methods(self):
        methods = []
        methods.append('Confidence Map estimation')
        return methods

    def _summary(self):
        summary = []
        if not self.isFinished():
            summary.append("Confidence map not ready yet.")
        else:
            summary.append("Confidence map estimated")

    def _validate(self):
        errors = []
        if not self.inputMap.get():
            errors.append("You must provide an input map")

        return errors
