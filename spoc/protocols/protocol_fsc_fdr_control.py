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
import os.path
from os.path import abspath
import numpy as np

from pwem.emlib.image import ImageHandler
from pwem.objects import FSC, Volume
from pwem.protocols import ProtAnalysis3D

from pyworkflow.protocol import PointerParam, BooleanParam, FloatParam, IntParam, StringParam
from pyworkflow import BETA
import pyworkflow.utils as pwutils

import spoc


class ProtFscFdrControl(ProtAnalysis3D):
    """
    Significance analysis of FSC curves
    """

    """
        FSC-FDR Control (ProtFscFdrControl) — User Manual

            Overview

            The FSC-FDR Control protocol performs significance analysis of
            Fourier Shell Correlation (FSC) curves using false discovery
            rate (FDR) statistics. Its main purpose is to provide a more
            statistically robust estimation of cryo-EM resolution by
            evaluating the significance of correlations between two half
            maps.

            Inputs and General Workflow

            The protocol requires two half maps representing independent
            reconstructions of the same structure. These maps can either
            be provided directly as separate volumes or extracted from a
            single input volume containing embedded half maps. Optionally,
            a mask can be included to restrict the analysis to biologically
            relevant regions during local resolution estimation.

            During execution, the input maps are converted into compatible
            formats and processed using the FSC_FDRcontrol program. The
            protocol computes statistically significant FSC curves while
            taking into account symmetry, sampling rate, and optional map
            sharpening parameters.

            Resolution Significance Analysis

            The protocol supports both global FSC significance analysis and
            local resolution estimation. In the standard mode, the method
            evaluates the FSC curve and determines the resolution threshold
            based on FDR statistics, providing a more reliable interpretation
            of map resolution than fixed FSC cutoffs alone.

            When local resolution analysis is enabled, the protocol computes
            voxel-dependent resolution estimates across the volume. Additional
            parameters such as step size, lowest resolution limit, symmetry,
            number of asymmetric units, and B-factor sharpening can be used
            to refine the analysis and adapt it to different biological
            datasets.

            From a biological perspective, careful masking and appropriate
            symmetry definition are important for avoiding misleading FSC
            estimations, especially in flexible regions or highly symmetric
            macromolecular assemblies.

            Outputs and Interpretation

            After execution, the protocol produces either a statistically
            evaluated FSC curve or a local resolution map, depending on the
            selected analysis mode. The FSC output can be used to determine
            the effective resolution of the reconstruction, while the local
            resolution map allows visualization of structural variability
            across different regions of the volume.

            The resulting outputs preserve the sampling information from the
            original maps and can be integrated into downstream cryo-EM
            interpretation and validation workflows.

            Final Perspective

            ProtFscFdrControl provides a statistically grounded framework
            for FSC analysis and local resolution estimation in cryo-EM
            reconstructions. By combining FDR-based significance evaluation
            with optional local analysis, the protocol helps users obtain
            more reliable resolution estimates and improve the interpretation
            of structural quality and heterogeneity in macromolecular maps.
        """
    _label = 'FSC-FDR control'
    _devStatus = BETA

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('halfWhere', BooleanParam,
                      label='Are the half maps stored with the input volume?',
                      default=False)
        form.addParam('halfOne', PointerParam, pointerClass="Volume",
                      label='First half map', important=True, condition="halfWhere==False")
        form.addParam('halfTwo', PointerParam, pointerClass="Volume",
                      label='Second half map', important=True, condition="halfWhere==False")
        form.addParam('inputVol', PointerParam, pointerClass="Volume",
                      label='Volume with half maps', important=True, condition="halfWhere==True")
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
        if self.halfWhere.get():
            file_halfOne, file_halfTwo = self.inputVol.get().getHalfMaps().split(",")
        else:
            file_halfOne = self.halfOne.get().getFileName()
            file_halfTwo = self.halfTwo.get().getFileName()
        if pwutils.getExt(file_halfOne) == '.mrc':
            file_halfOne += ':mrc'
        if pwutils.getExt(file_halfTwo) == '.mrc':
            file_halfTwo += ':mrc'
        ih.convert(file_halfOne, self._getTmpPath('halfone.mrc'))
        ih.convert(file_halfTwo, self._getTmpPath('halftwo.mrc'))

    def computeControlStep(self):
        path_half_one = os.path.abspath(self._getTmpPath('halfone.mrc'))
        path_half_two = os.path.abspath(self._getTmpPath('halftwo.mrc'))
        args = '--halfmap1 %s --halfmap2 %s' \
               ' --symmetry %s' % (path_half_one, path_half_two, self.sym.get().upper())

        if self.halfWhere.get():
            args += " --apix %f " % self.inputVol.get().getSamplingRate()
        else:
            args += " --apix %f " % self.halfOne.get().getSamplingRate()

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
        if self.localRes.get():
            _volume = Volume()
            _volume.setFileName(self._getExtraPath("halfone_localResolutions.mrc"))
            if self.halfWhere.get():
                _volume.setSamplingRate(self.inputVol.get().getSamplingRate())
                self._defineOutputs(outputLocalResMap=_volume)
                self._defineSourceRelation(self.inputVol, _volume)
            else:
                _volume.setSamplingRate(self.halfOne.get().getSamplingRate())
                self._defineOutputs(outputLocalResMap=_volume)
                self._defineSourceRelation(self.halfOne, _volume)
                self._defineSourceRelation(self.halfTwo, _volume)
        else:
            # Remove unused files
            pwutils.cleanPath(self._getExtraPath('FSC.pdf'))
            pwutils.cleanPath(self._getExtraPath('GuinierPlot.pdf'))
            pwutils.cleanPath(self._getExtraPath('postProcessed.mrc'))

            _fsc = FSC(objLabel='FSC')
            data = np.loadtxt(self._getExtraPath('FSC.txt'))
            _fsc.setData(data[0].tolist(), data[1].tolist())
            self._defineOutputs(outputFSC=_fsc)
            if self.halfWhere.get():
                self._defineSourceRelation(self.inputVol, _fsc)
            else:
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
            if self.localRes.get():
                summary.append("Local resolution computed from the half maps")
            else:
                stdout = self._getLogsPath('run.stdout')
                with open(stdout) as file:
                    for num, line in enumerate(file, 1):
                        if 'Resolution at 1 % FDR-FSC' in line:
                            res = [float(s) for s in line.split() if s.replace(".", "", 1).isdigit()][1]
                summary.append('Resolution at 1 %% FDR-FSC: %.2f Angstrom' % res)
        return summary
