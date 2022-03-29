# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     J.L. Vilas (jlvilas@cnb.csic.es)
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
import os

import matplotlib.pyplot as plt
from matplotlib import cm
from pwem import splitRange
from pwem.constants import COLOR_OTHER, AX_Z
from pwem.wizards import ColorScaleWizardBase
from pyworkflow.utils import replaceExt
from os.path import exists

from pyworkflow.protocol.params import (LabelParam, EnumParam, PointerParam,
                                        IntParam, LEVEL_ADVANCED)
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER

from pwem.viewers import (LocalResolutionViewer, EmPlotter, ChimeraView,
                          DataView)
from pwem.emlib.metadata import MetaData, MDL_X, MDL_COUNT

import pyworkflow.utils as pwutils
from pwem.emlib.image import ImageHandler

import chimera

from spoc.protocols.protocol_confidence_map import (ProtConfidenceMap, OUTPUT_MAP, INPUT_MAP)
from pyworkflow.gui import plotter


class XmippConfidenceMapViewer(LocalResolutionViewer):
    """Visualization tools for confidence maps results. """

    _environments = [DESKTOP_TKINTER]
    _targets = [ProtConfidenceMap]
    _label = 'viewer'

    @staticmethod
    def getColorMapChoices():
        return plt.colormaps()

    def __init__(self, *args, **kwargs):
        ProtocolViewer.__init__(self, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Visualization')

        form.addParam('doShowConfidenceMapSlices', LabelParam,
                      label="Show confidence map slices")

        form.addParam('doShowOriginalVolumeSlices', LabelParam,
                      label="Show original volume slices")

        form.addParam('doShowResHistogram', LabelParam,
                      label="Show confidence map histogram")

        group = form.addGroup('Colored resolution Slices and Volumes')
        group.addParam('sliceAxis', EnumParam, default=AX_Z,
                       choices=['x', 'y', 'z'],
                       display=EnumParam.DISPLAY_HLIST,
                       label='Slice axis')

        group.addParam('doShowVolumeColorSlices', LabelParam,
                       label="Show colored slices")

        group.addParam('doShowOneColorslice', LabelParam,
                       expertLevel=LEVEL_ADVANCED,
                       label='Show selected slice')
        group.addParam('sliceNumber', IntParam, default=-1,
                       expertLevel=LEVEL_ADVANCED,
                       label='Show slice number')

        group.addParam('doShowChimera', LabelParam,
                       label="Show Resolution map in Chimera")
        group.addParam('structure', PointerParam, pointerClass='AtomStruct',
                       allowsNull=True,
                       label="Structure to be colored",
                       help='Structure to be colored based on the confidence map. This parameter is '
                            'optional. If left empty, the previously provided input volume will be '
                            'colored instead.')

        ColorScaleWizardBase.defineColorScaleParams(group, defaultLowest=0, defaultHighest=1)

    def getOutputFileName(self, inputMap, outputMap):
        fnbase = pwutils.removeBaseExt(inputMap)
        return self._getExtraPath(fnbase + outputMap)

    def getImgData(self, imgFile):
        return LocalResolutionViewer.getImgData(self, imgFile)

    def _getVisualizeDict(self):
        return {
            'doShowConfidenceMapSlices': self._showConfidenceMapSlices,
            'doShowOriginalVolumeSlices': self._showOriginalVolumeSlices,
            'doShowVolumeColorSlices': self._showVolumeColorSlices,
            'doShowOneColorslice': self._showOneColorslice,
            'doShowResHistogram': self._plotHistogram,
            'doShowChimera': self._showChimera
        }

    def _showConfidenceMapSlices(self, param=None):
        cm = DataView(self.protocol.confidenceMap.getFileName())
        cm2 = DataView(self.protocol.confidenceMap_log10FDR.getFileName())
        return [cm, cm2]

    def _showOriginalVolumeSlices(self, param=None):
        cm = DataView(self.protocol.inputMap.get().getFileName())
        return [cm]

    def _showVolumeColorSlices(self, param=None):
        imageFile = self.protocol.confidenceMap.getFileName()
        imgData, _, _, _ = self.getImgData(imageFile)

        xplotter = EmPlotter(x=2, y=2, mainTitle="Confidence Map Slices "
                                                 "along %s-axis."
                                                 % self._getAxis())
        # The slices to be shown are close to the center. Volume size is divided
        # in segments, the fourth central ones are selected i.e. 3,4,5,6
        for i in list(range(3, 7)):
            sliceNumber = self.getSlice(i, imgData)
            a = xplotter.createSubPlot("Slice %s" % (sliceNumber + 1), '', '')
            matrix = self.getSliceImage(imgData, sliceNumber, self._getAxis())
            plot = xplotter.plotMatrix(a, matrix, self.lowest.get(), self.highest.get(),
                                       cmap=self.getColorMap(),
                                       interpolation="nearest")
        xplotter.getColorBar(plot)
        return [xplotter]

    @classmethod
    def getBackGroundValue(cls, data):
        return max(data) - 1

    def _showOneColorslice(self, param=None):
        imageFile = self.protocol.confidenceMap.getFileName()
        imgData, _, _, volDims = self.getImgData(imageFile)
        xplotter = EmPlotter(x=1, y=1, mainTitle="Confidence map Slices "
                                                 "along %s-axis."
                                                 % self._getAxis())
        sliceNumber = self.sliceNumber.get()
        if sliceNumber < 0:
            sliceNumber = volDims[0] / 2
        else:
            sliceNumber -= 1
        # sliceNumber has no sense to start in zero
        a = xplotter.createSubPlot("Slice %s" % (sliceNumber + 1), '', '')
        matrix = self.getSliceImage(imgData, sliceNumber, self._getAxis())
        plot = xplotter.plotMatrix(a, matrix, self.lowest.get(), self.highest.get(),
                                   cmap=self.getColorMap(),
                                   interpolation="nearest")
        xplotter.getColorBar(plot)
        return [xplotter]

    def _plotHistogram(self, param=None):
        imageFile = self.protocol.confidenceMap.getFileName()
        img = ImageHandler().read(imageFile)
        imgData = img.getData()
        imgList = imgData.flatten()
        # imgDataMax = self.getBackGroundValue(imgList)
        imgListNoZero = imgList  # filter(lambda x: 0 < x < imgDataMax, imgList)
        nbins = 10
        plotter = EmPlotter(x=1, y=1, mainTitle="  ")
        plotter.createSubPlot("Resolution histogram",
                              "Resolution (A)", "# of Counts")
        plotter.plotHist(imgListNoZero, nbins)
        return [plotter]

    def _getAxis(self):
        return self.getEnumText('sliceAxis')

    def _showChimera(self, param=None):
        views = []
        errors = []

        fnResVol = self.protocol.confidenceMap.getFileName()
        vol = self.protocol.inputMap.get()

        fnOrigMap = vol.getFileName()
        sampRate = vol.getSamplingRate()
        OPEN_FILE = "open %s\n"

        if self.structure.get() and '1.3' in chimera.Plugin.getHome():

            imageFile = os.path.abspath(fnResVol)
            _, minRes, maxRes, voldim = self.getImgData(imageFile)
            # Narrow the color range to the highest resolution range
            lowResLimit = min(maxRes, minRes + 5)
            highResLimit = minRes
            stepColors = splitRange(highResLimit, lowResLimit,
                                    splitNum=13)
            colorList = plotter.getHexColorList(len(stepColors), self._getColorName())
            colorStr = 'key{} fontSize 15 size 0.025,0.4 pos 0.01,0.3\n'
            keyStr = ''
            stepColors.reverse(), colorList.reverse()
            for step, color in zip(stepColors, colorList):
                keyStr += ' {}:{}'.format(color, step)
            colorStr = colorStr.format(keyStr)

            scriptFile = self.protocol._getExtraPath('show_confidence_map.cxc')
            fhCmd = open(scriptFile, 'w')

            fhCmd.write(OPEN_FILE % imageFile)
            fhCmd.write("volume #1 voxelSize %f\n" % sampRate)
            fhCmd.write("hide #1\n")
            fhCmd.write(OPEN_FILE % os.path.abspath(self.structure.get().getFileName()))
            # fhCmd.write("start Model Panel\n")
            fhCmd.write("show cartoons\n")
            fhCmd.write("cartoon style width 1.5 thick 1.5\n")
            fhCmd.write("style stick\n")
            fhCmd.write("measure mapvalues #1 atoms #2 attribute confidence\n")
            fhCmd.write("color byattribute confidence\n")
            fhCmd.write(colorStr)
            fhCmd.write("view\n")
            fhCmd.close()
            views.append(ChimeraView(scriptFile))
        elif self.structure.get() and '1.3' not in chimera.Plugin.getHome():
            errors.append('Atomic structure colouring by map values is only supported by ChimeraX-1.3. '
                          'Please, update the ChimeraX binaries if you would like to colour an atomic '
                          'structure with the confidence map values.')
            self.errorList(errors, views)
        else:
            scriptFile = self.protocol._getExtraPath('show_confidence_map.py')
            self.createChimeraScript(scriptFile, fnResVol, fnOrigMap, sampRate,
                                     numColors=self.intervals.get(),
                                     lowResLimit=self.highest.get(),
                                     highResLimit=self.lowest.get())
            views.append(ChimeraView(scriptFile))

        return views

    def getColorMap(self):
        cmap = cm.get_cmap(self.colorMap.get())
        if cmap is None:
            cmap = cm.jet
        return cmap
