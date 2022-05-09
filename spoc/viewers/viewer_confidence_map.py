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


import matplotlib.pyplot as plt
from matplotlib import cm
from pwem.wizards import ColorScaleWizardBase

import os
from pyworkflow.protocol.params import (LabelParam, EnumParam, PointerParam,
                                        IntParam, LEVEL_ADVANCED)
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER
from pwem.viewers import (EmPlotter, ChimeraView, DataView)
from pwem.viewers.viewer_localres import LocalResolutionViewer

from pwem.emlib.metadata import MetaData, MDL_X, MDL_COUNT
from pwem.constants import AX_Z

from pwem import splitRange
from pyworkflow.utils import replaceExt
from os.path import exists

import pyworkflow.utils as pwutils
from pwem.emlib.image import ImageHandler

import chimera

from spoc.protocols.protocol_confidence_map import (ProtConfidenceMap, OUTPUT_MAP, INPUT_MAP)
from pyworkflow.gui import plotter


class ProtConfidenceMapViewer(LocalResolutionViewer):
    """Visualization tools for confidence maps results. """

    _label = 'viewer confidence maps'
    _targets = [ProtConfidenceMap]
    _environments = [DESKTOP_TKINTER]

    @staticmethod
    def getColorMapChoices():
        return plt.colormaps()

    def __init__(self, *args, **kwargs):
        ProtocolViewer.__init__(self, *args, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Visualization')

        form.addParam('doShowConfidenceMapSlices', LabelParam,
                      label="Show confidence map slices")

        form.addParam('doShowOriginalVolumeSlices', LabelParam,
                      label="Show original volume slices")

        form.addParam('doShowResHistogram', LabelParam,
                      label="Show confidence map histogram")

        group = form.addGroup('Colored Slices and Volumes')
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
                       label="Show confidence map in Chimera")

        ColorScaleWizardBase.defineColorScaleParams(group, defaultLowest=0, defaultHighest=1)

        group.addParam('structure', PointerParam, pointerClass='AtomStruct',
                       allowsNull=True,
                       label="Structure to be colored",
                       help='Structure to be colored based on the confidence map. This parameter is '
                            'optional. If left empty, the previously provided input volume will be '
                            'colored instead.')


    def _getVisualizeDict(self):
        return {
            'doShowConfidenceMapSlices': self._showConfidenceMapSlices,
            'doShowOriginalVolumeSlices': self._showOriginalVolumeSlices,
            'doShowVolumeColorSlices': self._showVolumeColorSlices,
            'doShowOneColorslice': self._showOneColorslice,
            'doShowResHistogram': self._plotHistogram,
            'doShowChimera': self._showChimera,
        }

    def getOutputFileName(self, inputMap, outputMap):
        fnbase = pwutils.removeBaseExt(inputMap)
        return self._getExtraPath(fnbase + outputMap)

    def getImgData(self, imgFile):
        return LocalResolutionViewer.getImgData(self, imgFile)

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
        imgListNoZero = filter(lambda x: x > 0.1, imgList)
        nbins = 10
        plotter = EmPlotter(x=1, y=1, mainTitle="  ")
        plotter.createSubPlot("Confidence histogram",
                              "Significance", "# of Counts")
        plotter.plotHist(imgListNoZero, nbins)
        return [plotter]

    def _getAxis(self):
        return self.getEnumText('sliceAxis')

    def _showChimera(self, param=None):


        fnResVol = self.protocol.confidenceMap.getFileName()
        scriptFile = self.protocol._getExtraPath('show_confidence_map.cxc')
        imageFile = os.path.abspath(fnResVol)
        vol = self.protocol.inputMap.get()
        sampRate = vol.getSamplingRate()
        imgh = ImageHandler()
        x, y, z, n = imgh.getDimensions(imageFile)

        colorStr, scolorStr = self.defineColorMapandColorKey()

        OPEN_FILE = "open %s\n"
        fhCmd = open(scriptFile, 'w')
        fhCmd.write(OPEN_FILE % imageFile)
        fhCmd.write("volume #1 voxelSize %f\n" % sampRate)
        fhCmd.write("volume #1 origin %f,%f,%f\n" % (-x / 2 * sampRate, -y / 2 * sampRate, -z / 2 * sampRate))
        imageFileOrig = os.path.abspath(self.protocol.inputMap.get().getFileName())
        fhCmd.write("hide #1\n")

        views = []

        if self.structure.get():

            fhCmd.write(OPEN_FILE % os.path.abspath(self.structure.get().getFileName()))
            fhCmd.write("show cartoons\n")
            fhCmd.write("cartoon style width 1.5 thick 1.5\n")
            fhCmd.write("style stick\n")
            fhCmd.write("measure mapvalues #1 atoms #2 attribute confidence\n")
            fhCmd.write("color byattribute confidence palette %s\n" % scolorStr)
            fhCmd.write(colorStr+"\n")
            fhCmd.write("view\n")
            fhCmd.write(OPEN_FILE % imageFileOrig)
            fhCmd.write("volume #4 voxelSize %f\n" % sampRate)
            fhCmd.write("volume #4 origin %f,%f,%f\n" % (-x / 2 * sampRate, -y / 2 * sampRate, -z / 2 * sampRate))
            fhCmd.write("hide #4\n")
        else:
            fhCmd.write(OPEN_FILE % imageFileOrig)
            fhCmd.write("volume #2 voxelSize %f\n" % sampRate)
            fhCmd.write("volume #2 origin %f,%f,%f\n" % (-x / 2 * sampRate, -y / 2 * sampRate, -z / 2 * sampRate))
            fhCmd.write("color sample #2 map #1 palette %s\n" % scolorStr)
            fhCmd.write(colorStr + "\n")
            fhCmd.write("view\n")

        fhCmd.close()

        views.append(ChimeraView(scriptFile))

        return views

    def defineColorMapandColorKey(self):
        stepColors = splitRange(0, 1,
                                splitNum=13)
        colorList = plotter.getHexColorList(len(stepColors), self._getColorName())
        colorStr = 'key{} fontSize 15 size 0.025,0.4 pos 0.01,0.3\n'
        keyStr = ''
        scolorStr = ''
        stepColors.reverse(), colorList.reverse()
        for step, color in zip(stepColors, colorList):
            keyStr += ' %s:%s' % (color, step)
            scolorStr += '%s,%s:' % (step, color)
        colorStr = colorStr.format(keyStr)
        scolorStr = scolorStr[:-1]
        return colorStr, scolorStr

    def getColorMap(self):
        cmap = cm.get_cmap(self.colorMap.get())
        if cmap is None:
            cmap = cm.jet
        return cmap
