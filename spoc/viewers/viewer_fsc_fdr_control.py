# **************************************************************************
# *
# * Authors:  David Herreros Calero (dherreros@cnb.csic.es)
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

from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
import pyworkflow.protocol.params as params
from pwem.wizards import ColorScaleWizardBase
from pwem.constants import COLOR_OTHER, AX_Z
from pyworkflow.protocol.params import (LabelParam, EnumParam, PointerParam,
                                        IntParam, LEVEL_ADVANCED)
from pwem.viewers import ChimeraView, FscViewer
from pwem.emlib.image import ImageHandler
import matplotlib.pyplot as plt
from matplotlib import cm
from pwem.viewers import (LocalResolutionViewer, EmPlotter, ChimeraView,
                          DataView)

from spoc.protocols.protocol_fsc_fdr_control import ProtResolutionAnalysisFSCFDR

class ViewerFscFdrControl(LocalResolutionViewer):
    """ Visualize local resolution """
    _label = 'viewer local resolution'
    _targets = [ProtResolutionAnalysisFSCFDR]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    OPEN_FILE = "open %s\n"
    VOXEL_SIZE = "volume #%d voxelSize %s\n"
    VOL_HIDE = "vol #%d hide\n"
    VIEW = "view\n"

    def _defineParams(self, form):
        form.addSection(label='Show results')
        #if hasattr(self.protocol, "outputFSC"):
        form.addParam('doShowFSC', params.LabelParam,
                          label="Display the FSC curve")

        if (self.protocol.localRes.get()):
            form.addParam('doShowLocalRes', params.LabelParam,
                          label="Display the local resolution")

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
                           label="Show Resolution map in ChimeraX")

            group.addParam('sharpenedMap', PointerParam, pointerClass='Volume',
                           label="(Optional) Color a sharpen map by local resolution in ChimeraX",
                           allowsNull=True,
                           help='Local resolution should be estimated with the raw maps instead'
                                ' of sharpen maps. Information about this in (Vilas et al '
                                'Current Opinion in Structural Biology 2021). This entry parameter '
                                'allows to color the local resolution in'
                                'a different map')

            ColorScaleWizardBase.defineColorScaleParams(group, defaultLowest=10, defaultHighest=1)

    def getImgData(self, imgFile):
        return LocalResolutionViewer.getImgData(self, imgFile)

    def merge_two_dicts(self, d1, d2):
        z = d1.copy()  # start with keys and values of x
        z.update(d2)  # modifies z with keys and values of y
        return z

    def _getVisualizeDict(self):
        d1 = {'doShowFSC': self._doShowFSC}
        d2= {}
        if self.protocol.localRes.get():
            d2 = {'doShowFSC': self._doShowFSC,
                'doShowLocalRes': self._doShowLocalRes,
                'doShowVolumeColorSlices': self._showVolumeColorSlices,
                'doShowOneColorslice': self._showOneColorslice,
                'doShowChimera': self._showChimera
                }
        return self.merge_two_dicts(d1, d2)


    def _doShowLocalRes(self, param=None):
        cm = DataView(self.protocol.outputLocalResMap.getFileName())
        return [cm]


    def _showOriginalVolumeSlices(self, param=None):
        cm = DataView(self.protocol.inputMap.get().getFileName())
        return [cm]


    def _showVolumeColorSlices(self, param=None):
        imageFile = self.protocol.outputLocalResMap.getFileName()
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
            plot = xplotter.plotMatrix(a, matrix, self.highest.get(), self.lowest.get(),
                                       cmap=self.getColorMap(),
                                       interpolation="nearest")
        xplotter.getColorBar(plot)
        return [xplotter]

    @classmethod
    def getBackGroundValue(cls, data):
        return max(data) - 1

    def _showOneColorslice(self, param=None):
        imageFile = self.protocol.outputLocalResMap.getFileName()
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
        imageFile = self.protocol.outputLocalResMap.getFileName()
        img = ImageHandler().read(imageFile)
        imgData = img.getData()
        imgList = imgData.flatten()

        imgListNoZero = imgList
        nbins = 10
        plotter = EmPlotter(x=1,y=1,mainTitle="  ")
        plotter.createSubPlot("Resolution histogram",
                              "Resolution (A)", "# of Counts")
        plotter.plotHist(imgListNoZero, nbins)
        return [plotter]

    def _getAxis(self):
        return self.getEnumText('sliceAxis')

    def getColorMap(self):
        cmap = cm.get_cmap(self.colorMap.get())
        if cmap is None:
            cmap = cm.jet
        return cmap

    def _showChimera(self, param=None):
        views = []
        errors = []

        fnResVol = self.protocol.outputLocalResMap.getFileName()

        if self.protocol.halfWhere.get():
            fnOrigMap = abspath(self.protocol.inputVol.get().getFileName())
            sampRate = self.protocol.inputVol.get().getSamplingRate()
        else:
            fnOrigMap = abspath(self.protocol.halfOne.get().getFileName())
            sampRate = self.protocol.halfOne.get().getSamplingRate()
        OPEN_FILE = "open %s\n"

        if self.sharpenedMap.get():

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

            scriptFile = self.protocol._getExtraPath('show_chimerca.cxc')
            fhCmd = open(scriptFile, 'w')

            fhCmd.write(OPEN_FILE % imageFile)
            fhCmd.write("volume #1 voxelSize %f\n" % sampRate)
            fhCmd.write("hide #1\n")
            fhCmd.write(OPEN_FILE % os.path.abspath(self.sharpenedMap.get().getFileName()))
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
        else:
            scriptFile = self.protocol._getExtraPath('show_confidence_map.py')
            self.createChimeraScript(scriptFile, fnResVol, fnOrigMap, sampRate,
                                     numColors=self.intervals.get(),
                                     lowResLimit=self.highest.get(),
                                     highResLimit=self.lowest.get())
            views.append(ChimeraView(scriptFile))

        return views


    def _doShowFSC(self, param=None):
        import subprocess
        file = self.protocol._getExtraPath('FSC.pdf')
        subprocess.call(["xdg-open", file])
        #subprocess.Popen([path], shell=True)
        fscViewer = FscViewer(project=self.protocol.getProject(),
                              threshold=0.143,
                              protocol=self.protocol,
                              figure=None,
                              addButton=True)
        fscViewer.visualize(self.protocol.outputFSC)
        return [fscViewer]

