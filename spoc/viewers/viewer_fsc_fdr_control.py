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

from pwem.viewers import ChimeraView, FscViewer

from spoc.protocols import ProtFscFdrControl


class ViewerFscFdrControl(ProtocolViewer):
    """ Visualize local resolution """
    _label = 'viewer local resolution'
    _targets = [ProtFscFdrControl]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    OPEN_FILE = "open %s\n"
    VOXEL_SIZE = "volume #%d voxelSize %s\n"
    VOL_HIDE = "vol #%d hide\n"
    VIEW = "view\n"

    def _defineParams(self, form):
        form.addSection(label='Show results')
        if hasattr(self.protocol, "outputFSC"):
            form.addParam('doShowFSC', params.LabelParam,
                          label="Display the FSC curve")
        if hasattr(self.protocol, "outputLocalResMap"):
            form.addParam('doShowLocalRes', params.LabelParam,
                          label="Display the local resolution")

    def _getVisualizeDict(self):
        return {'doShowFSC': self._doShowFSC,
                'doShowLocalRes': self._doShowLocalRes
                }


    def _doShowLocalRes(self, param=None):
        scriptFile = self.protocol._getPath('localres_chimera.cxc')
        fhCmd = open(scriptFile, 'w')
        if self.protocol.halfWhere.get():
            fnVol = abspath(self.protocol.inputVol.get().getFileName())
            smprt = self.protocol.inputVol.get().getSamplingRate()
        else:
            fnVol = abspath(self.protocol._getExtraPath('halfone.mrc'))
            smprt = self.protocol.halfOne.get().getSamplingRate()
        fnResMap = abspath(self.protocol._getExtraPath("halfone_localResolutions.mrc"))

        fhCmd.write(self.OPEN_FILE % fnVol)
        fhCmd.write(self.OPEN_FILE % fnResMap)
        counter = 1
        fhCmd.write(self.VOXEL_SIZE % (counter, str(smprt)))
        counter += 1
        fhCmd.write(self.VOXEL_SIZE % (counter, str(smprt)))
        fhCmd.write(self.VOL_HIDE % counter)
        fhCmd.write('color sample #%d map #%d palette rainbow\n' % (counter - 1, counter))
        fhCmd.write(self.VIEW)
        fhCmd.close()

        view = ChimeraView(scriptFile)
        return [view]

    def _doShowFSC(self, param=None):
        fscViewer = FscViewer(project=self.protocol.getProject(),
                              threshold=0.143,
                              protocol=self.protocol,
                              figure=None,
                              addButton=True)
        fscViewer.visualize(self.protocol.outputFSC)
        return [fscViewer]

