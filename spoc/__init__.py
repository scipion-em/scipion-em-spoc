# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es) [1]
# *
# * [1] National Center for Biotechnology (CSIC)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

import pwem
import pyworkflow.utils as pwutils

import spoc.constants as spocConst

__version__ = "1.0.1"
_logo = "spoc_logo.png"
_references = ['Beckers2019']
_url = "https://github.com/scipion-em/scipion-em-spoc"


SCRATCHDIR = pwutils.getEnvVariable('SPOCSCRATCHDIR', default='/tmp/')


class Plugin(pwem.Plugin):
    _homeVar = spocConst.SPOC_HOME
    _pathVars = [spocConst.SPOC_HOME]
    _supportedVersions = [spocConst.V_CB]

    @classmethod
    def _defineVariables(cls, version=spocConst.V_CB):
        cls._defineEmVar(spocConst.SPOC_HOME, 'spoc-' + version)

    @classmethod
    def getEnviron(cls):
        pass

    @classmethod
    def getActiveVersion(cls, home=None, versions=None):
        home = os.path.basename(home or cls.getHome())
        versions = versions or cls.getSupportedVersions()
        currentVersion = home.split('-')[-1]

        for v in versions:
            if v == currentVersion:
                return v

        return ''

    @classmethod
    def isVersion(cls, version=spocConst.V_CB):
        return cls.getActiveVersion() == version

    @classmethod
    def getProgram(cls, program):
        """ Return the program binary that will be used. """
        program = os.path.join(cls.getHome('spoc-source'), program)
        return 'python %(program)s ' % locals()

    @classmethod
    def getSpocCommand(cls, program, args):
        return cls.getProgram(program) + args

    @classmethod
    def defineBinaries(cls, env):
        # For Spoc-CB
        spocCB_commands = []
        spocCB_commands.append(('wget -c https://github.com/MaximilianBeckers/SPOC/archive/%s.tar.gz' % spocConst.COMMIT,
                                "%s.tar.gz" % spocConst.COMMIT))
        spocCB_commands.append(("tar -xvf %s.tar.gz" % spocConst.COMMIT, []))
        spocCB_commands.append(("mv SPOC*/ spoc-source", "spoc-source"))

        env.addPackage('spoc', version=spocConst.V_CB,
                       commands=spocCB_commands,
                       tar='void.tgz',
                       default=True)
