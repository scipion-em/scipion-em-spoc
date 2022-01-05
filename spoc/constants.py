# **************************************************************************
# *
# *  Authors:     David Herreros Calero (dherreros@cnb.csic.es)
# *
# * National Center for Biotechnology (CSIC)
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
import spoc

SPOC_HOME = 'SPOC_HOME'

# Supported stable versions
V1_0 = '1.0'

# Supported continuous build versions
# !IMPORTANT: Whenever we update the continuous version, we need to update the commit to the one
# we have been testing and the version. Version format is the following:
#  'Previous_stable_spoc_version' + '_' + 'date_of_release_in_format_year+month+day'
# Example:
#  V_CB = '1.0_220105' (Stable version 1.0 and commit downloaded on the 05/01/2022)
COMMIT = '55b4f82'
V_CB = '1.0_220105'
