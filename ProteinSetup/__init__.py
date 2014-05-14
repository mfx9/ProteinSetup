#-------------------------------------------------------------------------------
# . File      : __init__.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
"""
The purpose of this module is to make easier setting up QC/MM models by defining
a helper class called ProteinSetup.

The ProteinSetup class contains lists of indices of quantum, fully flexible and
restrained atoms as well as a list with restraint values for all atoms.

It also contains the charge and multiplicity of the QC-region since these do
not change during the calculations.

The QC-charge, if not given, can be calculated automatically.

An instance of the ProteinSetup class should be pickled after the setup for further
use in QC/MM calculations.

The method UpdateSystem should be used anytime a ProteinSetup is unpickled.
"""

from ProteinSetup       import ProteinSetup, MakeScratch
