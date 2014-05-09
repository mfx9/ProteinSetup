#===============================================================================
# Preoptimize in MM
#===============================================================================
from pBabel           import XYZFile_FromSystem
from pCore            import Unpickle, logFile, Selection
from pMolecule        import NBModelABFS
from pMoleculeScripts import ConjugateGradientMinimize_SystemGeometry

from os               import uname


TolCrit = 0.04

logFile.Header ("Running job on %s" % uname ()[1])


mol      = Unpickle ("mol.pkl")

protein  = Unpickle ("protein.pkl")

nb_model = NBModelABFS ()

protein.UpdateSystem (mol, nbModel = nb_model)


#===============================================================================
XYZFile_FromSystem ("before.xyz", mol, label = "Starting geometry from CHARMM")

ConjugateGradientMinimize_SystemGeometry (mol, logFrequency = 1, maximumIterations = 9999, rmsGradientTolerance = TolCrit)

XYZFile_FromSystem ("after.xyz", mol, label = "MM-optimized geometry with soft constraints")


logFile.Footer ()
