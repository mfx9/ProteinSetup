#===============================================================================
# Optimize in QC/MM
#===============================================================================
from pBabel           import XYZFile_FromSystem, XYZFile_ToCoordinates3
from pCore            import Unpickle, logFile
from pMolecule        import NBModelORCA, QCModelORCA
from pMoleculeScripts import ConjugateGradientMinimize_SystemGeometry

from os               import uname

from ProteinSetup     import MakeScratch


TolCrit = 0.04

logFile.Header ("Running job on %s" % uname ()[1])


mol      = Unpickle ("../mol.pkl")

protein  = Unpickle ("../protein.pkl")

protein.Summary ()


#===============================================================================
scratch  = MakeScratch ()

qc_model = QCModelORCA ("B3LYP:6-31G*", "SCFCONV10", "PAL4", command = "/home/mikolaj/local/opt/orca_2_9_1_linux_x86-64/orca", scratch = scratch, job = "job", deleteJobFiles = False)

nb_model = NBModelORCA ()

protein.UpdateSystem (mol, qcModel = qc_model, nbModel = nb_model)

mol.coordinates3 = XYZFile_ToCoordinates3 ("../after.xyz")


#===============================================================================
ConjugateGradientMinimize_SystemGeometry (mol, logFrequency = 1, maximumIterations = 9999, rmsGradientTolerance = TolCrit)

XYZFile_FromSystem ("after_qcmm.xyz", mol, label = "QC/MM-optimized geometry with soft constraints")


logFile.Footer ()
