#===============================================================================
# Setup of the QC/MM system
#===============================================================================
from pBabel       import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, CHARMMCRDFile_ToCoordinates3
from pCore        import Pickle

from ProteinSetup import ProteinSetup


quantumRegion   = (
  ("GOLA", "GOL", 801, ("O1", "H1", "C1", "H11", "H12", "C2", "O2", "H2", "H21", "C3", "H31", "H32", "O3", "H3")),
  ("PRTA", "ASN", 156, ("HD21", "HD22", "ND2", "CG", "OD1", "CB", "HB1", "HB2")),
  ("PRTA", "CYS", 433, ("N", "H", "CA", "HA", "CB", "HB1", "HB2", "SG", "C", "O")),
  ("PRTA", "VAL", 434, ("N", "H", "CA", "HA")),
  ("PRTA", "ASP", 447, ("OD2", "HE2", "OD1", "CG", "CB", "HB1", "HB2")),
  ("PRTA", "GLU", 435, ("HE2", "OE2", "OE1", "CD", "CG", "HG1", "HG2")),
  ("PRTA", "HIS", 164, ),
  ("PRTA", "VAL", 162, ("O", "C")),
  ("PRTA", "GLY", 163, ("N", "H", "CA", "HA1", "HA2")),
  ("PRTA", "GLN", 549, ("HE21", "HE22", "NE2", "CD", "OE1")),
  ("PRTA", "GLY", 432, ("C", "O")),
  ("PRTA", "SER", 282, ),
  ("PRTA", "HSE", 281, ),
  ("WATA", "HOH", 71,  ("O", "H1", "H2")),
)

#===============================================================================
par = CHARMMParameterFiles_ToParameters (["charmm/toppar/par_all27_prot_na.inp"])

mol = CHARMMPSFFile_ToSystem ("charmm/monomer_xplor.psf", isXPLOR = True, parameters = par)

mol.coordinates3 = CHARMMCRDFile_ToCoordinates3 ("charmm/monomer.crd")

Pickle ("mol.pkl", mol)


protein = ProteinSetup (multiplicity = 2)

protein.Initialize (mol, quantumRegion)

protein.Summary ()


protein.WriteQuantum (mol)

protein.WritePDB (mol)

Pickle ("protein.pkl", protein)
