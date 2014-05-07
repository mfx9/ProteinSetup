#===============================================================================
# Setup of the QC/MM system
#===============================================================================
from pBabel       import CHARMMParameterFiles_ToParameters, CHARMMPSFFile_ToSystem, XYZFile_ToCoordinates3
from pCore        import Pickle

from ProteinSetup import ProteinSetup


chromophore = ("C1", "N2", "N3", "C2", "O2", "CA", "CB", "HB", "CG", "CD1", "HD1", "CD2", "HD2", "CE1", "HE1", "CE2", "HE2", "CZ", "OH", "CAG", "HA1", "HA2", "N", "H", "CAH", "HAH", "HH")

quantumRegion   = (
  ("PRTA", "CR8",   64, chromophore),
  ("OXYA", "O2",     1, ("O1", "O2")),
  ("PRTA", "MET",  159),
  ("PRTA", "SER",  142),
  ("WATA", "TP3M", 166),
  ("WATA", "TP3M",  17),
)

#===============================================================================
par = CHARMMParameterFiles_ToParameters (["initialize/toppar/par_all27_prot_na.inp", "initialize/toppar/par_cro_irisfp_anionic.inp"])

mol = CHARMMPSFFile_ToSystem ("initialize/monomer_waterbox.psf", isXPLOR = True, parameters = par)

mol.coordinates3 = XYZFile_ToCoordinates3 ("initialize/final_frame.xyz")

Pickle ("mol.pkl", mol)


protein = ProteinSetup ()

protein.Initialize (mol, quantumRegion)

protein.Summary ()


protein.WriteQuantum (mol)

protein.WritePDB (mol)

Pickle ("protein.pkl", protein)
