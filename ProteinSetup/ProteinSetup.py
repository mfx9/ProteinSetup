#-------------------------------------------------------------------------------
# . File      : ProteinSetup.py
# . Program   : pDynamo-1.8.0                           (http://www.pdynamo.org)
# . Copyright : CEA, CNRS, Martin  J. Field  (2007-2012),
#                          Mikolaj J. Feliks (2014)
# . License   : CeCILL French Free Software License     (http://www.cecill.info)
#-------------------------------------------------------------------------------
import exceptions, os, getpass, time

from pBabel            import PDBFile_FromSystem, XYZFile_FromSystem

from pCore             import logFile, LogFileActive, Selection, Clone

from pMolecule         import System, SoftConstraintContainer, SoftConstraintEnergyModelHarmonic, SoftConstraintTether, ElectronicState, QCModelORCA, NBModelORCA

from pMoleculeScripts  import PruneByAtom



_DefaultAtomLibraryCHARMM = {
"ALA" : ("CA", "HA", "CB", "HB1", "HB2", "HB3"),
\
"ARG" : ("CG", "HG1", "HG2", "CD", "HD1", "HD2", "NE", "HE", "CZ", "NH1", "HH11", "HH12", "NH2", "HH21", "HH22"),
\
"ASN" : ("CB", "HB1", "HB2", "CG", "OD1", "ND2", "HD21", "HD22"),
\
"ASP" : ("CB", "HB1", "HB2", "CG", "OD1", "OD2"),
\
"CYS" : ("CB", "HB1", "HB2", "SG", "HG1"),
\
"GLN" : ("CB", "HB1", "HB2", "CG", "HG1", "HG2", "CD", "OE1", "NE2", "HE21", "HE22"),
\
"GLU" : ("CB", "HB1", "HB2", "CG", "HG1", "HG2", "CD", "OE1", "OE2"),
\
"HIS" : ("CB", "HB1", "HB2", "CD2", "HD2", "CG", "NE2", "HE2", "ND1", "HD1", "CE1", "HE1"),
\
"ILE" : ("CB", "HB", "CG2", "HG21", "HG22", "HG23", "CG1", "HG11", "HG12", "CD1", "HD11", "HD12", "HD13"),
\
"LEU" : ("CB", "HB1", "HB2", "CG", "HG", "CD1", "HD11", "HD12", "HD13", "CD2", "HD21", "HD22", "HD23"),
\
"LYS" : ("CB", "HB1", "HB2", "CG", "HG1", "HG2", "CD", "HD1", "HD2", "CE", "HE1", "HE2", "NZ", "HZ1", "HZ2", "HZ3"),
\
"MET" : ("CB", "HB1", "HB2", "CG", "HG1", "HG2", "SD", "CE", "HE1", "HE2", "HE3"),
\
"PHE" : ("CB", "HB1", "HB2", "CG", "CD1", "HD1", "CE1", "HE1", "CZ", "HZ", "CD2", "HD2", "CE2", "HE2"),
\
"SER" : ("CB", "HB1", "HB2", "OG", "HG1"),
\
"THR" : ("CB", "HB", "OG1", "HG1", "CG2", "HG21", "HG22", "HG23"),
\
"TRP" : ("CB", "HB1", "HB2", "CG", "CD1", "HD1", "NE1", "HE1", "CE2", "CD2", "CE3", "HE3", "CZ3", "HZ3", "CZ2", "HZ2", "CH2", "HH2"),
\
"TYR" : ("CB", "HB1", "HB2", "CG", "CD1", "HD1", "CE1", "HE1", "CZ", "OH", "HH", "CD2", "HD2", "CE2", "HE2"),
\
"VAL" : ("CB", "HB", "CG1", "HG11", "HG12", "HG13", "CG2", "HG21", "HG22", "HG23"),
\
"TP3M" : ("OH2", "H1", "H2"),
}

_DefaultAtomLibraryCHARMM ["HSE"] = _DefaultAtomLibraryCHARMM ["HIS"]
_DefaultAtomLibraryCHARMM ["HSD"] = _DefaultAtomLibraryCHARMM ["HIS"]

_DefaultAtomLibraryCHARMM ["WAT"] = _DefaultAtomLibraryCHARMM ["TP3M"]
_DefaultAtomLibraryCHARMM ["HOH"] = _DefaultAtomLibraryCHARMM ["TP3M"]
_DefaultAtomLibraryCHARMM ["TIP"] = _DefaultAtomLibraryCHARMM ["TP3M"]


#-------------------------------------------------------------------------------
class ProteinSetupError (exceptions.StandardError):
  """A class for handling errors in ProteinSetup."""
  pass


#-------------------------------------------------------------------------------
class ProteinSetup (object):
  """A class for storing the system setup, i.e. indices of restrained atoms, reference coordinates, QC-region, etc."""

  defaultAttributes = {
                  "isInitialized"      : False  ,
                  "indicesQuantum"     : None   ,
                  "indicesFlexible"    : None   ,
                  "indicesRestrained"  : None   ,
                  "reference"          : None   ,
                  "restraints"         : None   ,
                  "forceConstant"      : 50.0   ,
                  "distanceStart"      :  8.0   ,
                  "distanceEnd"        : 16.0   ,
                  "charge"             : None   ,
                  "multiplicity"       :    1   ,
                      }

  def __init__ (self, *arguments, **keywordArguments):
    """Constructor."""
    for (key, value) in self.__class__.defaultAttributes.iteritems (): setattr (self, key, value)
    for (key, value) in                 keywordArguments.iteritems (): setattr (self, key, value)


  def Initialize (self, system, quantumResidues, log = logFile):
    """Calculates restraints and total charge of the QC-region."""
    if not self.isInitialized:

      # Determine indices of QC-atoms
      self.indicesQuantum = []

      for quantumResidue in quantumResidues:
        segName, resName, resNum = quantumResidue[:3]

        if len (quantumResidue) > 3:
          atomNames = quantumResidue[3]
        else:
          if not _DefaultAtomLibraryCHARMM.has_key (resName):
            raise ProteinSetupError ("The atoms of the residue %s must be specified explicitly." % resName)
          atomNames = _DefaultAtomLibraryCHARMM[resName]

        segment = system.sequence.childIndex[segName]
        residue = segment.childIndex["%s.%s" % (resName, resNum)]

        indices = []
        for atomName in atomNames:
          if residue.childIndex.has_key (atomName):
            indices.append (residue.childIndex[atomName].index)
          else:
            if LogFileActive (log):
              log.Text ("Warning: atom name %s not found in residue %s" % (atomName, resName))

        if indices:
          self.indicesQuantum.extend (indices)

      self.indicesQuantum.sort ()


      # Determine restrained atoms and calculate restraints
      self.restraints        = []
      self.indicesFlexible   = []
      self.indicesRestrained = []

      natoms = system.sequence.NumberOfAtoms ()

      for atomIndex in range (0, natoms):
        restr = 0.0

        if not atomIndex in self.indicesQuantum:
          distances = []

          for atomIndexQuantum in self.indicesQuantum:
            distance = system.coordinates3.Distance (atomIndex, atomIndexQuantum)
            distances.append (distance)

          shortestDist = min (distances)

          if shortestDist > self.distanceEnd:
            restr = self.forceConstant
          else:
            if shortestDist > self.distanceStart:
              restr = self.forceConstant / (self.distanceEnd - self.distanceStart) * (shortestDist - self.distanceStart)

        if restr > 0.0:
          self.indicesRestrained.append (atomIndex)
        else:
          self.indicesFlexible.append (atomIndex)
        self.restraints.append (restr)

      # self.reference = system.coordinates3[:]
      self.reference = Clone (system.coordinates3)


      # Calculate charge of the QC-region automatically if not specified
      if self.charge is None:
        charges     = system.AtomicCharges ()
        totalCharge = sum (map (lambda index: charges[index], self.indicesQuantum))

        roundCharge = round (totalCharge)
        if abs (roundCharge - totalCharge) > 0.001:
          raise ProteinSetupError ("The calculated total charge of the QC-region is not integral.")

        self.charge = int (roundCharge)

      self.isInitialized = True


  def Summary (self, log = logFile):
    """Summary."""
    if self.isInitialized:
      if LogFileActive (log):
        summary = log.GetSummary ()
        summary.Start ("Protein Setup Summary")

        nquantum    = len (self.indicesQuantum)
        nflexible   = len (self.indicesFlexible)
        nrestrained = len (self.indicesRestrained)

        summary.Entry ("QC-charge",              "%d" % self.charge)
        summary.Entry ("QC-multiplicity",        "%d" % self.multiplicity)
        summary.Entry ("QC-atoms",               "%d" % nquantum)
        summary.Entry ("Flexible atoms (QC+MM)", "%d" % nflexible)
        summary.Entry ("Flexible atoms (MM)",    "%d" % (nflexible - nquantum))
        summary.Entry ("Restrained atoms",       "%d" % nrestrained)
        summary.Entry ("Buffer start",           "%f" % self.distanceStart)
        summary.Entry ("Buffer stop",            "%f" % self.distanceEnd)
        summary.Entry ("Maximum force constant", "%f" % self.forceConstant)

        summary.Stop ()


  def _WriteList (self, filename, toWrite):
    """Write a list to a file.
    """
    output = open (filename, "w")
    output.write (" ".join (map (str, toWrite)))
    output.close ()


  def WriteQuantum (self, system, filestem = "qc"):
    """Write coordinates of QC-atoms.
    """
    if self.isInitialized:
      selection = PruneByAtom (system, Selection (self.indicesQuantum))

      XYZFile_FromSystem ("%s.xyz" % filestem, selection, label = "QC-region")

      self._WriteList (filestem, self.indicesQuantum)


  def WriteFlexible (self, system, filestem = "flex"):
    """Write coordinates of fully flexible atoms.
    """
    if self.isInitialized:
      selection = PruneByAtom (system, Selection (self.indicesFlexible))

      XYZFile_FromSystem ("%s.xyz" % filestem, selection, label = "Fully flexible atoms")

      self._WriteList (filestem, self.indicesFlexible)


  def WriteRestrained (self, system, filestem = "restr"):
    """Write coordinates of restrained atoms.
    """
    if self.isInitialized:
      selection = PruneByAtom (system, Selection (self.indicesRestrained))

      XYZFile_FromSystem ("%s.xyz" % filestem, selection, label = "Restrained atoms")

      self._WriteList (filestem, self.indicesRestrained)


  def WritePDB (self, system, filename = "see_restr.pdb"):
    """Write a PDB file containing restraints in the B-factor column.

    The restraints are scaled so that they run from 0 to 1."""
    if self.isInitialized:
      scale = 1.0 / self.forceConstant
      PDBFile_FromSystem (filename, system, occupancies = map (lambda r: r * scale, self.restraints))


  def GenerateConstraintContainer (self):
    """Return container of soft constraints."""
    if self.isInitialized:
      container = SoftConstraintContainer ()

      for atomIndex in self.indicesRestrained:
        restraint = self.restraints[atomIndex]
        model     = SoftConstraintEnergyModelHarmonic (0.0, restraint)
        container["restr%s" % atomIndex] = SoftConstraintTether (atomIndex, self.reference.GetRow (atomIndex), model)
      return container


  def UpdateSystem (self, system, qcModel = None, nbModel = None, container = None, log = logFile):
    """Apply restraints, electronic state and QC/NB models to the system."""
    if self.isInitialized:
      if container is None:
        container = self.GenerateConstraintContainer ()

      system.DefineSoftConstraints (container)

      # Electronic state, QC- and NB-models
      system.electronicState = ElectronicState (charge = self.charge, multiplicity = self.multiplicity)

      # if isinstance (qcModel, QCModelORCA):
      if qcModel is not None:
        system.DefineQCModel (qcModel, qcSelection = Selection (self.indicesQuantum))

      # if isinstance (nbModel, NBModelORCA):
      if nbModel is not None:
        system.DefineNBModel (nbModel)


#===============================================================================
# Helper functions
#===============================================================================
def MakeScratch (scratchRoot = "/tmp", log = logFile):
  user    = getpass.getuser ()
  pid     = os.getpid ()
  date    = time.strftime ("%y%m%d%H%M%S")
  scratch = "%s/pdynamo_%s_%s_%s" % (scratchRoot, user, date, pid)

  os.makedirs (scratch)
#  try:
#    os.makedirs (scratch)
#  except:
#    raise ProteinSetupError ("Cannot create scratch directory %s" % scratch)

  if LogFileActive (log):
    logFile.Text ("\nScratch directory: %s\n" % scratch)

  return scratch


#===============================================================================
if __name__ == "__main__" : pass
