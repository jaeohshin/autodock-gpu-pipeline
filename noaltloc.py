from Bio.PDB import PDBParser, PDBIO, Select

class NoAltLocSelect(Select):
    def accept_atom(self, atom):
        return atom.altloc in (' ', 'A')

parser = PDBParser(QUIET=True)
structure = parser.get_structure("receptor", "receptor.pdb")
io = PDBIO()
io.set_structure(structure)
io.save("receptor_noalt.pdb", select=NoAltLocSelect())

