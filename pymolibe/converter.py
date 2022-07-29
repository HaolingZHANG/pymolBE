from os import environ
from enum import Enum
from warnings import filterwarnings
from Bio.PDB.Atom import Atom


# Type to save a temp PDB file.
# 'average' is to save the average individual calculated from all the individual.
# 'individual' is to save files of each individual.
class TempFileType(Enum):
    average = 0
    individual = 1


def set_pdb_file(segment, structures, cluster_index=None,
                 temp_file_type=TempFileType.average, save_path=None):
    """
    Set the token file.

    :param segment: token segment.
    :type segment: str

    :param structures: structure groups.
    :type structures: list

    :param cluster_index: cluster index of a token.
    :type cluster_index: int

    :param temp_file_type: token type, average of cluster or all in cluster.
    :type temp_file_type: converter.TempFileType

    :param save_path: path to save protein file.
    :type save_path: str or None

    :return: file path(s) of protein.
    :rtype: str or list
    """
    protein_letters = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
                       "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
                       "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
                       "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"}

    filterwarnings("ignore")
    if save_path is None:
        save_path = environ["TMP"].replace("\\", "/") + "/"
    form = ["N", "CA", "C"]
    if temp_file_type == TempFileType.average:
        average_position_data = structures[0]
        for position_data in structures[1:]:
            average_position_data += position_data
        average_position_data = average_position_data / len(structures)
        file_path = save_path + str(len(segment)) + "_" + segment + "_A_"
        file_path += ("C_" + str(cluster_index + 1) if cluster_index else "") + ".pdb"
        with open(file_path, "w") as pdb_file:
            for position_index, position in enumerate(average_position_data):
                atom = Atom(name=form[position_index % len(form)], serial_number=position_index + 1,
                            fullname=" " + form[position_index % len(form)].ljust(3), altloc=" ",
                            bfactor=0.0, coord=position.tolist(), occupancy=1.0)
                residue = protein_letters[segment[position_index // len(form)]].upper()
                line = write_atom_line(atom=atom, atom_number=position_index + 1, residue=residue, chain_index="A")
                pdb_file.write(line)
            pdb_file.write("TER\nEND")
        return file_path
    elif temp_file_type == TempFileType.individual:
        file_paths = []
        for file_index, position_data in enumerate(structures):
            file_path = save_path + str(len(segment)) + "_" + segment + "_I_"
            file_path += ("C_" + str(cluster_index + 1) + "_") if cluster_index else ""
            file_path += str(file_index + 1) + ".pdb"
            with open(file_path, "w") as pdb_file:
                for position_index, position in enumerate(position_data):
                    atom = Atom(name=form[position_index % len(form)], serial_number=position_index + 1,
                                fullname=" " + form[position_index % len(form)].ljust(3), altloc=" ",
                                bfactor=0.0, coord=position.tolist(), occupancy=1.0)
                    residue = protein_letters[segment[position_index // len(form)]].upper()
                    line = write_atom_line(atom=atom, atom_number=position_index + 1, residue=residue, chain_index="A")
                    pdb_file.write(line)
                pdb_file.write("TER\nEND")
            file_paths.append(file_path)
        return file_paths
    else:
        raise ValueError("No such type of token file!")


def write_atom_line(atom, atom_number, residue, chain_index):
    """
    Write the atom line as PDB file.

    :param atom: atom data.
    :type atom: Bio.PDB.Atom.Atom

    :param atom_number: number of atom.
    :type atom_number: int

    :param residue: residue name.
    :type residue: str

    :param chain_index: chain index in the protein.
    :type chain_index: str

    :return: the atom line.
    :rtype str
    """
    record_type = "ATOM  "
    element = atom.element.strip().upper().rjust(2)
    x, y, z = atom.get_coord()
    line = "%s%5i %-4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%s%6.2f      %4s%2s%2s\n" % \
           (record_type, atom_number, atom.get_fullname(), " ", residue, chain_index, (atom_number - 1) // 3, " ",
            x, y, z, "  1.00", 0.0, " ", element, "  ")
    return line
