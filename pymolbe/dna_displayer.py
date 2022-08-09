from io import BytesIO
import enum
import matplotlib
from numpy import array
import pymol2
from PIL import Image
from matplotlib import pyplot, animation
from os import listdir
matplotlib.use('Agg')


# color to display each residue.
RESIDUE_COLOR_DICT = {"DA": "0xff0000", "DT": "0x0000ff", "DC": "0x00ff00", "DG": "0xffff00", "DI": "0xffff00"}


# type to display a DNA or a token.
class DisplayType(enum.Enum):
    picture = 0
    movie = 1


def draw_dna_structure(cif_path, cif_name, save_path, display_type):
    """
    Draw a dna structure.
    |cite| Warren L. DeLano (2002) CCP4 Newsletter on protein crystallography

    :param cif_path: path to load CIF file.
    :type cif_path: str

    :param cif_name: display name for saving.
    :type cif_name: str

    :param save_path: parent path to save file.
    :type save_path: str

    :param display_type: type to display.
    :type display_type: DNAStructure.displayer.DisplayType
    """
    with pymol2.PyMOL() as mol:
        mol.cmd.load(cif_path, quiet=1)
        set_initial_state(mol=mol)
        set_residue_colors(mol=mol)
        save_as_required(mol=mol, display_type=display_type, save_path=save_path, file_name=cif_name)


def set_initial_state(mol):
    """
    Set the initial state of dna.
    |cite| Warren L. DeLano (2002) CCP4 Newsletter on protein crystallography

    :param mol: PyMOL interface.
    :type mol: pymol2.PyMOL
    """
    mol.cmd.hide(representation="cartoon")
    mol.cmd.show(representation="sticks")
    mol.cmd.hide(selection="(r. HOH)")
    mol.cmd.hide(selection="(r. DOD)")
    mol.cmd.orient()
    mol.cmd.center()
    mol.cmd.zoom()


def set_residue_colors(mol):
    """
    Set residue color for the structure.
    |cite| Warren L. DeLano (2002) CCP4 Newsletter on protein crystallography

    :param mol: PyMOL interface.
    :type mol: pymol2.PyMOL
    """
    for residue, residue_color in RESIDUE_COLOR_DICT.items():
        mol.cmd.color(color=residue_color, selection="(r. " + residue + ")")


def save_as_required(mol, display_type, save_path, file_name):
    """
    Save structure file as the requirement.
    |cite| Warren L. DeLano (2002) CCP4 Newsletter on protein crystallography

    :param mol: PyMOL interface.
    :type mol: pymol2.PyMOL

    :param display_type: type to display.
    :type display_type: DNAStructure.displayer.DisplayType

    :param save_path: path to save file.
    :type save_path: str

    :param file_name: file name.
    :type file_name: str
    """
    if display_type == DisplayType.picture:
        mol.cmd.ray(quiet=1)
        mol.cmd.png(filename=save_path + file_name, dpi=600, quiet=1)
    elif display_type == DisplayType.movie:
        pyplot.axis("off")
        frames, degree = [], 15
        for axis in ["x", "y", "z"]:
            for time in range(360 // degree + 1):
                mol.cmd.rotate(axis=axis, angle=degree)
                mol.cmd.ray(quiet=1)
                image_bits = mol.cmd.png(filename=None, dpi=600, quiet=1)
                frames.append(array(Image.open(BytesIO(image_bits))))
        patch = pyplot.imshow(frames[0])
        worker = animation.FuncAnimation(fig=pyplot.gcf(), func=lambda frame_index: patch.set_data(frames[frame_index]),
                                         frames=len(frames), interval=1)
        worker.save(save_path + file_name + ".gif", writer="pillow", fps=6)
        pyplot.close('all')
    else:
        raise ValueError("No such save type!")


if __name__ == '__main__':

    # Draw DNA structures based on different bases.
    displaytype = DisplayType.movie
    pdbfile_path = 'E:/DNAStructure/PDB_data/'
    pdb_files = listdir(pdbfile_path)
    output_path = 'E:/DNAStructure/PDB_structure/'

    for pdb_file in pdb_files:
        filename = pdb_file.split('.')[0]
        draw_dna_structure(cif_path=pdbfile_path + pdb_file, cif_name=filename,
                           save_path=output_path, display_type=displaytype)

    # draw_dna_structure(cif_path='E:/DNAStructure/PDB_data/1ana.cif', cif_name='1ana',
    #                    save_path='E:/DNAStructure/PDB_structure/', display_type=displaytype)
