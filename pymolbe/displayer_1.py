from io import BytesIO
from numpy import array
from pymol2 import PyMOL
from PIL import Image
from matplotlib import pyplot, animation


def draw_an_individual(segment, load_path, save_path, is_movie=False,
                       cluster_index=None, cluster_total=None, member_number=None):
    """
    Draw the structures of token one by one (or only the average structure).

    :param segment: segment.
    :type segment: str

    :param cluster_index: cluster index in the cluster group of this token.
    :type cluster_index: int

    :param cluster_total: the cluster number of the cluster group.
    :type cluster_total: int

    :param member_number: member number in this cluster.
    :type member_number: int

    :param load_path: path to load structure file.
    :type load_path: str

    :param save_path: parent path to save file.
    :type save_path: str

    :param is_movie: display through animation.
    :type is_movie: bool

    .. note::
        reference: Warren L. DeLano (2002) CCP4 Newsletter on protein crystallography
    """
    token_info = "[" + str(len(segment)) + "]" + segment
    if cluster_index is not None and cluster_total is not None:
        if member_number is not None:
            cluster_info = "[" + str(cluster_index + 1) + "-" + str(cluster_total) + "(" + str(member_number) + ")]"
        else:
            cluster_info = "[" + str(cluster_index + 1) + "-" + str(cluster_total) + "]"

    else:
        cluster_info = ""

    with PyMOL() as mol:
        mol.cmd.load(load_path, quiet=1)
        set_initial_state(mol=mol)
        set_residue_colors(mol=mol)
        save_as_required(mol=mol, is_movie=is_movie, save_path=save_path, file_name=token_info + cluster_info)


def draw_a_population(token_segment, cluster_index, cluster_total, load_paths, save_path,
                      is_movie=False, display_limit=1):
    """
    Draw the structures of token as the cluster.

    :param token_segment: token segment.
    :type token_segment: str

    :param cluster_index: cluster index in the cluster group of this token.
    :type cluster_index: int

    :param cluster_total: the cluster number of the cluster group.
    :type cluster_total: int

    :param load_paths: paths to load different structure files.
    :type load_paths: list

    :param save_path: parent path to save file.
    :type save_path: str

    :param is_movie: display through animation.
    :type is_movie: bool

    :param display_limit: maximum number (limit) of individual needs to be display.
    :type display_limit: int

    .. note::
        reference: Warren L. DeLano (2002) CCP4 Newsletter on protein crystallography
    """
    token_info = "[" + str(len(token_segment)) + "]" + token_segment
    cluster_info = "[" + str(cluster_index + 1) + "-" + str(cluster_total) + "(" + str(len(load_paths)) + ")]"

    with PyMOL() as mol:
        if display_limit is None or display_limit >= len(load_paths):
            for load_path in load_paths:
                mol.cmd.load(load_path, quiet=1)
        else:
            for display_index in range(0, len(load_paths), len(load_paths) // display_limit):
                mol.cmd.load(load_paths[display_index], quiet=1)
        set_initial_state(mol=mol)
        set_residue_colors(mol=mol)
        save_as_required(mol=mol, is_movie=is_movie, save_path=save_path, file_name=token_info + cluster_info)


def draw_motifs_in_a_protein(motif_colors, load_path, structure_name, save_path, is_movie=False):
    """
    Draw a protein by the motif vision.

    :param motif_colors: pair of the motif and the color.
    :type motif_colors: dict

    :param load_path: path to load structure file.
    :type load_path: str

    :param structure_name: display name for saving.
    :type structure_name: str

    :param save_path: parent path to save file.
    :type save_path: str

    :param is_movie: display through animation.
    :type is_movie: bool

    .. note::
        reference: Warren L. DeLano (2002) CCP4 Newsletter on protein crystallography
    """
    with PyMOL() as mol:
        mol.cmd.load(load_path, quiet=1)
        set_initial_state(mol=mol)
        set_motif_colors(mol=mol, motif_colors=motif_colors)
        save_as_required(mol=mol, is_movie=is_movie, save_path=save_path, file_name=structure_name)


def draw_bases_in_a_dna(residue_colors, load_path, structure_name, save_path, is_movie=False):
    """
    Draw a dna by the base type.

    :param residue_colors: pair of the residue and the color.
    :type residue_colors: dict

    :param load_path: path to load structure file.
    :type load_path: str

    :param structure_name: display name for saving.
    :type structure_name: str

    :param save_path: parent path to save file.
    :type save_path: str

    :param is_movie: display through animation.
    :type is_movie: bool

    .. note::
        reference: Warren L. DeLano (2002) CCP4 Newsletter on protein crystallography
    """
    with PyMOL() as mol:
        mol.cmd.load(load_path, quiet=1)
        set_initial_state(mol=mol)
        set_residue_colors(mol=mol, residue_colors=residue_colors)
        save_as_required(mol=mol, is_movie=is_movie, save_path=save_path, file_name=structure_name)


def set_initial_state(mol, representation="cartoon"):
    """
    Set the initial state of a structure.

    :param mol: PyMOL interface.
    :type mol: pymol2.PyMOL

    :param representation: representation type.
    :type representation: str
    """
    if representation is not "cartoon":
        mol.cmd.hide(representation="cartoon")
        mol.cmd.show(representation=representation)
    mol.cmd.hide(selection="(r. HOH)")
    mol.cmd.orient()
    mol.cmd.center()
    mol.cmd.zoom()


def set_motif_colors(mol, motif_colors, neglected_color=None):
    """
    Set the colors of token.

    :param mol: PyMOL interface.
    :type mol: pymol2.PyMOL

    :param motif_colors: pair of the motif and the color.
    :type motif_colors: dict

    :param neglected_color: neglected color.
    :type neglected_color: str
    """
    mol.cmd.color(color=neglected_color if neglected_color is not None else "0xFFFFCC", selection="(all)")
    for motif, concerned_color in motif_colors.items():
        mol.cmd.color(color=concerned_color, selection="(ps. " + motif + ")")


def set_residue_colors(mol, residue_colors=None, neglected_color=None):
    """
    Set residue color for the structure.

    :param mol: PyMOL interface.
    :type mol: pymol2.PyMOL

    :param residue_colors: pair of the residue and the color.
    :type residue_colors: dict

    :param neglected_color: neglected color.
    :type neglected_color: str
    """
    if residue_colors is None:
        residue_colors = {"ALA": "0x8A685C", "ARG": "0x402F42", "ASN": "0x7fA9C2", "ASP": "0x632D3B", "CYS": "0x7D779D",
                          "GLN": "0x853D2F", "GLU": "0x88191F", "GLY": "0x945A4F", "HIS": "0xA29DB3", "ILE": "0x645D87",
                          "LEU": "0x82A293", "LYS": "0xA58121", "MET": "0x6273A1", "PHE": "0xB33C24", "PRO": "0x73584D",
                          "SER": "0x4F698A", "THR": "0xB9B9BB", "TRP": "0x686C47", "TYR": "0x674E3A", "VAL": "0x9D491B",
                          "DA": "0xff0000", "DT": "0x0000ff", "DC": "0x00ff00", "DG": "0xffff00"}

    mol.cmd.color(color=neglected_color if neglected_color is not None else "0xFFFFCC", selection="(all)")
    for residue, hex_color in residue_colors.items():
        mol.cmd.color(color=hex_color, selection="(r. " + residue + ")")


def save_as_required(mol, save_path, file_name, dpi=600, is_movie=False, degree=30, fps=20):
    """
    Save structure file as the requirement.

    :param mol: PyMOL interface.
    :type mol: pymol2.PyMOL

    :param save_path: path to save file.
    :type save_path: str

    :param file_name: file name.
    :type file_name: str

    :param dpi: dots per inch.
    :type dpi: int

    :param is_movie: display through animation.
    :type is_movie: bool

    :param degree: single rotation angle (available at is_movie=True).
    :type degree: int

    :param fps: frames per second (available at is_movie=True).
    :type fps: int
    """
    if is_movie:
        pyplot.axis("off")
        frames = []
        for axis in ["x", "y", "z"]:
            for time in range(360 // degree + 1):
                mol.cmd.rotate(axis=axis, angle=degree)
                mol.cmd.ray(quiet=1)
                image_bits = mol.cmd.png(filename=None, dpi=dpi, quiet=1)
                frames.append(array(Image.open(BytesIO(image_bits))))
        patch = pyplot.imshow(frames[0])
        worker = animation.FuncAnimation(fig=pyplot.gcf(), frames=len(frames), interval=1,
                                         func=lambda frame_index: patch.set_data(frames[frame_index]))
        worker.save(save_path + file_name + ".gif", writer="pillow", fps=fps)
        pyplot.close()
    else:
        mol.cmd.ray(quiet=1)
        mol.cmd.png(filename=save_path + file_name, dpi=dpi, quiet=1)
