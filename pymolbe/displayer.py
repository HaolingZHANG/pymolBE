from io import BytesIO
from numpy import array
from os import environ, remove
from pymol2 import PyMOL
from PIL import Image, ImageSequence
from matplotlib import pyplot, animation


def draw_colorfully(load_path, save_path, information,
                    representation=None, residue_colors=None, neglected_color=None,
                    dpi=400, is_movie=False, shafts="xyz", degree=10, fps=10):
    """
    Draw a molecular structure colorfully.

    :param load_path: path to load structure file.
    :type load_path: str

    :param save_path: path to save display file.
    :type save_path: str

    :param information: segment information.
    :type information: str

    :param representation: representation type.
    :type representation: str

    :param residue_colors: pair of the residue and the color.
    :type residue_colors: dict or None

    :param neglected_color: neglected color.
    :type neglected_color: str or None

    :param dpi: dots per inch.
    :type dpi: int

    :param is_movie: display through animation.
    :type is_movie: bool

    :param shafts: rotating shafts.
    :type shafts: str

    :param is_movie: display through animation.
    :type is_movie: bool

    :param degree: single rotation angle (available at is_movie=True).
    :type degree: int

    :param fps: frames per second (available at is_movie=True).
    :type fps: int
    """
    with PyMOL() as mol:
        mol.cmd.load(load_path, quiet=1)
        set_initial_state(mol=mol, representation=representation)
        set_residue_colors(mol=mol, residue_colors=residue_colors, neglected_color=neglected_color)
        save_figure(mol=mol, is_movie=is_movie, save_path=save_path, file_name=information,
                    shafts=shafts, dpi=dpi, degree=degree, fps=fps)


def draw_specially(load_path, save_path, information,
                   representation=None, motif_colors=None, neglected_color=None,
                   dpi=400, is_movie=False, shafts="xyz", degree=10, fps=10):
    """
    Draw a molecular structure with special motifs.

    :param load_path: path to load structure file.
    :type load_path: str

    :param save_path: path to save display file.
    :type save_path: str

    :param information: segment information.
    :type information: str

    :param representation: representation type.
    :type representation: str

    :param motif_colors: pair of the motif and the color.
    :type motif_colors: dict

    :param neglected_color: neglected color.
    :type neglected_color: str or None

    :param dpi: dots per inch.
    :type dpi: int

    :param is_movie: display through animation.
    :type is_movie: bool

    :param shafts: rotating shafts.
    :type shafts: str

    :param is_movie: display through animation.
    :type is_movie: bool

    :param degree: single rotation angle (available at is_movie=True).
    :type degree: int

    :param fps: frames per second (available at is_movie=True).
    :type fps: int
    """
    with PyMOL() as mol:
        mol.cmd.load(load_path, quiet=1)
        set_initial_state(mol=mol, representation=representation)
        set_motif_colors(mol=mol, motif_colors=motif_colors, neglected_color=neglected_color)
        save_figure(mol=mol, is_movie=is_movie, save_path=save_path, file_name=information,
                    shafts=shafts, dpi=dpi, degree=degree, fps=fps)


def merge_pictures(load_paths, save_path, row_number, column_number, figure_size, titles=None, dpi=200):
    """
    Merge multiple structure pictures into a picture.

    :param load_paths: paths to load structure picture.
    :type load_paths: list

    :param save_path: path to save display file.
    :type save_path: str

    :param row_number: picture number in a row.
    :type row_number: int

    :param column_number: picture number in a column.
    :type column_number: int

    :param figure_size: size of saved figure.
    :type figure_size: tuple

    :param titles: titles of each load paths if required.
    :type titles: list or None

    :param dpi: dots per inch.
    :type dpi: int
    """
    assert len(load_paths) <= row_number * column_number
    if titles is not None:
        assert len(titles) == len(load_paths)

    pyplot.figure(figsize=figure_size)
    for location in range(len(load_paths)):
        pyplot.subplot(row_number, column_number, location + 1)
        if titles is not None:
            pyplot.title(titles[location])
        pyplot.imshow(Image.open(load_paths[location]))
        pyplot.xticks([])
        pyplot.yticks([])
        pyplot.axis('off')

    pyplot.savefig(save_path, dpi=dpi)
    pyplot.close()


def merge_animations(load_paths, save_path, row_number, column_number, figure_size, titles=None, fps=10, dpi=200):
    """
    Merge multiple structure animations into a animation.

    :param load_paths: paths to load structure animation.
    :type load_paths: list

    :param save_path: path to save display file.
    :type save_path: str

    :param row_number: picture number in a row.
    :type row_number: int

    :param column_number: picture number in a column.
    :type column_number: int

    :param figure_size: size of saved figure.
    :type figure_size: tuple

    :param titles: titles of each load paths if required.
    :type titles: list or None

    :param fps: frames per second (available at is_movie=True).
    :type fps: int

    :param dpi: dots per inch.
    :type dpi: int
    """
    assert len(load_paths) <= row_number * column_number
    if titles is not None:
        assert len(titles) == len(load_paths)

    group, number = None, None
    for load_path in load_paths:
        with Image.open(load_path) as image:
            frames = ImageSequence.all_frames(image)
            if number is not None:
                assert number == len(frames)
            else:
                number = len(frames)
                group = [[] for _ in range(number)]
            for frame_index, frame in enumerate(frames):
                group[frame_index].append(frame)

    new_frames = []
    for frame_index, images in enumerate(group):
        pyplot.figure(figsize=figure_size)
        for location, image in enumerate(images):
            pyplot.subplot(row_number, column_number, location + 1)
            if titles is not None:
                pyplot.title(titles[location])
            pyplot.imshow(Image.open(load_paths[location]))
            pyplot.xticks([])
            pyplot.yticks([])
            pyplot.axis('off')

        temp_path = environ["TEMP"] + str(frame_index) + ".png"
        pyplot.savefig(temp_path, dpi=dpi)
        pyplot.close()

        new_frames.append(array(Image.open(temp_path)))
        remove(path=temp_path)

    patch = pyplot.imshow(frames[0])
    worker = animation.FuncAnimation(fig=pyplot.gcf(), frames=len(frames), interval=1,
                                     func=lambda index: patch.set_data(frames[index]))
    worker.save(save_path, writer="pillow", fps=fps)
    pyplot.close()


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
                          "DA": "0xf2521b", "DT": "0xfabc09", "DC": "0x81cc28", "DG": "0x00aef0"}

    mol.cmd.color(color=neglected_color if neglected_color is not None else "0xFFFFCC", selection="(all)")
    for residue, hex_color in residue_colors.items():
        mol.cmd.color(color=hex_color, selection="(r. " + residue + ")")


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


def save_figure(mol, save_path, file_name, shafts="xyz", dpi=400, is_movie=False, degree=10, fps=10):
    """
    Save structure file as the requirement.

    :param mol: PyMOL interface.
    :type mol: pymol2.PyMOL

    :param save_path: path to save file.
    :type save_path: str

    :param file_name: file name.
    :type file_name: str

    :param shafts: rotating shafts.
    :type shafts: str

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
        for axis in shafts:
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
