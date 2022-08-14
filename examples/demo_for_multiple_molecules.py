from os import listdir, remove
from zipfile import ZipFile

from pymolbe import draw_colorfully, draw_specially, merge_pictures, merge_animations


if __name__ == "__main__":
    ZipFile("./data/m.zip", "r").extractall(path="./data/")
    child_paths = listdir("./data/m/")

    temp_paths, titles = [], []
    for index, child_path in enumerate(child_paths):
        temp_path = "./figures/temp" + str(index + 1).zfill(2) + ".png"
        draw_colorfully(load_path="./data/m/" + child_path, save_path=temp_path)
        titles.append(child_path[:4])
        temp_paths.append(temp_path)

    merge_pictures(load_paths=temp_paths, save_path="./figures/m_p.png", row_number=5, column_number=8,
                   figure_size=(10, 6), titles=titles)

    for used_path in temp_paths:
        remove(used_path)

    temp_paths, titles = [], []
    for index, child_path in enumerate(listdir("./data/m/")):
        temp_path = "./figures/temp" + str(index + 1).zfill(2) + ".gif"
        draw_specially(load_path="./data/m/" + child_path, save_path=temp_path, motif_colors={"HHHHH": "0xff0000"},
                       is_movie=True, shafts="y")
        titles.append(child_path[:4])
        temp_paths.append(temp_path)

    merge_animations(load_paths=temp_paths, save_path="./figures/m_ps.gif",
                     row_number=5, column_number=8, figure_size=(10, 6), titles=titles,
                     simultaneous=True)

    merge_animations(load_paths=temp_paths, save_path="./figures/m_pi.gif",
                     titles=titles, simultaneous=False, fps=20)

    for used_path in temp_paths:
        remove(used_path)

    for child_path in child_paths:
        remove("./data/m/" + child_path)
    remove("./data/m/")
