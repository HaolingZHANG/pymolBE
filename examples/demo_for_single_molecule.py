from pymolbe import draw_colorfully, draw_specially

if __name__ == "__main__":
    draw_colorfully(load_path="./data/p.pdb", save_path="./figures/ps_c.png",
                    hides=["(r. HOH)"])
    draw_colorfully(load_path="./data/p.pdb", save_path="./figures/ps_c.gif",
                    hides=["(r. HOH)"], is_movie=True, shafts="y")
    draw_colorfully(load_path="./data/d.cif", save_path="./figures/ds_c.gif",
                    hides=["(r. HOH)"], is_movie=True, shafts="y")
    draw_colorfully(load_path="./data/d.cif", save_path="./figures/ds_s.gif",
                    representation="sticks",  hides=["(r. HOH)"], is_movie=True, shafts="y")
    draw_specially(load_path="./data/p.pdb", save_path="./figures/ps_s.gif", motif_colors={"RG": "0xff0000"},
                   hides=["(r. HOH)"], is_movie=True, shafts="y")
    draw_specially(load_path="./data/p.pdb", save_path="./figures/ps_xyz.gif", motif_colors={"RG": "0xff0000"},
                   hides=["(r. HOH)"], is_movie=True, shafts="xyz")
