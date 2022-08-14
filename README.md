# pymolBE
Customized batch exhibition of three-dimensional molecules based on the pymol framework

## Installation
You can install this package using pip:
```sh
pip install pymolbe
```
The packages requires a Python version >=3.7, 
as well as some basic libraries 
[PyMOL 2.5.0](https://pymol.org/2/),
[biopython 1.78](https://pypi.org/project/biopython/),
[matplotlib 3.1.1](https://pypi.org/project/matplotlib/),
[Pillow 8.2.0](https://pypi.org/project/Pillow/), and
[numpy](https://pypi.org/project/numpy/).


## Repository Structure
The structure of this tool is shown below:
```html
├── examples                                // Examples of pymolBE
├── pymolbe                                 // Source codes of pymolBE
│    ├── __init__.py                        // Calls of this tool
│    ├── converter.py                       // Conversion between different information
│    │    ├── protein_structure_to_pdb_file // Save 3D structure into pdb file
│    │    ├── pdb_file_to_protein_structure // Load 3D structure and its information from pdb file
│    │    ├── align                         // Rotate candidate structure unto reference structure using Kabsch algorithm based on each position
│    ├── displayer.py                       // Exhibition the structure data
│    │    ├── draw_colorfully               // Draw a molecular structure colorfully
│    │    ├── draw_specially                // Draw a molecular structure with special motifs
│    │    ├── merge_pictures                // Merge multiple structure pictures into a picture
│    │    ├── merge_animations              // Merge multiple structure animations into a animation
│    │    ├── set_initial_state             // Set the initial state of a structure with specific representation
│    │    ├── set_residue_colors            // Set residue color for the structure
│    │    ├── set_motif_colors              // Set the color of each motif
├── README.md                               // Description document of library
```

## Usage & Parameters
Please check each interface annotation for more information.
- [protein_structure_to_pdb_file](https://github.com/HaolingZHANG/pymolBE/blob/main/pymolbe/converter.py#L10)
- [pdb_file_to_protein_structure](https://github.com/HaolingZHANG/pymolBE/blob/main/pymolbe/converter.py#L100)
- [align](https://github.com/HaolingZHANG/pymolBE/blob/main/pymolbe/converter.py#L258)
- [draw_colorfully](https://github.com/HaolingZHANG/pymolBE/blob/main/pymolbe/displayer.py#L9)
- [draw_specially](https://github.com/HaolingZHANG/pymolBE/blob/main/pymolbe/displayer.py#L58)
- [merge_pictures](https://github.com/HaolingZHANG/pymolBE/blob/main/pymolbe/displayer.py#L107)
- [merge_animations](https://github.com/HaolingZHANG/pymolBE/blob/main/pymolbe/displayer.py#L150)
- [set_initial_state](https://github.com/HaolingZHANG/pymolBE/blob/main/pymolbe/displayer.py#L220)
- [set_residue_colors](https://github.com/HaolingZHANG/pymolBE/blob/main/pymolbe/displayer.py#L247)
- [set_motif_colors](https://github.com/HaolingZHANG/pymolBE/blob/main/pymolbe/displayer.py#L272)


## Exhibition Examples
### Draw a static picture for a PDB file.
run the following statement(s)

```python
from pymolbe import draw_colorfully

draw_colorfully(load_path="./data/p.pdb", 
                save_path="./figures/ps_c.png",
                hides=["(r. HOH)"])
```
and obtain

<p align="center">
<img src="./examples/figures/ps_c.png" alt="ps_c.png" title="ps_c.png", width=50%/>
</p>


### Draw an animation for a PDB file that revolved around the y-axis.
run the following statement(s)

```python
from pymolbe import draw_colorfully

draw_colorfully(load_path="./data/p.pdb", 
                save_path="./figures/ps_c.gif",
                hides=["(r. HOH)"], 
                is_movie=True, 
                shafts="y")
```
and obtain

<p align="center">
<img src="./examples/figures/ps_c.gif" alt="ps_c.png" title="ps_c.png", width=50%/>
</p>

### Draw an animation for a mmCIF file that revolved around the y-axis.
run the following statement(s)

```python
from pymolbe import draw_colorfully

draw_colorfully(load_path="./data/d.cif", 
                save_path="./figures/ds_c.gif",
                hides=["(r. HOH)"], 
                is_movie=True, 
                shafts="y")
```
and obtain

<p align="center">
<img src="./examples/figures/ds_c.gif" alt="ds_c.gif" title="ds_c.gif", width=50%/>
</p>


### Draw an animation for a mmCIF file that revolved around the y-axis with sticks format.
run the following statement(s)

```python
from pymolbe import draw_colorfully

draw_colorfully(load_path="./data/d.cif", 
                save_path="./figures/ds_s.gif",
                representation="sticks",  
                hides=["(r. HOH)"], 
                is_movie=True, 
                shafts="y")
```

and obtain

<p align="center">
<img src="./examples/figures/ds_s.gif" alt="ds_s.gif" title="ds_s.gif", width=50%/>
</p>

### Draw an animation for a PDB file that revolved around the y-axis with motif "RG".
run the following statement(s)

```python
from pymolbe import draw_specially

draw_specially(load_path="./data/p.pdb", 
               save_path="./figures/ps_s.gif", 
               motif_colors={"RG": "0xff0000"},
               hides=["(r. HOH)"], 
               is_movie=True, 
               shafts="y")
```

and obtain

<p align="center">
<img src="./examples/figures/ps_s.gif" alt="ps_s.gif" title="ps_s.gif", width=50%/>
</p>

### Draw an animation for a PDB file that revolved around the xyz-axis with motif "RG".
run the following statement(s)

```python
from pymolbe import draw_specially

draw_specially(load_path="./data/p.pdb", 
               save_path="./figures/ps_xyz.gif", 
               motif_colors={"RG": "0xff0000"},
               hides=["(r. HOH)"], 
               is_movie=True, 
               shafts="xyz")
```

and obtain

<p align="center">
<img src="./examples/figures/ps_xyz.gif" alt="ps_xyz.gif" title="ps_xyz.gif", width=50%/>
</p>
