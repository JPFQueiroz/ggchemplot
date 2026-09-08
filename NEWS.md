# ggchemplot 0.3.4 (September 8, 2026)

* Added color column to collapsed H labels data.  

# ggchemplot 0.3.3 (September 7, 2026)

* Added `show_atom_note()` to selectively set a label nearby an atom.  
* Added `flip_double_bond_side()` to flip a double bond between a pair of specified atoms.  

# ggchemplot 0.3.2 (September 7, 2026)

* Added `place_collapsed_h()` to selectively set the placement of collapsed H atom labels.  

# ggchemplot 0.3.1 (September 2, 2026)

* Added `show_atom_label()` to selectively restore explicit atom labels.  

# ggchemplot 0.3.0 (September 1, 2026)

* Improved logic for geometrical placement of collapsed hydrogen labels.
* Added options to flip structures horizontally and vertically in `ggchemplot1()`.
* Added option to read and plot stereo chemical bonds from SDF files.
* Added `uncollapse_hydrogens()` to selectively restore explicit hydrogens.
* Improved double and triple bond rendering following IUPAC recommendations
  (offset segment shortening based on adjoining bonds and double bond sideness).
* Added `save_chemplot()` with recommended SVG support and automatic
  background rectangle removal. The SVG file is useful for preparation of
  publication-quality figures using external tools (e.g. Inkscape).
* Fixed parameter resolution issues with `NULL` sizes and ggplot2 ≥ 4.0.
* Various robustness improvements.
* Added `normalize_structure()`, internally used by `ggchemplot1()` and 
  `ggchemplot2()` to normalize bond lengths from an SDF file to a 
  target value in ggplot2 coordinates system.
