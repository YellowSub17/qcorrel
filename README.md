# Qcorrel

Repo for diffraction intensity correlations in Qspace.

## File Outline

#### AM_inten_correl
Original qspace correlation script by Andrew Martin

#### correlate
My version on correlation script

#### edit_cif
Util module for editing cif files

#### padf_plug
Util module for connecting to padfpy (Complied only on Linux OS)

#### plot_and_process_utils
Util module for plotting correlations

#### Rfactors
Util module for determining errors and differences between correlations

#### symmetry
Util module for applying symmetry operations of hkl reflections


## .gitignore info

Witin my local repo, I have a few folders that aren't commited. This is to save space on commits. These are:

#### cifs
A directory for storing cif data

#### dbins
A directory for storing 3D arrays saved as dbins. These are usually Q-space correlation volumes are read by padf.

#### saved_plots
A directory where I save matplotlib figures

#### tiffs
A directory where I save tiffs
