# STANCE
IMPORTANT: Requires SPM8 to run.

SPM8 installation: http://www.fil.ion.ucl.ac.uk/spm/software/spm8/

If installing to the filepath e.g. /C:/spm/spm8 it is suggested to install STANCE in the parent dir

The first time running STANCE in the installation folder run the script:

STANCE_initialize_STANCE.m

NOTE: STANCE will prompt a new user only once for the SPM8 path, but it must exists for STANCE to initialize properly. STANCE will save SPMpath to STANCE.mat. If your spm8 directory gets moved you may need to edit STANCE.mat to change the path.

STANCE is useful for developing scripts for fMRI data simulation.

See the User Guide, Simulation Pipeline and the demos included in the "scripts_for_demos" folder or their PDF reports in the "reports" folder for more details.

7 July 2017, Xiangyu Liu and Jason E. Hill.
