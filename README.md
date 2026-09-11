# MITgcmCoastalAntarctica
Code for three MITgcm models of the Antarctic regional circulation with a thermodynamically active ice shelf, sea ice parameterizations, and KPP mixing. Idealized coastal geometries are intended to produce general results.

ISEF: One year simulations of the on-shelf circulation. Code used to produce results in https://journals.ametsoc.org/view/journals/phoc/55/11/JPO-D-24-0217.1.xml. Takes ~12 hours to run on 100 processes. Zenodo DOI https://doi.org/10.5281/zenodo.20722285

CDWonshelf: Ten year simulations of the on-shelf circulation designed for longer runtimes. Takes ~16 hours to run on 144 processes. Zenodo DOI https://doi.org/10.5281/zenodo.16955044

ShelfBreak: One year simulations of the coastal circulation with model geometry including a continental shelf break. Takes ~12 hours to run on 288 processes. Not on Zenodo yet, too many bots there!

The steps to run these codes are all essentially the same. You will need access to an HPC with MATLAB, the MITgcm source code, and all necessary packages. The basic workflow is

(1) Compile MITgcm into a fresh build directory using the modifications in the compile directory and the appropriate options file for your system, eg
    mkdir build
    cd build
    /pathtoroot/MITgcm_20200601/tools/genmake2 -mpi -mods=/pathtomodel/compile_code -of=/pathtooptfile -rootdir=/pathtoroot/MITgcm_20200601
    make depend
    make
    ls mitgcmuv
    cd ../run

(2) The executable mitgcmuv should have appeared in your build directory. Then from the run directory, generate the initial conditions as bin files by running setup_all_run.m in MATLAB. 

(3) Run the MITgcm job with the specified number of processes. If you use a different number of processes, you will need to edit compile_code/SIZE.h. My sample scheduling script script_submit is provided for SLURM. The run command will look something like this depending on your system
  mpirun --mca btl_vader_fbox_max 0 -np NumProc /pathtomodel/build/mitgcmuv

(4) When the run is finished, glue the output together in the run directory. The script glue_all.sh needs to be modified first to target the right output (set NumProc and the desired output directory).
  nohup glue_all.sh > logfile_glue_all.log 2>&1 &

(5) Now you have run the model and created a series of netCDF files that can be processed in the programming language of your choice. Some sample output figures are provided under each model branch.
