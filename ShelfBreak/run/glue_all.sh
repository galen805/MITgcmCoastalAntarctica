#!/bin/sh

# glues desired mnc_test*/whatever files together
# Make sure you are in the run directory
# Creates output directory if it doesn't already exist

directory="Output"
range="{0001..NumProc}"

mkdir -p "$directory"

target_file="$directory/sice_glued.nc"
mnc_files="mnc_test_$range/sice.*.nc"
filegroup=$(eval ls "$mnc_files")
gluemncbig --many -o "$target_file" $filegroup

target_file="$directory/ptracers_glued.nc"
mnc_files="mnc_test_$range/ptracers.*.nc"
filegroup=$(eval ls "$mnc_files")
gluemncbig --many -o "$target_file" $filegroup

target_file="$directory/state_glued.nc"
mnc_files="mnc_test_$range/state.*.nc"
filegroup=$(eval ls "$mnc_files")
gluemncbig --many -o "$target_file" $filegroup

target_file="$directory/DensityAnomaly_glued.nc"
mnc_files="mnc_test_$range/DensityAnomaly.*.nc"
filegroup=$(eval ls "$mnc_files")
gluemncbig --many -o "$target_file" $filegroup

target_file="$directory/VorticityAvg_glued.nc"
mnc_files="mnc_test_$range/VorticityAvg.*.nc"
filegroup=$(eval ls "$mnc_files")
gluemncbig --many -o "$target_file" $filegroup

target_file="$directory/kpp_state_glued.nc"
mnc_files="mnc_test_$range/kpp_state.*.nc"
filegroup=$(eval ls "$mnc_files")
gluemncbig --many -o "$target_file" $filegroup


target_file="$directory/SHIFlxAvg_glued.nc"
mnc_files="mnc_test_$range/SHIFlxAvg.*.nc"
filegroup=$(eval ls "$mnc_files")
gluemncbig --many -o "$target_file" $filegroup


target_file="$directory/grid_glued.nc"
mnc_files="mnc_test_$range/grid.*.nc"
filegroup=$(eval ls "$mnc_files")
gluemncbig --many -o "$target_file" $filegroup

target_file="$directory/seaiceFluxAvg_glued.nc"
mnc_files="mnc_test_$range/seaiceFluxAvg.*.nc"
filegroup=$(eval ls "$mnc_files")
gluemncbig --many -o "$target_file" $filegroup


