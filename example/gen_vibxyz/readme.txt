transfer the normal modes writed out by DL-FIND to vibmodes.xyz format - this file will be used to visualise the modes pointing arrows using VESTA

to run:
python3 gen_PCSvibxyz.py normal_modes.dat 45 15
input: normal_modes.dat  45(number of modes) 15(number of atoms)
output: vibmodes.xyz
