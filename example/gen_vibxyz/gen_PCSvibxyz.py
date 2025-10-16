from sys import argv
import numpy as np
def gen_vib_xyz():
    # file name of modes - columns are modes;
    # rows are atoms xyz/disp
    fn = argv[1]
    n_modes = int(argv[2])
    n_atoms = int(argv[3])
    # vib modes array used to write out modes.xyz
    vibxyz  = np.zeros((n_atoms*n_modes, 3), dtype=float)
    read_ncol = False
    arow = None
    freq_file  = "freqs_wavenumber.dat"
    outputfile = "vibmodes.xyz"
    with open(fn, 'r') as f:
        blockno = 0
        for line in f:
            line = line.strip().split()
            print(line)
            # read num of cols once
            if not read_ncol:
                ncol = len(line)
                print("tot num of cols: {}".format(ncol))
                read_ncol = True
                arow = np.zeros(ncol, dtype=float)
            # read line no
            if len(line) != 0:
                lineno = int(line[0])
            else:
                print("empty line")
                print(line)
                lineno = 0
            arow = line[1:] 
            print("a row")
            print(arow)
            for colno, modedisp in enumerate(arow):
                print("dbg only: modedisp={}".format(modedisp))
                #if lineno%3 == 0:
                print("(lineno-1)%3={}".format((lineno-1)%3))
                # ncol - 1 as the 0th is the line no.
                # for example nh3, each 12 modes is a block
                blocknum = blockno*(ncol-1)*n_atoms
                # in each block there are 12 modes
                modenum  = colno * n_atoms
                # in each mode, there are 15 atoms - x,y,z
                atomnum  = int((lineno-1)/3)
                print("atom number: {}".format(atomnum))
                vibxyz[blocknum+modenum+atomnum, (lineno-1)%3] = modedisp
                #else:
                #    vibxyz[blockno*ncol+colno, colno%3-1] = modedisp
            print("current vib xyz array:")
            print(vibxyz)
            #quit()
            print("line no.={}".format(lineno))
            if lineno == 3*n_atoms:
                blockno += 1
                print ("block no.{}".format(blockno))
        print("vib xyz array:")
        print(vibxyz)
        np.savetxt("vibmodes.xyz", vibxyz)

    with open(freq_file, "r") as f_freq:
        with open(outputfile,"w") as fout:
            for modeid in range(n_modes):
                freq = float(f_freq.readline().strip().split()[0])
                fout.write("{}\n".format(n_atoms))
                fout.write("{}  {}\n".format(modeid+1, freq))
                for atomno in range(n_atoms):
                    print(vibxyz[modeid*n_atoms + atomno])
                    for xyz in range(3):
                        fout.write("{}    ".format(vibxyz[modeid*n_atoms + atomno][xyz]))
                    fout.write("\n")
                fout.write("END MODE \n")
                    


gen_vib_xyz()
