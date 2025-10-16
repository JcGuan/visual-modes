# This script reads info from PCS results and generates a .vesta file to visualise qm or qm/mm structure
# plus associted vibrational modes

# CELLP lattice vectors lengths set to 1
# charge the coords in STRUC to xyz coords
# SITE change the colours consistent with Aten
# VECTR add vibrational modes here
# VECTT radius and colours of arrows, last digit = 2 arrows start at atom

"""
Author: J.Guan@UCL
jingcheng.guan@ucl.ac.uk
use upon request, pls respect the author's coding efforts
"""

from sys import argv
from colour_table import table_colours 
#print("dbg only: table_colours")
#for i in range(len(table_colours)):
#    print(table_colours[i][8])
#print(table_colours[i] for i in range(len(table_colours)))
class modesvisual():

    def visual_setup(self):
        """
        user set up the lib for visualisation
        """
        self.atoms_colourslib={'Si': "200 200 200 0 0 0 0",
                               'O' : "200 200 200 0 0 0 0",
                               'Cu': "200 200 200 0 0 0 0",
                               'Al': "200 200 200 0 0 0 0",
                               'H' : "200 200 200 0 0 0 0",
                               'N' : "200 200 200 0 0 0 0",
                               'C' : "0 0 255 240 255 255 255",
                               'P' : "200 200 200 0 0 0 0",
                         }
        self.atoms_sizescale={'H': 0.4,
                              'N': 0.6}
        
        self.rgba_table = {}
        self.rgba_spec  = 255
        self.coval_rad  = {}
        for i in range(len(table_colours)):
            print(table_colours[i][8])
            # specie - rgba
            self.rgba_table[table_colours[i][4]] = table_colours[i][8]

        self.atom_radii = {} # radii of atoms themselves
        for i in range(len(table_colours)):
            # specie - atomic radius
            self.atom_radii[table_colours[i][4]] = table_colours[i][7] 

        print("RGBA table:\n {} \n   ".format(self.rgba_table))
        print("atomic radii: \n {} \n".format(self.atom_radii))
        # read eledata and get the covalence radii of elements
        with open("eledata", 'r') as f:
            self.eledata = f.readlines()
            for i, line in enumerate(self.eledata):
                if i == 0:
                    continue
                else:
                    specie = line.strip().split()[0]
                    self.coval_rad[specie] = float(line.strip().split()[3])
        print("covalence radii: \n {} \n".format(self.coval_rad))
                

    def __init__(self):
        try:
            self.struc_filename = argv[1]
        except:
            exit("no strucutre filename (.xyz) given")
        try:
            self.vibmodes_filename = argv[2]
        except:
            self.vibmodes_filename = "vibmodes.xyz"
        try:
            self.systemname = argv[3]
        except:
            self.systemname = "XXX"
        try:
            self.qmmm_config = argv[4]
        except:
            self.qmmm_config = False
        if self.qmmm_config:
            print("PLEASE NOTE:")
            print("MAKE SURE the vib active atoms are at the very top of the xyz file")

        self.visual_setup()

        self.vec_lscale     = 8.0 # vib eigenvectors length scaling control
        self.modevec_radius = 0.2
        self.modevec_colour = '0 0 0' # r|g|b (green: 0 255 0; black: 0 0 0) 
        self.modevec_start  = 2 # 2: mode vec start from the centre of atom
        self.natoms         = 0
        self.natoms1        = 0 # number of vib active atoms
        self.nmodes         = 0
        self.coord_scale    = 1.0
        self.bond_tolerance = 2.5 # angstrom, within this radius there's a bond, default=2.5 (A)
        self.bond_scale     = 1.1 # default = 1.0, i.e. use the default value, could be 1.1/.2
        self.bond_radius    = 0.2 
        self.species        = []
        self.pairs_species  = []
        self.vesta_filename = "{}.mode{}.{}.vesta"
        self.vesta_filestr  = ""
        self.vesta_front    = """#VESTA_FORMAT_VERSION 3.3.0
MOLECULE

TITLE
{} | MODE-{} | {} cm^-1                                                                          

GROUP
1 1 P 1
SYMOP
 0.000000  0.000000  0.000000  1  0  0   0  1  0   0  0  1   1
 -1.0 -1.0 -1.0  0 0 0  0 0 0  0 0 0
TRANM 0
 0.000000  0.000000  0.000000  1  0  0   0  1  0   0  0  1
LTRANSL
 -1
 0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
LORIENT
 -1   0   0   0   0
 1.000000  0.000000  0.000000  1.000000  0.000000  0.000000
 0.000000  0.000000  1.000000  0.000000  0.000000  1.000000
LMATRIX
 1.000000  0.000000  0.000000  0.000000
 0.000000  1.000000  0.000000  0.000000
 0.000000  0.000000  1.000000  0.000000
 0.000000  0.000000  0.000000  1.000000
 0.000000  0.000000  0.000000
CELLP
  1.000000   1.000000   1.000000  90.000000  90.000000  90.000000
  0.000000   0.000000   0.000000   0.000000   0.000000   0.000000
"""
        
        self.vesta_struc = ""
        self.vesta_bound = """BOUND
       0        1         0        1         0        1   
  0   0   0   0  0
"""
        self.vesta_sbond = ""
        self.vesta_sitet = ""
        self.vesta_vectr = ""        
        self.vesta_vectt = ""

    def genvesta_sbond(self):
        """
        gen vesta sbond
        """
        self.vesta_sbond  = "SBOND \n"
        for p in range(len(self.pairs_species)):
            specie1 = self.pairs_species[p][0]
            specie2 = self.pairs_species[p][1]
            self.bond_tolerance = self.coval_rad[specie1]+self.coval_rad[specie2]
            self.bond_tolerance = self.bond_scale * float(self.bond_tolerance)
            print("specie1: {}, specie2: {}, bond tolerance: {}".format(specie1, specie2, self.bond_tolerance))
            # self.bond_radius is the size of the bonding sticks in visualisation
            # self.bond_tolerance is the covalence radius
            self.vesta_sbond += "   1     {}     {}    0.00000    {}  0  1  0  0  1  {}  2.000 127 127 127 \n".format(specie1, specie2, self.bond_tolerance, self.bond_radius) 
        self.vesta_sbond += "0 0 0 0 \n" 

    def readpcs_genvesta_struc(self):
        """
        read atomic structure (.xyz file) and      
        for qm/mm it will read the whole qm region
        """
        with open(self.struc_filename) as xyzf:
            self.xyzlines = xyzf.readlines()

        try:
            self.natoms = int(self.xyzlines[0].strip())
        except:
            exit("cannot read num of atoms from the 1st line of xyz file")

        self.vesta_struc = "STRUC \n"

        #print(self.xyzlines)
        for i in range(self.natoms):
            #print("dbg only")
            #print(self.xyzlines[i+2].strip('\n').split()) 
            coord_str = ""  
            for xyz in self.xyzlines[i+2].strip('\n').split()[1:4]: #skipo natoms and string
                coord_str += "{}  ".format(xyz)
            #coord_str += "\n"
            #print(coord_str)
            specie = self.xyzlines[i+2].strip().split()[0]

            if specie not in self.species:
                self.species.append(specie)

            #print(specie)
            
            #Z specie specie scale dx dy dz 1a 1
            #0.0 0.0 0.0 0.0
            self.vesta_struc += "{} {} {} {} {} {} \n {} \n".format(i+1, specie, specie, self.coord_scale, coord_str, '1a 1', '0.0 0.0 0.0 0.0')

        self.vesta_struc += "  0 0 0 0 0 0 0 \n"  

        self.gen_speciespairs()

    def gen_speciespairs(self):
        """
        to gen bond info between pairs of atoms
        """
        for i in range(len(self.species)):
            for j in range(len(self.species)): 
                sp1 = self.species[i]
                sp2 = self.species[j]
                #if sp1 != sp2:
                self.pairs_species.append([sp1, sp2])
                #else:
                #    continue
        #print("species pairs:")
        #print(self.pairs_species)
    

    def genvesta_sitet(self):
        """
        fixme: not sure what is the use of this section
        """
        self.vesta_sitet = "SITET \n"
        for i in range(self.natoms):
            specie = self.xyzlines[i+2].strip().split()[0] 
            self.vesta_sitet += "{} {} {} {} \n".format(i+1, specie, self.atoms_sizescale[specie], self.atoms_colourslib[specie])
        self.vesta_sitet += "  0 0 0 0 0 0 \n"

#    def readpcs_vibmodes_genvesta_vectr(self):
#        """read pcs vib modes (similar format as .xyz file) and generate the vesta VECTR"""
#
#        self.vesta_vectr = "VECTR \n"
# 
#        for m in range(self.nmodes):
#            
#            self.modeid, self.wavenumber = self.modeslines[m*(3+self.natoms) + 1].strip().split()
#
#            for i in range(self.natoms):                    
#            
#                modevector_peratom = "{} {} \n".format(i+1, self.modeslines[m*(3+self.natoms) + i+2])
#                modevector_peratom += "{}  0   0    0    0 \n".format(i+1)
#                modevector_peratom += " 0 0 0 0 0 \n"
#
#            self.vesta_vectr += modevector_peratom 
#            self.vesta_vectr += " 0 0 0 0 0"

    def genvesta_vectt(self):
        """
        write out a line in the VECTT section
        VECTT controls the colours and radius of the modes vectors
        """ 
        self.vesta_vectt = "VECTT \n" 
        for i in range(self.natoms1):
            self.vesta_vectt += "{} {} {} {} \n".format(i+1, self.modevec_radius, self.modevec_colour, self.modevec_start)
        self.vesta_vectt += '0 0 0 0 0 \n'

    def genvesta_atomt(self):
        """
        generate the ATOMT section
        ATOMT controls the colours of atoms
        """
        self.vesta_atomt = "ATOMT \n"
        for i in range(self.natoms):
            specie = self.xyzlines[i+2].strip().split()[0]
            r = int(self.rgba_table[specie][0]*self.rgba_spec)
            g = int(self.rgba_table[specie][1]*self.rgba_spec)
            b = int(self.rgba_table[specie][2]*self.rgba_spec)
            print("dbg only: rgba: {}".format(r))
            print("dbg only: rgba: {}".format(g))
            print("dbg only: rgba: {}".format(b))
            self.vesta_atomt += "{} {} {}  {} {} {}  0 0 0 0 \n".format(i+1, specie, self.atom_radii[specie], r, g, b) 
        self.vesta_atomt += '0 0 0 0 0 0 \n'

    def read_vibmodesxyz(self):
        with open(self.vibmodes_filename) as modesf: 
            self.modeslines = modesf.readlines()
        try:
            self.natoms1 = int(self.modeslines[0])
        except:
            exit("cannot read num of atoms from vib modes file")
        if not self.qmmm_config and self.natoms != self.natoms1:
            exit("ERROR: number of atoms in vib modes WRONG - unless this is a qmmm config")
        # get the total number of modes
        for line in self.modeslines:
            if "END MODE" in line:
                self.nmodes += 1
        print("total number of modes = {}".format(self.nmodes))

    def pcs2vesta_combine(self):
        """
        generate vesta file from combination of pcs results
        """
        # gen vesta struc section
        self.readpcs_genvesta_struc()
        # read and gen the vib active atoms
        self.read_vibmodesxyz()
        # gen vesta sitet section
        #self.genvesta_sitet()
        # gen vesta vectt section
        self.genvesta_vectt()
        # gen vesta sbond section
        self.genvesta_sbond()
        # gen vesta atomt section (colours of atoms)
        self.genvesta_atomt()

        # generate vesta vectr section for each vib mode (looped)
        # VECTR controls the directions of mode vectors
        # gen vesta file per mode
        for m in range(self.nmodes):
            self.vesta_vectr = "VECTR \n"
            
            self.modeid, self.wavenumber = self.modeslines[m*(3+self.natoms1) + 1].strip().split()
            print("check wave number read from file:{}".format(self.wavenumber))

            modevector_peratom = ""
            for i in range(self.natoms1): 
                modevec_str = ""
                #print(self.modeslines[m*(3+self.natoms) + i+2].strip('\n').split())
                for dxyz in self.modeslines[m*(3+self.natoms1) + i+2].strip('\n').split():
                    # do a nil vib vector check
                    if abs(float(dxyz)) <= 0.001:
                        dxyz = 0.0
                    modevec_str += "{} ".format(self.vec_lscale * float(dxyz))
                print("dbg only1: mode vector str")
                print(modevec_str)
                modevector_peratom += "{} {} \n".format(i+1, modevec_str)
                modevector_peratom += "{}  0   0    0    0 \n".format(i+1)
                modevector_peratom += " 0 0 0 0 0 \n"

            self.vesta_vectr += modevector_peratom 
            self.vesta_vectr += " 0 0 0 0 0 \n"
            # gen the beginning of the vesta file with title section
            # til the struc section
            self.vesta_front_valued = self.vesta_front.format(self.systemname, self.modeid, self.wavenumber)
            print("mode no.{}".format(self.modeid))
            print("wavenumber={}".format(self.wavenumber))
            vesta_filename  = self.vesta_filename.format(self.systemname, self.modeid, int(float(self.wavenumber)))
            print("vesta filename:")
            print(vesta_filename)
            self.vesta_filestr += self.vesta_front_valued
            self.vesta_filestr += self.vesta_struc
            self.vesta_filestr += self.vesta_bound
            self.vesta_filestr += self.vesta_sbond
            #self.vesta_filestr += self.vesta_sitet
            self.vesta_filestr += self.vesta_vectr
            self.vesta_filestr += self.vesta_vectt
            self.vesta_filestr += self.vesta_atomt
 
            with open(vesta_filename, 'w') as vestaf:
                vestaf.write(self.vesta_filestr)
                self.vesta_filestr = ""



mv = modesvisual()
mv.pcs2vesta_combine()
