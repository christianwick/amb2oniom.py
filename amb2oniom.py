#!/usr/bin/env python3
"""
use parmed to read in topology and write out a gaussian input file for ONIOM
calculations with 2 or 3 layers.

Author: Christian R. Wick
"""


from parmed import __version__ as pmd_version
from parmed import amber
from parmed.periodic_table import Element

from collections import OrderedDict
from itertools import permutations
import sys
import re



#-----------------------------------------------------------------------------#

class amb2oniom():
    """ amb2oniom base class

    Arguments:
    -----------
    top : :class:`parmed.amber.AmberParm`
        Amber topology
    xyz : str or TextIOWrapper
        amber coordinate / restart / ncrst file
    params : :class:`parmed.amber.AmberParamSet`
    layermask : list
        list of ambmask strings to select layers (hi,med)
    cm : string
        gaussian Charge + Multiplicity line e.g. "0 1 0 1 0 1"
    link0 : :class:`amb2oniom.link0`
        gaussian link0 commands
    route : string
        gaussian route commands
    linkatoms : list
        list of all link atoms (as strings)
    optflag : list
        list of all gaussian optimization flags (-1 or 0)
    comment : str
        gaussian comment line
    """
    def __init__(self,top,xyz,layermask=["*"],cm="0 1 0 1 0 1",
            comment="amb2oniom",a2g=False,keep_types=False):
        # initialize all Arguments
        self.top = amber.AmberParm(top,xyz=xyz)
        self.params = amber.AmberParameterSet.from_structure(self.top)
        self.link0 = link0()
        route="#P ONIOM(bp86/def2SVP empiricalDispersion=GD3:amber=Softonly) "
        route+="IOp(2/15=3) nosymm "
        route+="geom=connect "
        self.route = route
        self.cm = cm
        self.set_layers(*layermask)
        self.linkatoms = ["" for x in range(self.top.ptr("NATOM"))]
        self.set_optflag("*",symbol1=0)
        self.comment=comment
        # mainly debug options and testing:
        self._a2g = a2g
        self._keep_types = keep_types

    def __repr__(self):
        retstr=["<{} top='{}'; ".format(type(self).__name__,self.top)]
        retstr.append("route='{}'; ".format(self.route))
        retstr.append("link0={}; ".format(self.link0))
        if len(self.layers) > 5:
            layers="['" + "','".join(self.layers[:3]) + "' ... '" \
                    + "','".join(self.layers[-2:]) + "']"
        else: layers="['"+ "','".join(self.layers) + "']"
        retstr.append("layers={}; ".format(layers))
        if len(self.linkatoms) > 5:
            linkatoms="['" + "','".join(self.linkatoms[:3]) + "' ... '" \
                    + "','".join(self.linkatoms[-2:]) + "']"
        else: linkatoms="['"+ "','".join(self.linkatoms) + "']"
        retstr.append("linkatoms={}; ".format(linkatoms))
        if len(self.optflag) > 5:
            optflag="[" + ",".join([str(x) for x in self.optflag[:3]]) \
                + " ... " + ",".join([str(x) for x in self.optflag[-2:]])\
                + "]"
        else: optflag="["+ ",".join([str(x) for x in self.optflag]) + "]"
        retstr.append("optflag={}; ".format(optflag))
        retstr.append("comment='{}'; ".format(self.comment))
        retstr.append("a2g={}; keep_types={}>".format(self._a2g,
            self._keep_types))
        return(''.join(retstr))

    #---------------------------------#

    def find_link_atoms(self,low="L"):
        """ search all link atoms and print bond partners

        Parameters:
        -----------
        low : low level identifier (default: "L")

        """
        print("\nSearching link atoms..")
        #self.linkatoms = ["" for x in range(self.top.ptr("NATOM"))]
        for n,label1 in enumerate(self.layers):
            if label1 == low: continue
            for partner in self.top.atoms[n].bond_partners:
                if self.layers[partner.idx] == low:
                    self.linkatoms[partner.idx] = "H " + str(n+1)
                    print("\n  found linkatom: LAH {}.{}.{}.{} is bonded to LAC {}.{}.{}.{}".format(
                        partner.idx+1,partner.name,partner.element,partner.type,n+1,self.top.atoms[n].name,
                        self.top.atoms[n].element,self.top.atoms[n].type))
                    self._print_partners(n)
        print("\nfound {} link atoms.\n".format(len(
            [x for x in self.linkatoms if x != "" ])))

    def set_layers(self,mask1="*",mask2=""):
        """ set self.layers flags for all atoms

        Parameters:
        -----------
        mask1 : str ; ambmask string
        mask2 : str, optional ; ambmask string

        """
        if mask1:
            layers1 = self.select_atoms(self.top,mask1,symbol1="H",symbol2="L")
        if mask2:
            layers2 = self.select_atoms(self.top,mask2,symbol1="M",symbol2="L")
            for n in range(len(layers1)):
                if layers2[n] == "M" and layers1[n] == "L": layers1[n] = "M"
        self.layers=layers1

    def set_optflag(self,mask="*",symbol1=0,symbol2=-1):
        """ set opflag for gaussian for all atoms

        Parameters:
        -----------
        mask : str ; ambmask string
        symbol1, symbol2 : str
            used to set matching atoms to symbol1 and not matching atoms to
            symbol2

        """
        self.optflag=self.select_atoms(self.top,mask,symbol1,symbol2)

    def strip(self,mask="! *", write_parm=True, parm_suffix="amb2oniom"):
        """ strip atoms from topology

        Parameters:
        -----------
        mask : str ; ambmaskstring

        """
        self.top.strip(mask)
        self.set_layers()
        self.set_optflag()
        self.linkatoms = ["" for x in range(self.top.ptr("NATOM"))]
        if write_parm:
            self.top.write_parm(parm_suffix+".parm7")
            self.top.write_rst7(parm_suffix+".rst7")

    def print_net_charge_per_layer(self):
        """ print the net charge for all layers and the total Charge
        """
        items=sorted(set(self.layers))
        total_charge=0.0
        total_num_electrons=0.0
        for item in items:
            charge=0.0
            num_electrons=0.0
            for idx,val in enumerate(self.layers):
                if val == item:
                    charge+=self.top.atoms[idx].charge
                    num_electrons+=self.top.atoms[idx].atomic_number
            total_charge+=charge
            total_num_electrons+=num_electrons
            print("Netcharge for layer {}: {:7.4f} ({:7d} electrons)".format(item,charge,round(num_electrons)))
        print("Total charge: {:7.4f} ({:7d} electrons)".format(total_charge, round(total_num_electrons)))


    def write_com(self,outf=sys.stdout):
        """ write gaussian com file

        Parameters:
        -----------
        outf : TextIOWrapper (Default: sys.stdout)

        """
        # gaussian link0 commands
        self.link0._write(outf)
        # gaussian route
        outf.write(self.route)
        outf.write("\n\n")
        # gaussian Title
        outf.write(self.comment)
        outf.write("\n\n")
        #charge and multiplicity
        outf.write(self.cm)
        outf.write("\n")
        # write coordinates
        self._write_coordinates(outf)
        outf.write("\n")
        self._write_connectivity(outf)
        outf.write("\n")
        # amber parameters
        # gaussian force field setup
        self._write_magic_line(outf)
        # vdw terms
        self._write_vdw(outf)
        # bonds
        self._write_bonds(outf)
        # angles
        self._write_angles(outf)
        # dihedrals
        self._write_dihedrals(outf)
        # impropers
        # pre parmed v3.0.0: impropers may have wrong order.
        # we try to read the parmed version and decide if we
        # need a fix or not:
        if _check_parmed_version(): self._write_impropers(outf)
        else: self._write_impropers_fix(outf)
        outf.write("\n")

    #---------------------------------#

    def _adjust_atom_types(self,at):
        """ change amber atom types for gaussian com

        Parameters:
        -----------
        at : str or list
            amber atom types to change

        Returns:
        --------
        at : str or list
            modified atom type

        Note:
        -----
        Gaussian shows strange behaviour when original amber atom types are
        used directly, even with "Softonly". Therefore, we change the atom types
        by adding "i" or "j" to gaff or Amber atom types. This is the only
        change in a2g mode. Otherwise, we also remove leading numbers and
        replace potential wild cards.

        """
        if isinstance(at,list):
            for n in range(len(at)):
                at[n]=self._adjust_atom_types(at[n])
        else:
            if not self._a2g:
                at=self._wrap_leading_numbers(at)
                if "*" in at: at=at.replace("*","star")
                if "+" in at: at=at.replace("+","plus")
            if not self._keep_types:
                if at.isupper(): at+="j"
                else: at+="i"
        return(at)

    def _wrap_leading_numbers(self,at):
        """
        wrap leading numbers in atom types such
        as 2C, 3C, ... and return C_2, C_3,...

        Parameters:
        -----------
        at : str; amber atom type
        """
        if at[0].isdigit():
            at=at[1:]+"_"+at[0]
            at=self._wrap_leading_numbers(at)
        return(at)

    def _print_partners(self,idx):
        """ print bond partners for a given atom

        Parameters:
        -----------
        idx : int; atom index
        """
        for partner in self.top.atoms[idx].bond_partners:
            print("   atom {}.{}.{}.{} has a bond to {}.{}.{}.{}".format(
                idx+1,self.top.atoms[idx].name,Element[self.top.atoms[idx].element],self.top.atoms[idx].type,
                partner.idx+1,partner.name,Element[partner.element],partner.type))

    def _write_coordinates(self,outf=sys.stdout):
        """ write com coordinate section

        Parameters:
        -----------
        outf : TextIOWrapper (Default: sys.stdout)
        """
        for n,atom in enumerate(self.top.atoms):
            if self._a2g:
                atom_name="{:s}-{:s}-{:<7.4f}".format(
                    Element[atom.element],self._adjust_atom_types(atom.type),
                    atom.charge)

                outf.write("{:13s}          {:3d}   {:12.7f}{:12.7f}{:12.7f} {:s} {:s}\n".format(
                      atom_name,self.optflag[n],atom.xx,atom.xy,atom.xz,self.layers[n],
                      self.linkatoms[n]))
            else:
                atom_name="{:s}-{:s}-{:<.7f}".format(
                    Element[atom.element],self._adjust_atom_types(atom.type),
                    atom.charge)
                outf.write("{:17s}  {:d} {:12.6f} {:12.6f} {:12.6f} {:s} {:s}\n".format(
                      atom_name,self.optflag[n],atom.xx,atom.xy,atom.xz,self.layers[n],
                      self.linkatoms[n]))

    def _write_connectivity(self,outf=sys.stdout):
        """ write com connectivity section

        Parameters:
        -----------
        outf : TextIOWrapper (Default: sys.stdout)
        """
        for atom in self.top.atoms:
            outf.write("{:6d} ".format(atom.idx+1))
            for partner in atom.bond_partners:
                outf.write("{:6d} 1.0 ".format(partner.idx+1))
            outf.write("\n")

    def _write_magic_line(self,outf=sys.stdout):
        """ write com nonbond line

        Parameters:
        -----------
        outf : TextIOWrapper (Default: sys.stdout)
        """
        if self._a2g:
            outf.write("NonBon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 0.8333333333333 \n")
        else:
            outf.write("NonBon 3 1 0 0 0.0 0.0 -2.0 0.0 0.0 -1.2 \n")

    def _write_vdw(self,outf=sys.stdout):
        """ write com vdw section

        Parameters
        ----------
        outf : TextIOWrapper (Default: sys.stdout)
        """
        for atom,at in self.params.atom_types.items():
            if self._a2g:
                outf.write("VDW {:s}    {:7.4f} {:7.4f}\n".format(
                   self._adjust_atom_types(atom),at.rmin,at.epsilon))
            else:
                outf.write("VDW {:s}    {:7.4f} {:7.4f}\n".format(
                   self._adjust_atom_types(atom),at.rmin,at.epsilon))

    def _write_bonds(self,outf=sys.stdout):
        """ write com bond section

        Parameters
        ----------
        outf : TextIOWrapper (Default: sys.stdout)
        """
        #inspired by parmed's frcmod function
        past=set()
        for bond,bt in self.params.bond_types.items():
             if id(bt) in past : continue
             past.add(id(bt))
             if self._a2g:
                 outf.write("HrmStr1  {:3s} {:3s} {:8.2f} {:6.3f}\n".format(
                     *self._adjust_atom_types(list(bond)),bt.k,bt.req))
             else:
                 outf.write("HrmStr1 {:3s} {:3s} {:8.3f} {:9.4f}\n".format(
                     *self._adjust_atom_types(list(bond)),bt.k,bt.req))

    def _write_angles(self,outf=sys.stdout):
        """ write com angle section

        Parameters
        ----------
        outf : TextIOWrapper (Default: sys.stdout)
        """
        past=set()
        for angle,at in self.params.angle_types.items():
            if id(at) in past: continue
            past.add(id(at))
            # same as a2g, precision ok.
            outf.write("HrmBnd1  {:3s} {:3s} {:3s} {:8.2f} {:9.3f}\n".format(
                     *self._adjust_atom_types(list(angle)),at.k,at.theteq))
        # tip3p water has no angles but a bond between HW-HW
        # we have to add the parameters if water is present!
        if ("HW","HW") in self.params.bond_types.keys() and \
          not ("HW","OW","HW") in self.params.angle_types.keys():
            outf.write("HrmBnd1  {:s} {:s} {:s} {:8.2f} {:9.3f}\n".format(
                 *self._adjust_atom_types(["HW","OW","HW"]),0.0,0.0))
            outf.write("HrmBnd1  {:s} {:s} {:s} {:8.2f} {:9.3f}\n".format(
                 *self._adjust_atom_types(["HW","HW","OW"]),0.0,0.0))

    def _write_dihedrals(self,outf=sys.stdout):
        """ write com dihedral section

        Parameters
        ----------
        outf : TextIOWrapper (Default: sys.stdout)
        """
        past=set()
        for dihedral,d_types in self.params.dihedral_types.items():
            if id(d_types) in past: continue
            past.add(id(d_types))
            outf.write("AmbTrs  {:3s} {:3s} {:3s} {:3s} ".format(
                *self._adjust_atom_types(list(dihedral))))
            # gaussian allows a series of up to 4 dyhedrals
            # we stick to gaussian convention
            i=[1,2,3,4]
            Mag_i=[0,0,0,0]
            PO_i=[0,0,0,0]
            NPaths=1.0
            for dt in d_types:
                try:
                    idx = i.index(dt.per)
                except ValueError():
                    sys.stderr.write("!!Attention, periodicity of {}".format(
                          dt.per))
                    sys.stderr.write("is not supported by gaussian!!\n")
                else:
                    Mag_i[idx] = dt.phi_k
                    PO_i[idx] = int(dt.phase)
            outf.write(" {:3d} {:3d} {:3d} {:3d} {:10.7f} {:10.7f} {:10.7f} {:10.7f} {:3.1f}\n".format(
                  *PO_i,*Mag_i,NPaths))

    def _write_impropers(self,outf=sys.stdout):
        """ write com improper section

        Parameters
        ----------
        outf : TextIOWrapper (Default: sys.stdout)t

        WARNING: atom order of dihedrals is wrong in parmed versions pre 3.X.X!
        the third atom is always the central atom, permutations of a1 a2 a4 are
        all referring to the same improper!
        """
        print("  using parmed 3. \n")
        past=set()
        for (a1,a2,a3,a4),d_type in self.params.improper_periodic_types.items():
            a1,a2,a4 = sorted([a1,a2,a4])
            if (a1,a2,a3,a4) in past: continue
            past.add((a1,a2,a3,a4))
            # we stick to gaussian convention
            Mag=d_type.phi_k
            PO=d_type.phase
            Period=d_type.per
            outf.write("ImpTrs {:3s} {:3s} {:3s} {:3s} {:11.7f} {:7.2f} {:4.1f} \n".format(
                *self._adjust_atom_types([a1,a2,a3,a4]),Mag,PO,Period))

    def _write_impropers_fix(self,outf=sys.stdout):
        """ write com dihedral section

        Parameters
        ----------
        outf : TextIOWrapper (Default: sys.stdout)

        Note: We do not use a AmberParamSet because some impropers
        have wrong atom order! (for example hn-ca-nh-hn became
        hn-nh-hn-ca in test cases (anniline//gaff)
        """
        print("  using parmed 2. workaround..\n")
        past=set()
        for dihedral in self.top.dihedrals:
            # impropers are also not always correctly tagged
            # with the improper flag. therefore we have to
            # cycle through all of them...
            a1 = dihedral.atom1
            a2 = dihedral.atom2
            a3 = dihedral.atom3
            a4 = dihedral.atom4
            # the central atom in amber impropers has to
            # be at the 3rd position. otherwise the atom
            # order is reversed. if that is also not true,
            # there is something fishy with the dihedral
            for atom in [a1,a2,a3,a4]:
                for oatom in [a1,a2,a3,a4]:
                    if oatom is atom: 
                        continue
                    if oatom not in atom.bond_partners:
                        break
                # this will execute if we did not break out of
                # the last for loop -> we found the central atom
                else:
                    at1,at2,at4 = sorted([x.type for x in [a1,a2,a3,a4] if x is not atom])
                    at3 = atom.type
                    if (at1,at2,at3,at4) in past: continue
                    past.add((at1,at2,at3,at4))
                    # by now we have an unique improper with correct atom order.
                    Mag = dihedral.type.phi_k
                    PO = dihedral.type.phase
                    Period = dihedral.type.per
                    outf.write("ImpTrs {:3s} {:3s} {:3s} {:3s} {:11.7f} {:7.2f} {:4.1f} \n".format(
                        *self._adjust_atom_types([at1,at2,at3,at4]),Mag,PO,Period))

    #---------------------------------#

    @staticmethod
    def select_atoms(top,mask,symbol1="H",symbol2="L"):
        """ label atoms in mask

        Parameters:
        -----------
        top : :class:`parmed.amber.AmberParm`
        mask : str ; ambmask string
        symbol1, symbol2 : str, optional

        Returns:
        --------
        selection : list
            list with symbol1 for matching atoms and symbol2 for non-matching
            atoms
        """
        selection=amber.AmberMask(top,mask).Selection()
        for idx in range(len(selection)):
            if selection[idx] == 1:
                selection[idx] = symbol1
            elif selection[idx] == 0:
                selection[idx] = symbol2
        return (selection)


#------------------------------------------------------------------------------#

class link0():
    """ class with link0 commands

    Arguments:
    ----------
    link0 : OrderedDict

    Methods:
    --------
    add(key,value) : add key and value to link0
    deledte(key) : delete key from link0
    _write(outf) : write com link0 section

    """
    def __init__(self):
        self.link0=OrderedDict()
        self.add("MEM","24GB")
        self.add("NPROC","24")

    def __repr__(self):
        retstr=["<{} link0=".format(type(self).__name__)]
        retstr.append("{")
        for k,val in self.link0.items():
            retstr.append("%{} : '{}', ".format(k,val))
        retstr.append("}>")
        return(''.join(retstr))

    def delete(self,a):
        """ delete key from link0

        Parameters
        ----------
        a :: string

        """
        try:
            self.link0.pop(a)
        except KeyError:
            print("KeyError: {}".format(a))

    def add(self,a,b):
        """ add link 0 values

        Parameters:
        -----------
        a,b : str
            key value pair to add to OrderedDict
        """
        self.link0.update({a.upper():b.upper()})

    def _write(self,outf=sys.stdout):
        """ write link0 com selection

        Parameter:
        ----------
        outf : TextIOWrapper (Default: sys.stdout)
        """
        for k,val in self.link0.items():
            outf.write("%{}={}\n".format(k,val))
        outf.write("\n")


#------------------------------------------------------------------------------#

def _check_parmed_version():
    """  try to read version information from parmed

    parmed versions pre 3.0 have (as far as I know) the improper
    bug. therefore we check if the version is at least 3.0

    Returns:
    --------
    retval: boolean; True : version 3.0+
                     False : otherwise
    """
    ver =  pmd_version.split(".")
    # untagged versions have a different format, tagged versions have a
    # string of type "X.X.X". if we don't understand the version,
    # we return false to be save.
    if len(ver[0]) == 1: major = ver[0]
    else: major = ver[1][0]
    if major.isdigit():
        if int(major) >= 3: retval = True
        else: retval = False
    else: retval = False
    return(retval)

def parser(args=None,namespace=None):
    """ argparse Argument parser

    Parameters
    ----------
    args: arguments to parse (Default: None = sys.argv)
    namespace: argparse Namespace object (Default: None)

    Returns
    -------
    args : Namespace object
    """
    parser=argparse.ArgumentParser(
        description='Create ONIOM input from amber top and crd files',
        epilog='report bugs to christian.wick@fau.de')
    grp_files=parser.add_argument_group("Files")
    grp_files.add_argument("-f", metavar="FILE", type=argparse.FileType("r"),
        help="input file")
    grp_files.add_argument("-p", "--parm", type=argparse.FileType("r"),
        help="Amber Topology")
    grp_files.add_argument("-c", "--crd", type=str,
        help="Amber Coordinate/Restart/netcdf file")
    grp_files.add_argument("-o", metavar="COM", type=argparse.FileType("w"),
        default=sys.stdout,help="oniom com file")
    grp_lay=parser.add_argument_group("Layers and active atoms")
    grp_lay.add_argument("-hi", type=str, help="high layer mask")
    grp_lay.add_argument("-med", type=str, help="medium layer mask (ambmask)")
    grp_lay.add_argument("-act", type=str, help="active atoms mask (ambmask)")
    grp_top=parser.add_argument_group("Topology changes")
    grp_top.add_argument("-strip", type=str, help="strip mask from top (ambmask)")
    grp_top.add_argument("-suffix", type=str, default="amb2oniom",
        help="suffix of stripped topology and restart files")
    grp_top.add_argument("-no_parm", action="store_true", default=False,
        help="do not write stripped parameter and restart file after stripping")
    grp_com=parser.add_argument_group("Gaussian input")
    grp_com.add_argument("-cm", type=str, default="0 1 0 1 0 1",
        help="charge and multiplicity line")
    grp_com.add_argument("-route", help="Gaussian route",type=str)
    #parser.add_argument("-i", help="run in interactive mode",
    #      action="store_true",default=False)
    grp_test=parser.add_argument_group("Compatibility and Testing")
    grp_test.add_argument("-a2g", action="store_true", default=False,
        help="reproduce a2g output format")
    grp_test.add_argument("-keep_types", action="store_true", default=False,
        help="keep original atom types. use only for comparisions!")
    args = parser.parse_args(args,namespace)
    return(args)

#------------------------------------------------------------------------------#

if __name__ == "__main__":
    """ start main """

    import argparse

    args = parser()
    print("\nWelcome to amb2oniom.py!\n")
    print("using parmed {} \n".format(pmd_version))
    if args.f:
        # read in input file and add args to args.
        # will overwrite command line arguments.
        content=args.f.readlines()
        file_args=[]
        for line in content:
            temp = line.strip().split(maxsplit=1)
            for val in temp: file_args.append(val)
        args=parser(file_args,args)
    if not args.parm or not args.crd:
        sys.stderr.write("Parm and crd files are mandatory in batch mode!\n\n")
        sys.exit(1)
    top=amb2oniom(args.parm,xyz=args.crd,a2g=args.a2g,keep_types=args.keep_types)
    if args.strip:
        if args.no_parm: top.strip(args.strip,write_parm=False)
        else: top.strip(args.strip,write_parm=True,parm_suffix=args.suffix)
    if args.hi and args.med: top.set_layers(args.hi,args.med)
    elif args.hi: top.set_layers(args.hi)
    elif args.med: top.set_layers(args.med)
    top.print_net_charge_per_layer()
    if args.act: top.set_optflag(args.act)
    top.find_link_atoms()
    if args.med: top.find_link_atoms("M")
    if args.cm: top.cm = args.cm
    if args.route: top.route = args.route
    print("\nWriting com file.. \n \n")
    top.write_com(args.o)
