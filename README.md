Introduction:
We developed a supporting code of Gaussian 09 to study the charge transfer mechanism between two fragments of a molecule. This code follows a simple yet powerful approach for analyzing a chemical bond in terms of the electron density rearrangement occurring upon formation of an adduct from two constituting fragments. This approach relies on the charge displacement (CD) function, measuring the electron charge displacement along a chosen axis occurring upon chemical bond formation. This method has been successfully used for quantifying the charge transfer (CT) component of the interaction in weakly interacting adducts and more recently in the characterization of electronic transitions. 
We use Fortran 90 programming language to produce these code.

Explanation of the Code charge_transfar.f90:

module commondat:  module commondat is the location of the common data used in the program.  Common data means all the data (can be a constant value, a name, or a variable value which updated in different subroutines) used in one or more subroutines. Common data can be generated in any subroutine and instantly update the value in the Module commondata.  In this program the covalent radius of the different atoms corresponding to their names are listed in this module. 

subroutine prepname: Here program prepare the various file-names required for the different subroutines. 

subroutine shiftmol: this subroutine shift the centre point of the donor-acceptor connected line of the molecule to the origin. 

subroutine rottoz: this subroutine rotate the molecule towards the z axis. The rotation is done by operating a rotational matrix operator on each of the atomic vectors. The rotational matrix which rotate a point around a line. Here the line is perpendicular to both the atomic axis and z axis. That means this line is a perpendicular vector of the plane containing atomic axis and z axis. So I take the cross product of them to get it. Now apply the rotational matrix operator on the each atom. 

Subroutine output1: here the atoms of the molecule is shorted out as monomer1 and monomer2. There are three options to select the atoms. They are “distance”, “radius” and “user”
In the first one, Atoms are selected according to their distances from acceptor and donor. 
