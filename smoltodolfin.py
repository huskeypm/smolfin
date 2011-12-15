def write_dolfin_files(filename, mesh, u,markers):
    from dolfin import File

    File(filename+"_mesh.xml.gz") << mesh
    File(filename+"_values.xml.gz") << u
    File(filename+"_markers.xml.gz") << markers

def do_read_write(mcsffilename,apbsfilename):
    import numpy
    from apbstodolfin import *
    from mcsftodolfin import *

    # read apbs
    acoordinates, acells, avalues = read_apbs_file(apbsfilename)
    #read gamer
    mcoordinates, mcells, mmarkers= read_mcsf_file(mcsffilename)


    # interpolate apbs values onto mcoordinates grid
    from scipy.interpolate import griddata
    mvalues = griddata(acoordinates, avalues,(mcoordinates),method='linear')


    # create mesh, etc
    mesh = generate_dolfin_mesh(mcoordinates, mcells)
    values = generate_dolfin_function(mesh, mvalues)


    # save
    write_dolfin_files(mcsffilename.replace(".m", ""), mesh, values,mmarkers)


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        raise RuntimeError("expected an 1) mcsf and 2) apbs file")
    mfilename = sys.argv[1]
    afilename = sys.argv[2]

    # ~/localTemp/NBCR/smol/apbs_fe/potential-0.dx
    #apbsfilename = "example/molecule/potential-0.dx"
    # ~/localTemp/NBCR/smol/gamer/p.pqr.output.out.m
    #mcsffilename = "example/molecule/p.pqr.output.out.m"

    do_read_write(mcsffilename,apbsfilename)

