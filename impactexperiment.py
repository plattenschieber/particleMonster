#!/usr/bin/env python
# vim: et sts=4 ts=4 sw=4 ai

import sys, getopt, csv, os, math, re, tempfile, shutil, exceptions

def usage():
    print "Usage:"
    print sys.argv[0], " experimentname\n"
    print "Options:\n\t-s,--sigma\tsigma [=1]"
    print "\t-e,--epsilon\tepsilon [=1]"
    print "\t-v,--velocity\tvelocity of the object1 [=-10]"
    print "\t-d,--deltat\tdelta t [=0.001]"
    print "\t--end\tend_t[=3]"
    print "\t-n1 i   \tobject1 of i*i particles [i=10]"
    print "\t-n2 j   \tobject2 of i*j particles [j=90]"
    print "\t-m1     \tmass of particles in object1 [=1]"
    print "\t-m2     \tmass of particles in object2 [=1]"

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "he:v:d:s:", ["help","epsilon=","eps=","deltat=","delta_t=","n1=","n2=","m1=","m2=","end="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    if len(args) != 1:
        print "error: no experiment name given", args
        usage()
        sys.exit(2)
    # default values
    velocity = 10
    sigma = 1
    epsilon = 1
    delta_t = 0.001
    t_end = 3
    n1 = 10
    n2 = 90
    m1 = m2 = 1
    for o, a in opts:
        print o
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-s", "--sigma"):
            sigma = float(a)
        if o in ("-e", "--epsilon", "--eps"):
            epsilon = float(a)
        if o in ("-d", "--deltat", "--delta_t"):
            delta_t = float(a)
        if o == "--n1":
            n1 = int(a)
        if o == "--n2":
            n2 = int(a)
        if o == "--m1":
            m1 = float(a)
        if o == "--m2":
            m2 = float(a)
        if o == "--end":
            t_end = float(a)

    # first argument is the experiment name
    name = args[0]

    # calc derived values
    rcut = sigma*2.5
    gridstep = pow(2.,1./6.)*sigma
    length_object1 = n1 * gridstep
    length_object2 = n2 * gridstep
    lx = math.ceil(length_object2 * 1.5)
    ly = math.ceil(length_object1 * 5)
    lz = rcut*2
    z = rcut

    # write parameter file
    f = open(name+'.parameter', 'w+')
    f.write("name %s\n" % name)
    f.write("delta_t %g\n" % delta_t)
    f.write("t_end %g\n" % t_end)
    f.write("epsilon %g\n" % epsilon)
    f.write("sigma %g\n" % sigma)
    f.write("cell_r_cut %g\n" % rcut)
    f.write("length %g %g %g\n" % (lx,ly,lz) )
    f.write("upper_border leaving leaving leaving\n")
    f.write("lower_border leaving leaving leaving\n")
    f.close()
    print name+'.parameter done.'
    # write particle file
    count = 0
    # open particle file
    f = open(name+'.particle', 'w+')
    # write object1
    for i in range(0,n1):
        x = (lx/2.)-(length_object1/2.) + i * gridstep;
        for j in range(0,n1):
            y = (ly/4.)-(length_object1) + j * gridstep;
            f.write("%i %g %g %g %g %g %g %g\n" % (count,m1,x,y,z,0,velocity,0))
            count+=1
    print name+'.particle done.'
    # write object2
    for i in range(0,n2):
        x = (lx/2.)-(length_object2/2.) + i * gridstep;
        for j in range(0,n1):
            y = 2*(ly/4.)-(length_object1) + j * gridstep;
            f.write("%i %g %g %g %g %g %g %g\n" % (count,m1,x,y,z,0,0,0))
            count+=1
    # close particle file
    f.close()

if __name__ == "__main__":
    main(sys.argv[1:])
