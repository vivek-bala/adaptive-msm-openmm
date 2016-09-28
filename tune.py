import argparse

def tune(confFile, tprFile):
    """Set max. run based on configuration file."""
    # TODO: fix this. For now, only count the number of particles and
    # the system size.
    # read the system size 
    inf=open(confFile, 'r')
    i=0
    for line in inf:
        lastsplit=line.split()
        if i==1:
            N=int(line)
        i+=1
    sx=float(lastsplit[0])
    sy=float(lastsplit[0])
    sz=float(lastsplit[0])
    # as a rough estimate, the max. number of cells is 1 per rounded nm
    mincellsize=1.2

    print sx, sy, sz
    '''
    Nsize = int(sx/mincellsize)*int(sy/mincellsize)*int(sz/mincellsize)
    # and the max. number of processors should be N/250
    NN = int(N/250)
    if Nmax is None:
        # the max number of processors to use should be the minimum of these
        Nmax = min(Nsize, NN)
    Nmax = max(1, Nmax)

    while True:
        # make sure we return a sane number:
        # It's either 4 or smaller, 6, or has at least 3 prime factors. 
        if (Nmax < 5 or Nmax == 6 or
            (Nmax < 32 and len(primefactors(Nmax)) > 2) or
            (Nmax < 32 and len(primefactors(Nmax)) > 3) ):
            canRun, stdo=tryRun(tprFile, testRunDir, Nmax)
            if canRun:
                break
        Nmax -= 1
        if Nmax < 1:
            raise GromacsError("Can't run simulation: %s"%stdo)
    '''


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--conf", help="configuration filename")
    parser.add_argument("--tpr", help="tpr filename")
    args = parser.parse_args()

    tune(args.conf, args.tpr)
