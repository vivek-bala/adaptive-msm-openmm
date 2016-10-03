import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--inp", help="input filename")
    parser.add_argument("--out", help="output filename")
    args = parser.parse_args()

    inf=open(args.inp,'r')
    outf=open(args.out,'w')
    for line in inf:
        sp=line.split('=')
        if len(sp) == 2:
            key=sp[0].strip().replace('-', '_').lower()
            if key == 'gen_vel':
                outf.write('%s = %s\n'%('gen_vel', 'yes'))
            else:
                outf.write(line)

    # and write our remaining options 
    #for key, value in repl.iteritems():
    #    outf.write('%s = %s\n'%(key, str(value)))

    inf.close()
    outf.close()