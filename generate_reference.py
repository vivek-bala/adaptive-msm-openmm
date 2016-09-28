import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--inp", help="input filename")
    parser.add_argument("--out", help="output filename")
    args = parser.parse_args()

    inf=open(args.inp,'r')
    outf=open(args.out,'w')
    for line in inf:
        # exclude any lines that could lead to errors
        if line.startswith( ('TITLE', 'MODEL', 'ATOM', 'TER', 'ENDMOL') ):
            outf.write(line)
    outf.close()
    inf.close()
