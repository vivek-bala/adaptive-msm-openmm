import re
import sys
import fileinput

if __name__ == '__main__':
    
    ret1 = sys.stdin.readlines()
    f = ''
    for item in ret1:
        f += item
    step = re.compile('^Step\s*([0-9.]*)\s*([0-9.]*)', re.MULTILINE)

    match = step.search(f)
    frames = int(match.group(1))
    dtstring = match.group(2)

    if dtstring.strip() == "":
        dt=0
    else:
        dt=float(match.group(2))

    ns=(frames-1)*dt/1000

    with open('traj_info.txt','w') as f:
        f.write('ns = {0}\n'.format(ns))
        f.write('dt = {0}\n'.format(dt))
        f.write('frames = {0}\n'.format(frames))
