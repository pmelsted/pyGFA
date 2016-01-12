import sys,collections

class GFAGraph:

    def __init__(self):
        self.segments = {}
        self.links = {}

def parseGFA(fn):
    if fn == '-':
        f = sys.stdin
    else:
        f = open(fn,'r')

    G = GFAGraph()

    for line in f:
        line = line.strip()
        
        
    if fn != '-':
        f.close()

    return G


if __name__ == '__python__':
    cmds = {}
    v = sys.argv[1:]
    if len(v) < 2:
        print 'Usage: python pyGFA.py <command> <in.GFA> [options]'
    else:        
        cmd = v[0]
        if cmd in cmds:
            G = parseGFA(v[1])
            cmds[cmd](G,v[2:])
        
