import sys,collections

class GFAGraph:

    def __init__(self):
        self.segments = {} # id -> seq
        self.info = {} # id -> info for segment
        self.links = {} # (id,bool) -> {(id,bool)-> (cig,opts)}
        self.paths = {} # ?!?

    def addSegment(self,name,seq,opts={}):
        self.segments[name] = seq
        self.info[name] = opts

    def addLink(self,id1,o1,id2,o2,cig,opts={}):
        if (id1,o1) in self.links:
            if (id2, o2) in self.links[(id1,o1)]:
                oldcig,oldopts = self.links[(id1,o1)][(id2,o2)]
                print(oldcig,oldopts,cig,opts)
                assert rev_cig(cig) == oldcig and "Repeated cigar string does not match!"
                assert oldopts == opts and "Repeated options don't match"
                return

        for x in (id1,id2):
            if x not in self.links:
                self.links[(x,True)] = {}
                self.links[(x,False)] = {}

        self.links[(id1,o1)][(id2,o2)] = (cig,opts)
        self.links[(id2,not o2)][(id1, not o1)] = (rev_cig(cig),opts)



def rev_cig(cig):
    return cig # not implemented yet


__constr = {'A':str,'Z':str,'f':float,'i':int,'H':lambda val: [int(x,16) for x in val]}

def parseOpts(l):
    d = {}
    for x in l:
        t = x.split(':')
        assert len(t) == 3, "Incorrect format for option: " + x
        assert t[1] in __constr, "Unknown format: " + ':'.join(t)
        d[t[0]] = __constr[t[1]](t[2])

    return d

def printOpts(d):
    r = []
    for key in d:
        val = d[key]
        typ = 'Z'
        if type(val) == float:
            typ='f'
        elif type(val) == int:
            typ='i'
        r.append(':'.join([key,typ,str(val)]))
    return ' '.join(r)

def parseGFA(fn,skipSeq=False):
    f = None
    if fn == '-':
        f = sys.stdin
    else:
        f = open(fn,'r')

    G = GFAGraph()

    for line in f:
        line = line.strip()

        if len(line) == 0:
            pass
        elif line[0] == 'H':
            pass
        elif line[0] == 'S':
            l = line.split()
            assert len(l) >= 3, "Missing data in line " + line
            assert l[0] == 'S', "Incorrect line specifier " + line
            name = l[1]
            assert name not in G.segments, "Segments cannot be repeated in GFA"
            seq = l[2] if not skipSeq else '*'
            opts = parseOpts(l[3:])
            G.addSegment(name,seq,opts)
        elif line[0] == 'L':
            l = line.split()
            assert len(l) >= 6, "Missing data in line " + line
            assert l[0] == 'L', "Incorrect line specifier " + line
            id1,id2 = l[1],l[3]
            o1,o2 = None,None
            if l[2] == '+':
                o1 = True
            elif l[2] == '-':
                o1 = False
            else:
                assert "Incorrect direction for line " + line

            if l[4] == '+':
                o2 = True
            elif l[4] == '-':
                o2 = False
            else:
                assert "Incorrect direction for line " + line

            cigar = l[5]
            opts = parseOpts(l[6:])

            G.addLink(id1,o1,id2,o2,cigar,opts)

    if fn != '-':
        f.close()

    return G


def printGFA(G,opt):
    # ignores opts, prints to stdout
    f = sys.stdout

    # header
    f.write('H\tVN:Z:1.0\n')

    ids = list(G.segments.keys())
    allint = True
    try:
        for x in ids:
            int(x)
    except ValueError:
        allint = False

    keyfn = id
    if allint:
        keyfn = int

    ids.sort(key=keyfn)

    for x in ids:
        f.write('S\t%s\t%s%s\n'%(x,G.segments[x],'\t'+printOpts(G.info[x]) if x in G.info else ''))


    lkeyfn = lambda x: (keyfn(x[0]),x[1])
    links = list(G.links.keys())
    links.sort(key=lkeyfn)

    ori = {True:'+',False:'-'}
    for x,o in links:
        for xx,oo in G.links[(x,o)]:
            cig,opts = G.links[(x,o)][(xx,oo)]
            f.write('L\t%s\t%s\t%s\t%s\t%s%s\n'%(x,ori[o],xx,ori[oo],cig,'\n'+printOpts(opts) if not len(opt)==0 else ''))


if __name__ == '__main__':
    cmds = {'print':printGFA}
    v = sys.argv[1:]
    if len(v) < 2:
        print('Usage: python pyGFA.py <command> <in.GFA> [options]')
    else:
        cmd = v[0]
        if cmd in cmds:
            G = parseGFA(v[1])
            cmds[cmd](G,v[2:])
