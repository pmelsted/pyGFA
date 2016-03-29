import sys,collections

class GFAGraph:

    def __init__(self):
        self.segments = {} # id -> seq
        self.info = {} # id -> info for segment
        self.links = {} # (id,bool) -> {(id,bool)-> (cig,opts)}
        self.paths = {} # id -> (ids,oris,cigs,opts)

    def addSegment(self,name,seq,opts={}):
        self.segments[name] = seq
        self.info[name] = opts
        for o in True,False:
            if (name,o) not in self.links:
                self.links[(name,o)] = {}

    def addPath(self,name,ids,oris,cigs,opts={}):
        assert len(ids) == len(oris) and "Error: orienations and ids don't match"
        assert (cigs==['*'] or len(cigs) == len(ids)-1 or len(cigs) == len(ids)) and "Error incorrect number of cigar strings"
        if cigs == ['*']:
            cigs = ['*' for x in range(len(ids)-1)]
        assert name not in self.paths and "Error: cannot repeat paths in GFA"

        self.paths[name] = (ids,oris,cigs,opts)


    def addLink(self,id1,o1,id2,o2,cig,opts={}):
        if (id1,o1) in self.links:
            if (id2, o2) in self.links[(id1,o1)]:
                oldcig,oldopts = self.links[(id1,o1)][(id2,o2)]
                assert rev_cig(cig) == oldcig and "Repeated cigar string does not match!"
                assert oldopts == opts and "Repeated options don't match"
                return

        for x in (id1,id2):
            if (x,True) not in self.links:
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
        elif line[0] == 'P':
            l = line.split()
            assert len(l) >= 4 and "Error: too few columns in path"
            assert l[0] == 'P' and "Incorrect path specifier"
            name = l[1]
            segnames = l[2].split(',')
            orifn = {'+':True,'-':False}
            oris = list(orifn[x[-1]] for x in segnames)
            segs = list(x[:-1] for x in segnames)
            cigs = l[3].split(',')
            opts = parseOpts(l[4:])
            G.addPath(name,segs,oris,cigs,opts)

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
            if o != oo and keyfn(x) > keyfn(xx):
                continue
            elif o == False and oo == False:
                continue

            cig,opts = G.links[(x,o)][(xx,oo)]
            f.write('L\t%s\t%s\t%s\t%s\t%s%s\n'%(x,ori[o],xx,ori[oo],cig,'\t'+printOpts(opts) if not len(opts)==0 else ''))

    paths = list(G.paths.keys())
    allpathint = True
    try:
        for x in paths:
                int(x)
    except ValueError:
        allpathint = False
    pathkeyfn = id
    if allpathint:
        pathkeyfn = int
    paths.sort(key=pathkeyfn)

    for path in paths:
        segs,oris,cigs,opts = G.paths[path]
        cigstr = None
        if len(cigs) == len(segs)-1 and all(x=='*' for x in cigs):
            cigstr = '*'
        else:
            cigstr = ','.join(cigs)

        pstr = ','.join(a+ori[b] for a,b in zip(segs,oris))
        f.write('P\t%s\t%s\t%s%s\n'%(path,pstr,cigstr,'\t'+printOpts(opts) if not len(opts)==0 else ''))

def connected(G,opts):
    anchor = opts[0]
    assert anchor in G.segments and "Error: specified id = %s is not in the graph"%(anchor)

    subG = GFAGraph()
    stack = [anchor]
    while stack:
        x = stack.pop()
        if x not in subG.segments:
            subG.addSegment(x,G.segments[x],G.info[x])
            for o in True,False:
                for xx,oo in G.links[(x,o)]:
                    if xx in subG.segments:
                        ## both x and xx are in the graph, add the links
                        cig,opts = G.links[(x,o)][(xx,oo)]
                        subG.addLink(x,o,xx,oo,cig,opts)
                    else:
                        ## xx not yet in graph, put it on the list
                        stack.append(xx)
    return subG

def subgraph(G,opts):
    ids = set(opts[0].split(','))
    for x in ids:
        assert x in G.segments and "Error, segment not found in graph"

    subG = GFAGraph()
    # add segments
    for x in ids:
        subG.addSegment(x,G.segments[x], G.info[x])

    # add links
    for x in ids:
        for o in True,False:
            if (x,o) in G.links:
                for xx,oo in G.links[(x,o)]:
                    if xx in ids:
                        cig,opts = G.links[(x,o)][(xx,oo)]
                        subG.addLink(x,o,xx,oo,cig,opts)



    # add paths
    # TODO: handle logic for subgraphing circular paths
    for path in G.paths:
        segs,oris,cigs,opts = G.paths[path]
        circular = len(segs) == len(cigs)
        i = 1
        a,b = 0,0
        firstPath = None
        lastPath = None
        while a < len(segs):
            b = a + 1
            if segs[a] in ids:
                while b < len(segs) and segs[b] in ids:
                    b += 1
                # now segs[a:b] are all in the subgraph
                # copy the path
                if circular and a==0 and b < len(segs):
                    firstPath = b
                else:
                    if a==0 and b == len(segs):
                        subG.addPath(path,segs,oris,cigs,opts)
                        break # special case copy the original path
                    elif b-a > 1:
                        if circular and b == len(segs) and firstPath is not None:
                            #special case, original path is circular
                            lastPath = a
                        else:
                            cigB = b if not circular else b+1
                            subG.addPath("%s_%s"%(path,i), segs[a:b],oris[a:b],cigs[a:cigB],opts)
                            i += 1

            if firstPath is not None:
                if lastPath is not None:
                    subG.addPath("%s_%s"%(path,i), segs[lastPath:] + segs[:firstPath], oris[lastPath:] + oris[:firstPath], cigs[lastPath:] + cigs[:firstPath] ,opts)
                    i += 1
                else:
                    subG.addPath("%s_%s"%(path,i), segs[:firstPath], oris[:firstPath], cigs[:firstPath], opts)
                    i += 1
            elif lastPath is not None:
                subG.addPath("%s_%s"%(path,i), segs[lastPath:], oris[lastPath:], cigs[lastPath:], opts)
            a = b # next to check

    return subG



if __name__ == '__main__':
    cmds = {'print':printGFA,
            'connected': lambda G,opts: printGFA(connected(G,opts),opts),
            'subgraph': lambda G,opts : printGFA(subgraph(G,opts),opts)

            }
    v = sys.argv[1:]
    if len(v) < 2:
        print('Usage: python pyGFA.py <command> <in.GFA> [options]')
    else:
        cmd = v[0]
        if cmd in cmds:
            G = parseGFA(v[1])
            cmds[cmd](G,v[2:])
