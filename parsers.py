import sys,collections,re,os.path


def parse_spades(v):
    prefix,output,k = v[0:3]
    k = int(k)


    of = sys.stdout
    if output != '-':
        ofn = os.path.join(prefix,output)
        of = open(ofn,'w')
    of.write('H\tVN:Z:1.0\n')

    node_re = re.compile('NODE_(\d+)_length_(\d+)_cov_(.*?)_ID_(\d+)')
    edge_re = re.compile("EDGE_(\d+)_length_(\d+)_cov_(\d+(.\d+)*)(')?")
    ori_re = re.compile('(\d+)([+-])')

    def parse_e(s):
        x = edge_re.match(s)
        if x is None or len(x.groups()) != 5:
            return None
        else:
            return x.group(1),True if x.group(5) is None else False




    # get all edges
    E = {} # (id,ori) -> [(id,ori)]
    fn = os.path.join(prefix,'assembly_graph.fastg')
    with open(fn) as f:
        for line in f:
            line = line.strip()
            if len(line) > 0 and line[0] == '>' and line[-1] == ';':
                l = line[1:-1].split(':')
                e_from = parse_e(l[0])
                E[e_from] = []
                if len(l) > 1:
                    l2 = l[1].split(',')
                    for x2 in l2:
                        e_to = parse_e(x2)
                        E[e_from].append(e_to)

    fn = os.path.join(prefix,'contigs.fasta')
    names = []
    name,length,cov,seq = None,None,None,[]
    with open(fn) as f:
        for line in f:
            line = line.strip()
            if line[0] == '>':
                if name is not None:
                    of.write('\t'.join(['S',name,''.join(seq),'LN:i:'+length,'KC:i:'+str(int((int(length)-k+1)*cov)),'\n']))
                x = node_re.match(line[1:])
                name,length,cov = x.group(1,2,3)
                cov = float(cov)
                names.append(name)
                seq = []
            else:
                seq.append(line)
        if name is not None:
            of.write('\t'.join(['S',name,''.join(seq),'LN:i:'+length,'KC:i:'+str(int((int(length)-k+1)*cov)),'\n']))

    toends = {} # maps (node,ori) -> (edge,ori)
    fromends = {} # maps (edge,ori) -> (node,ori)
    fn = os.path.join(prefix, 'contigs.paths')
    name,ori = None,None
    with open(fn) as f:
        for line in f:
            line = line.strip()
            if line[:4] == 'NODE':
                ori = True
                if line[-1] == "'":
                    ori = False
                    line = line[:-1]
                name = node_re.match(line).group(1)
            else:
                if name is None:
                    continue
                if line[-1] == ';':
                    line = line[:-1]
                x = ori_re.match(line.split(',')[0])
                eid,eo = x.group(1),True if x.group(2) == '+' else False
                toends[(name,not ori)] = (eid, not eo)
                fromends[(eid, not eo)] = (name,not ori)
                name,ori = None,None

    omap = {True:'+',False:'-'}

    for name in names:
        for no in True,False:
            if (name,no) in toends:
                eid,eo = toends[(name,no)]
                if (eid,eo) in E:
                    for eeid,eeo in E[(eid,eo)]:
                        if (eeid,not eeo) in fromends:
                            nname,nno = fromends[(eeid,not eeo)]
                            nno = not nno
                            of.write('\t'.join(['L',name,omap[no],nname,omap[nno],str(k)+'M\n']))


    if of != sys.stdout:
        of.close()

if __name__ == "__main__":

    parsers = {'spades': (parse_spades,"""
spades:
    python parsers.py spades PREFIX OUTPUT K

      PREFIX = directory where spades stores output
      OUTPUT = GFA file to write
      K = largest k used for assembly, look for directory name k55 in output directory
    """)}
    v = sys.argv[1:]
    if len(v) < 2:
        print('Usage: python parsers.py <type> PREFIX OUTPUT [options]')
        print('')
        print('       available parsers:')
        for p in parsers:
            print('')
            print(parsers[p][1])
    else:
        name = v[0]
        if name in parsers:
            parsers[name][0](v[1:])
