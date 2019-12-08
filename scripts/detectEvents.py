#!/usr/bin/python3

import sys, os, itertools, operator
import argparse
import math

from Bio import SeqIO
import gffutils

from BitVector import BitVector

# Computes edit distance between two strings
def editDistance(s1, s2):
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]


def levenshtein(seq1, seq2):
    size_x = len(seq1) + 1
    size_y = len(seq2) + 1
    #matrix = np.zeros ((size_x, size_y))
    matrix = [[0 for y in range(size_y)] for x in range(size_x)]
    for x in range(size_x):
        matrix[x][0] = x
    for y in range(size_y):
        matrix[0][y] = y

    for x in range(1, size_x):
        for y in range(1, size_y):
            if seq1[x-1] == seq2[y-1]:
                matrix[x][y] = min(
                    matrix[x-1][y] + 1,
                    matrix[x-1][y-1],
                    matrix[x][y-1] + 1
                )
            else:
                matrix[x][y] = min(
                    matrix[x-1][y] + 1,
                    matrix[x-1][y-1] + 1,
                    matrix[x][y-1] + 1
                )
    return (matrix[size_x - 1][size_y - 1])

# L -> (l0,l1), (l1,l2), (l2, l3), ...
def pairwise(L):
    L0, L1 = itertools.tee(L)
    next(L1, None)
    return zip(L0,L1)

# Opens gtf (gffutils wrapper)
def openGTF(gtfPath):
    try:
        gtf = gffutils.FeatureDB("{}.db".format(gtfPath),
                                 keep_order=True)
    except ValueError:
        gtf = gffutils.create_db(gtfPath,
                                 dbfn="{}.db".format(gtfPath),
                                 force=True, keep_order=True,
                                 disable_infer_genes=True,
                                 disable_infer_transcripts=True,
                                 merge_strategy='merge',
                                 sort_attribute_values=True)
        gtf = gffutils.FeatureDB("{}.db".format(gtfPath), keep_order=True)
    return gtf

# Extracts strand, transcripts and introns from gtf
def extractFromGTF(gtf):
    strand = "+"
    introns = set()
    transcripts = {}
    for g in gtf.features_of_type('gene'):
        strand = g.strand
        for tr in gtf.children(g, featuretype='transcript', order_by='start'):
            exons = list(gtf.children(tr, featuretype='exon', order_by='start'))

            transcript = [(ex.start, ex.end) for ex in exons]
            transcripts.update({tr.id:transcript})

            introns_ = set(zip([ex.end+1 for ex in exons[:-1]], [ex.start-1 for ex in exons[1:]]))
            introns = introns | introns_
    return strand, transcripts, introns

# TODO: extract directly from FASTA + GTF
# Extracts text and exon positions from "index" file
def extractFromInfoFile(infoPath):
    lines = open(infoPath).readlines()
    text = lines[1].strip("\n")
    exons = [(int(p[0]), int(p[1])) for p in [pos.split(",") for pos in lines[4].strip("\n").split()]]
    return text, exons

# Extracts elements from one line of the spliced graph-alignments file
def readLine(line):
    # 0: strand
    # 1: ID
    # 2: errors
    # 3 to -1: mems
    # -1: read
    line = line.strip("\n").strip(" ").split(" ")
    strand = line[0]
    readID = line[1]
    err = int(line[2])
    mems = line[3:-1]
    read = line[-1]
    return strand, readID, err, mems, read

# Removes annotated introns
def filterAnnotated(newIntrons, annIntrons):
    introns = {}
    annFoundIntrons = {}
    for (p1,p2),w in newIntrons.items():
        if (p1,p2) not in annIntrons:
            introns.update({(p1,p2):w})
        else:
            annFoundIntrons.update({(p1,p2):w})
    return introns, annFoundIntrons

# Filters new introns that are not sufficiently covered
def filterLowCovered(introns, tresh):
    filtIntrons = {}
    for (p1,p2),w in introns.items():
        if w >= tresh:
            filtIntrons.update({(p1,p2):w})
        #else:
        #    print("# W {} {}".format(p1,p2))
    return filtIntrons

# Reconciliate a given intron (start/end) position with respect to the input pattern
def reconciliate(pos, ref, patt, isStart):
    off = 1
    MaxOff = 3 # TODO: this one could be a parameter
    while off<=MaxOff:
        newPos = pos-off
        if isStart:
            newPatt = str(ref[newPos-1:newPos+1].seq)
        else:
            newPatt = str(ref[newPos-2:newPos].seq)
        if newPatt == patt:
            return True,newPos

        newPos = pos+off
        if isStart:
            newPatt = str(ref[newPos-1:newPos+1].seq)
        else:
            newPatt = str(ref[newPos-2:newPos].seq)
        if newPatt == patt:
            return True,newPos
        off+=1
    return False,-1

# Cleans introns basing on canonical patterns (reconciliation)
def reconciliateIntrons(introns, ref, strand):
    recIntrons = {}
    for (p1,p2),w in introns.items():
        intronPref = str(ref[p1-1:p1+1].seq)
        intronSuff = str(ref[p2-2:p2].seq)
        if strand == '+':
            canIP = {"GT":["AG"], "GC":["AG"]}
            canIPrev = {"AG":["GT", "GC"]}
        else:
            canIP = {"CT":["AC", "GC"]}
            canIPrev = {"AC":["CT"], "GC":["CT"]}
        if intronPref in canIP.keys() and intronSuff in canIP[intronPref]:
            key = (p1,p2)
            recIntrons[key] = recIntrons[key]+w if key in recIntrons else w
        elif intronPref in canIP.keys():
            newPos = float('inf')
            for acceptedIntronicPattern in canIP[intronPref]:
                found,pos = reconciliate(p2, ref, acceptedIntronicPattern, False)
                if found:
                    if abs(pos-p2) < abs(newPos-p2):
                        newPos = pos
            if newPos == float('inf'):
                newPos = p2
            key = (p1,newPos)
            recIntrons[key] = recIntrons[key]+w if key in recIntrons else w
        elif intronSuff in canIPrev:
            newPos = float('inf')
            for acceptedIntronicPattern in canIPrev[intronSuff]:
                found,pos = reconciliate(p1, ref, acceptedIntronicPattern, True)
                if found:
                    if abs(pos-p1) < abs(newPos-p1):
                        newPos = pos
            if newPos == float('inf'):
                newPos = p1
            key = (newPos,p2)
            recIntrons[key] = recIntrons[key]+w if key in recIntrons else w
        else:
            off = 1
            MaxOff = 3
            while off <= MaxOff:
                newP1, newP2 = p1-1, p2-1
                intronPref = str(ref[newP1-1:newP1+1].seq)
                intronSuff = str(ref[newP2-2:newP2].seq)
                if intronPref in canIP.keys() and intronSuff in canIP[intronPref]:
                    key = (newP1,newP2)
                    recIntrons[key] = recIntrons[key]+w if key in recIntrons else w
                    break

                newP1, newP2 = p1+1, p2+1
                intronPref = str(ref[newP1-1:newP1+1].seq)
                intronSuff = str(ref[newP2-2:newP2].seq)
                if intronPref in canIP.keys() and intronSuff in canIP[intronPref]:
                    key = (newP1,newP2)
                    recIntrons[key] = recIntrons[key]+w if key in recIntrons else w
                    break
                off+=1
    return recIntrons

# Returns all the exons starting and ending close to a given position
def getExonsCloseTo(Es,p):
    exsS = set()
    exsE = set()
    for (s,e) in Es:
        if 0 <= s-p <= 100-30: # TODO: set 100-30 to readLen - 2*L
            exsS.add((s,e))
        elif 0 <= p-e <= 100-30: # TODO: set 100-30 to readLen - 2*L
            exsE.add((s,e))
    return (list(exsS),list(exsE))

# Returns all the exons containing a position
def getExonsContaining(Es,p):
    exs = set()
    for (s,e) in Es:
        if s <= p <= e:
            exs.add((s,e))
    return list(exs)

# Returns True if there exists an intron ending at the given position, False otherwise
def existsIntronEndingAt(introns, p):
    for (s,e) in introns:
        if e == p:
            return True
    return False

# Returns True if there exists an intron starting at the given position, False otherwise
def existsIntronStartingAt(introns, p):
    for (s,e) in introns:
        if s == p:
            return True
    return False

# Returns all the (possibly overlapping) introns that follow the given intron
def getSuccIntrons(Introns, I):
    Introns.sort()
    i = Introns.index(I)

    Succ = set()
    if i<len(Introns)-1:
        found = False
        while not found and i<len(Introns)-1:
            i+=1
            (s,e) = Introns[i]
            if I[1]<s:
                found = True
                Succ.add((s, e))

        flag = True
        while flag and i<len(Introns)-1:
            i+=1
            (s_,e_) = Introns[i]
            if s<=s_<=e:
                Succ.add((s_, e_))
            else:
                flag = False
    return list(Succ)

# Returns all the (possibly overlapping) introns that precede the given intron
def getPrecIntrons(Introns, I):
    Introns.sort(key=operator.itemgetter(1))
    i = Introns.index(I)

    Prec = set()
    if i>0:
        found = False
        while not found and i>0:
            i-=1
            (s,e) = Introns[i]
            if e<I[0]:
                found = True
                Prec.add((s, e))

        flag = True
        while flag and i>0:
            i-=1
            (s_,e_) = Introns[i]
            if s<=e_<=e:
                Prec.add((s_, e_))
            else:
                flag = False
    return list(Prec)

# Extracts events from introns
def checkNewIntrons(newIntrons, allIntrons, strand, transcripts):
    allIntrons = list(allIntrons)
    allIntrons.sort()
    events = {'ES': {}, 'A3': {}, 'A5': {}, 'IR': {}}
    for (p1,p2),w in newIntrons.items():
        # Getting introns preceding and following the considered intron
        precIntrons = getPrecIntrons(allIntrons, (p1,p2))
        succIntrons = getSuccIntrons(allIntrons, (p1,p2))

        # For each transcript...
        for trID,exons in transcripts.items():
            tranSt,tranEnd = exons[0][0], exons[-1][1]
            intronsEnds = [ex[0]-1 for ex in exons]
            intronsStarts = [ex[1]+1 for ex in exons]

            # Checking exon skippings
            if p1 in intronsStarts and p2 in intronsEnds:
                i1 = intronsStarts.index(p1)
                i2 = intronsEnds.index(p2)
                if i1 != i2-1:
                    key = (p1,p2,w)
                    if key not in events['ES']:
                        events['ES'][key] = []
                    events['ES'][key].append(trID)

            # Checking intron retentions
            for (s,e) in exons:
                if s < p1 < p2 < e:
                    if (len(precIntrons) == 0 or s-1 in [i[1] for i in precIntrons]) and (len(succIntrons) == 0 or e+1 in [i[0] for i in succIntrons]):
                        key = (p1,p2,w)
                        if key not in events['IR']:
                            events['IR'][key] = []
                        events['IR'][key].append(trID)

            # Checking alternative splice sites
            if p1 in intronsStarts and p1 != tranEnd and p2 not in intronsEnds and p2 < tranEnd:
                exonsContaining = getExonsContaining(exons,p2)
                exonsCloseS,_ = getExonsCloseTo(exons,p2)
                if len(exonsContaining)>0 or len(exonsCloseS)>0:
                    if len(succIntrons) == 0 or any([i[0] in intronsStarts for i in succIntrons]) or tranEnd in [e[1] for e in exonsCloseS + exonsContaining]:
                        if strand == '+':
                            t = 'A3'
                        else:
                            t = 'A5'
                        key = (p1,p2,w)
                        if key not in events[t]:
                            events[t][key] = []
                        events[t][key].append(trID)

            if p1 not in intronsStarts and p2 in intronsEnds and p2 != tranSt and p1 > tranSt:
                exonsContaining = getExonsContaining(exons,p1)
                _,exonsCloseE = getExonsCloseTo(exons,p1)
                if len(exonsContaining)>0 or len(exonsCloseE)>0:
                    if len(precIntrons) == 0 or any([i[1] in intronsEnds for i in precIntrons]) or tranSt in [e[0] for e in exonsCloseE + exonsContaining]:
                        if strand == '+':
                            t = 'A5'
                        else:
                            t = 'A3'
                        key = (p1,p2,w)
                        if key not in events[t]:
                            events[t][key] = []
                        events[t][key].append(trID)
    return events

# Printing events (TODO: they can be printed when found, but maybe the dict could be useful for some analysis)
def printEvents(events, outPath):
    out = open(outPath, 'w')
    out.write("Type,Start,End,Support,Transcripts\n")
    for t,evs in events.items():
        for (p1,p2,w),trs in evs.items():
            out.write("{},{},{},{},{}\n".format(t,p1,p2,w,"/".join(trs)))

def extractIntrons(memsPath, Ref, exons, BitV, errRate, onlyPrimary):
    introns = {}
    lastID = ""
    for line in open(memsPath, 'r').readlines():
        alStrand, readID, err, mems, read = readLine(line)
        if onlyPrimary:
            if readID == lastID:
                continue
            lastID = readID

        if len(mems) > 1:
            for mem1,mem2 in pairwise(mems):
                # Remove ( and ) from mem and cast to int
                mem1 = [int(x) for x in mem1[1:-1].split(",")]
                id1 = BitV.rank(mem1[0] - 1)

                mem2 = [int(x) for x in mem2[1:-1].split(",")]
                id2 = BitV.rank(mem2[0] - 1)

                Poverlap = mem2[1]-mem1[1]-mem1[2]
                if id1 == id2: #MEMs inside the same exon
                    Toverlap = mem2[0]-mem1[0]-mem1[2]
                    # Intron Retention
                    if Poverlap <= 0 and Toverlap > 0:
                        gap = Toverlap+abs(Poverlap)-1
                        if gap > 0:
                            pos1 = exons[id1-1][0] + mem1[0] + mem1[2] - BitV.select(id1) + Poverlap
                            pos2 = pos1 + gap
                            key = (pos1, pos2)
                            introns[key] = introns[key]+1 if key in introns else 1
                else: #MEMs on different exons
                    offset1 = BitV.select(id1 + 1) - (mem1[0] + mem1[2])
                    offset2 = mem2[0] - BitV.select(id2)-1
                    if Poverlap <= 0:
                        Poverlap = abs(Poverlap)
                        # No gap on P: possible Int.Alt.S.S.
                        if offset1 == 0:
                            offset2 += Poverlap
                        else: #anyway, not only if offset2 == 0 !!! maybe this is wrong
                            offset1 += Poverlap
                        pos1 = exons[id1-1][1] - offset1 + 1
                        pos2 = exons[id2-1][0] + offset2 - 1
                        key = (pos1, pos2)
                        introns[key] = introns[key]+1 if key in introns else 1
                    else:
                        #Gap on P
                        if offset1 == 0 and offset2 == 0:
                            #No gap on T -> possible Ext.Alt.S.S.
                            intronStart, intronEnd  = exons[id1-1][1] + 1, exons[id2-1][0] - 1
                            intronString = Ref[intronStart-1:intronEnd-1] #-1 since strings are indexed starting from 0, gtf from 1
                            readGapString = read[mem1[1]+mem1[2]-1:mem2[1]-1]
                            maxErr = round(len(read)*errRate/100)
                            err1 = editDistance(readGapString, intronString[:len(readGapString)])
                            err2 = editDistance(readGapString, intronString[len(intronString)-len(readGapString):])

                            pos1, pos2 = -1, -1
                            if err1 <= err2:
                                # We check the start of the intron
                                if err1 + err <= maxErr:
                                    pos1 = intronStart + Poverlap
                                    pos2 = intronEnd
                            else:
                                #We check the end of the intron
                                if err2 + err <= maxErr:
                                    pos1 = intronStart
                                    pos2 = intronEnd - Poverlap

                            if pos1 != -1 and pos2 != -1:
                                key = (pos1, pos2)
                                introns[key] = introns[key]+1 if key in introns else 1
    return introns













class Vertex:
    def __init__(self, node, mtype):
        self.id = node
        self.adjacent = {}
        self.type = mtype

    def __str__(self):
        return str(self.id) + ' adjacent: ' + str([(x.id, x.get_weight(self)) for x in self.adjacent])

    def add_neighbor(self, neighbor, weight=0):
        self.adjacent[neighbor] = weight

    def get_connections(self):
        return self.adjacent.keys()

    def get_id(self):
        return self.id

    def get_adjacent(self):
        return self.adjacent

    def get_weight(self, neighbor):
        return self.adjacent[neighbor]

class Graph:
    def __init__(self):
        self.vert_dict = {}
        self.num_vertices = 0

    def __iter__(self):
        return iter(self.vert_dict.values())

    def has_vertex(self, node):
        return node in self.vert_dict

    def add_vertex(self, node, mtype):
        self.num_vertices = self.num_vertices + 1
        new_vertex = Vertex(node, mtype)
        self.vert_dict[node] = new_vertex
        return new_vertex

    def get_vertex(self, n):
        if n in self.vert_dict:
            return self.vert_dict[n]
        else:
            return None

    def add_edge(self, frm, to, cost = 0, mtype1 = "E", mtype2 = "E"):
        if frm not in self.vert_dict:
            self.add_vertex(frm, mtype1)
        if to not in self.vert_dict:
            self.add_vertex(to, mtype2)

        self.vert_dict[frm].add_neighbor(self.vert_dict[to], cost)
        self.vert_dict[to].add_neighbor(self.vert_dict[frm], cost)

    def get_vertices(self):
        return self.vert_dict.keys()

    def get_num_vertices(self):
        return self.num_vertices

    def print_graph(self):
        for v in self.vert_dict:
            print(self.vert_dict[v])



readSize = 100
L = 1
eps = 15
t1 = eps * readSize * L
t2 = eps * readSize
K0 = math.ceil(t1/100)
K1 = math.ceil(t1/100) - L + 1
K2 = math.ceil(t2/100)
# List of mems
intronInsertion = []

# Builds a mems graph including intronic mems 
def buildMEMsGraph(refPath, e_memsPath, i_memsPath, gtfPath):
    gtf = openGTF(gtfPath)

    # Read the reference genome
    genome = ""
    text = "|"
    intrText = "|"
    with open(refPath) as f:
        for line in f:
            if line[0] != ">":
                genome = genome + line[:-1]

    # Init Adjacency Matrix
    adjm = [[0, 0], [0, 0]]

    exons = set()
    introns = set()
    currentID = 1

    # Get the Exons and Introns from the GTF
    for g in gtf.features_of_type('gene'):
        for tr in gtf.children(g, featuretype='transcript', order_by='start'):
            exonsList = list(gtf.children(tr, featuretype='exon', order_by='start'))
            exons_ = set([(int(ex.start), int(ex.end)) for ex in exonsList])
            introns_ = set(zip([ex.end+1 for ex in exonsList[:-1]], [ex.start-1 for ex in exonsList[1:]]))

            # TODO check strand +/-

            # Add the exon text (if it's not already there)
            for ex in sorted(exons_):
                if ex not in exons:
                    text = text + genome[ex[0] : ex[1]] +  "|"
                    if currentID > 1:
                        # Expand Adjacency matrix
                        adjm = expandMatrix(adjm)

                        # Add arc e1 -> e2
                        adjm[currentID-1][currentID] = 1

                    currentID = currentID + 1
            
            # Same for the introns
            for intr in sorted(introns_):
                if intr not in introns:
                    intrText = intrText + genome[intr[0] : intr[1]] +  "|"

            exons = exons | exons_
            introns = introns | introns_

    # Is sorted needed?
    exons = sorted(exons)
    introns = sorted(introns)

    # Transitive Closure for adjm (TODO needed?)
    h = 1
    for e1 in exons:
        k = 1
        for e2 in exons:
            if e1 != e2 and e1[1] < e1[0]:
                if adjm[h][k] == 0:
                    adjm[h][k] = 2
            k = k + 1

        h = h + 1
    

    # Build the BitVector
    BitV = BitVector(text)
    IntrBitV = BitVector(intrText)

    # Build parents and sons
    parents, sons = buildParentsAndSons(adjm)

    # Empty the genome var (TODO needed?)
    genome = ""

    # Remove intron duplicates
    for i in introns:
        for o in introns:
            if i != o and i[0] >= o[0] and i[1] <= o[1]:
                introns.remove(i)
                break


    # Exon-Intron relationship
    # ex_in = [(e1_indx, e2_indx, i1_indx), (e, e, i)]
    ex_in = []
    if len(exons) > 1:
        exID = 1
        for e1, e2 in pairwise(exons):
            link = ()
            inID = 1

            for i in introns:
                if e1[1]+1 == i[0] and e2[0]-1 == i[1]:
                    link = (exID, exID+1, inID)
                    break
                inID = inID + 1

            if link:
                ex_in.append(link)

            exID = exID + 1



    # Init MEMs graph
    g = Graph()
    g.add_vertex("start", "Start")
    g.add_vertex("end", "End")
    
    
    # Store all MEMs in a single list [(mems, read), (ms, r), ...]
    emems = []
    for line in open(e_memsPath, 'r').readlines():
        line = line.replace("\t", ",")
        mem = strToMem(line)

        exon_text = getTextFromID(BitV, text, BitV.rank(mem[0] - 1))
        read = exon_text[mem[0] : mem[0] + mem[2]]
        emems.append((mem, read))

    # Store Intronic MEMs
    imems = []
    for line in open(i_memsPath, 'r').readlines():
        line = line.replace("\t", ",")
        mem = strToMem(line)

        intron_text = getTextFromID(IntrBitV, intrText, IntrBitV.rank(mem[0] - 1))
        read = intron_text[mem[0] : mem[0] + mem[2]]
        imems.append((mem, read))

    # Convert the MEMs list to a list where 
    # list[p] = [((t1, p, l1), r1), ...]
    e_mems = [[]] * (readSize + 1)
    for mr in emems:
        p = mr[0][1]
        if e_mems[p] == []:
            e_mems[p] = [mr]
        else:
            e_mems[p].append(mr)

    # Same for intron MEMs
    i_mems = [[]] * (readSize + 1)
    for mr in imems:
        p = mr[0][1]
        if i_mems[p] == []:
            i_mems[p] = [mr]
        else:
            i_mems[p].append(mr)
    

    
    # Insert the MEMs into the graph
    for mlist in e_mems:
        for mr in mlist:
            if mr:
                m = mr[0]
                read = mr[1]
                
                startInfo = isValidStart(m, parents, read, BitV, text)
                err = startInfo[1]

                # START
                if startInfo[0]:
                    if isNew(g, m):
                        # (If it doesn't have a father?)
                        if True:
                            # Add the new node and link it to the start
                            g.add_vertex(str(m), "E")       # mtype = E for Exons, = I for Introns
                            g.add_edge("start", str(m), err)
                    
                else:
                    if isNew(g, m):
                        continue

                # EXTEND
                isExt = False
                isNov = False
                maxRead = readSize - L + 1

                # Scan every MEM possibly connected to m (until maxRead)
                nextP = m[1] + 1

                while nextP <= maxRead:
                    for mr2 in e_mems[nextP]:
                        m2 = mr2[0]
                        if isNew(g, m2):
                            # True if m, m2 from the same exon, False otherwise
                            linkageInfo = checkMEMs(adjm, m, m2, BitV, text)
                            flag = linkageInfo[0]
                            err = linkageInfo[1]

                            if err >= 0:
                                if isNew(g, m2):
                                    g.add_vertex(str(m2), "E")

                                if flag:
                                    g.add_edge(str(m), str(m2), err)
                                    isExt = True
                                    isNov = True
                                else:
                                    # Add an arc in the Novel graph (not needed?)
                                    isNov = True

                    nextP += 1


                # END
                if not isExt and not isNov: 
                    endInfo = isValidEnd(m, sons, read, BitV, text)
                    endFlag = endInfo[0]
                    err = endInfo[1]

                    # If the mem doesn't have a son, link it to the end
                    if endFlag:
                        g.add_vertex(str(m), "E")
                        g.add_edge(str(m), "end", err)

                    
            
    print("> Built MG")
    #g.print_graph()
    print("")


    print("> Possible intr-fillable gaps:")

    # Try to fill gaps with Intron MEMs
    for (mem1, mem2) in intronInsertion:
        m1t = mem1[0]
        m1p = mem1[1]
        m1l = mem1[2]
        m2t = mem2[0]
        m2p = mem2[1]
        m2l = mem2[2]

        id1 = BitV.rank(m1t - 1)
        id2 = BitV.rank(m2t - 1)

        # Ignore same-exon mems for now? (TODO)
        if id1 == id2:
            continue

        # Find the intron id between the exons
        for (eID1, eID2, iID) in ex_in:
            #if (id1 == eID1 and id2 == eID2) or (id2 == eID1 and id1 == eID2): TODO?
            if id1 == eID1 and id2 == eID2:
                # There might be an intronic MEM that better connects m1 -> intr -> m2
                # Update the graph if we can insert an iMEM between m1 and m2
                print(mem1, "-->", mem2)
                exon1_text = getTextFromID(BitV, text, id1)
                exon2_text = getTextFromID(BitV, text, id2)
                intron_text = getTextFromID(IntrBitV, intrText, iID) # Needed? TODO

                exon1_read = exon1_text[m1t : m1t + m1l]
                exon2_read = exon2_text[m2t : m2t + m2l]

                if m2p + m2l > m1p + m1l:
                    gapP = m2p - m1p - m1l
                    print("Gap:", gapP)
                    gapE1 = 0
                    gapE2 = 0

                for imem in i_mems:
                    for imr in imem:
                        if imr:
                            im = imr[0]
                            iread = imr[1]

                            # Ignore imems from other introns, we only check the one in between
                            if iID == IntrBitV.rank(im[0] - 1):
                                # Check overlap mem1 -> intron
                                leftOverlap = 0
                                rightOverlap = 0
                                for k in range(len(exon1_read)):
                                    if overlap(exon1_read, iread, k):
                                        thisOverlap = len(exon1_read) - k
                                        if leftOverlap < thisOverlap:
                                            leftOverlap = thisOverlap
                                
                                for h in range(len(iread)):
                                    if overlap(iread, exon2_text, h):
                                        thisOverlap = len(iread) - h
                                        if leftOverlap < thisOverlap:
                                            leftOverlap = thisOverlap    
                                
                                newGapP = gapP - im[2] + leftOverlap + rightOverlap
                                # Found linking intron mem?
                                if newGapP <= K2:
                                    print("  >", mem1, im, mem2)

                                # TODO update graph
                break


# Computes the overlap between 2 strings, the first starting at s1Start
def overlap(s1, s2, s1Start):
    j = 0
    s1 = s1[s1Start:]

    for i in range(len(s1)):
        if s1[i] != s2[j]:
            return False
        j = j + 1

    return True




# Returns (true, error) if a mem is a valid starting one
def isValidStart(mem, parents, read, BitV, text):
    mt = mem[0]
    mp = mem[1]
    ml = mem[2]

    if mp <= K0:
        err = K2 + 1

        if mp == 1:
            err = 0
        else:
            id = BitV.rank(mt - 1)
            exon_text = getTextFromID(BitV, text, id)

            l = mp - 1 # Left in R

            subP = ""
            subE = ""
            if l != 0:
                subP = read[0 : l]

            exon_pref_len = mt - BitV.select(id)-1 - 1

            if exon_pref_len < l:
                shared_pref_len = l - exon_pref_len
                exon_pref = exon_text[0 : exon_pref_len]
                err = l
                
                # Exclude the fisrt parent = []
                for par in parents[id][1:]:
                    par_text = getTextFromID(BitV, text, par)

                    par_sel = BitV.select(par)
                    par_next = BitV.select(par + 1)
                    
                    # Check if the father is long enough to get its suffix
                    if par_next - shared_pref_len - par_sel - 1 >= 0:
                        subE = par_text[par_next - shared_pref_len - par_sel - 1 : shared_pref_len] + exon_pref
                    # ELse, get the entire label
                    else:
                        subE = par_text + exon_pref

                    newErr = editDistance(subP, subE)

                    if newErr < err:
                        err = newErr
            else:
                subE = exon_text[mt - l - BitV.select(id) - 2 : l]
                err = editDistance(subP, subE)

        if err <= K2:
            return (True, err)   

    return (False, K2 + 1)

# Returns (true, err) if the mem is a valid end
def isValidEnd(m, sons, read, BitV, text):
    mt = m[0]
    mp = m[1]
    ml = m[2]

    if mp + ml >= readSize - K0:
        err = K2 + 1

        if mp + ml == readSize + 1:
            err = 0
        else:
            id = BitV.rank(mt - 1)            
            exon_text = getTextFromID(BitV, text, id)

            l = readSize - (mp + ml) + 1
            subP = read

            subP = read[mp + ml - 1 : mp + ml - 1 + l]
            subE = ""

            exon_suff_len = BitV.select(id + 1) - (mt + ml) + 1

            if exon_suff_len < l:
                shared_suff_len = l - exon_suff_len
                exon_suff = ""
                if exon_suff_len != 0:
                    exon_suff = exon_text[mt + ml - s - 2 : mt + ml - s - 2 + exon_suff_len]

                err = l

                # Exclude the fisrt son = []
                for son in sons[id][1:]:
                    son_text = getTextFromID(BitV, text, son)
                    subE = exon_suff + son_text[0 : shared_suff_len]

                    curr_err = editDistance(subP, subE)

                    if curr_err < err:
                        err = curr_err

                
            else:
                subE = exon_text[mt + ml - s - 2 : mt + ml - s - 2 + l]
                err = editDistance(subP, subE)


        if err <= K2:
            return (True, err)
        
    return (False, K2 + 1)
    

def isNew(graph, mem):
    return not(graph.has_vertex(str(mem)))

# Expands an nxn matrix with 0s
def expandMatrix(m):
    nr = len(m[0]) + 1
    for r in m:
        r.append(0)
    m.append([0] * nr)

    return m


def buildParentsAndSons(sgm):
    parents = [[]] * (len(sgm))
    sons = [[]] * (len(sgm))

    for i in range(len(sgm)):
        for j in range(len(sgm[i])):
            if sgm[i][j] > 0:
                if parents[j] == []:
                    parents[j] = [i]
                else:
                    parents[j].extend([i])
                if sons[i] == []:
                    sons[i] = [j]
                else:
                    sons[i].extend([j])

    return (parents, sons)


# Returns true if 2 exons are linked
def contains(adjm, id1, id2):
    return adjm[id1][id2] >= 1

# Returns true if the link between the exons is new in sg
def isNewLink(adjm, id1, id2):
    return adjm[id1][id2] > 1

# Converts a mem string into a tuple of ints
def strToMem(str):
    mem = tuple([int(x) for x in str[0:-1].split(",")])
    return mem


def getTextFromID(BitV, text, eid):
    s = BitV.select(eid)
    e = BitV.select(eid+1)-1

    return text[s : e]

# Checks if 2 MEMs are connected and calcualtes the weight of the arc
def checkMEMs(adjm, mem1, mem2, BitV, text):
    m1t = mem1[0]
    m1p = mem1[1]
    m1l = mem1[2]
    m2t = mem2[0]
    m2p = mem2[1]
    m2l = mem2[2]

    # Get the exon IDs
    id1 = BitV.rank(m1t - 1)
    id2 = BitV.rank(m2t - 1)

    exon1_text = text[m1t : m1t + m1l]
    exon2_text = text[m2t : m2t + m2l]

    # Default error = -1
    err = -1
    resType = True

    # If mem1 and mem2 are of the same exon
    if id1 == id2:
        if m2p + m2l > m1p + m1l and m1t < m2t and m1t + m1l < m2t + m2l:
            gapP = m2p - m1p - m1l # Gap in R
            gapE = m2t - m1t - m1l # Gap in Z

            if gapP >= 0 and gapE >= 0:
                if gapP == 0: # If the gap is only in Z
                    if gapE > K2: # Possible Intron
                        err = 0
                        resType = False
                    else:         # Nothing
                        err = gapE
                        resType = True

                elif abs(gapP - gapE) <= K2: # Nothing / SNV Single Nucleotide Variation
                    subP = read[m1p + m1l - 1 : m1p + m1l - 1 + gapP]
                    subE = exon1_text[m1t + m1l - BitV.select(id1) - 2 : m1t + m1l - BitV.select(id1) - 2 + gapE]
                    
                    err = editDistance(subP, subE)
                    resType = True

            elif gapP <= 0 and gapE <= 0:

                err = abs(gapP - gapE)
                resType = True
                                                                            # TODO check intronInsertion here as well
            elif gapP <= 0 and gapE > K2: # Possible Intron with Overlap

                err = 0
                resType = False

            else:
                err = abs(gapP) + abs(gapE)
                resType = True

    # Else, if mem1 and mem2 are from different exons
    else:
        if contains(adjm, id1, id2): # If there is an edge id1 -> id2
            if m2p + m2l > m1p + m1l:
                gapP = m2p - m1p - m1l
                gapE1 = BitV.select(id1 + 1) + 1 - m1t - m1l
                gapE2 = m2t - BitV.select(id2) - 2

                if gapP <= 0:
                    err = 0
                    if not isNewLink(adjm, id1, id2) and gapE1 == 0 and gapE2 == 0:
                        resType = True
                   
                    elif err <= K2:
                        resType = False
                    else:
                        err = -1
                else:
                    if gapE1 == 0 and gapE2 == 0:
                        # Gap only in R, possible insertion
                        if not isNewLink(adjm, id1, id2):
                            err = 0
                            resType = True #FALSE  <-----------------------------------------
                            # Link the mems anyway but
                            # remember to check for introns later
                            intronInsertion.append((mem1, mem2))

                    else:
                        if abs(gapP - (gapE1 + gapE2)) <= K2: # SNV
                            subP = read[m1p + m1l - 1 : m1p + m1l - 1 + gapP]
                            subE1 = exon1_text[m1t + m1l - BitV.select(id1) - 2 : m1t + m1l - BitV.select(id1) - 2 + gapE1]
                            subE2 = exon2_text[0 : gapE2]
                            subE = subE1 + subE2

                            err = editDistance(subP, subE)

                            if not isNewLink(adjm, id1, id2):
                                resType = True
                            else:
                                resType = False
    
    # Final check on err
    if err > K2:
        err = -1

    return (resType, err)




















def main(memsPath, refPath, gtfPath, errRate, tresh, outPath, allevents):
    #TODO: add this as cmd line parameter
    onlyPrimary = False

    # Reading reference genome
    Ref = list(SeqIO.parse(refPath, "fasta"))[0]

    # !
    buildMEMsGraph(refPath, "input/case_2.ENST00000624081_e_9593_21_5012664.exons.mem", "input/case_2.ENST00000624081_e_9593_21_5012664.introns.mem", gtfPath)
    # Exit after the graph
    return 0

    # Reading annotation
    gtf = openGTF(gtfPath)
    strand, transcripts, annIntrons = extractFromGTF(gtf)

    # Extracting text and exons from "index" file
    text, exons = extractFromInfoFile(gtfPath + ".sg")
    BitV = BitVector(text)

    # Extracting introns from spliced graph-alignments
    introns = extractIntrons(memsPath, Ref, exons, BitV, errRate, onlyPrimary)

    # Cleaning introns
    if not(allevents):
        newIntrons, annFoundIntrons = filterAnnotated(introns, annIntrons)
        newIntrons = reconciliateIntrons(newIntrons, Ref, strand)
        newIntrons, annFoundIntrons_ = filterAnnotated(newIntrons, annIntrons)
        newIntrons = filterLowCovered(newIntrons, tresh)
        annFoundIntrons = filterLowCovered(annFoundIntrons, tresh)
        annFoundIntrons_ = filterLowCovered(annFoundIntrons_, tresh)
        allIntronsKey = set(newIntrons.keys()) | set(annFoundIntrons.keys())
    else:
        newIntrons = reconciliateIntrons(introns, Ref, strand)
        newIntrons = filterLowCovered(newIntrons, tresh)
        allIntronsKey = newIntrons.keys()

    # Extracting events from introns
    events = checkNewIntrons(newIntrons, allIntronsKey, strand, transcripts)
    printEvents(events, outPath)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Detects alternative splicing events from splice-aware alignments to a splicing graph")
    parser.add_argument('-g', '--genome', required=True, help='FASTA input file containing the reference')
    parser.add_argument('-a', '--annotation', required=True, help='GTF input file containing the gene annotation')
    parser.add_argument('-m', '--mems', required=True, help='input file containing the alignments to the splicing graph')
    parser.add_argument('-o', '--output', required=True, help='SAM output file')
    parser.add_argument('-e', '--erate', required=False, default=3, type=int, help='error rate of alignments (from 0 to 100, default: 3)')
    parser.add_argument('-w', '--support', required=False, default=3, type=int, help='minimum number of reads needed to confirm an event (default: 3)')
    parser.add_argument('--allevents', required=False, action='store_true', help='output all events, not only the novel ones')

    args = parser.parse_args()

    memsPath = args.mems
    refPath = args.genome
    gtfPath = args.annotation
    errRate = args.erate
    tresh = args.support
    outPath = args.output
    allevents = args.allevents

    main(memsPath, refPath, gtfPath, errRate, tresh, outPath, allevents)



# python3 ./scripts/detectEvents.py -g input/Homo_sapiens.GRCh38.dna.chromosome.21.fa -a input/ENSG00000279493.no2.gtf -m out/case_2.ENST00000624081_e_9593_21_5012664.exons.mem -o out/output.events.csv
