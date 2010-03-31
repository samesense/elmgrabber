#---------------------------------------
#
# Author     Perry Evans
#            evansjp@mail.med.upenn.edu
# 2008
#---------------------------------------
""" Functions for reading the my protein 
    annotation file, which is organized:

    geneid st stp annotation_name seq annotation_tool
    
    and making annotated multiple alignment plots. """
import utils_graph, inspect
import math, time, sys
import PyMozilla, ClientCookie, ClientForm, urllib2,urllib
import table_parser

cdict = {
'red'  :  ((0., 0., 0.), (0.5, 0.25, 0.25), (1., 1., 1.)),
'green':  ((0., 1., 1.), (0.7, 0.0, 0.5), (1., 1., 1.)),
'blue' :  ((0., 1., 1.), (0.5, 0.0, 0.0), (1., 1., 1.))
}

def annotation2protein(afile, tool_dict):
    d = dict()
    f = open(afile)
    for line in f.xreadlines():
        #if len(line.split('\t')) != 6:
        #    print 'ERROR', afile, line
        #    sys.exit(0)
        [geneid, start, stop,
         domain, seq_desc, tool] = [x.strip() for x in line.split('\t')]
        start = int(start)
        stop = int(stop)
        if tool_dict.has_key( tool ):
            if not d.has_key( domain ): d[domain] = {}
            if not d[domain].has_key(geneid): d[domain][geneid] = []
            add = True
            for [st, stp, old_desc] in d[domain][geneid]:
                if start == st and stp == stop:
                    add = False
            if add:
                d[domain][geneid].append([start, stop, seq_desc])
    f.close()
    return d

def annotation2protein_forProteins(afile, tool_dict, proteins):
    d = dict()
    f = open(afile)
    for line in f.xreadlines():
        #if len(line.split('\t')) != 6:
        #    print 'ERROR', afile, line
        #    sys.exit(0)
        [geneid, start, stop,
         domain, seq_desc, tool] = [x.strip() for x in line.split('\t')]
        if proteins.has_key(geneid):
            start = int(start)
            stop = int(stop)
            if tool_dict.has_key( tool ):
                if not d.has_key( domain ): d[domain] = {}
                if not d[domain].has_key(geneid): d[domain][geneid] = []
                add = True
                for [st, stp, old_desc] in d[domain][geneid]:
                    if start == st and stp == stop:
                        add = False
                if add:
                    d[domain][geneid].append([start, stop, seq_desc])
    f.close()
    return d

def protein2annotation(afile, tool_dict):
    """ Parse the given file and return {}
        mapping proteins to their annotations.

    @param afile: annotation file; format
    geneid st stp annotation seq tool
    """

    d = dict()
    f = open(afile)
    for line in f.xreadlines():
        [geneid, start, stop,
         domain, seq_desc, tool] = [x.strip() for x in line.split('\t')]
        start = int(start)
        stop = int(stop)
        if tool_dict.has_key( tool ):
            if not d.has_key( geneid ): d[geneid] = {}
            if not d[geneid].has_key(domain):
                d[ geneid ][ domain ] = []
            add = True
            for [st, stp, old_desc] in d[geneid][domain]:
                if start == st and stp == stop:
                    add = False
            if add:
                d[geneid][domain].append([start, stop, seq_desc])
    f.close()
    return d

def protein2annotation_forProteins(afile, tool_dict, proteins):
    """ Parse the given file and return {}
        mapping proteins to their annotations.

    @param afile: annotation file; format
    geneid st stp annotation seq tool
    """

    d = dict()
    f = open(afile)
    for line in f.xreadlines():
        [geneid, start, stop,
         domain, seq_desc, tool] = [x.strip() for x in line.split('\t')]
        if proteins.has_key(geneid):
            start = int(start)
            stop = int(stop)
            if tool_dict.has_key( tool ):
                if not d.has_key( geneid ): d[geneid] = {}
                if not d[geneid].has_key(domain):
                    d[ geneid ][ domain ] = []
                add = True
                for [st, stp, old_desc] in d[geneid][domain]:
                    if start == st and stp == stop:
                        add = False
                if add:
                    d[geneid][domain].append([start, stop, seq_desc])
    f.close()
    return d

def protein2annotation_forMotifs(afile, motifs):
    d = dict()
    f = open(afile)
    for line in f.xreadlines():
        [geneid, start, stop,
         domain, seq_desc, tool] = [x.strip() for x in line.split('\t')]
        start = int(start)
        stop = int(stop)
        if motifs.has_key(domain):
            if not d.has_key( geneid ): d[ geneid ] = {}
            if not d[geneid].has_key(domain):
                d[ geneid ][ domain ] = []
            add = True
            for [st, stp, old_desc] in d[geneid][domain]:
                if start == st and stp == stop:
                    add = False
            if add:
                d[geneid][domain].append([start, stop, seq_desc])
    f.close()
    return d

def annotation2protein_forMotifs(afile, motifs):
    d = dict()
    f = open(afile)
    for line in f.xreadlines():
        [geneid, start, stop,
         domain, seq_desc, tool] = [x.strip() for x in line.split('\t')]
        start = int(start)
        stop = int(stop)
        if motifs.has_key(domain):
            if not d.has_key( domain ): d[ domain ] = {}
            if not d[domain].has_key(geneid):
                d[ domain ][ geneid ] = []
            add = True
            for [st, stp, old_desc] in d[domain][geneid]:
                if start == st and stp == stop:
                    add = False
            if add:
                d[domain][geneid].append([start, stop, seq_desc])
    f.close()
    return d

def motif2protein_forProteinLs(afile, ls):
    d = dict()
    f = open(afile)
    for line in f.xreadlines():
        sp = line.split('\t')
        if ls.has_key( sp[0] ):
            if not d.has_key( sp[3] ): d[ sp[3] ] = dict()
            d[ sp[3] ][ sp[0] ] = True
    f.close()
    return d

def mergeLs(ls):
    correctOrder = True
    for i in xrange(len(ls)-1):
        if ls[i][1] >= ls[i+1][0]:
            correctOrder = False
            new_st = ls[i][0]
            new_stp = ls[i+1][1]
            del ls[i+1]
            ls[i] = [new_st, new_stp]
            break
    return correctOrder

def getOverlap_runner(st, stp, ls):
    overlap = 0
    for [ls_st, ls_stp] in ls:
        if st >= ls_st and stp <= ls_stp:
            # search region is entirely in this region
            overlap = stp - st + 1
            break
        elif st <= ls_st and stp >= ls_stp:
            # search region contains this region
            overlap += ls_stp - ls_st + 1
        elif st <= ls_st and stp >= ls_st and stp <= ls_stp:
            # search region begins before this region and ends inside it
            overlap += stp - ls_st + 1
        elif st >= ls_st and st <= ls_stp and stp >= ls_stp:
            # search region starts inside this resion and ends outside it
            overlap += ls_stp - st + 1
        elif stp < ls_st:
            break
    return overlap

def getOverlap(ls1, ls2):
    """ Given 2 sets of annotations, find the base/residue overlap between them. """

    overlap = 0
    for [st, stp] in ls1:
        overlap += getOverlap_runner(st, stp, ls2)
    return overlap

def mergeMotifs(start, end, ls):
    """ Merge based on overlaping motifs and adjacent motifs. """

    merged = False
    for i in xrange(len(ls)):
        oldStart = ls[i][0]
        oldEnd = ls[i][1]
        if start >= oldStart and end <= oldEnd:
            # search region is entirely in this region
            merged = True
        elif start <= oldStart and end >= oldEnd:
            # search region contains this region
            ls[i] = [start, end]
            merged = True
        elif start <= oldStart and end >= oldStart-1 and end <= oldEnd:
            # search region begins before this region and ends inside it
            ls[i] = [start, oldEnd]
            merged = True
        elif start >= oldStart and start <= oldEnd+1 and end >= oldEnd:
            # search region starts inside this resion and ends outside it
            ls[i] = [oldStart, end]
            merged = True
    if not merged:
        ls.append( [start, end] )
    ls.sort()
    correctOrder = mergeLs(ls)
    while correctOrder == False:
        correctOrder = mergeLs(ls)

def printAnnotation_runner(annotationLs, fasta):
    toColor = []
    for [start, stop, seq] in annotationLs:
        mergeMotifs(start, stop, toColor)
    line = ''
    toColor = sorted(toColor)    
    previousStop = 0
    for i in xrange( len(toColor) -1 ):        
        [start, stop] = toColor[i]
        [nextStart, nextStop] = toColor[i+1]
        line = line + fasta[previousStop:start-1] + '\033[1;31m' + fasta[start-1:stop] + '\033[0m'
        previousStop = stop
    [start, stop] = toColor[-1]
    line = line + fasta[previousStop:start-1] + '\033[1;31m' + fasta[start-1:stop] + '\033[0m' + fasta[stop:]
    return line

def printAnnotation(annotationLs, fasta):
    print printAnnotation_runner(annotationLs, fasta)

def printAnnotation_tofile(annotationLs, fasta, afile):
    afile.write(printAnnotation_runner(annotationLs, fasta) + '\n')

def printAnnotation01_runner(annotationLs, fasta):
    toColor = []
    for [start, stop, seq] in annotationLs:
        mergeMotifs(start, stop, toColor)
    line = ''
    toColor = sorted(toColor)
    for [st,stp] in toColor:
        if stp<st:
            print 'error'
    for i in xrange(len(toColor)-1):
        if toColor[i][1] >= toColor[i+1][0]:
            print 'error in merge or order', toColor[i], toColor[i+1]
    previousStop = 0
    for i in xrange( len(toColor) -1 ):        
        [start, stop] = toColor[i]
        [nextStart, nextStop] = toColor[i+1]
        for x in xrange(previousStop,start-1):
            line = line + '0 '
        for x in xrange(start-1,stop):
            line = line + '1 '
        #line = line + fasta[previousStop:start-1] + '\033[1;31m' + fasta[start-1:stop] + '\033[0m'
        previousStop = stop
    [start, stop] = toColor[-1]
    #line = line + fasta[previousStop:start-1] + '\033[1;31m' + fasta[start-1:stop] + '\033[0m' + fasta[stop:]
    for x in xrange(previousStop, start-1):
        line = line + '0 '
    for x in xrange(start-1, stop):
        line = line + '1 '
    for x in xrange(stop, len(fasta)):
        line = line + '0 '
    return line.strip()

def printAnnotation0X_runner(annotationLs, fasta, X):
    toColor = []
    for [start, stop, seq] in annotationLs:
        mergeMotifs(start, stop, toColor)
    line = ''
    toColor = sorted(toColor)
    for [st,stp] in toColor:
        if stp<st:
            print 'error'
    for i in xrange(len(toColor)-1):
        if toColor[i][1] >= toColor[i+1][0]:
            print 'error in merge or order', toColor[i], toColor[i+1]
    previousStop = 0
    for i in xrange( len(toColor) -1 ):        
        [start, stop] = toColor[i]
        [nextStart, nextStop] = toColor[i+1]
        for x in xrange(previousStop,start-1):
            line = line + str(X) + ' '
        for x in xrange(start-1,stop):
            line = line + '0 '
        #line = line + fasta[previousStop:start-1] + '\033[1;31m' + fasta[start-1:stop] + '\033[0m'
        previousStop = stop
    [start, stop] = toColor[-1]
    #line = line + fasta[previousStop:start-1] + '\033[1;31m' + fasta[start-1:stop] + '\033[0m' + fasta[stop:]
    for x in xrange(previousStop, start-1):
        line = line + str(X) + ' '
    for x in xrange(start-1, stop):
        line = line + '0 '
    for x in xrange(stop, len(fasta)):
        line = line + str(X) + ' '
    return line.strip()

def printAnnotation01_runner_2annotations(annotationLs1,
                                          annotationLs2, fasta):
    toColor = []
    for [start, stop, seq] in annotationLs1:
        mergeMotifs(start, stop, toColor)
    line = ''
    toColor = sorted(toColor)    
    previousStop = 0
    for i in xrange( len(toColor) -1 ):        
        [start, stop] = toColor[i]
        [nextStart, nextStop] = toColor[i+1]
        for x in xrange(previousStop,start-1):
            line = line + '0 '
        for x in xrange(start-1,stop):
            line = line + '1 '

        previousStop = stop
    
    [start, stop] = toColor[-1]

    for x in xrange(previousStop, start-1):
        line = line + '0 '
    for x in xrange(start-1, stop):
        line = line + '1 '
    for x in xrange(stop, len(fasta)):
        line = line + '0 '
    return line.strip()

def printAnnotation01(annotationLs, fasta):
    print printAnnotation01_runner(annotationLs, fasta)

def printAnnotation01_tofile(annotationLs, fasta, afile):
    afile.write(printAnnotation01_runner(annotationLs, fasta) + '\n')

def printAnnotation0X_tofile(annotationLs, fasta, afile, X):
    afile.write(printAnnotation0X_runner(annotationLs, fasta, X) + '\n')

def printAnnotation01_tofile_2annotations(annotationLs1, 
                                          annotationLs2, fasta, afile):
    afile.write(printAnnotation01_runner_2annotations(annotationLs1,
                                                      annotaitonLs2, 
                                                      fasta) + '\n')
        
def extendSeq(max_len, seq):
    """ All sequnces in the alignment must
        be the same length.  This fills in
        the seq with 0 until max_len is reached.
    
    @param max_len: desired seq len
    @param seq as a list of 1 and 0 for motif presence/absense
    @return: seq w/ filled in 0's
    """

    while len(seq) != max_len:
        seq.append(0)
    return seq

def loadMatrix(afile):
    alignment_matrix = []
    longest_len = 0
    matrix_file = open(afile)
    for line in matrix_file:
        alignment_matrix.append([])
        for item in line.strip().split():
            alignment_matrix[-1].append(int(item))
        seq_len = len(alignment_matrix[-1])
        if seq_len > longest_len:
            longest_len = seq_len    
    return [[extendSeq(longest_len, seq) for seq in alignment_matrix],
            longest_len]

def mkBlankMatrix(matrix):
    new_matrix = []
    for row in matrix:
        new_matrix.append([])
        for item in row:
            new_matrix[-1].append(0)
    return new_matrix

def getShortestAlignmentLength(motif_name_file_ls, motif_matrix_dir):
    shortest_alignment_len = 5000
    for motif_name_file in motif_name_file_ls:
        afile = open(motif_name_file)
        for line in afile:
            [location, name] = line.strip().split('\t')
            [matrix, ll] = loadMatrix(motif_matrix_dir + location)
            alignment_len = len(matrix)
            if alignment_len < shortest_alignment_len:
                shortest_alignment_len = alignment_len
        afile.close()
    return shortest_alignment_len

def concatMatrix(motif, motifs_to_add, shortest_alignment_len):
    cat_matrix = [[]]
    for [name, motifs, blank] in motifs_to_add:
        if motifs.has_key(motif):
            for index in xrange(shortest_alignment_len):
                if len(cat_matrix) - 1 < index:
                    cat_matrix.append([])
                for item in motifs[motif][0][index]:
                    cat_matrix[index].append(item)
        else:
            for index in xrange(shortest_alignment_len):
                if len(cat_matrix) - 1 < index:
                    cat_matrix.append([])
                #else: print 'no', len(cat_matrix), index
                #print index, cat_matrix
                #print len(cat_matrix)
                for item in blank[index]:
                    cat_matrix[index].append(item)
    print len(cat_matrix)
    return cat_matrix

def concatCols(matrix_ls):
    new_matrix = []
    for matrix in matrix_ls:
        for row in matrix:
            new_matrix.append(row)
    return new_matrix        

def myParse(page):
    suffix = page.split("REFRESH")[1].split('URL=')[1].split("\"")[0]
    if suffix.find('http') == -1:
        return 'http://elm.eu.org/basicELM/' + suffix.replace(';','&')
    return suffix.replace(';', '&')

def getELMpage(protein, sequence, dump_file):
    print 'getting page', protein
    request=urllib2.Request('http://elm.eu.org/')
    response = urllib2.urlopen(request)
    forms = ClientForm.ParseResponse(response)
    form = forms[0]
    response.close()
    form['sequence'] = '>' + protein + '\n' + sequence
    request2 = form.click()
    opener = ClientCookie.build_opener(ClientCookie.HTTPEquivProcessor,
                                       ClientCookie.HTTPRefreshProcessor)
    ClientCookie.install_opener(opener)
    page = ClientCookie.urlopen(request2).read()
    newURL = myParse(page)
    res = ClientCookie.urlopen(newURL)
    response = res.read()
    res.close()
    print 'made it'
    currentTime = time.clock()
    oldTime = currentTime
    #with open(protein + '.init.html', 'w') as f:
    #    f.write(response)
    print 'writing file'
    while response.find('atient') != -1:
        if currentTime-oldTime > float(20):
            print 'reread'
            newURL = myParse(response)
            try:
                print newURL
                res = ClientCookie.urlopen(newURL)
                response = res.read()
                res.close()
                #with open(protein + '.' + str(currentTime) + '.html', 'w') as f:
                #    f.write(response)
                #print response
            except:
                print 'FAIL', protein
            oldTime = currentTime
        currentTime = time.clock()
    f = open(dump_file, 'w')
    f.write(response)
    f.close()

def findFieldIndex(field_name, ls):
    for i in xrange(ls):
        if ls[i].strip().lstrip() == field_name:
            return i

def fixPositionLs(position_ls):
    new_ls = []
    for item in position_ls:
        if item[0] != '':
            new_ls.append(item)
    return new_ls

def parseELMpage(gene_name, elm_file):
    f = urllib.urlopen(elm_file)
    p = table_parser.TableParser()
    p.feed(f.read())
    ret = []    
    if len(p.doc) < 4: print 'problem w/ file', gene_name
    table_index = 0
    for table in p.doc:
        if len(table) > 1:
            if len(table[0]) > 1:
                if table[0][0].strip().lstrip() == 'Elm Name':
                    if table[0][1].strip().lstrip() == 'Instances\n(Matched Sequence)':
                        for row in p.doc[table_index]:
                            if row[0].strip().lstrip() != 'Elm Name':
                                #print row, elm_file
                                if type(row[2]) == type(table_parser.Table()):
                                    row[2] = fixPositionLs(row[2])
                                    for i in xrange(len(row[1])):
                                        [start, stop] = row[2][i][0].split('-')
                                        if len(row[2][i]) > 1 or len(row[1][i]) > 1:
                                            print 'problem'
                                        # gene start stop motif_type seq pattern                                        
                                        ret.append( [gene_name, start, stop, row[0].strip().lstrip(),
                                                     row[1][i][0].strip().lstrip(), 'ELM'] ) 
                                else:
                                    positions = str(row[2]).split('\n')
                                    for i in xrange(len(row[1])):            
                                        [start, stop] = positions[i].split('-')
                                        seq = row[1][i]
                                        if type(seq) == type(table_parser.Row()):
                                            seq = seq[0]
                                        ret.append( [gene_name, start, stop,
                                                     str(row[0]).strip().lstrip(),
                                                     seq.strip().lstrip(), 'ELM'] )
                            
        table_index += 1       
    return ret

def parseELMpage_v2(gene_name, elm_file):
    #print 'CALLED'
    f = urllib.urlopen(elm_file)
    p = table_parser.TableParser()
    p.feed(f.read())
    ret = []    
    if len(p.doc) < 5: print 'problem w/ file', gene_name, len(p.doc)
    for table_index in xrange(3, len(p.doc)-1):
        for row in p.doc[table_index][1:]:
            if type(row[2]) == type(table_parser.Table()):
                row[2] = fixPositionLs(row[2])
                for i in xrange(len(row[1])):
                    [start, stop] = row[2][i][0].split('-')
                    if len(row[2][i]) > 1 or len(row[1][i]) > 1:
                        print 'problem'
                # gene start stop motif_type seq pattern                                        
                    ret.append( [gene_name, start, stop, row[0].strip().lstrip(),
                                 row[1][i][0].strip().lstrip(), 'ELM'] ) 
            else:
                positions = str(row[2]).split('\n')
                for i in xrange(len(row[1])):            
                    [start, stop] = positions[i].split('-')
                    seq = row[1][i]
                    if type(seq) == type(table_parser.Row()):
                        seq = seq[0]
                    ret.append( [gene_name, start, stop,
                                 str(row[0]).strip().lstrip(),
                                 seq.strip().lstrip(), 'ELM'] )

    return ret

def parseExcludedELMs(gene_name, elm_file, protein_seq):
    """ Rip out the excluded ELMs for the ELM
        resource HTML.  They are in the last
        table. """

    f = urllib.urlopen(elm_file)
    p = table_parser.TableParser()
    p.feed(f.read())
    ret = []    
    if len(p.doc) < 4: print 'problem w/ file', gene_name
    table_index = 0
    for table in p.doc:
        if len(table) > 1:
            if len(table[0]) > 1:
                if table[0][0].strip().lstrip() == 'Elm Name':
                    if table[0][1].strip().lstrip() == 'Positions':
                        for row in p.doc[table_index]:
                            if row[0].strip().lstrip() != 'Elm Name':
                                if type(row[1]) == type(table_parser.Table()):
                                    row[1] = fixPositionLs(row[1])
                                    for i in xrange(len(row[1])):
                                        [start, stop] = row[1][i][0].split('-')
                                        if len(row[1][i]) > 1:
                                            sys.stderr.write('problem\n')
                                        # gene start stop motif_type seq pattern
                                        seq = protein_seq[int(start)-1:int(stop)]
                                        ret.append( [gene_name, start, stop, row[0].strip().lstrip(),
                                                     seq, 'ELM'] ) 
                                else:
                                    positions = str(row[1]).split('\n')
                                    for i in xrange(len(row[1])):                                        
                                        [start, stop] = positions[i].split('-')
                                        seq = protein_seq[int(start)-1:int(stop)]
                                        ret.append( [gene_name, start, stop,
                                                     row[0].strip().lstrip(),
                                                     seq, 'ExcludedELM'] )
        table_index += 1       
    return ret

#def renderDiagram():
#    parser = Bio.GenBank.FeatureParser()
#    fhandle = open('NC_005213.gbk')
#    genbank_entry = parser.parse(fhandle)
    
#    feature_set = GenomeDiagram.GDFeatureSet('test')
#    for feature in genbank_entry.features:
#        if feature.type == 'CDS':
#            feature_set.add_feature(feature)
#    gdgs = GenomeDiagram.GDGraphSet('GC Content')
#    graphdata = GenomeDiagram.GDUtilities.gc_content(genbank_entry.seq, 100)
#    gdgs.new_graph(graphdata, 'GC content', style='line')
#    gdt1 = GenomeDiagram.GDTrack('CDS features', greytrack=1)
#    gdt2 = GenomeDiagram.GDTrack('GC Content', greytrack=1)
#    gdt1.add_set(feature_set)
#    gdt2.add_set(feature_set)
#    gdd = GenomeDiagram.GDDiagram('test')
#    gdd.add_track(gdt1, 2)
#    gdd.add_track(gdt2, 4)
#    gdd.draw(format='linear', orientation='landscape',
#             tracklines=0, pagesize='A5', fragments=5, circular=0)
#    gdd.write('test.ps', 'PS')

def annotation2frac(annotation2proteins, bg_genes):
    fracs = {}
    bg_len = float(len(bg_genes.keys()))
    for annotation in annotation2proteins.keys():
        proteins_w_annotation = float(len(annotation2proteins[annotation].keys()))
        fracs[annotation] = proteins_w_annotation/bg_len
    return fracs

def all_perms(str):
    if len(str) <=1:
        yield str
    else:
        for perm in all_perms(str[1:]):
            for i in range(len(perm)+1):
                yield perm[:i] + str[0:1] + perm[i:]

def all_perms_ELM(elm_reg_exp):
    # the | operator must be split 
    # ex. .RK|RR[^KR]
    # will be permuted .RK then RR[^KR]
    # and combined
    pass

                
