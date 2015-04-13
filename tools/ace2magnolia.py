#!/usr/local/bin/python
"""Get readcounts from an .ACE file per sample

Reads are coupled to a samples by using regular expressions.
Provide a set of regular expressions s.t. each read matches to one and only one expression.

USAGE:  ace2magnolya -a 454Contigs.ace -r "original;perturbed" -t newbler -o counts.txt

PARAMETERS:
-a    ACE file
-r    Regular expressions, each sample separated by semicolon
-o    Output file with read counts
-t    Assembler specific .ACE format [newbler, abyss or default] (optional)
-x    hits files prefix. List all alignment positions of the reads (optional) 

"""

import re
import getopt
import pdb
import sys
import numpy as np

class Sample():
    """Sequencing sample"""
    def __init__(self, regexp):
        self.pattern = re.compile(regexp)
    
    def match(self,str):
        """Check if readname belongs to sample"""
        m = self.pattern.match(str)
        return m


class Samples():
  def __init__(self, regexps):
        self.sample     = []
        for si in regexps:
            self.sample.append(Sample(si))

  def getSample(self,readname):

      return 0
      
      
class AF():
  def __init__(self,name, orient, offset,skip):
    self.name   = name
    self.orient = orient
    self.offset = offset
    self.skip   = skip

class NewblerAF(AF):
  def __init__(self,afline,pairtags):
    af,name,orient,offset = afline.split(" ")
    name, skip = self.process(name,pairtags)
    AF.__init__(self,name,orient,int(offset),skip)

  def process(self,str,pairtags):
  
    vals = {}
    
    #if str == "SRR096469.135931.164-278.fm11239" or str == "SRR096469.135931.1-163.to14808":
    #    pdb.set_trace()
    
    # Eat pair information
    pr   = re.compile('(.*)\.pr([0-9]+)$')
    m    = pr.match(str)
    if m is not None:
        str  = m.group(1)
        vals["pr"] = m.group(2)

    # Eat to edge read
    to  = re.compile('(.*)\.to([0-9]+)$')
    m   = to.match(str)
    if m is not None:
        str = m.group(1)
        vals["to"] = m.group(2)
    
    # Eat from edge read
    fm   = re.compile('(.*)\.fm([0-9]+)$')
    m    = fm.match(str)
    if m is not None:
        str  = m.group(1)
        vals["fm"] = m.group(2)
    
    # Eat range
    ra  = re.compile('(.*)\.([0-9]+)-([0-9]+)$')
    m   = ra.match(str)
    if m is not None:
        str = m.group(1)
        vals["ra"] = (int(m.group(2)),int(m.group(3)))
    
    # Eat pair information?            
    for pairtag in pairtags:                
        pt  = re.compile('(.*)(' + pairtag + ')')
        m   = pt.match(str)
        if m is not None:
            str = m.group(1)
            vals["pt"] = m.group(2)
            break
    
    # TODO Check with template, if not the read name has been eaten
    #vals['template']     = self.contig.reads[ai].ds.template
    vals['clippedname']  = str
    #vals['padded_start'] = self.contig.af[ai].padded_start
    
    if vals.has_key("fm"): skip = True
    else:                  skip = False
    
    return vals['clippedname'],skip
    
class AbyssAF(AF):
  def __init__(self,afline):
    af,name,orient,offset = afline.split(" ")
    name, skip = self.process(name)
    AF.__init__(self,name,orient,int(offset),skip)
    return
    
  def process(self, name):
    
    skip = False
    
    # Check if read started in other contig
    overlap = re.compile('(.*)_([0-9]+)\/([0-9]+)$')
    m       = overlap.match(name)
    if m is not None:
      str    = m.group(1)
      number = m.group(2)
      total  = m.group(3)
      if number != 1:
        skip = True
    else:
      str = name

    # Check if this is really a read
    contig = re.compile('CONTIG_[0-9]+$')
    m      = contig.match(name)
    if m is not None:
      skip = True
    else:
      str = name
      
    return str, skip


class run():
  
  def __init__(self,acefile,regexp,pairtags,type,outfile, hitsprefix=None):
    
    samples   = Samples(regexp)
    counts    = [0]*len(regexp)
    allcounts = {}
    contigID = None
    
    # Open one file for each sample to write the hits, i.e. "contigID \t position" for each read
    if hitsprefix is not None:
      hitshandles = []
      for i in range(1,len(regexp)+1):
        hitshandles.append(open(hitsprefix + str(i) +".txt","w"))
    
    print "type: ", type

    inseq = False
    infile = open(acefile,"r")
    while 1:

      line = infile.readline()
      if not line:
        break

      # New contig
      if line[0:3] == "CO ":
        inseq  = True
        offset = 0
        if contigID is not None:
          allcounts[contigID] = [counts, int(contigLen)]
        counts = [0]*len(regexp)
          
        vals = line.split(" ")
        if len(vals) != 6:
          raise Exception, ("Unexpected contig identifier line (not length 6): " + line)
        contigID  = vals[1]
        contigLen = vals[2]
        if hitsprefix is not None:
          nrgaps    = np.zeros(int(contigLen))
        print "Processing contig: ", contigID
      elif line[0:2] == "BQ":
        inseq = False
      elif inseq:
        if hitsprefix is not None:
          line   = line.rstrip()
          nrgaps = self.getGapPositions(line,offset,nrgaps)
          offset += len(line)
        


      # New read line
      if line[0:3] == "AF ":
        if   type == "default": af = AF(line)
        elif type == "newbler": af = NewblerAF(line,pairtags)
        elif type == "abyss":   af = AbyssAF(line)
        
        if af.skip is False:
          si = samples.getSample(af.name)
          counts[si] += 1
          #Print to hits files
          if hitsprefix is not None:
            if af.offset > 0:
              gaplesspos = af.offset - nrgaps[af.offset-1]
              print >>hitshandles[si], '%s\t%d' % (contigID, gaplesspos)
        
    infile.close() 
    
    if hitsprefix is not None:    
      for h in hitshandles:
        h.close()
    if contigID is not None:
       allcounts[contigID] = [counts, int(contigLen)]
    self.printCounts(allcounts,outfile)
  
  def getGapPositions(self,line,offset,nrgaps):
    """Cumalative number of gaps in sequence"""
    for i in range(len(line)):
      if offset == 0 and i == 0:
        if line[i] == "*":
          nrgaps[0] = 1
        else:
          nrgaps[0] = 0
      else:
        if i+offset >= len(nrgaps):
          pdb.set_trace()
        if line[i] == "*":
          nrgaps[i+offset] = nrgaps[i+offset-1]+1
        else:
          nrgaps[i+offset] = nrgaps[i+offset-1]
    return nrgaps
        
  def printCounts(self,allcounts,outfile):
    out = open(outfile,"w")
    ids = allcounts.keys()
    ids.sort()
    for id in ids:
      str = "%s\t%d" % (id, allcounts[id][1])
      for i in range(len(allcounts[id][0])):
        str += "\t%d" % (allcounts[id][0][i])
      print >>out, str
    out.close()
    

def main():
    
    # Set default
    acefile    = None
    regexp     = None
    outfile    = None
    type       = "default"
    pairtags   = ['/1','/2']
    hitsprefix = None
    
    # parse command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "x:ha:r:p:t:o:", ["help"])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        sys.exit(2)
    # process options
    for o, a in opts:
        if o in ("-h", "--help"):
            print __doc__
            sys.exit(0)
        if o == "-a":
            acefile = a
        if o == "-r":
            regexp   = a.split(";")
        if o =="-p":
            pairtags = a.split(";")
        if o =="-t":
            type = a
        if o =="-o":
            outfile = a
        if o =="-x":
            hitsprefix = a
    if acefile is None:
        print "Specifiy the folder with the acefile"
        print "for help use --help"
        sys.exit(2)
    if regexp is None:
        print "Specifiy a regular expressions per sample for read identification"
        print "for help use --help"
        sys.exit(2)
    if outfile is None:
        print "Specifiy an output file"
        print "for help use --help"
        sys.exit(2)
    run(acefile,regexp,pairtags,type,outfile, hitsprefix)
    
    
if __name__ == "__main__":
    main()
            
