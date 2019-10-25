############################################################################
#
#  bin.py
#  cmd-line sets parameters as well as positions of data in records
#  data is read into Record objects,
#
#  author: eightbitastronomy (eightbitastronomy@protonmail.com)
#  license: none. 
#           This source code may be used freely and without constraint.
#  Last updated: 25 oct 2019 for commenting.
#
############################################################################
#
#  File format:
#  Commented lines are accepted (see global specifiers, below) 
#  Columns (fields) should be space-separated, and column headings should
#    be specified by the "spec" line.
#  For example, for # (comment) and ## (spec)...
#
#  # someprog output 1/1/1999
#  ##var1 var2 var3 output1 r-squared
#  .5 -1.1 10.5 100.0 .996
#  .5 -1.2 13.6 101.0 .993
#  .4 -.3  9.8  1.2e+10 .544
#  ...etc etc
#
#  Use:
#  User must tell the script the column/field which will be used
#    for binning (--binposition). If from each bin a min or max value should
#    be calculated from another field and output, then --minimum/--maximum
#    should be used (along with --valposition). If not, use --histogram.
#  Predetermined bin sizes are specified with --size, and lower and upper
#    bounds on bins for the data as a whole are specified with --lower and
#    --upper.
#  Due to the nature of floating point comparison, relative error in
#    comparisons may be specified with --relerror (see specifiers for the
#    default.
#  Script accepts pipes (e.g., cat mydata.out | bin.py --blablabla). If
#    no piping is done, then the final cmdline argument must be the
#    name/path of the input data file.
#  Output is to stdout; user must redirect it howsoever they wish.
#  Yes, unix-style options flags are available. Use -h or --help.
#
############################################################################


'''Binning script for numerical/floating point data. 
   Possibly more useful as a skeleton for future tools
   than actually as a bin utility.'''


import argparse
import sys
import os




### global specifiers ###


SPEC    = '##'
REMARK  = '#'
RELERROR = .0001



### functions ###


def getLen(items, dummy):
    '''Why is this useful? Good question.'''
    return len(items)


def getMax(itemlist, pos):
    '''A maximum-finding function'''
    if len(itemlist) == 0 :
        return None
    runningitem = itemlist.pop(0)
    runningmax = runningitem[pos]
    for item in itemlist:
        if item[pos] > runningmax :
            runningmax = item[pos]
            runningitem = item
    return runningitem


def getMin(itemlist, pos):
    '''A minimum-finding function'''
    if len(itemlist) == 0 :
        return None
    runningitem = itemlist.pop(0)
    runningmin = runningitem[pos]
    for item in itemlist:
        if item[pos] < runningmin :
            runningmin = item[pos]
            runningitem = item
    return runningitem


def prepareBins(lower, upper, size):
    '''Create bins to be used for histogram'''
    counter = 0
    current = lower
    binlist = []
    emptylists = []
    while current <= upper:
        counter += 1
        nextbin = lower + counter*size
        binlist.append((current, nextbin,))
        emptylists.append([])
        current = nextbin
    return (binlist, emptylists, )


def joinString(prefix, item):
    '''Preparation of values for printing'''
    strbuffer = prefix
    for v in item.values():
        strbuffer = strbuffer + "{:1.8e} ".format(v)
    return strbuffer

def placeInBin(item, bins, bintemplate, params):
    '''Place a record line/item into a bin'''
    guess = int((item[params.binposition]-params.lower)/params.size)
    if guess > len(bintemplate):
        return
    for i in range(guess-1, guess+2):
        if (item[params.binposition] < bintemplate[i][1]) and (item[params.binposition] >= bintemplate[i][0]):
            bins[i].append(item)
            return
    return


### Command-line argument parser ###


def parseCommandLine(filetest):
    '''Parse arguments: argument indicates whether to expect a filename'''
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument("-b","--binposition",required=True,help="position 0..n-1, or name, of bin-key",metavar="")
    parser.add_argument("-v","--valposition",default=None,help="position 0..n-1, or name, of comparison-key",metavar="")
    parser.add_argument("-s","--size",required=True,type=float,help="bin size",metavar="")
    parser.add_argument("-l","--lower",required=True,type=float,help="lowest bin lower limit",metavar="")
    parser.add_argument("-u","--upper",required=True,type=float,help="highest bin upper limit",metavar="")
    parser.add_argument("-r","--relerror",type=float,default=RELERROR,help="relative error for comparisons",metavar="")
    group.add_argument("-H","--histogram",action="store_true",help="histogram mode: store bin counts")
    group.add_argument("-M","--maximum",action="store_true",help="maximum mode: store maximum values of bins")
    group.add_argument("-m","--minimum",action="store_true",help="minimum mode: store minimum values of bins")
    if filetest:
        parser.add_argument("filename")
    args = parser.parse_args()
    if ( ( args.maximum is True ) or ( args.minimum is True ) ) and ( args.valposition is None ) :
        print("Comparison binning (-M or -m) requires value-key, -v")
        return None
    if ( args.valposition is not None ) and ( ( args.maximum is False ) and ( args.minimum is False ) ) :
        print("Value-key (-v) found but comparison binning (-M or -m) is not requested")
        return None
    return args






############################# Begin ################################


### Check for a pipe ###

if os.isatty(0): # this only works in linux as far as I know!
    # no pipe:
    guide = parseCommandLine(True)
    if guide is not None :
        data = open(guide.filename,"r",-1)
    else:
        sys.exit(0)
else:
    # yes pipe:
    guide = parseCommandLine(False)
    if guide is not None :
        data = sys.stdin
    else:
        sys.exit(0)

        
### declare variables for processing data ###


template = []
process = None


### set functions according to cmdline args ###


if guide.histogram is True:
    process = getLen
elif guide.maximum is True:
    process = getMax
else:
    process = getMin
    

### Sanity checks on user imperatives (i.e., cmdline args) ###


if guide.binposition.isdigit():
    Record = type("Record",(dict,object),{ '__init__': lambda self,l: dict.__init__(self,[(str(i), float(l[i])) for i in range(len(l))]) }    )
    if ( guide.valposition is not None ) and ( not guide.valposition.isdigit() ):
        print("Conflicting field reference types: int and str")
        data.close()
        sys.exit(0)
else:
    if guide.valposition.isdigit():
        print("Conflicting field reference types: str and int")
        data.close()
        sys.exit(0)
    #read in the file definition
    for line in data:
        if line[0:len(SPEC)] == SPEC:
            print(line[0:len(line)-1])
            template = ( line[2:len(line)-1] ).split()
            Record = type("Record",(dict,object),{ '__init__': lambda self,l: dict.__init__(self,[(template[i], float(l[i])) for i in range(len(template))]) }    )
            break
        if line[0:len(REMARK)] == REMARK:
            print(line[0:len(line)-1])
            continue


### Prepare bins for histogram ###


(mastertemplate,listoflists) = prepareBins(guide.lower,guide.upper,guide.size)


### Begin reading in records, culling comments ###


for line in data:
    if line[0:len(REMARK)] == REMARK:
        print(line[0:len(line)-1])
        continue
    lbuffer = Record( (line[0:len(line)-1]).split() )
    placeInBin( lbuffer, listoflists, mastertemplate, guide )


### Output ###


for i,sublist in enumerate(listoflists):
    strbuffer = "{:0.3f} {:0.3f} ".format(mastertemplate[i][0],mastertemplate[i][1])
    slbuffer = process(sublist,guide.valposition)
    if slbuffer is None :
        print(strbuffer)
        continue
    if type(slbuffer) is int:
        print( strbuffer + str(slbuffer) )
    else:
        print( joinString(strbuffer,slbuffer) )


### Finish up ###


data.close()
sys.exit(0)
