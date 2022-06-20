import string
import sys
import re
import os
from pyspark import SparkConf, SparkContext
import time
from datetime import timedelta

conf = SparkConf()
conf.setAppName("proteinscoring")

_matchscore = 1
_missscore   = 0
_firstgappenalty = -2
_furthergappentaly = -1
_topN = 25
_minScoreToPrint = 1

_gap = "-"
_totalPairsEvaluated =0

def _GetScoreOnePair (s1, s2 ):
    global _totalPairsEvaluated

    thisscore = 0
    firstgap = False
    
    m = min(len(s1), len(s2))
    for pos in range (0, m):
        if (s1[pos] == s2[pos]) and (s2[pos] != _gap) :
            thisscore += _matchscore
            firstgap = False
        elif (s1[pos] == _gap) and ( s2[pos] != _gap ) :
            if firstgap == True:
                thisscore += _furthergappentaly
            else:
                firstgap = True
                thisscore += _firstgappenalty
        elif (s2[pos] == _gap) and ( s1[pos] != _gap ) :
            if firstgap == True:
                thisscore += _furthergappentaly
            else:
                firstgap = True
                thisscore += _firstgappenalty
        else :
            firstgap = False
            thisscore += _missscore

    _totalPairsEvaluated = _totalPairsEvaluated + 1
    if ( _totalPairsEvaluated%100000) == 0 :
        print (_totalPairsEvaluated)
        
    return thisscore
                

def  _OutputProteinPairWithScore (_sl1, _sl2, score):
    str1 = _sl1.rstrip()
    str2 = _sl2.rstrip()
    m = max(len(str1), len(str2));


    if score >= _minScoreToPrint:
        print ("pos", end = " " )
        for pos in range (0, m):
            print(pos, end = " ")
        print("  score = ", score)

        print ("   ", end = " " )
        for pos in range (0, len(str1)):
            print (str1[pos], end = " " )
        print("")

        print ("   ", end = " " )
        for pos in range (0, len(str2)):
            print (str2[pos], end = " " )
        print("")
        print("")

def _CreateProteinVariantsStartingSpaces ( _s1, _diff_len ):
    result = []
    result.append(_s1)
    
    for i in range (1, _diff_len):
        s = _s1.rjust(len(_s1)+i)
        result.append(s)
    return result

def _CreateProteinVariantsWithGap (_sl1, _num_gaps, _max_num_gaps, _min_index, _max_index ):
    strings = []
    if  _num_gaps == _max_num_gaps:
        strings.append(_sl1)
        
    strinlen = len(_sl1)
    if strinlen >  _max_index:
        strinlen =  _max_index
        
    #do not add a gap at the beginning and the end    
    for i in range (_min_index+1, strinlen-1):
        stl = _sl1[0:i] + _gap + _sl1[i:];
        strings.append(stl)
        #print (_min_index, _max_index, stl)
        if  _num_gaps > 1 :
            stl2 = _CreateProteinVariantsWithGap(stl, _num_gaps-1, _max_num_gaps,  _min_index, _max_index )
            for s in stl2:
                strings.append(s)
            
    return strings

def _GetScoresWithVariants ( _sl1, _sl2, max_gaps ):
    currentmin = -1000 # Current minimum value to be added to the vec
    currentmax = -1000 # Maximum sore observed so far
    scorevec = []
    s1vec    = []
    
    s1 = _sl1.upper()
    s2 = _sl2.upper()
    #print("printing s1",s1)
    #print ("printing s2",s2)

    diff_len = len(s2) - len(s1)
    if diff_len < 0 :
        diff_len = 0
    s1x = _CreateProteinVariantsStartingSpaces (s1, diff_len )
    for s in s1x:
        max_len = min (len(s), len(s2)) + max_gaps
        min_index = 0
        for j in range (0, len(s) ):
            if s[j] == " ":
                min_index = min_index + 1
            else:
                break

        s1var = _CreateProteinVariantsWithGap(s, max_gaps, max_gaps, min_index, max_len )
        for sx in s1var:
            s1vec.append(sx)
            
    #String s2 is assumed to be from a database, so we do not calculate variants on it

    #remove duplicates
    s1variants = list(dict.fromkeys(s1vec))
    
    i = 0
    for s1var in s1variants:
        score = _GetScoreOnePair (s1var, s2 )
        tup = (score, s1var, s2)
        
        if i < _topN:
            scorevec.append(tup)
            if  i == 0:
                currentmin = score
            elif  score < currentmin:
                currentmin = score
                if i == 0 :
                    currentmax = score
                elif score > currentmax:
                    currentmax = score
        else:
            if  score > currentmin:
                scorevec.append(tup)
        i = i+1
        scorevec.sort(key = lambda x : -x[0]);
        scorevec = scorevec[0:_topN]
        currentmin = min (scorevec)[0]
        currentmax = max (scorevec)[0]
            
    return scorevec

if __name__ == "__main__":

    print("Usage : python3 proteinscoring-with-test.py proteinfile.txt dbfile.txt max_gaps")
    start_time = time.time()
    allvec = []
    num_entries_in_db = 0
    proteinsequence = open(sys.argv[1]).readline()
    
    dbfile = sys.argv[2]
    
    max_gaps = int(sys.argv[3])
    sc = SparkContext.getOrCreate()
    
    text=sc.textFile(dbfile)
    num_entries_in_db = text.count()
    variantTuple = text.flatMap(lambda s2: _GetScoresWithVariants(proteinsequence, s2,max_gaps))
    #sorted_score = variantTuple.sortBy(lambda x:x[0], ascending=False) # part 2
    sorted_score = variantTuple.sortBy(lambda x:x[0], ascending=False, numPartitions=4) # part 3 and 4
    RDD_topN = sorted_score.take(_topN)
    f = open('FinalOutputFile.txt', 'w')
    for t in RDD_topN:
        f.writelines(str(t) + "\n")
    f.close()
    print("File created")
    print("Number of entries in database " , num_entries_in_db)    
    print("Total Pairs evaluated " , _totalPairsEvaluated)
    
    elapsed_time_secs = time.time() - start_time
    print("Execution took: %s secs " % elapsed_time_secs)
    #sc.exit()