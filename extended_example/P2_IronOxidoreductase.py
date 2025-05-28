import os,sys,glob,random
from functools import reduce
#from fasta_parser import fasta_iter
#from itertools import combinations
from collections import Counter
from statistics import mean 
#import pickle

def isnonzero(value):
    if value == '1': return "Y"
    elif value == '0': return "N"
    else: raise ValueError("Unexpected value detected")


def IronOxidation():
    fegenieSummaryFN, keyCheckFN = sys.argv[1], sys.argv[2]
    name1,name2,name3 = ["Neutrophilic","Acidophilic","Netrophilic Probable"]
    summaryFP = open(fegenieSummaryFN,'r')
    genomeL = summaryFP.readline().strip().split('\t')[2:]
    keywordL = [("iron_oxidation","Cyc1"),("iron_oxidation","Cyc2_repCluster1"), ("iron_oxidation","Cyc2_repCluster2"), ("iron_oxidation","Cyc2_repCluster3"),
                 ("possible_iron_oxidation_and_possible_iron_reduction","MtoA"),("possible_iron_oxidation_and_possible_iron_reduction","MtrB_TIGR03509")]
    includingD = {}
    for line in summaryFP:
        for keywords in keywordL:
            ## Criteria for selecting line
            if all([keyword in line for keyword in keywords]):
                values = line.strip().split('\t')[2:]
                includingGenomes = [genomeName for genomeName, value in zip(genomeL,values) if value == "TRUE" or int(value) > 0]
                includingD[keywords] = set(includingGenomes)
                break ## Required because either (Rev) or ?
    for keywords in keywordL:
        if keywords not in includingD:   includingD[keywords] = set()
    keywords1, keywords2, keywords3 =  (("iron_oxidation","Cyc2_repCluster1"),), (("iron_oxidation","Cyc1"),("iron_oxidation","Cyc2_repCluster2"),("iron_oxidation","Cyc2_repCluster3")), (("possible_iron_oxidation_and_possible_iron_reduction","MtoA"),("possible_iron_oxidation_and_possible_iron_reduction","MtrB_TIGR03509"))

    genomeL_new = open(keyCheckFN,'r').readline().strip().split('\t')[5:] ##!
    if set(genomeL) != set(genomeL_new): raise ValueError(genomeL, genomeL_new)
    wFP = open(keyCheckFN,'a')
    for name, keywords in zip([name1,name2,name3],[keywords1,keywords2,keywords3]):
        anyOfEnzyme = reduce(lambda x,y:x.union(y),[includingD[keyword] for keyword in keywords])
        anyOfEnzyme = dict.fromkeys(anyOfEnzyme,"1")

        wFP.write(f"Iron\tReaction.Half\t{name} iron oxidation\t{len(keywords)}\t\t")
        wFP.write('\t'.join([isnonzero(anyOfEnzyme.get(genome,"0")) for genome in genomeL]) + '\n')
    wFP.close()


def IronReduction():
    fegenieSummaryFN, keyCheckFN = sys.argv[1], sys.argv[2]
    name1,name2 = ["iron_reduction","probable_iron_reduction"]
    summaryFP = open(fegenieSummaryFN,'r')
    genomeL = summaryFP.readline().strip().split('\t')[2:]
    keywordL = [("iron_reduction","DFE_0461"),("probable_iron_reduction","MtrA"),("iron_reduction","MtrA")]
    includingD = {}
    for line in summaryFP:
        for keywords in keywordL:
            ## Criteria for selecting line
            if all([keyword in line for keyword in keywords]):
                values = line.strip().split('\t')[2:]
                includingGenomes = [genomeName for genomeName, value in zip(genomeL,values) if value == "TRUE" or int(value) > 0]
                includingD[keywords] = set(includingGenomes)
                break ## Required because either (Rev) or ?

    for keywords in keywordL:
        if keywords not in includingD:   includingD[keywords] = set()

    keywords1, keywords2 =  (("iron_reduction","DFE_0461"),("iron_reduction","MtrA")),(("probable_iron_reduction","MtrA"),)

    genomeL_new = open(keyCheckFN,'r').readline().strip().split('\t')[5:] ##!
    if set(genomeL) != set(genomeL_new): raise ValueError(genomeL, genomeL_new)
    wFP = open(keyCheckFN,'a')
    for name, keywords in zip(["iron_reduction",],[(("iron_reduction","DFE_0461"),("iron_reduction","MtrA"),("probable_iron_reduction","MtrA")),]):
        anyOfEnzyme = reduce(lambda x,y:x.union(y),[includingD[keyword] for keyword in keywords])
        anyOfEnzyme = dict.fromkeys(anyOfEnzyme,"1")

        wFP.write(f"Iron\tReaction.Half\t{name}\t{len(keywords)}\t\t")
        wFP.write('\t'.join([isnonzero(anyOfEnzyme.get(genome,"0")) for genome in genomeL]) + '\n')

    wFP.close()

    

if __name__ == "__main__": 
    IronOxidation()
    IronReduction()
