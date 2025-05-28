import os,sys,glob,random
#from fasta_parser import fasta_iter
#from itertools import combinations
from collections import Counter
#from statistics import mean 
class SAG:
    def __init__(self,name):
        self.genomeAcc = name
        pass

def readKeyEnzymes(listFN):
    listFP = open(listFN)
    keyEnzymeD = {}
    reactionNames = []
    for line in listFP:
        if line.startswith("#"): continue
        reactionName, genesetName, geneAccLStr = line.strip().split('\t')
        keyEnzymeD.setdefault(reactionName, [])
        geneAccL = []
        for geneAccStr in geneAccLStr.replace(' ','').split(','): ## K17222-7,K22622
            if geneAccStr.startswith("BLAST"):
                geneAccL.append(f"BLAST.{genesetName}")
            elif geneAccStr.startswith("K") and '-' in geneAccStr: ## K17222-7
                koStart, koEnd = geneAccStr.split('-')
                koEnd = koStart[:-len(koEnd)] + koEnd
                geneAccL.extend([f"K{koNum:05d}" for koNum in range(int(koStart[1:]),int(koEnd[1:])+1)])
            else:
                geneAccL.append(geneAccStr)
        keyEnzymeD[reactionName].append((genesetName,geneAccL))
        if reactionName not in reactionNames: reactionNames.append(reactionName)
    return reactionNames, keyEnzymeD            

def countNonZero(values):
    return countUp(values,0)

def countUp(values, threshold):
    if threshold != 0:
        return len([value for value in values if value >= threshold])
    else:
        return len([value for value in values if value > threshold])

def replaceNumber(value):
    if value >= 1: return 'Y'
    else: return 'N'

def checkKeyEnzymes(koDir,outFN,keyFile,blastSulfurFN,idx):
    keyword = '.'.join(os.path.basename(keyFile).split('.')[:1])
    keyEnzymeD = {}
    reactionList = []
    pReactionList, pkeyEnzymeD = readKeyEnzymes(keyFile)
    reactionList.extend(pReactionList)
    keyEnzymeD.update(pkeyEnzymeD)
#    print(keyEnzymeD)
    koFNs = sorted(glob.glob(f"{koDir}/*.ko"))
    qcPassed = []
    for koFN in koFNs:
        genomeAcc = os.path.basename(koFN).replace('.ko','')
        sag = SAG(genomeAcc)    
        sag.koCnt = Counter([line.rstrip().split('\t')[1] for line in open(koFN,'r') if len(line.rstrip().split('\t')) > 1])
        qcPassed.append(sag)
    sagD = {sag.genomeAcc:sag for sag in qcPassed}
    for sag in qcPassed: sag.blastAnnotation = {}
    for line in open(blastSulfurFN,'r'):
        locusTag, geneName = line.strip().split('\t')
        genomeAcc = locusTag[:locusTag.rfind('_')]
        if genomeAcc not in sagD: continue
        sagD[genomeAcc].blastAnnotation.setdefault(geneName,[])
        sagD[genomeAcc].blastAnnotation[geneName].append(locusTag)

    blastCntD = {}
    for sag in qcPassed:
        blastCntD[sag.genomeAcc] = Counter({key.lower():len(values) for key, values in sagD[sag.genomeAcc].blastAnnotation.items()})
    
    for sag in qcPassed: sag.genesets, sag.reactions = {}, {}    

    for reactionName, genesetList in sorted(keyEnzymeD.items(),key=lambda x:reactionList.index(x[0])): ## Reaction : aprAB or aprM or qmoABC
        for sag in qcPassed: sag.reactions[reactionName] = []
        for genesetName, geneAccL in genesetList: ## Gene set : aprAB
            for sag in qcPassed: sag.genesets[genesetName] = []
            for geneAcc in geneAccL: ## Single gene : aprA
                for sag in qcPassed:
                    if geneAcc.startswith("BLAST"):
                        geneCnt = blastCntD[sag.genomeAcc][geneAcc[6:].lower()] ## BLAST.{genesetName}
                    elif geneAcc.startswith("K") and len(geneAcc) == 6:
                        geneCnt = sag.koCnt[geneAcc]
                    else: raise ValueError(f"Unexpected geneAcc, {geneAcc}")
                    sag.genesets[genesetName].append(geneCnt)
            for sag in qcPassed: sag.reactions[reactionName].append(countNonZero(sag.genesets[genesetName])/len(geneAccL))

    if idx == 0:
        summaryFP = open(outFN,'w')
        summaryFP.write("Category\tRank\tReaction\tGeneset/#Geneset\tGeneName/#Gene\t" + '\t'.join([sag.genomeAcc for sag in qcPassed]) + '\n')
    else:
        summaryFP = open(outFN,'r')
        genomeIDs = summaryFP.readline().rstrip('\n').split('\t')[5:]
        if genomeIDs != [sag.genomeAcc for sag in qcPassed]: raise ValueError(f"Mismatch found from sample order in column : {genomeIDs}, {[sag.genomeAcc for sag in qcPassed]}")
        summaryFP.close()
        summaryFP = open(outFN,'a')

    for reactionName, genesetList in sorted(keyEnzymeD.items(),key=lambda x:reactionList.index(x[0])): ## Reaction : aprAB or aprM or qmoABC
        summaryFP.write(f"{keyword}\tReaction.Half\t{reactionName}\t{len(genesetList)}\t\t" + '\t'.join([replaceNumber(countUp(sag.reactions[reactionName],0.5)) for sag in qcPassed]) + '\n')
    summaryFP.close()                
                

def main():
    koDir,outFN = sys.argv[1],sys.argv[2]
    keyFiles = sys.argv[3].split(',')
    blastSulfurFN = sys.argv[4]
    for idx, keyFile in enumerate(keyFiles):
        checkKeyEnzymes(koDir,outFN,keyFile,blastSulfurFN,idx)

if __name__ == "__main__": main()
