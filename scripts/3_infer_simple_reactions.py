import os,sys,glob,random,argparse
from collections import Counter

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
            if geneAccStr.startswith("K") and '-' in geneAccStr: ## K17222-7
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



def checkKeyEnzymes(koDir,outFN,keyFile,idx,debug=False):
    keyword = '.'.join(os.path.basename(keyFile).split('.')[:1])
    keyEnzymeD = {}
    reactionList = []
    pReactionList, pkeyEnzymeD = readKeyEnzymes(keyFile)
    reactionList.extend(pReactionList)
    keyEnzymeD.update(pkeyEnzymeD)
    koFNs = sorted(glob.glob(f"{koDir}/*.ko"))
    qcPassed = []
    for koFN in koFNs:
        genomeAcc = os.path.basename(koFN).replace('.ko','')
        sag = SAG(genomeAcc)    
        sag.koCnt = Counter([line.rstrip().split('\t')[1] for line in open(koFN,'r') if len(line.rstrip().split('\t')) > 1])
        qcPassed.append(sag)
    
    for sag in qcPassed: sag.genesets, sag.reactions = {}, {}    

    for reactionName, genesetList in sorted(keyEnzymeD.items(),key=lambda x:reactionList.index(x[0])): ## Reaction : aprAB or aprM or qmoABC
        for sag in qcPassed: sag.reactions[reactionName] = []
        for genesetName, geneAccL in genesetList: ## Gene set : aprAB
            for sag in qcPassed: sag.genesets[genesetName] = []
            for geneAcc in geneAccL: ## Single gene : aprA
                for sag in qcPassed:
                    if geneAcc.startswith("K") and len(geneAcc) == 6:
                        geneCnt = sag.koCnt[geneAcc]
                    else: raise ValueError(f"Unexpected geneAcc, {geneAcc}")
                    sag.genesets[genesetName].append(geneCnt)
            for sag in qcPassed: sag.reactions[reactionName].append(countNonZero(sag.genesets[genesetName])/len(geneAccL))

    if idx == 0:
        summaryFP = open(outFN,'w')
        summaryFP.write("Category\tReaction\tGeneset/#Geneset\tGeneName/#Gene\t" + '\t'.join([sag.genomeAcc for sag in qcPassed]) + '\n')
    else:
        summaryFP = open(outFN,'r')
        genomeIDs = summaryFP.readline().rstrip('\n').split('\t')[4:]
        if genomeIDs != [sag.genomeAcc for sag in qcPassed]: raise ValueError(f"Mismatch found from sample order in column : {genomeIDs}, {[sag.genomeAcc for sag in qcPassed]}")
        summaryFP.close()
        summaryFP = open(outFN,'a')

    for reactionName, genesetList in sorted(keyEnzymeD.items(),key=lambda x:reactionList.index(x[0])): ## Reaction : aprAB or aprM or qmoABC
        summaryFP.write(f"{keyword}\t{reactionName}\t{len(genesetList)}\t\t" + '\t'.join([replaceNumber(countUp(sag.reactions[reactionName],0.5)) for sag in qcPassed]) + '\n')
        if debug:
            for idx, genesetInfo in enumerate(genesetList): ## Gene set : aprAB
                genesetName, geneAccL = genesetInfo 
                summaryFP.write(f"{keyword}\t{reactionName}\t{genesetName}\t{len(geneAccL)}\t" + '\t'.join([f"{countNonZero(sag.genesets[genesetName])/len(geneAccL):.2f}" for sag in qcPassed]) + '\n')
                for geneidx, geneAcc in enumerate(geneAccL): ## Single gene : aprA
                    summaryFP.write(f"{keyword}\t{reactionName}\t{genesetName}\t{geneAcc}\t" + '\t'.join([str(sag.genesets[genesetName][geneidx]) for sag in qcPassed]) + '\n')                
    summaryFP.close()                
                

def main():
    parser = argparse.ArgumentParser(description="Search key enzymes for single-step reactions")
    parser.add_argument("-k", "--ko_dir", required=True, help="Directory containing .ko files for each genome")
    parser.add_argument("-g", "--keyenzyme_dir", required=True, help="Directory containing KOs for simple reactions")
    parser.add_argument("-o", "--output", required=True, help="Output file for step coverage results")
    parser.add_argument("--debug", action="store_true", help="Enable debug mode")
    args = parser.parse_args()

    # 예시 print
    print(f"Reading KO files from: {args.ko_dir}")
    print(f"Check key enzymes written in : {args.keyenzyme_dir}")
    print(f"Saving step coverage result to: {args.output}")    

    koDir,outFN = args.ko_dir, args.output
    debug = args.debug
    keyFiles = glob.glob(f"{args.keyenzyme_dir}/*.txt")

    if len(keyFiles) == 0: raise ValueError("No key enzyme files are detected")
    for idx, keyFile in enumerate(keyFiles):
        checkKeyEnzymes(koDir,outFN,keyFile,idx,debug)

if __name__ == "__main__": main()
