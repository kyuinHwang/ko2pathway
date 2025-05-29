import os,sys,glob,random,gzip,argparse
import re
from collections import Counter

def moduleIter(iFP):
    lines = []
    for line in iFP:
        if line.rstrip() == "///":
            yield lines
            lines = []
        else: lines.append(line)
    if len(lines) > 0:  yield lines


def findChars(targetStr,chars,start,direction):
    ## Return the position of one of target chars which are nearest with start position
    if direction == -1:
        targetStr = targetStr[:start]
        return max([targetStr.rfind(char_s) for char_s in chars])
    else:
        idxs = [targetStr.find(char_s,start) for char_s in chars]
        idxs = list(filter(lambda x:x > -1,idxs))
        if len(idxs) == 0: return len(targetStr)+1
        else: return min(idxs)


def splitLevelone(moduleStr):
    level = 0
    idxs = []
    for idx, s_chr in enumerate(moduleStr):
        if s_chr == "(":
            level += 1
        elif s_chr == ")":
            level -= 1
        elif level == 0 and s_chr == "*":
            idxs.append(idx)
    part_strs = []
    for idx_start, idx_end in zip([0,]+idxs,idxs+[len(moduleStr)+1]):
        if idx_start == 0: part_strs.append(moduleStr[idx_start:idx_end])
        else: part_strs.append(moduleStr[idx_start+1:idx_end])
    return part_strs


def parseModuleDefinition(moduleDefinition,moduleID,debug):
    ## Reform
    if debug: print(f"Original {moduleID}:",moduleDefinition)
    moduleStr = moduleDefinition[:]

    ## replace " " => "*",  "," => "+"
    ## M00375 "  " ==> "*"
    moduleStr = moduleStr.replace(" -- "," ").replace(" --","").replace("-- ","").replace("  ","*").replace(' ',"*")
    if debug: print("Reformed:",moduleStr)

    ## minus sign denotes a non-essential component in the complex # 2023.08
    ## remove from euqation
    while moduleStr.find("-") > 0:
        if debug: print(f"Warning : Minus is included in {moduleID}")
        idxMinus = moduleStr.find("-")
        if moduleStr[idxMinus+1] == "(": ## -(K02134+K02314) ## Remove
            posEndR = findChars(moduleStr,")",idxMinus+1,1)
            moduleStr = moduleStr[:idxMinus] + moduleStr[posEndR+1:]
        elif moduleStr[idxMinus+1] == "K" or "M": ## -K12304 ## Remove
            posEndR = findChars(moduleStr,"+-() ",idxMinus+1,1)
            moduleStr = moduleStr[:idxMinus] + moduleStr[posEndR:]
        else:
            raise ValueError("Unexpected minus sign is found")
        if debug: print("Reformed:",moduleStr)

    ## plus signs are used to represent a complex
    ## replace "A+B" => "(A*B)"
    while moduleStr.find("+") > 0:
        idxPlus = moduleStr.find("+")
        posInsL = findChars(moduleStr,"(, *", idxPlus,-1) + 1
        posInsR = findChars(moduleStr,"), *", idxPlus+1,1)
        moduleStr = moduleStr[:posInsL] + '('  + moduleStr[posInsL:posInsR].replace('+','*')  + ')'  + moduleStr[posInsR:]
        if debug: print("Reformed:",moduleStr)

    moduleStr = moduleStr.replace(',',"+") ## comma separated K numbers indicate alternative
    if debug: print("Reformed:",moduleStr)
    splitedModule = splitLevelone(moduleStr)
    if debug: print("Splited:",splitedModule)
    return splitedModule


def parseModule(moduleFN,debug=False,submodule=False):
    print(f"Parsing {moduleFN}")
    moduleFP = gzip.open(moduleFN,'rt')
    moduleDict = {}
    moduleNameDict = {}
    cntMulti, cntSubmodule = 0, 0
    for linesModule in moduleIter(moduleFP):
        moduleID = linesModule[0].split()[1]
        if not moduleID.startswith("M"): break ## e.g. aor_M00001 ## module for each organisms
        moduleName = " ".join(linesModule[1].split()[1:])
        moduleNameDict[moduleID] = moduleName
        moduleDefinition = linesModule[2]
        if not moduleDefinition.startswith("DEFINITION"): raise ValueError(f"3rd line of module {moduleID} is not startswith DEFINITION")
        if debug: print(f"###################################,{moduleID}")
        ## 2020.07.07. Multiline definition exists for several modules (e.g. M00127 )
        ## option 1 (submodule=True) : divide each line as individual module (i.e. submodule)
        ## option 2 (else) : concatenate lines and regarded as one single module
        if submodule:
            splitedModule = parseModuleDefinition(moduleDefinition[12:].strip(),moduleID,debug)
            moduleDict[moduleID] = splitedModule    
            for idx, defLine in enumerate(linesModule[3:]):
                if defLine.startswith("ORTHOLOGY"): break
                splitedModule = parseModuleDefinition(defLine.strip().strip(),moduleID,debug)
                moduleDict[moduleID+f'_{idx+1}'] = splitedModule
                moduleNameDict[moduleID+f'_{idx+1}'] = moduleName
                cntSubmodule += 1
            if idx > 1: cntMulti += 1
        else:
            defLines = []
            for idx, defLine in enumerate(linesModule[2:]):
                if defLine.startswith("ORTHOLOGY"): break
                defLines.append(defLine)
            defLines[0] = defLines[0][12:] ## trimout "DEFINITION "
            moduleDef = ' '.join([line.strip() for line in defLines])
            splitedModule = parseModuleDefinition(moduleDef,moduleID,debug)
            moduleDict[moduleID] = splitedModule    
            if len(defLines) > 1: cntMulti += 1
            
    print(f"{len(moduleDict)} modules detected and {cntMulti} modules contained multiple line definition")
    if submodule: print(cntSubmodule, "submodules defined")
    return moduleDict, moduleNameDict

        
        


def main():
    parser = argparse.ArgumentParser(description="Track module step coverage based on KO profiles.")
    parser.add_argument("-k", "--ko_dir", required=True, help="Directory containing .ko files for each genome")
    parser.add_argument("-m", "--module_file", required=True, help="KEGG module definition file")
    parser.add_argument("-o", "--output", required=True, help="Output file for step coverage results")
    parser.add_argument("--debug", action="store_true", help="Enable debug mode")

    args = parser.parse_args()

    # 예시 print
    print(f"Reading KO files from: {args.ko_dir}")
    print(f"Using module definitions from: {args.module_file}")
    print(f"Saving step coverage result to: {args.output}")

    koDir, outFN = args.ko_dir, args.output
    debug = args.debug
    debug_module = args.debug
    usesubmodule = False
    moduleFN = args.module_file #"kegg/module/module.gz"
    moduleDict, moduleNameDict = parseModule(moduleFN, debug_module, submodule=usesubmodule)

    ## Update this
    ## ['(K00844+K12407+K00845+K00886+K08074+K00918)', '(K01810+K06859+K13810+K15916)', '(K00850+K16370+K00918)', '(K01623+K01624+K11645+K16305+K16306)', 'K01803', '((K00134+K00150)*K00927+K11389)', '(K01834+K15633+K15634+K15635)', 'K01689', '(K00873+K12406)']
    moduleDict["M00001"][2]='(K00850+K16370+K00918+K21071)' ## Update

    summaryDict = {}#{moduleID:{} for moduleID in moduleDict}
    #targetModules = {line.strip():0 for line in open("targetModule.txt")}#dict.fromkeys(["M00001","M00003","M00004","M00008"],0)
    targetModules = {moduleID:0 for moduleID in "M00001 M00003 M00009 M00012 M00309 M00086 M00087".split()} ## Target Modules ## Only M00009 include minus (-)
    for moduleID, moduleContents in sorted(moduleDict.items()):
        for stepIDX, modulePart in enumerate(moduleContents): summaryDict[(moduleID,stepIDX)] = {"Def":modulePart}

    genomeAccL = []
    for koFN in sorted(glob.glob(f"{koDir}/*.ko")):
        koCnt = Counter([line.rstrip().split('\t')[1] for line in open(koFN,'r') if len(line.rstrip().split('\t')) > 1])
        genomeAcc = os.path.basename(koFN).replace('.ko','')
        genomeAccL.append(genomeAcc)
        for moduleID, moduleContents in sorted(moduleDict.items()):
            if debug:
                print("====",moduleID,"====")
                print(moduleContents)
            partSucess = []
            if moduleID not in targetModules: continue
            for stepIDX, modulePart in enumerate(moduleContents):
                ## See M00611
                transModulePart = modulePart[:]
                moIDs = re.findall("M\d\d\d\d\d",modulePart)                
                while len(moIDs) > 0:
                    for moID in moIDs: transModulePart = transModulePart.replace(moID,'*'.join(moduleDict[moID]))                    
                    moIDs = re.findall("M\d\d\d\d\d",transModulePart)

                koIDs = re.findall("K\d\d\d\d\d",transModulePart)
                for koID in koIDs: transModulePart = transModulePart.replace(koID,str(koCnt[koID]))
                if debug:
                    print(modulePart)
                    print(transModulePart)
                    print(eval(transModulePart))
#                if eval(transModulePart) > 0:  partSucess.append(True)
#                else: partSucess.append(False)
                summaryDict[(moduleID,stepIDX)][genomeAcc] = eval(transModulePart) > 0
#            summaryDict[moduleID][genomeID] = sum(partSucess)

    ### Write result
    wFP = open(outFN,'w')
    wFP.write("Module\t#Steps(Lv.1)\tStepDef\t" + "\t".join(genomeAccL) + '\tName\n')
    for moduleID,stepIDX in sorted(summaryDict):
        if moduleID not in targetModules: continue
        moduleName = moduleNameDict[moduleID]
        stepDef = summaryDict[(moduleID,stepIDX)]["Def"]
        Ncontains = [summaryDict[(moduleID,stepIDX)][genomeAcc] for genomeAcc in genomeAccL]
        #if sum(Ncontains) == 0: continue
        strNcontains = '\t'.join(list(map(lambda x:str(x)[0],Ncontains)))
        wFP.write(f"{moduleID}\t{stepIDX+1}\t{stepDef}\t{strNcontains}\t{moduleName}\n") ##!! 0 based index to 1 based indexx
    wFP.close()

        
if __name__ == "__main__": main()
