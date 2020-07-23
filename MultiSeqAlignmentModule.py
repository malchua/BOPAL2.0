from Event import Event
import numpy as np
import globals

#############################################
### Multiple Sequence Alignment Functions ###
#############################################

######################################################
# computeMultiSeqAlignmentMatrix
# Parameters: Three orthologous operons which we are aligning
# Description: Creates a score matrix for determining the optimal alignment
######################################################
def computeMultiSeqAlignmentMatrix(strain1, strain2, strain3):
    if globals.printToConsole:
        print('Computing global alignment matrix for: {%s, %s}...' % (strain1.name, strain2.name))

    #initialize the matrix to store the global alignment scores
    #alignmentMatrix = [[[ 0.0 for x in range(0, len(strain2.genomeFragments))] for y in range(0, len(strain1.genomeFragments))] for z in range(0, len(strain3.genomeFragments))]
    #eventMatrix = [[None for x in range(0, len(strain2.genomeFragments))] for y in range(0, len(strain1.genomeFragments))]


def performMultiSequenceAlignment(operon1, operon2, operon3):

    event = Event(0)
    event.setScore(1.0)
    event.setDistance(abs(5 - 10))
    event.setFragmentDetails1(operon1)
    event.setFragmentDetails2(operon2)
    event.setGenome1Name("strain1")
    event.setGenome2Name("strain2")
    event.setTechnique('Multiple Sequence Alignment')

    #initialize the distance matrix
    scoreMatrix = np.zeros((len(operon1)+1, len(operon2)+1, len(operon3)+1))
    dirMatrix = np.zeros((len(operon1)+1, len(operon2)+1, len(operon3)+1), dtype=int)

    for a in range(1, len(operon1)+1):
        scoreMatrix[a][0][0] = scoreMatrix[a-1][0][0] + calculateSumOfPairsScore(operon1[a-1], '-', '-')
        dirMatrix[a][0][0] = 7

    for a in range(1, len(operon2)+1):
        scoreMatrix[0][a][0] = scoreMatrix[0][a-1][0] + calculateSumOfPairsScore('-', operon2[a-1], '-')
        dirMatrix[0][a][0] = 6

    for a in range(1, len(operon3)+1):
        scoreMatrix[0][0][a] = scoreMatrix[0][0][a-1] + calculateSumOfPairsScore('-', '-', operon3[a-1])
        dirMatrix[0][0][a] = 5

    for a in range(1, len(operon1)+1):
        for b in range(1, len(operon2)+1):
            spScore = calculateSumOfPairsScore(operon1[a-1], operon2[b-1], '-')
            scoreMatrix[a][b][0], dirMatrix[a][b][0] = findMax(scoreMatrix, a, b, 0, operon1, operon2, operon3)

    for a in range(1, len(operon1)+1):
        for c in range(1, len(operon3)+1):
            spScore = calculateSumOfPairsScore(operon1[a-1], '-', operon3[c-1])
            scoreMatrix[a][0][c], dirMatrix[a][0][c] = findMax(scoreMatrix, a, 0, c, operon1, operon2, operon3)

    for b in range(1, len(operon2)+1):
        for c in range(1, len(operon3)+1):
            spScore = calculateSumOfPairsScore('-', operon2[b-1], operon3[c-1])
            scoreMatrix[0][b][c], dirMatrix[0][b][c] = findMax(scoreMatrix, 0, b, c, operon1, operon2, operon3)

    #perform the Global Alignment
    print ""
    for a in range(1, len(operon1)+1):
        for b in range(1, len(operon2)+1):
            for c in range(1, len(operon3)+1):
                # spScore = calculateSumOfPairsScore(operon1[a-1], operon2[b-1], operon3[c-1])
                print str(spScore) + " ",
                scoreMatrix[a][b][c], dirMatrix[a][b][c] = findMax(scoreMatrix, a, b, c, operon1, operon2, operon3)
                # max(scoreMatrix[a-1][b-1][c-1] + spScore, scoreMatrix[a-1][b-1][c] + spScore, scoreMatrix[a-1][b][c-1] + spScore, scoreMatrix[a][b-1][c-1] + spScore, scoreMatrix[a][b][c-1] + spScore, scoreMatrix[a][b-1][c] + spScore, scoreMatrix[a-1][b][c] + spScore)
            print ""
        print "\n"

    #Compute the number of events that occured between the operons
    # event = alignmentTraceback(scoreMatrix, dirMatrix, operon1, operon2, operon3, event)

    np.set_printoptions(precision=3, linewidth=np.inf)
    print scoreMatrix
    print dirMatrix
    print ""

    print operon1
    print operon2
    print operon3
    print ""

    traceback(scoreMatrix, dirMatrix, operon1, operon2, operon3)

directions = {'blu': 1, 'l': 5, 'u': 6, 'b': 7}

def findMax(scoreMatrix, a, b, c, operon1, operon2, operon3):
    maxScore = -999.0
    # Match in gene1,gene2,gene3
    if a != 0 and b != 0 and c != 0:
        spScore = calculateSumOfPairsScore(operon1[a-1], operon2[b-1], operon3[c-1])
        if maxScore < scoreMatrix[a-1][b-1][c-1] + spScore:
            maxScore = scoreMatrix[a-1][b-1][c-1] + spScore
            direction = 1
    # Match in gene1,gene2. Gap in gene 3
    if a != 0 and b != 0:
        spScore = calculateSumOfPairsScore(operon1[a-1], operon2[b-1], '-')
        if maxScore < scoreMatrix[a-1][b-1][c] + spScore:
            maxScore = scoreMatrix[a-1][b-1][c] + spScore
            direction = 2
    # Match in gene1,gene3. Gap in gene 2
    if a != 0 and c != 0:
        spScore = calculateSumOfPairsScore(operon1[a-1], '-', operon3[c-1])
        if maxScore < scoreMatrix[a-1][b][c-1] + spScore:
            maxScore = scoreMatrix[a-1][b][c-1] + spScore
            direction = 3
    # Match in gene2,gene3. Gap in gene 1
    if b != 0 and c != 0:
        spScore = calculateSumOfPairsScore('-', operon2[b-1], operon3[c-1])
        if maxScore < scoreMatrix[a][b-1][c-1] + spScore:
            maxScore = scoreMatrix[a][b-1][c-1] + spScore
            direction = 4
    # Gap in gene1, gene2
    if c != 0:
        spScore = calculateSumOfPairsScore('-', '-', operon3[c-1])
        if maxScore < scoreMatrix[a][b][c-1] + spScore:
            print str(scoreMatrix[a][b][c-1]) + "+"  + str(spScore) + " ",
            maxScore = scoreMatrix[a][b][c-1] + spScore
            direction = 5
    # Gap in gene1, gene3
    if b != 0:
        spScore = calculateSumOfPairsScore('-', operon2[b-1], '-')
        if maxScore < scoreMatrix[a][b-1][c] + spScore:
            maxScore = scoreMatrix[a][b-1][c] + spScore
            direction = 6
    # Gap in gene2, gene3
    if a != 0:
        spScore = calculateSumOfPairsScore(operon1[a-1], '-', '-')
        if maxScore < scoreMatrix[a-1][b][c] + spScore:
            maxScore = scoreMatrix[a-1][b][c] + spScore
            direction = 7
    if maxScore == -999.0:
        print "ERROR: No max score found. Direction was not chosen. Exiting..."
        sys.exit(0)

    return maxScore, direction

# def findMax2(scoreMatrix, a, b, c, spScore):
#     maxScore = -999.0
#     # Match in gene1,gene2,gene3
#     if a != 0 and b != 0 and c != 0 and maxScore < scoreMatrix[a-1][b-1][c-1] + spScore:
#         maxScore = scoreMatrix[a-1][b-1][c-1] + spScore
#         direction = 1
#     # Match in gene1,gene2. Gap in gene 3
#     if a != 0 and b != 0 and maxScore < scoreMatrix[a-1][b-1][c] + spScore:
#         maxScore = scoreMatrix[a-1][b-1][c] + spScore
#         direction = 2
#     # Match in gene1,gene3. Gap in gene 2
#     if a != 0 and c != 0 and maxScore < scoreMatrix[a-1][b][c-1] + spScore:
#         maxScore = scoreMatrix[a-1][b][c-1] + spScore
#         direction = 3
#     # Match in gene2,gene3. Gap in gene 1
#     if b != 0 and c != 0 and maxScore < scoreMatrix[a][b-1][c-1] + spScore:
#         maxScore = scoreMatrix[a][b-1][c-1] + spScore
#         direction = 4
#     # Gap in gene1, gene2
#     if c != 0 and maxScore < scoreMatrix[a][b][c-1] + spScore:
#         print str(scoreMatrix[a][b][c-1]) + "+"  + str(spScore) + " ",
#         maxScore = scoreMatrix[a][b][c-1] + spScore
#         direction = 5
#     # Gap in gene1, gene3
#     if b != 0 and maxScore < scoreMatrix[a][b-1][c] + spScore:
#         maxScore = scoreMatrix[a][b-1][c] + spScore
#         direction = 6
#     # Gap in gene2, gene3
#     if a != 0 and maxScore < scoreMatrix[a-1][b][c] + spScore:
#         maxScore = scoreMatrix[a-1][b][c] + spScore
#         direction = 7
#     if maxScore == -999.0:
#         print "ERROR: No max score found. Direction was not chosen. Exiting..."
#         sys.exit(0)

#     return maxScore, direction

def traceback(matrix, dirMatrix, operon1, operon2, operon3):
    i = len(operon1)
    j = len(operon2)
    k = len(operon3)

    #Track the alignment (two in the event we have substitutions)
    alignmentSequence1 = []
    alignmentSequence2 = []
    alignmentSequence3 = []

    while i > 0 or j > 0 and k > 0:
        print dirMatrix[i][j][k]
        #Case 1: Perfect match
        if i > 0 and j > 0 and (dirMatrix[i][j][k] == 1 or dirMatrix[i][j][k] == 2) and operon1[i-1] == operon2[j-1]:
            #Self global alignment
            print "Case1"
            alignmentSequence1.insert(0, operon1[i-1])
            alignmentSequence2.insert(0, operon2[j-1])

            if dirMatrix[i][j][k] == 1:
                alignmentSequence3.insert(0, operon3[k-1])
                k -= 1
            else:
                alignmentSequence3.insert(0, "-")

            i -= 1
            j -= 1
            
        #Case 2: Codon mismatch
        elif i > 0 and j > 0 and (dirMatrix[i][j][k] == 1 or dirMatrix[i][j][k] == 2) and operon1[i-1].split('_')[0].strip() == operon2[j-1].split('_')[0].strip():
            print "Case2"
            alignmentSequence1.insert(0, operon1[i-1])
            alignmentSequence2.insert(0, operon2[j-1])

            if dirMatrix[i][j][k] == 1:
                alignmentSequence3.insert(0, operon3[k-1])
                k -= 1
            else:
                alignmentSequence3.insert(0, "-")

            i -= 1
            j -= 1
            
        #Case 3: Substitution
        elif i > 0 and j > 0 and (dirMatrix[i][j][k] == 1 or dirMatrix[i][j][k] == 2):
            print "Case3"
            alignmentSequence1.insert(0, operon1[i-1])
            alignmentSequence2.insert(0, operon2[j-1])

            if dirMatrix[i][j][k] == 1:
                alignmentSequence3.insert(0, operon3[k-1])
                k -= 1
            else:
                alignmentSequence3.insert(0, "-")

            i -= 1
            j -= 1
            
        #Case 4: Mismatch- Gap in operon 2
        elif i > 0 and (dirMatrix[i][j][k] == 3 or dirMatrix[i][j][k] == 7):
            print "Case4"
            alignmentSequence1.insert(0, operon1[i-1])
            alignmentSequence2.insert(0, "-")

            if dirMatrix[i][j][k] == 3:
                alignmentSequence3.insert(0, operon3[k-1])
                k -= 1
            else:
                alignmentSequence3.insert(0, "-")

            i -= 1

        #Case 5: Mismatch - Gap in operon 1
        elif j > 0 and (dirMatrix[i][j][k] == 4 or dirMatrix[i][j][k] == 6):
            print "Case5"
            alignmentSequence2.insert(0, operon2[j-1])
            alignmentSequence1.insert(0, "-")

            if dirMatrix[i][j][k] == 4:
                alignmentSequence3.insert(0, operon3[k-1])
                k -= 1
            else:
                alignmentSequence3.insert(0, "-")

            j -= 1

        #Case 6: Mismatch - Gap in operon 1 and operon 2 (need to move operon3 index)
        else:
            print "Case6"
            if dirMatrix[i][j][k] == 5:
                alignmentSequence3.insert(0, operon3[k-1])
                k -= 1
                alignmentSequence1.insert(0, "-")
                alignmentSequence2.insert(0, "-")
            else:
                print "Error: Case 6 encountered but direction was not 5. Direction was " + str(dirMatrix[i][j][k])

    print [gene.center(7) for gene in alignmentSequence1]
    print [gene.center(7) for gene in alignmentSequence2]
    print [gene.center(7) for gene in alignmentSequence3]

######################################################
# globalAlignmentTraceback
# Parameters:
# Description: Performs a traceback on a given matrix
######################################################
def alignmentTraceback(matrix, dirMatrix, operon1, operon2, operon3, event):

    i = len(operon1)
    j = len(operon2)
    k = len(operon3)

    match = 0
    codonMismatch = 0
    mismatch = 0
    substitution = 0

    #Tracks index and genes for codon mismatches in both strains
    codonMismatchIndexesStrain1 = []
    codonMismatchIndexesStrain2 = []
    codonMismatchGenesStrain1 = []
    codonMismatchGenesStrain2 = []

    #Tracks substitution indexes and genes for both strains
    substitutionIndexesStrain1 = []
    substitutionIndexesStrain2 = []
    substitutionGenesStrain1 = []
    substitutionGenesStrain2 = []

    #Tracks the genes in a gap and the index of those genes in both strains note: details stored in an array of arrays
    operon1Gaps = []
    operon1Gap = []
    operon1GapIndexes = [] #This is used to determine the position of the genes with respect to the genome
    operon1GapIndex = []
    operon1ConsecutiveGap = False #Tracks consecutive gaps

    operon2Gaps = []
    operon2Gap = []
    operon2GapIndexes = []
    operon2GapIndex = []
    operon2ConsecutiveGap = False #Tracks consecutive gaps

    #Track the alignment (two in the event we have substitutions)
    alignmentSequence1 = []
    alignmentSequence2 = []

    #Tracks where the extra genes are from
    gap1Indexes = [] #This is used to determine where to insert the genes into the alignment
    gap2Indexes = []
    
    selfDuplication = '' #Used only for self global alignment
    selfPosition = event.fragmentDetails1.startPositionInGenome

    while i > 0 or j > 0:
        #Case 1: Perfect match
        if i > 0 and j > 0 and (dirMatrix[i][j][k] == 1 or dirMatrix[i][j][k] == 2) and operon1[i-1] == operon2[j-1]:
            #Self global alignment
            if event.fragmentDetails1.isNegativeOrientation == False:
                selfDuplication = operon2[j-1] + ' ' + str((i-1) + selfPosition) + ', ' + selfDuplication
            else:
                selfDuplication = selfDuplication + operon2[j-1] + ' ' + str(len(operon1) - (i-1) + selfPosition) + ', '
                
            match += 1
            alignmentSequence1.insert(0, operon1[i-1])
            alignmentSequence2.insert(0, operon2[j-1])
            i -= 1
            j -= 1
            operon1ConsecutiveGap = False
            operon2ConsecutiveGap = False

            if dirMatrix[i][j][k] == 1:
                k -= 1
            
        #Case 2: Codon mismatch
        elif i > 0 and j > 0 and (dirMatrix[i][j][k] == 1 or dirMatrix[i][j][k] == 2) and operon1[i-1].split('_')[0].strip() == operon2[j-1].split('_')[0].strip():
            #Self global alignment
            if event.fragmentDetails1.isNegativeOrientation == False:
                selfDuplication = '!' + operon2[j-1] + ' ' + str(-1) + ', ' + selfDuplication
            else:
                selfDuplication = selfDuplication + '!' + operon2[j-1] + ' ' + str(-1) + ', '
                
            #Increment the Id counter to ensure Id id unique
            globals.codonMismatchId += 1
            
            codonMismatch += 1

            alignmentSequence1.insert(0, operon1[i-1] + '-#' + str(globals.codonMismatchId) + '#')
            alignmentSequence2.insert(0, operon2[j-1] + '-#' + str(globals.codonMismatchId) + '#')

            codonMismatchIndexesStrain1.append(i-1)
            codonMismatchGenesStrain1.append(operon1[i-1] + '-#' + str(globals.codonMismatchId) + '#')

            codonMismatchIndexesStrain2.append(j-1)
            codonMismatchGenesStrain2.append(operon2[j-1] + '-#' + str(globals.codonMismatchId) + '#')

            i -= 1
            j -= 1
            operon1ConsecutiveGap = False
            operon2ConsecutiveGap = False

            if dirMatrix[i][j][k] == 1:
                k -= 1
            
        #Case 3: Substitution
        elif i > 0 and j > 0 and (dirMatrix[i][j][k] == 1 or dirMatrix[i][j][k] == 2):
            #Self global alignment
            if event.fragmentDetails1.isNegativeOrientation == False:
                selfDuplication = '!' + operon2[j-1] + ' ' + str(-1) + ', ' + selfDuplication
            else:
                selfDuplication = selfDuplication + '!' + operon2[j-1] + ' ' + str(-1) + ', '
                
            #Increment the Id counter to ensure the ID is unique
            globals.substitutionId += 1
            
            substitution += 1
            
            alignmentSequence1.insert(0, operon1[i-1] + '-@' + str(globals.substitutionId) + '@')
            alignmentSequence2.insert(0, operon2[j-1] + '-@' + str(globals.substitutionId) + '@')

            substitutionIndexesStrain1.append(i-1)
            substitutionGenesStrain1.append(operon1[i-1] + '-@' + str(globals.substitutionId) + '@')

            substitutionIndexesStrain2.append(j-1)
            substitutionGenesStrain2.append(operon2[j-1] + '-@' + str(globals.substitutionId) + '@')

            i -= 1
            j -= 1
            operon1ConsecutiveGap = False
            operon2ConsecutiveGap = False

            if dirMatrix[i][j][k] == 1:
                k -= 1
            
        #Case 4: Mismatch- Gap in operon 2
        elif i > 0 and (dirMatrix[i][j][k] == 3 or dirMatrix[i][j][k] == 7):
            index = i-1
            mismatch += 1
            i -= 1
            operon1ConsecutiveGap = False
            #Check if this is a consecutive gap, if it is then append to the gap list if not then append to the list of gaps and start a new gap
            if operon2ConsecutiveGap:
                operon2Gap.insert(0, operon1[index])
                operon2GapIndex.insert(0, index)

                operon2ConsecutiveGap = True
            else:
                if len(operon2Gap) > 0:
                    operon2Gaps.insert(0, operon2Gap)
                    operon2GapIndexes.insert(0, operon2GapIndex)
                operon2Gap = []
                operon2GapIndex = []
                operon2Gap.insert(0, operon1[index])
                operon2GapIndex.insert(0, index)
                gap2Indexes.insert(0, len(alignmentSequence2))
                operon2ConsecutiveGap = True

            if dirMatrix[i][j][k] == 3:
                k -= 1

        #Case 5: Mismatch - Gap in operon 1
        elif i > 0 and (dirMatrix[i][j][k] == 4 or dirMatrix[i][j][k] == 6):
            #Self global alignment
            if event.fragmentDetails1.isNegativeOrientation == False:
                selfDuplication = '!' + operon2[j-1] + ' ' + str(-1) + ', ' + selfDuplication
            else:
                selfDuplication = selfDuplication + operon2[j-1] + ' ' + str(-1) + ', '
                
            index = j - 1
            mismatch += 1
            j -= 1
            operon2ConsecutiveGap = False
            #Check if this is a consecutive gap, if it is then append to the gap list if not then append to the list of gaps and start a new gap
            if operon1ConsecutiveGap:
                operon1Gap.insert(0, operon2[index])
                operon1GapIndex.insert(0, index)

                operon1ConsecutiveGap = True
            else:
                if len(operon1Gap) > 0:
                    operon1Gaps.insert(0, operon1Gap)
                    operon1GapIndexes.insert(0, operon1GapIndex)
                operon1Gap = []
                operon1GapIndex = []
                operon1Gap.insert(0, operon2[index])
                operon1GapIndex.insert(0, index)
                gap1Indexes.insert(0, len(alignmentSequence1))
                operon1ConsecutiveGap = True

            if dirMatrix[i][j][k] == 4:
                k -= 1

        #Case 6: Mismatch - Gap in operon 1 and operon 2 (need to move operon3 index)
        else:
            print "Last case"
            if dirMatrix[i][j][k] == 5:
                print "Direction is 5"
                k -= 1
    
    event.selfDuplication = selfDuplication[0:(len(selfDuplication) - 2)] + ';' #Remove the last comma and space and add a semicolon 
    #Empty any remaining gaps
    if len(operon1Gap) > 0:
        operon1Gaps.insert(0, operon1Gap)
        operon1GapIndexes.insert(0, operon1GapIndex)
        operon1Gap = []
        operon1GapIndex = []

    if len(operon2Gap) > 0:
        operon2Gaps.insert(0, operon2Gap)
        operon2GapIndexes.insert(0, operon2GapIndex)
        operon2Gap = []
        operon2GapIndex = []

    #The indexes values need to be flipped b/c right now they're oriented from right to left
    if len(gap1Indexes) > 0:
        for x in range(0, len(gap1Indexes)):
            gap1Indexes[x] = len(alignmentSequence1) - gap1Indexes[x]
    if len(gap2Indexes) > 0:
        for x in range(0, len(gap2Indexes)):
            gap2Indexes[x] = len(alignmentSequence2) - gap2Indexes[x]

    #Need to swap the gap lists since the gaps refer to extra genes
    temp = operon1Gaps
    operon1Gaps = operon2Gaps
    operon2Gaps = temp

    temp = operon1GapIndexes
    operon1GapIndexes = operon2GapIndexes
    operon2GapIndexes = temp

    temp = gap1Indexes
    gap1Indexes = gap2Indexes
    gap2Indexes = temp

    #Match, mismatch details
    event.setNumMatches(match)
    event.setNumMismatches(mismatch)

    #Codon Mismatch details
    event.setNumCodonMismatches(codonMismatch)
    event.setCodonMismatchGenesStrain1(codonMismatchGenesStrain1)
    event.setCodonMismatchGenesStrain2(codonMismatchGenesStrain2)
    event.setCodonMismatchIndexesStrain1(codonMismatchIndexesStrain1)
    event.setCodonMismatchIndexesStrain2(codonMismatchIndexesStrain2)

    #Substitution details
    event.setNumSubstitutions(substitution)
    event.setSubstitutionIndexesStrain1(substitutionIndexesStrain1)
    event.setSubstitutionIndexesStrain2(substitutionIndexesStrain2)
    event.setSubstitutionGenesStrain1(substitutionGenesStrain1)
    event.setSubstitutionGenesStrain2(substitutionGenesStrain2)

    #Alignment details
    event.setOperon1Alignment(alignmentSequence1)
    event.setOperon2Alignment(alignmentSequence2)

    #Gap details
    event.setOperon1Gaps(operon1Gaps)
    event.setOperon2Gaps(operon2Gaps)
    event.setOperon1GapPositions(operon1GapIndexes)
    event.setOperon2GapPositions(operon2GapIndexes)
    event.setOperon1GapIndexes(gap1Indexes)
    event.setOperon2GapIndexes(gap2Indexes)

    #Used for debugging
    #print('These are the operons being compared: %s, %s' %(operon1, operon2))
    #print('This is the resulting alignment: %s, %s' %(alignmentSequence1, alignmentSequence2))
    #print('These are the extra genes for operon 1: %s' %(operon1Gaps))
    #print('These are the indexes for extra genes in operon 1: %s' %(gap1Indexes))
    #print('These are the extra genes for operon 2: %s' %(operon2Gaps))
    #print('These are the indexes for extra genes in operon 2: %s' %(gap2Indexes))

    return event

def calculateSumOfPairsScore(gene1, gene2, gene3):
    sopScore = 0.0
    numPairings = 3.0

    sopScore += compareGenes(gene1, gene2)
    sopScore += compareGenes(gene1, gene3)
    sopScore += compareGenes(gene2, gene3)

    return sopScore/numPairings

def compareGenes(gene1, gene2):
    score = 0.0

    if gene1 != '-' and gene2 != '-':
        if gene1.split('_')[0].strip() == gene2.split('_')[0].strip():
            if gene1.strip() == gene2.strip():
                score = globals.match
            else:
                score = globals.codonCost
        else:
            score = globals.substitutionCost
    elif gene1 == '-' and gene2 == '-':
        # CHANGE THIS TO PROPER SCORE
        score = 0
    else:
        score = globals.deletionCost

    return score

def main():
    globals.initialize() #Initialize the globals file

    op1 = ["16S", "Ile_AUC", "Ala_GCA", "23S", "5S"]
    op2 = ["16S", "Ala_GCA", "23S", "5S"]
    op3 = ["16S", "Gln_CAA", "Ala_GCA", "23S", "5S", "16S"]

    performMultiSequenceAlignment(op1, op2, op3)

main()
