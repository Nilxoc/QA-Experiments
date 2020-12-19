"""Skript zur Abbildung des N-Damen Problems auf ein QUBO-Problem 
sowie Ausführung und Interpretation der Ergebnisse mit 
D-Wave Quanten Annealer oder Softwarelösungen
"""
import dimod
import math
import pickle
import minorminer

from dimod import BinaryQuadraticModel
from dwave.system import DWaveSampler
from dwave.system.composites import FixedEmbeddingComposite

#Indizes starten bei 1(nicht 0)
def inCol(i,j,n):
    """Gibt an, ob sich die Felder i und j auf einem n*n Schachbrett in der 
    selben Spalte befinden.
    """
    return (i%n)==(j%n)

def inRow(i,j,n):
    """Gibt an, ob sich die Felder i und j auf einem n*n Schachbrett in der 
    selben Zeile befinden.
    """
    return math.floor((i-1)/n) == math.floor((j-1)/n)

def inDiag(i,j,n):
    """Gibt an, ob die Felder i und j auf einem n*n Schachbrett sich eine
    Diagonale teilen
    """
    return abs((math.floor((i-1)/n)+1)-(math.floor((j-1)/n)+1)) == abs(((i-1)%n)-((j-1)%n))

def initMat(n):
    """Generiert die Hamiltonian-Matrix für n Damen
    """
    nsquared = n**2
    result = BinaryQuadraticModel(dimod.Vartype.BINARY)
    for i in range(1,nsquared+1):
        result.linear['x'+str(i)] = -1.0
        for j in range(i+1,nsquared+1):
            if inCol(i,j,n)\
            or inRow(i,j,n)\
            or inDiag(i,j,n):
                result.quadratic[('x'+str(i), 'x'+str(j))] = 1.0
    return result

def printRowDivider(n):
    """Hilfsfunktion für printField, die eine Trennung zwischen Zeilen ausgibt.
    """
    print('\nI', end='')
    for j in range(0, n):
        print('---I', end='')

def printHamiltonian(bqm,n):
    """Gibt das gegebene dimod.BinaryQuadraticModel als Diagonalmatrix aus
    """
    for i in range(1, n**2+1):
        for j in range(0, i-1):
            print(str(0.0).rjust(4),end=' ')
        for j in range(i, n**2+1):
            if i != j:
                key = ('x'+str(i),'x'+str(j))
                if key in bqm.quadratic:
                    print(str(bqm.quadratic[key]).rjust(4),end=' ')
                else:
                    print(str(0.0).rjust(4),end = ' ')
            else:
                print(str(bqm.linear['x'+str(i)]).rjust(4), end=' ')
        print('\n')

def printField(sample):
    """Gibt eine visuelle Darstellung des Spielfeldes wie es im
    gegebenen sample beschrieben ist auf dem terminal aus
    """
    n = math.floor(math.sqrt(len(sample)))
    printRowDivider(n)
    for i in range(0, n):
        print('\nI ',end='')
        for j in range(0, n):
            val = ' '
            if sample['x'+str(i*n+j+1)] == 1:
                val = 'D'
            print(val+' I ', end='')
        printRowDivider(n)
    print('\n',end='')

def qAnnealNQueens(n, num_reads):
    """Löst das n-Damen Problem auf dem in der Umgebung konfigurierten
    Quanten Annealer
    """
    bqm = initMat(n)
    Q = bqm.to_qubo()[0]

    solver = DWaveSampler()
    #Einbettung auf Chimera Graph(Struktur der Qubits)
    __, target_edges, target_adjacency = solver.structure
    embedding = minorminer.find_embedding(Q, target_edges) #Heuristisch
    sampler = FixedEmbeddingComposite(solver, embedding) #Composite übernimmt embedding und unembedding

    return sampler.sample_qubo(Q,num_reads=num_reads)

if __name__ == '__main__': #nur wenn als skript aufgerufen(nicht als import)
    print("WARNING: This will use computation time on the configured Quantum Annealer")

    print("Number of Queens: ",end='')
    nQueens = 0
    try:
        nQueens = int(input())
    except ValueError:
        print('ERROR: Number of Queens must be an Integer')
        exit()

    print('Number of Reads: ',end='')
    nReads = 0
    try:
        nReads = int(input())
    except ValueError:
        print('ERROR: Number of Reads must be an Integer')
        exit()

    sampleset = qAnnealNQueens(nQueens, nReads);
    
    print('Best Solution:')
    printField(sampleset.first.sample)
    if(sampleset.first.energy == -nQueens):
        print('Solution is correct')
    else:
        print('Solution is NOT correct')
    
    doneSaving = False;
    while(not doneSaving):
        print('Enter file name for sample set(Leave blank to skip): ', end='')
        filename = input()
        if(filename == ''):
            print('Skipping save')
        else:
            try:
                of = open(filename, mode='wb')
                pickle.dump(sampleset, of)
                doneSaving = True
            except OSError:
                print('Unable to open file ' + filename + ' for writing')
