# Submitted by - Mukund Aggarwal and Sarmistha Subhadarshini

from fractions import Fraction;

class GomoryCut:
    def __init__(self, A, b, c):
        self.tableau = self.makeTableau(A, b, c)
        self.rows = len(self.tableau)
        self.origRows = len(self.tableau)
        self.origCols = 0;
        self.cols = len(self.tableau[0])
        self.currentBFS = [-1] * (self.rows - 1)  # x_i corresponding to the rows of the equation
        self.finalRes = [0] * (self.rows)

    def makeTableau(self, A, b, c):
        tabs = [[0] + c] + A
        for i in range(1, len(tabs)):
            tabs[i] = [b[i - 1]] + tabs[i]
        return tabs

    def DualSimplex(self):
        while any(self.tableau[m][0] < 0 for m in range(1, self.rows)):
            for i in range(1, self.rows):
                if self.tableau[i][0] < 0:
                    minIndx = 1
                    for j in range(1, self.cols):
                        minIndx = j if ((self.tableau[i][j] < self.tableau[i][minIndx] and self.tableau[i][j] < 0) or
                                        self.tableau[i][minIndx] == 0) else minIndx
                    self.currentBFS[i - 1] = minIndx
                    for j in range(0, self.rows):
                        if i == j:
                            ref = self.tableau[j][minIndx]
                            for k in range(0, self.cols):
                                self.tableau[j][k] =  Fraction(self.tableau[j][k], ref);
                        else:
                            ref = Fraction(self.tableau[j][minIndx], self.tableau[i][minIndx]);
                            for k in range(0, self.cols):
                                self.tableau[j][k] -= (ref) * self.tableau[i][k]
                    break

    def pivotrowcol(self):
        cp = -1
        for i in range(1, self.cols):
            if self.tableau[0][i] < 0:
                if cp == -1:
                    cp = i
                elif self.tableau[0][cp] > self.tableau[0][i]:
                    cp = i;
        a = []
        for i in range(1, self.rows):
            u = self.tableau[i][cp]
            if u > 0:
                a.append([Fraction(self.tableau[i][0], u), i])
        if len(a) == 0:
            return
        p = min(a)
        rp = p[1]  # p[1] is pivot row
        self.currentBFS[p[1] - 1] = cp
        factor = self.tableau[rp][cp]
        for i in range(0, self.cols):
            self.tableau[rp][i] = Fraction(self.tableau[rp][i], factor)

        for i in range(0, self.rows):
            factor1 = self.tableau[i][cp]
            for j in range(0, self.cols):
                factor2 = self.tableau[rp][j]
                if i != rp:
                    self.tableau[i][j] = self.tableau[i][j] - factor2 * factor1
        self.tableau[rp][cp] = 1

    def Simplex(self):
        while any(self.tableau[0][m] < 0 for m in range(1, self.cols)):
            self.pivotrowcol()

    def GomoryCut(self):
        def simplexSolver():
            while (any(self.tableau[0][m] < 0 for m in range(1, self.cols)) or any(
                    self.tableau[m][0] < 0 for m in range(1, self.rows))):
                self.Simplex()
                self.DualSimplex()

        simplexSolver()
        # First we arrive at the Reduced Simplex Tableau, and then we add the fractional constraint(s)
        while any(self.tableau[i][0] != int(self.tableau[i][0]) for i in range(1, self.origRows)):
            minFracIndx = -1
            for i in range(1, self.rows):
                try:
                    if self.tableau[i][0].denominator != 1: #self.tableau[i][0] != int(self.tableau[i][0]):
                        minFracIndx = i
                        break
                except:
                    if self.tableau[i][0] != int(self.tableau[i][0]):
                        minFracIndx = i
                        break

            if minFracIndx == -1: break

            # Append new rows/constraints to the tableau and current BFS solns
            self.tableau.append([Fraction(0)] * self.cols)
            self.rows += 1
            for i in range(0, self.rows):
                self.tableau[i].append(Fraction(0))
            self.cols += 1
            self.currentBFS.append(Fraction(0))

            for i in range(0, self.cols):
                if self.tableau[minFracIndx][i] >= 0:
                    self.tableau[-1][i] = -1 * (self.tableau[minFracIndx][i] - int(self.tableau[minFracIndx][i]))
                else:
                    self.tableau[-1][i] = -1 * (self.tableau[minFracIndx][i] - int(self.tableau[minFracIndx][i] - 1))

            self.tableau[-1][-1] = 1

            simplexSolver()

        for i in range(len(self.currentBFS)):
            if self.currentBFS[i] <= len(self.finalRes) and self.currentBFS[i]!=-1:
                self.finalRes[self.currentBFS[i] - 1] = float(self.tableau[i + 1][0])
        self.finalRes = self.finalRes[0 : self.origCols];


def add_slacks(A, b, c):
    n = len(A)
    for i in range(len(c)):
        c[i] = Fraction(-c[i])

    c = c + [Fraction(0)] * n  # adding the coefficients of slack variable to the c
    I = []
    for i in range(n):
        e = [Fraction(0) for j in range(n)]
        I.append(e)
    for i in range(n):
        I[i][i] = Fraction(1)
    for i in range(n):
        A[i] = A[i] + I[i]  # adding the identity matrix
    return A, b, c


def gomory(filename: str):
    with open(filename, 'r') as fl:
        cols, rows = [int(i) for i in fl.readline().split(" ")]
        b = [Fraction(i) for i in fl.readline().split(" ")]
        c = [Fraction(i) for i in fl.readline().split(" ")]
        A = []
        for i in range(rows):
            A.append([Fraction(i) for i in fl.readline().split()])
        fl.close()
    A, b, c = add_slacks(A, b, c)
    obj = GomoryCut(A, b, c)
    obj.origCols = cols;
    obj.GomoryCut()
    k = obj.finalRes
    if len(k) < cols:
        for i in range(0, cols - len(k)): k.append(0)
    return k