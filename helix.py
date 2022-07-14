monomers_start = 20
monomers_step = 1
monomers_count = 1
hydrogenBetaE_start = 1.7
hydrogenBetaE_step = 1.0
hydrogenBetaE_count = 1
peptydeBetaE_start = 15.0
peptydeBetaE_step = 1.0
peptydeBetaE_count = 1
covalentBetaE_start = 100.0
covalentBetaE_step = 1.0
covalentBetaE_count = 1
helixFrom_start = 3
helixFrom_step = 1
helixFrom_count = 1
helixPeriod_start = 3
helixPeriod_step = 1
helixPeriod_count = 1
helixCycleLength_start = 10
helixCycleLength_step = 1
helixCycleLength_count = 1

from random import randint
from math import log, exp

def linspace(start, stop, n):
    if n == 1:
        yield stop
        return
    h = (stop - start) / (n - 1)
    for i in range(n):
        yield start + h * i

eps = 0.00001
INFINITY = 1000000000000000

class Matrix:

    def __init__(self, n, y):
        self.n = n
        self.cntE = len(y)
        self.visited = [False] * n
        self.yalls = y
        self.m = [[0] * n for _ in range(n)]
        self.yalha = [[] for _ in range(n)]
        self.targety = [[] for _ in range(n)]
        self.lastI = -1
        self.q = set()
        for i in range(self.cntE):
            self.targety[int(y[i] / n)] += [y[i] % n]
            self.targety[int(y[i] % n)] += [y[i] / n]
            self.swit(y[i])
        self.E = self.cntE
        self.V = n
        self.A = 0
        self.calcL()
        self.calcB()
        self.totalB = self.B
        self.makeRandomPos()

    def isBridge(self, x):
        return self.hasE(x) and not self.connected(x)

    def calcB(self):
        cnt = 0
        for i in range(self.n):
            for j in self.yalha[i]:
                if (i < self.getNode(j, i) and self.isBridge(j)):
                    cnt+=1
        self.B = cnt + self.calcL()

    def connected(self, num):
        a = int(num / self.n)
        b = int(num % self.n)
        xx = self.V
        cntA = len(self.yalha[a]) - 1
        cntB = len(self.yalha[b]) - 1
        self.V -= (1 if cntA == 0 else 0) + (1 if cntB == 0 else 0)
        if self.V < 4:
            self.V = xx
            return True
        if self.E < self.V:
            self.V = xx
            return False
        self.visited = [False] * self.n
        cnt = 0
        for i in range(self.n):
            if len(self.yalha[i]) - (1 if i == a or i == b else 0) > 0:
                cnt+=1
                self.visited[i] = True
                self.q.add(i)
                break
        while self.q:
            x = self.q.pop()
            for o in self.yalha[x]:
                if o != num:
                    z = int(self.getNode(o, x))
                    if not self.visited[z]:
                        cnt+=1
                        self.visited[z] = True
                        self.q.add(z)
        y = cnt == self.V
        self.V = xx
        return y

    def getNode(self, x, i):
        return (x - i * self.n) if x - i * self.n > i else ((x - i) / self.n)

    def calcL(self):
        cnt = 0
        for i in range(self.n):
            if len(self.yalha[i]) == 1:
                cnt+=1
        self.L = 1 if self.E == 1 else cnt
        return self.L

    def shift(self, i, j):
        self.swit(i * self.n + j)
        self.calcIt(i, j)
        self.calcB()

    def calcIt(self, a, b):
        self.E += 1 if self.m[a][b] else -1
        cntA = len(self.yalha[a])
        cntB = len(self.yalha[b])
        cntA1 = len(self.targety[a])
        cntB1 = len(self.targety[b])
        x = 0
        if self.m[a][b]:
            x += (1 if cntA == 1 else 0) + (1 if cntB == 1 else 0)
        else:
            x -= (0 if cntA else 1) + (0 if cntB else 1)
        self.V += x
        if x == 0:
            if self.m[a][b] == 1:
                self.A -= 1
            else:
                self.A += 1
        elif x == 1:
            if cntA == 1:
                self.A += self.calcA(a) - 1
            else:
                self.A += self.calcA(b) - 1
        elif x == 2:
            self.A += cntA1 + cntB1 - 2
        elif x == -1:
            if cntA:
                self.A -= self.calcA(b) - 1
            else:
                self.A -= self.calcA(a) - 1
        elif x == -2:
            self.A = 0

    def calcA(self, a):
        cnt = 0
        for i in self.targety[a]:
            cnt += 0 if len(self.yalha[int(i)]) else 1
        return cnt

    def makeRandomPos(self):
        self.makeEmpty()
        for _ in range(self.cntE ** 2):
            self.step()

    def makeEmpty(self):
        for i in range(self.n):
            while self.yalha[i]:
                self.swit(self.yalha[i][0])
        self.A = self.B = self.L = self.V = self.E = 0

    def swit(self, x):
        if self.hasE(x):
            self.yalha[int(x / self.n)].remove(x)
            self.yalha[int(x % self.n)].remove(x)
        else:
            self.yalha[int(x / self.n)]+=[x]
            self.yalha[int(x % self.n)]+=[x]
        self.m[int(x / self.n)][int(x % self.n)] = 1 - self.m[int(x / self.n)][int(x % self.n)]
        self.m[int(x % self.n)][int(x / self.n)] = 1 - self.m[int(x % self.n)][int(x / self.n)]

    def getStrongE(self):
        cnt = 0
        for i in range(self.n-1):
            if self.m[i][i+1]:
                cnt+=1
        return cnt

    def getCovPepEs(self):
        cov = 0
        pep = 0
        for i in range(self.n - 1):
            if self.m[i][i+1]:
                if i%3 == 0:
                    cov+=1
                else:
                    pep+=1
        return cov, pep

    def getE(self):
        return self.E

    def getnp(self):
        return self.A

    def getnn(self):
        return self.E - self.B + self.L

    def step(self):
        self.lastI = self.getRand()
        self.shift(int(self.lastI / self.n), int(self.lastI % self.n))

    def stepBack(self):
        self.shift(int(self.lastI / self.n), int(self.lastI % self.n))

    def hasE(self, i):
        return self.m[int(i / self.n)][int(i % self.n)] != 0

    def canPut(self, x):
        if not self.E:
            return True
        i = int(x / self.n)
        j = int(x % self.n)
        return self.yalha[i] or self.yalha[j]

    def getRand(self):
        z = randint(0, self.cntE - 1)
        x = self.yalls[z]
        while (self.hasE(x) and self.getE() == 1 or not self.hasE(x) and self.getE() == self.cntE - 1 or self.isBridge(x) or not self.hasE(x) and not self.canPut(x)):
            z = randint(0, self.cntE - 1)
            x = self.yalls[z]
        return x

    def getQ(self):
        return self.E * 2.0 / self.V

    def getV(self):
        return self.V

    def getB(self):
        return self.B

    def getTotalB(self):
        return self.totalB

def getNum(a, b, n):
    a-=1
    b-=1
    if a > b:
        b += a
        a = b - a
        b -= a
    return a * n + b

def helix(_from, period, cycle, n):
    anslist = []
    f = [0] * n
    for i in range(_from-1, n):
        f[i] = cycle * (0 if (i+1-_from)%period else 1)
    print("%" + str(f)[1:-1])
    for i in range(n-1):
        index = i + 1
        anslist+=[getNum(index, index + 1, n)]
        if i + f[i] < n and f[i] > 0:
            anslist+=[getNum(index, index + f[i], n)]
    return anslist

def customGraph():
    anslist = []
    n=13
    anslist+=[getNum(1, 2, n)]
    anslist+=[getNum(1, 3, n)]
    anslist+=[getNum(1, 4, n)]
    anslist+=[getNum(1, 5, n)]
    anslist+=[getNum(1, 7, n)]
    anslist+=[getNum(2, 6, n)]
    anslist+=[getNum(2, 7, n)]
    anslist+=[getNum(3, 6, n)]
    anslist+=[getNum(3, 7, n)]
    anslist+=[getNum(4, 6, n)]
    anslist+=[getNum(4, 7, n)]
    anslist+=[getNum(5, 6, n)]
    anslist+=[getNum(5, 7, n)]
    anslist+=[getNum(6, 7, n)]
    anslist+=[getNum(7, 8, n)]
    anslist+=[getNum(7, 9, n)]
    anslist+=[getNum(7, 10, n)]
    anslist+=[getNum(7, 11, n)]
    anslist+=[getNum(7, 12, n)]
    anslist+=[getNum(7, 13, n)]
    anslist+=[getNum(8, 9, n)]
    anslist+=[getNum(8, 10, n)]
    anslist+=[getNum(8, 11, n)]
    anslist+=[getNum(8, 12, n)]
    anslist+=[getNum(9, 13, n)]
    anslist+=[getNum(10, 13, n)]
    anslist+=[getNum(11, 13, n)]
    anslist+=[getNum(12, 13, n)]
    return anslist


class MonteCarlo:

    def __init__(self, n, btmph, btmpp, btmpc, t):
        self.n = n
        self.btmph = btmph
        self.btmpp = btmpp
        self.btmpc = btmpc
        self.t = t
        self.e = len(t)
        self.lnhEV = [[0.0] * (self.n+1) for _ in range(self.e+1)]
        self.bEV = [[0.0] * (self.n+1) for _ in range(self.e+1)]
        self.qEV = [[0.0] * (self.n+1) for _ in range(self.e+1)]
        self.ehEV = [[0.0] * (self.n+1) for _ in range(self.e+1)]
        self.epEV = [[0.0] * (self.n+1) for _ in range(self.e+1)]
        self.ecEV = [[0.0] * (self.n+1) for _ in range(self.e+1)]


    def run(self):
        MCsteps = 100
        skiptohist = 10
        lnf = 1
        lnfmin = 0.0001
        mat = Matrix(self.n, self.t)
        hist2 = [[0] * (self.n+1) for _ in range(self.e+1)]
        totalnum = [[0] * (self.n+1) for _ in range(self.e+1)]
        state0 = [0] * 4
        state1 = [0] * 4
        while lnf > lnfmin:
            MCcounter = 1
            skipcounter = skiptohist
            while MCcounter:
                sigmahist = 0
                for i in range(self.n ** 2):
                    state0[0] = mat.getE()
                    state0[1] = mat.getV()
                    state0[2] = mat.getnp()
                    state0[3] = mat.getnn()
                    mat.step()
                    state1[0] = mat.getE()
                    state1[1] = mat.getV()
                    state1[2] = mat.getnp()
                    state1[3] = mat.getnn()
                    r = randint(0, 1000) / 1000.0
                    lnre = exp(self.lnhEV[state0[0]][state0[1]] - self.lnhEV[state1[0]][state1[1]])
                    p1 = (state0[3] / state1[2]) if state1[2] else 1
                    p2 = state0[2] / state1[3]
                    p = p1 if p1 < p2 else p2
                    if (r >= lnre * p):
                        mat.stepBack()
                    self.lnhEV[mat.getE()][mat.getV()] += lnf
                    self.bEV[mat.getE()][mat.getV()] += mat.getB()
                    self.qEV[mat.getE()][mat.getV()] += mat.getQ()
                    cov, pep = mat.getCovPepEs()
                    self.epEV[mat.getE()][mat.getV()] += pep
                    self.ecEV[mat.getE()][mat.getV()] += cov
                    self.ehEV[mat.getE()][mat.getV()] += mat.getE() - mat.getStrongE()
                    hist2[mat.getE()][mat.getV()]+=1
                    totalnum[mat.getE()][mat.getV()]+=1
                    sigmahist+=1
                MCcounter+=1
                skipcounter+=1
                if MCcounter > MCsteps and skipcounter > skiptohist:
                    hmin2 = INFINITY
                    lnhmin2 = INFINITY
                    Nstates2 = 0
                    for i in range(self.e + 1):
                        for j in range(self.n + 1):
                            if hist2[i][j]:
                                Nstates2+=1
                            if hist2[i][j] and hist2[i][j] < hmin2:
                                hmin2 = hist2[i][j]
                            if self.lnhEV[i][j] > eps and self.lnhEV[i][j] < lnhmin2:
                                lnhmin2 = self.lnhEV[i][j]
                    havg2 = sigmahist / Nstates2
                    if hmin2 / havg2 >= 0.8:
                        for i in range(self.e + 1):
                            for j in range(self.n + 1):
                                hist2[i][j] = 0
                                if self.lnhEV[i][j] > eps:
                                    self.lnhEV[i][j] = self.lnhEV[i][j] - lnhmin2 + log(2)
                        lnf /= 2
                        MCcounter = 0
                    else:
                        skipcounter = 1
        self.bEV[self.e][self.n] = mat.getTotalB()
        totalnum[self.e][self.n] = 1
        for i in range(self.e + 1):
            for j in range(self.n + 1):
                if totalnum[i][j]:
                    self.bEV[i][j] /= totalnum[i][j]
                    self.qEV[i][j] /= totalnum[i][j]
                    self.ehEV[i][j] /= totalnum[i][j]
                    self.epEV[i][j] /= totalnum[i][j]
                    self.ecEV[i][j] /= totalnum[i][j]
                else:
                    self.bEV[i][j] = None
                    self.qEV[i][j] = None
                    self.ehEV[i][j] = None
                    self.epEV[i][j] = None
                    self.ecEV[i][j] = None

        ans = [0.0] * (self.n + 1)
        for j in range(2,self.n+1):
            ans[j] = self.getBetaF(j) - self.getBetaF(2)
            print(str(j) + "\t\t" + str(ans[j]))

    def getLnZ(self, e, v):
        lnqc = log(self.qEV[e][v])
        lnqd = lnqc
        lnro = log(62.0 ** -3)
        return (self.ecEV[e][v] * self.btmpc + self.ehEV[e][v] * self.btmph + self.epEV[e][v] * self.btmpp + v * lnro - (v - 1) * lnqc - (v - self.bEV[e][v] - 1) * lnqd)

    def getBetaF(self, v):
        ans = 0
        for i in range(self.e+1):
            if ((self.lnhEV[i][v] < eps and (i != self.e or v != self.n)) or self.bEV[i][v] == None or self.qEV[i][v] == 0):
                continue
            lnZ = self.getLnZ(i, v)
            x = ans
            y = self.lnhEV[i][v] + lnZ
            x = log(exp(y) + (exp(x) if (x > eps or x < -eps) else 0))
            ans = x
        return -ans


graph = "mamooli"
if graph == "gerafe gerd":
    MonteCarlo(13, 7.8, 7.8, 7.8, customGraph()).run()
if graph == "mamooli":
    print("data={")
    for monomers_i in range(monomers_start, monomers_start + monomers_step * monomers_count, monomers_step):
        for helixFrom_i in range(helixFrom_start, helixFrom_start + helixFrom_step * helixFrom_count, helixFrom_step):
            for helixPeriod_i in range(helixPeriod_start ,helixPeriod_start + helixPeriod_step * helixPeriod_count, helixPeriod_step):
                for helixCycleLength_i in range(helixCycleLength_start, helixCycleLength_start + helixCycleLength_step * helixCycleLength_count, helixCycleLength_step):
                    for peptydeBetaE_i in linspace(peptydeBetaE_start, peptydeBetaE_start + peptydeBetaE_step * (peptydeBetaE_count - 1), peptydeBetaE_count):
                        for hydrogenBetaE_i in linspace(hydrogenBetaE_start, hydrogenBetaE_start + hydrogenBetaE_step * (hydrogenBetaE_count - 1), hydrogenBetaE_count):
                            for covalentBetaE_i in linspace(covalentBetaE_start, covalentBetaE_start + covalentBetaE_step * (covalentBetaE_count - 1), covalentBetaE_count):
                                alg = MonteCarlo(monomers_i, hydrogenBetaE_i, peptydeBetaE_i, covalentBetaE_i, helix(helixFrom_i, helixPeriod_i, helixCycleLength_i, monomers_i))
                                print("% n: " + str(monomers_i) + ", betaE_h: " + str(hydrogenBetaE_i) + ", betaE_p: " + str(peptydeBetaE_i) + ", betaE_c: " + str(covalentBetaE_i))
                                print("[")
                                alg.run()
                                print("],")
    print("};")