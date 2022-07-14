import java.util.*;

import static java.lang.Math.*;

public class Main {

    private static int monomers_start = 10;
    private static int monomers_step = 1;
    private static int monomers_count = 1;
    private static double hydrogenBetaE_start = 1.7;
    private static double hydrogenBetaE_step = 1.0;
    private static int hydrogenBetaE_count = 1;
    private static double peptydeBetaE_start = 15.0;
    private static double peptydeBetaE_step = 1.0;
    private static int peptydeBetaE_count = 1;
    private static double thirdBetaE_start = 0.0;
    private static double thirdBetaE_step = 1.0;
    private static int thirdBetaE_count = 1;
    private static int helixFrom_start = 3;
    private static int helixFrom_step = 1;
    private static int helixFrom_count = 1;
    private static int helixPeriod_start = 3;
    private static int helixPeriod_step = 1;
    private static int helixPeriod_count = 1;
    private static int helixCycleLength_start = 3;
    private static int helixCycleLength_step = 1;
    private static int helixCycleLength_count = 1;

    public static void main(String[] args) {
        MonteCarlo alg;
        switch ("mamooli") {
            case "gerafe gerd":
                alg = new MonteCarlo();
                alg.run();
                break;
            case "mamooli":
                System.out.println("\n\ndata={");
                for(int monomers_i = monomers_start; monomers_i < monomers_start + monomers_step * monomers_count; monomers_i += monomers_step){
                for(int helixFrom_i = helixFrom_start; helixFrom_i < helixFrom_start + helixFrom_step * helixFrom_count; helixFrom_i += helixFrom_step){
                for(int helixPeriod_i = helixPeriod_start; helixPeriod_i < helixPeriod_start + helixPeriod_step * helixPeriod_count; helixPeriod_i += helixPeriod_step){
                for(int helixCycleLength_i = helixCycleLength_start; helixCycleLength_i < helixCycleLength_start + helixCycleLength_step * helixCycleLength_count; helixCycleLength_i += helixCycleLength_step){
                for(double peptydeBetaE_i = peptydeBetaE_start; peptydeBetaE_i < peptydeBetaE_start + peptydeBetaE_step * peptydeBetaE_count; peptydeBetaE_i += peptydeBetaE_step){
                for(double hydrogenBetaE_i = hydrogenBetaE_start; hydrogenBetaE_i < hydrogenBetaE_start + hydrogenBetaE_step * hydrogenBetaE_count; hydrogenBetaE_i += hydrogenBetaE_step){
                for(double thirdBetaE_i = thirdBetaE_start; thirdBetaE_i < thirdBetaE_start + thirdBetaE_step * thirdBetaE_count; thirdBetaE_i += thirdBetaE_step){
                    alg = new MonteCarlo(monomers_i, hydrogenBetaE_i, peptydeBetaE_i, thirdBetaE_i, helixFrom_i, helixPeriod_i, helixCycleLength_i);
                    System.out.println("% n: " + monomers_i + ", betaE_h: " + hydrogenBetaE_i + ", betaE_p: " + peptydeBetaE_i);
                    System.out.println("[");
                    alg.run();
                    System.out.println("],");
                }}}}}}}
                System.out.println("};\n");
                break;
            default:
                break;
        }
    }
}

class MonteCarlo {
    private double bEV[][];
    private int e;
    private double eps = 0.00001;
    private double lnhEV[][];
    private double qEV[][];
    private double ehEV[][];
    private double epEV[][];
    private int n;
    private double btmph = 1.7;
    private double btmpp = 15.0;
    private Integer t[];

    MonteCarlo(int n, double btmph, double btmpp, double btmpc, int from, int period, int cycle){ // gerafe mamooli ba parameter haa
        this.n = n;
        this.btmph = btmph;
        this.btmpp = btmpp;
        // btmpc
        t = Matrix.helix(from, period, cycle, n);
    }

    MonteCarlo(){ // gerafe gerde 13 tayi
        n = 13;
        t = Matrix.customGraph(13);
        btmph = 7.8;
        btmpp = btmph;
    }

    public void run() {
        int MCsteps = (int) 1e+2;
        int skiptohist = 10;
        double lnf = 1;
        double lnfmin = 1.0e-4;
        e = t.length;
        int cmax = e - n + 1;
        Matrix mat = new Matrix(n, t);
        lnhEV = new double[e + 1][n + 1];
        int hist2[][] = new int[e + 1][n + 1];
        int totalnum[][] = new int[e + 1][n + 1];
        bEV = new double[e + 1][n + 1];
        qEV = new double[e + 1][n + 1];
        ehEV = new double[e + 1][n + 1];
        epEV = new double[e + 1][n + 1];
        int state0[] = new int[4];
        int state1[] = new int[4];
        for (int i = 0; i < e + 1; i++) {
            for (int j = 0; j < n + 1; j++) {
                lnhEV[i][j] = bEV[i][j] = qEV[i][j] = epEV[i][j] = ehEV[i][j] = hist2[i][j] = totalnum[i][j] = 0;
            }
        }
        Random rand = new Random();
        int sigmahist;
        while (lnf > lnfmin) {
            int MCcounter = 1;
            int skipcounter = skiptohist;
            while (MCcounter != 0) {
                sigmahist = 0;
                for (int i = 0; i < n * n; i++) {
                    state0[0] = mat.getE();
                    state0[1] = mat.getV();
                    state0[2] = mat.getnp();
                    state0[3] = mat.getnn();
                    mat.step();
                    state1[0] = mat.getE();
                    state1[1] = mat.getV();
                    state1[2] = mat.getnp();
                    state1[3] = mat.getnn();
                    double r = ((double) rand.nextInt(1000)) / 1000;
                    double lnre = exp(lnhEV[state0[0]][state0[1]] - lnhEV[state1[0]][state1[1]]);
                    double p1 = (state1[2] == 0) ? 1 : ((double) state0[3]) / state1[2];
                    double p2 = ((double) state0[2]) / state1[3];
                    double p = (p1 < p2) ? p1 : p2;
                    if (r >= lnre * p)
                        mat.stepBack();
                    lnhEV[mat.getE()][mat.getV()] += lnf;
                    bEV[mat.getE()][mat.getV()] += mat.getB();
                    qEV[mat.getE()][mat.getV()] += mat.getQ();
                    epEV[mat.getE()][mat.getV()] += mat.getStrongE();
                    ehEV[mat.getE()][mat.getV()] += mat.getE() - mat.getStrongE();
                    hist2[mat.getE()][mat.getV()]++;
                    totalnum[mat.getE()][mat.getV()]++;
                    sigmahist++;
                }
                MCcounter++;
                skipcounter++;
                if (MCcounter > MCsteps && skipcounter > skiptohist) {
                    int hmin2 = (int) 1.0e+10;
                    double lnhmin2 = 1.0e+10;
                    int Nstates2 = 0;
                    for (int i = 0; i < e + 1; i++) {
                        for (int j = 0; j < n + 1; j++) {
                            if (hist2[i][j] != 0)
                                Nstates2++;
                            if (hist2[i][j] != 0 && hist2[i][j] < hmin2)
                                hmin2 = hist2[i][j];
                            if (lnhEV[i][j] > eps && lnhEV[i][j] < lnhmin2)
                                lnhmin2 = lnhEV[i][j];
                        }
                    }
                    double havg2 = ((double) sigmahist) / Nstates2;
                    if (hmin2 / havg2 >= 0.8) {
                        for (int i = 0; i < e + 1; i++) {
                            for (int j = 0; j < n + 1; j++) {
                                hist2[i][j] = 0;
                                if (lnhEV[i][j] > eps)
                                    lnhEV[i][j] = lnhEV[i][j] - lnhmin2 + log(2);
                            }
                        }
                        lnf /= 2;
                        MCcounter = 0;
                    } else {
                        skipcounter = 1;
                    }
                }
            }
        }
        bEV[e][n] = mat.getTotalB();
        totalnum[e][n] = 1;
        double x[][] = new double[cmax + 1][e + 1];
        for (int i = 0; i < cmax + 1; i++) {
            for (int j = 0; j < e + 1; j++) {
                x[i][j] = -1;
            }
        }
        for (int i = 0; i < e + 1; i++) {
            for (int j = 0; j < n + 1; j++) {
    
                bEV[i][j] /= ((double) totalnum[i][j]);
                qEV[i][j] /= ((double) totalnum[i][j]);
                ehEV[i][j] /= ((double) totalnum[i][j]);
                epEV[i][j] /= ((double) totalnum[i][j]);
                int c = i - j + 1;
                if (c >= 0 && c < cmax + 1)
                    x[c][i] = log(exp(lnhEV[i][j]) + ((x[c][i] < -eps) ? 0 : exp(x[c][i])));
                
                    System.out.print(bEV[i][j] + " ");
                }
            System.out.println();
        }
        System.out.println();
        double[] ans = new double[n + 1];
        for (int j = 2; j < n + 1; j++) {
            ans[j] = (getBetaF(j) - getBetaF(2));
            System.out.println(j + "\t\t" + (Double.isInfinite(ans[j]) ? "-Inf" : ans[j]));
        }
    }

    private double getLnZ(int e, int v) {
        double lnqc = log(qEV[e][v]);
        double lnqd = log(qEV[e][v]);
        double lnro = log(pow(62.0, -3));
        return (ehEV[e][v] * btmph + epEV[e][v] * btmpp + v * lnro - (v - 1) * lnqc - (v - bEV[e][v] - 1) * lnqd);
    }

    private double getBetaF(int v) {
        double ans = 0;
        for (int i = 0; i < e + 1; i++) {
            double lnZ = getLnZ(i, v);
            // System.out.println(log(qEV[i][v])+" "+lnZ);
            if ((lnhEV[i][v] < eps && (i != e || v != n)) || Double.isNaN(lnZ))
                continue;
            double x = ans;
            double y = lnhEV[i][v] + lnZ;
            x = log(exp(y) + ((x > eps || x < -eps) ? exp(x) : 0));
            if (!Double.isInfinite(x)) {
                ans = x;
            }
        }
        return -ans;
    }

   
}


class Matrix {
    private int n;
    private int m[][];
    private int yalls[];
    private int E;
    private int V;
    private int B;
    private int L;
    private int A;
    private HashSet<Integer> q = new HashSet<>();
    private boolean visited[];
    private int cntE;
    private int lastI = -1;
    private int totalB;
    private LinkedList<Integer>[] yalha;
    private LinkedList<Integer>[] targety;

    Matrix(int n, Integer[] y) {
        this.n = n;
        this.cntE = y.length;
        visited = new boolean[n];
        yalls = new int[cntE + 1];
        yalls[0] = cntE;
        m = new int[n][n];
        yalha = new LinkedList[n];
        targety = new LinkedList[n];
        for (int i = 0; i < n; i++) {
            yalha[i] = new LinkedList<>();
            targety[i] = new LinkedList<>();
        }
        for (int i = 0; i < cntE; i++) {
            yalls[i + 1] = y[i];
            targety[y[i] / n].add(y[i] % n);
            targety[y[i] % n].add(y[i] / n);
            swit(y[i]);
        }
        E = cntE;
        V = n;
        A = 0;
        calcL();
        calcB();
        totalB = B;
        makeRandomPos();
    }

    private boolean isBridge(int x) {
        return hasE(x) && !connected(x);
    }

    private void calcB() {
        int cnt = 0;
        for (int i = 0; i < n; i++)
            for (int j : yalha[i])
                if (i < getNode(j, i) && isBridge(j))
                    cnt++;
        B = cnt + calcL();
    }

    private boolean connected(int num) {
        int a = num / n, b = num % n, xx = V;
        int cntA = yalha[a].size() - 1, cntB = yalha[b].size() - 1;
        V -= ((cntA == 0) ? 1 : 0) + ((cntB == 0) ? 1 : 0);
        if (V < 4) {
            V = xx;
            return true;
        }
        if (E < V) {
            V = xx;
            return false;
        }
        for (int i = 0; i < n; i++)
            visited[i] = false;
        int cnt = 0;
        for (int i = 0; i < n; i++) {
            if (yalha[i].size() - ((i == a || i == b) ? 1 : 0) > 0) {
                cnt++;
                visited[i] = true;
                q.add(i);
                break;
            }
        }
        while (!q.isEmpty()) {
            int x = q.iterator().next();
            q.remove(x);
            for (int o : yalha[x]) {
                if (o == num)
                    continue;
                int z = getNode(o, x);
                if (!visited[z]) {
                    cnt++;
                    visited[z] = true;
                    q.add(z);
                }
            }
        }
        boolean y = (cnt == V);
        V = xx;
        return y;
    }

    private int getNode(int x, int i) {
        return (x - i * n > i) ? x - i * n : (x - i) / n;
    }

    private int calcL() {
        int cnt = 0;
        for (int i = 0; i < n; i++)
            if (yalha[i].size() == 1)
                cnt++;
        L = (E == 1) ? 1 : cnt;
        return L;
    }

    private void shift(int i, int j) {
        swit(i * n + j);
        calcIt(i, j);
        calcB();
    }

    private void calcIt(int a, int b) {
        E += (m[a][b] == 0) ? -1 : 1;
        int cntA = yalha[a].size(), cntB = yalha[b].size();
        int cntA1 = targety[a].size(), cntB1 = targety[b].size();
        int x = 0;
        if (m[a][b] == 1)
            x += ((cntA == 1) ? 1 : 0) + ((cntB == 1) ? 1 : 0);
        else
            x -= ((cntA == 0) ? 1 : 0) + ((cntB == 0) ? 1 : 0);
        V += x;
        switch (x) {
        case 0:
            if (m[a][b] == 1)
                A -= 1;
            else
                A += 1;
            break;
        case 1:
            if (cntA == 1)
                A += calcA(a) - 1;
            else
                A += calcA(b) - 1;
            break;
        case 2:
            A += cntA1 + cntB1 - 2;
            break;
        case -1:
            if (cntA == 0)
                A -= calcA(a) - 1;
            else
                A -= calcA(b) - 1;
            break;
        case -2:
            A = 0;
        }
    }

    private int calcA(int a) {
        ArrayList<Integer> ints = new ArrayList<>(targety[a]);
        int cnt = 0;
        for (int i : ints)
            cnt += (yalha[i].size() == 0) ? 1 : 0;
        return cnt;
    }

    private void makeRandomPos() {
        makeEmpty();
        for (int i = 0; i < cntE * cntE; i++)
            step();
    }

    private void makeEmpty() {
        for (int i = 0; i < n; i++)
            while (yalha[i].size() > 0)
                swit(yalha[i].getFirst());
        A = B = L = V = E = 0;
    }

    private void swit(int x) {
        if (!hasE(x)) {
            yalha[x / n].add(x);
            yalha[x % n].add(x);
        } else {
            yalha[x / n].removeFirstOccurrence(x);
            yalha[x % n].removeFirstOccurrence(x);
        }
        m[x / n][x % n] = 1 - m[x / n][x % n];
        m[x % n][x / n] = 1 - m[x % n][x / n];
    }

    int getStrongE() {
        int cnt = 0;
        for(int i = 0; i < n-1; i++){
            if(m[i][i+1] == 1){
                cnt++;
            }
        }
        return cnt;
    }

    int getE() {
        return E;
    }

    int getnp() {
        return A;
    }

    int getnn() {
        return E - B + L;
    }

    void step() {
        lastI = getRand();
        shift(lastI / n, lastI % n);
    }

    void stepBack() {
        shift(lastI / n, lastI % n);
    }

    private boolean hasE(int i) {
        return m[i / n][i % n] != 0;
    }

    private boolean canPut(int x) {
        if (E == 0)
            return true;
        int i = x / n;
        int j = x % n;
        return yalha[i].size() > 0 || yalha[j].size() > 0;
    }

    private int getRand() {
        Random r = new Random();
        int y = r.nextInt(99999999);
        int z = y % yalls[0] + 1;
        int x = yalls[z];
        while (hasE(x) && getE() == 1 || !hasE(x) && getE() == cntE - 1 || isBridge(x) || !hasE(x) && !canPut(x)) {
            y = r.nextInt(99999999);
            z = y % yalls[0] + 1;
            x = yalls[z];
        }
        return x;
    }

    private static int getNum(int a, int b, int n) {
        a--;
        b--;
        if (a > b) {
            b = a + b;
            a = b - a;
            b = b - a;
        }
        return a * n + b;
    }

    static Integer[] helix(int from, int period, int cycle, int n) {
        ArrayList<Integer> anslist = new ArrayList<>();
        int[] f = new int[n];
        for (int i = from-1; i < n; i++) f[i] = cycle * ((i+1-from)%period==0?1:0);
        /* printing */ System.out.print("% "); for (int aX : f) System.out.print(aX+","); System.out.println();
        for (int i = 0; i < n - 1; i++) {
            int index = i + 1;
            anslist.add(getNum(index, index + 1, n));
            if (i + f[i] < n && f[i] > 0) anslist.add(getNum(index, index + f[i], n));
        }
        return anslist.toArray(new Integer[0]);
    }

    static Integer[] customGraph(int n) {
        ArrayList<Integer> anslist = new ArrayList<>();
        n=13;
        anslist.add(getNum(1, 2, n));
        anslist.add(getNum(1, 3, n));
        anslist.add(getNum(1, 4, n));
        anslist.add(getNum(1, 5, n));
        anslist.add(getNum(1, 7, n));
        anslist.add(getNum(2, 6, n));
        anslist.add(getNum(2, 7, n));
        anslist.add(getNum(3, 6, n));
        anslist.add(getNum(3, 7, n));
        anslist.add(getNum(4, 6, n));
        anslist.add(getNum(4, 7, n));
        anslist.add(getNum(5, 6, n));
        anslist.add(getNum(5, 7, n));
        anslist.add(getNum(6, 7, n));
        anslist.add(getNum(7, 8, n));
        anslist.add(getNum(7, 9, n));
        anslist.add(getNum(7, 10, n));
        anslist.add(getNum(7, 11, n));
        anslist.add(getNum(7, 12, n));
        anslist.add(getNum(7, 13, n));
        anslist.add(getNum(8, 9, n));
        anslist.add(getNum(8, 10, n));
        anslist.add(getNum(8, 11, n));
        anslist.add(getNum(8, 12, n));
        anslist.add(getNum(9, 13, n));
        anslist.add(getNum(10, 13, n));
        anslist.add(getNum(11, 13, n));
        anslist.add(getNum(12, 13, n));
        return anslist.toArray(new Integer[0]);
    }

    double getQ() {
        return E * 2.0 / V;
    }

    int getV() {
        return V;
    }

    int getB() {
        return B;
    }

    int getTotalB() {
        return totalB;
    }
}