#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include "mpi.h"

using namespace std;

long long Combination(int n, int k)
{
    long long C = 1;
    for (int i = 1; i <= k; i++)
        C = C * (n - k + i) / i;
    return C;
}

void DFS(int st, int** graph, int n, int* visited)
{
    visited[st] = 1;
    for (int i = 0; i < n; i++)
        if (graph[st][i] != 0 && !visited[i])
            DFS(i, graph, n, visited);
}

int Check(int n, int k, int rank, int size, int** a)
{
    int res = 0, n1, f, k1 = 0, k2 = 0, mod, t, K = 0;
    bool ret = false;
    do {
        K++;
        n1 = n - K;
        int** a1 = new int* [n1];
        for (int i = 0; i < n1; i++)
            a1[i] = new int[n1];
        int* visited = new int[n1];
        int* b = new int[K + 1];
        for (int i = 0; i <= K; i++)
            b[i] = i;
        long long C = Combination(n, K);
        mod = C % size;
        long long S = C / size, s;
        if (rank + 1 == size)
            S += mod;
        long long I = S * rank + 1;
        for (long long i = 0; i < S; i++, I++)
        {
            s = 0;
            for (int j, t = 1; t <= K; t++)
            {
                j = b[t - 1] + 1;
                for (; j < n1 + t && s + (C = Combination(n - j, K - t)) < I; j++)
                    s += C;
                b[t] = j;
            }
            for (int y = 0; y < n; y++)
            {
                k2 = 0;
                f = 0;
                for (int x = 1; x <= K; x++)
                    if (b[x] == y + 1)
                    {
                        f = 1;
                        k1++;
                    }
                if (f == 1)
                    continue;
                for (int j1 = 0; j1 < n; j1++)
                {
                    f = 0;
                    for (int x = 1; x <= K; x++)
                        if (b[x] == j1 + 1)
                        {
                            f = 1;
                            k2++;
                        }
                    if (f == 1)
                        continue;
                    a1[y - k1][j1 - k2] = a[y][j1];
                }
            }
            f = 0;
            for (int i2 = 0; i2 < n1; i2++)
                visited[i2] = 0;
            DFS(0, a1, n1, visited);
            for (int x = 0; x < n1; x++)
                if (visited[x] == 0)
                {
                    cout << "The process under the number " << rank << " found the following combination of vertices №" << I << ", which leads to a rupture: " << '{';
                    for (int i1 = 1; i1 <= K; i1++)
                    {
                        cout << b[i1];
                        if (i1 < K) cout << ", ";
                    }
                    cout << '}' << endl;
                    ret = true;
                    break;
                }
            k1 = 0;
            for (int i1 = 0; i1 <= K; i1++)
                b[i1] = i1;
        }
        for (int i = 0; i < n1; i++)
        {
            delete[] a1[i];
        }
        delete[] a1;
        delete[] b;
        MPI_Allreduce(MPI_IN_PLACE, &ret, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
        if (ret)
            return K;
    } while (K != k);
    return k + 1;
}

int main(int argc, char** argv)
{
    setlocale(LC_ALL, "Russian");
    MPI_Init(&argc, &argv);
    double start_time = MPI_Wtime();
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int n, m, k;
    ifstream inputFile("input.txt");
    if (!inputFile)
    {
        cout << "Error! Failed to open input.txt file" << endl;
        double end_time = MPI_Wtime();//останавливам таймер
        double search_time = (end_time - start_time);
        cout << "Processor №" << rank << " has finished its work. Running time : " << search_time << "\n" << endl;
        MPI_Finalize();
        return 1;
    }
    inputFile >> n >> m >> k;
    int** a = new int* [n];
    for (int i = 0; i < n; i++)
        a[i] = new int[n];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            a[i][j] = 0;
    for (int i = 0; i < m; i++)
    {
        int u, v;
        inputFile >> u >> v;
        a[u - 1][v - 1] = 1;
        a[v - 1][u - 1] = 1;
    }
    if (k <= 0)
        cout << "The value entered is wrong! k must be > 0." << endl;
    else if (k == n - 1)
        cout << "In this case k = n - 1 and if the graph is complete then it will be vertex " << k << "-connected, and if not, try taking k < " << k << endl;
    else
    {
        int K = Check(n, k, rank, size, a);
        if (k >= K)
            cout << "Process №" << rank << ": the graph is vertex " << K << " - connected" << endl;
        else
            cout << "Process №" << rank << " did not find a single combination to break the graph. ";
    }
    double end_time = MPI_Wtime();//останавливам таймер
    double search_time = (end_time - start_time);
    cout << "Processor №" << rank << " has finished its work. Running time : " << search_time << "\n" << endl;
    for (int i = 0; i < n; i++)
    {
        delete[] a[i];
    }
    delete[] a;
    MPI_Finalize();
    return 0;
}