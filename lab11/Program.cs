using System.Diagnostics;
using System.Globalization;

int n = 6;
int m = 20;
int valbuilt = 103;
const double x = 0.1;
double[,] A = new double[n, m];
double[] B = new double[m];
double[] C = new double[m];
double[] D = new double[m];

for (int i = 1; i <= n; i++)
{
    for (int j = 1; j <= m; j++)
    {
        A[i - 1, j - 1] = Math.Pow(Math.Tan(Math.Pow(j + 1, 3)), 2) / (i - i + 3) - Math.Pow(i, 2 / x);
    }
}
for (int i = 1; i <= m; i++)
    B[i - 1] = Math.Pow(Math.Tan(Math.Pow(7 + i, 3 / i)), 2) / (i + x - 3);

int count = 0;
Stopwatch stopwatch1, stopwatch2;

while (count < 10)
{
    stopwatch1 = new();
    stopwatch2 = new();
    stopwatch1.Start();
    AlgorithmBlock(A, B, C, n, m);
    stopwatch1.Stop();
    Console.WriteLine("Execution time of the block algorithm with n = {0} and m = {1}: {2}", n, m, stopwatch1.Elapsed.TotalMilliseconds);
    stopwatch2.Start();
    AlgorithmRows(A, B, D, n, m);
    stopwatch2.Stop();
    Console.WriteLine("Execution time of the rows algorithm with n = {0} and m = {1}: {2}", n, m, stopwatch2.Elapsed.TotalMilliseconds);
    if (count < 4)
    {
        n += valbuilt;
    }
    else
    {
        m += valbuilt;
        B = new double[m];
        C = new double[m];
        D = new double[m];
        for (int i = 1; i <= m; i++)
            B[i - 1] = Math.Pow(Math.Tan(Math.Pow(7 + i, 3 / i)), 2) / (i + x - 3);
    }
    A = new double[n, m];
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= m; j++)
        {
            A[i - 1, j - 1] = Math.Pow(Math.Tan(Math.Pow(j + 1, 3)), 2) / (i - i + 3) - Math.Pow(i, 2 / x);
        }
    }
    count++;
}

void AlgorithmBlock(double[,] Matrix, double[] Vector, double[] Result, int rows, int cols)
{
    int ThreadsNum = 2;
    int BlockSize = rows / ThreadsNum;
    int BlockSizeC = cols / ThreadsNum;

    Parallel.For(0, ThreadsNum, ThreadIndex =>
    {
        // Рядки
        int StartIndexI = ThreadIndex * BlockSize;
        int EndIndexI = StartIndexI + BlockSize;
        
        if (ThreadIndex == ThreadsNum - 1 && EndIndexI < rows)
        {
            EndIndexI = rows;
        }
        Parallel.For(0, ThreadsNum, ThreadIndexC =>
        {
            //Стовпці
            int StartIndexJ = ThreadIndexC * BlockSizeC;
            int EndIndexJ = StartIndexJ + BlockSizeC;
            
            if (ThreadIndexC == ThreadsNum - 1 && EndIndexJ < cols)
            {
                EndIndexJ = cols;
            }
            for (int i = StartIndexI; i < EndIndexI; i++)
            {
                for (int j = StartIndexJ; j < EndIndexJ; j++)
                {
                    Result[j] += Matrix[i, j] * Vector[j];
                }
            }
        });
    });
}

void AlgorithmRows(double[,] Matrix, double[] Vector, double[] Result, int rows, int cols)
{
    Parallel.For(0, rows, ThreadIndex =>
    {
        for (int j = 0; j < cols; j++)
        {
            Result[j] += Matrix[ThreadIndex, j] * Vector[j];
        }
    });
}