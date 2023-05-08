#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#define pi 3.1415926535

using namespace std;
void Sample(double* x,int M,double fs)
{
    for (int i = 0; i <=M; i++) {
        x[i] = 0.8 * sin((double)(2 * pi * 103 * i) / fs) + sin((double)(2 * pi * 107 * i) / fs) + 0.1 * sin((double)(2 * pi * 115 * i) / fs);
    }
   
}
void Windowing(double* x,int M,int N,int window_type)
{
    switch (window_type) {
    case 1://矩形窗
        
        break;
    case 2://三角窗
        for (int i = 0; i <= (int)(M / 2); i++) {
            x[i] = x[i] * 2 * i / M;
        }
        for (int i = (int)(M / 2) + 1; i <= M; i++) {
            x[i] = 2*x[i]-x[i] * 2 * i / M;
        }
        break;
    case 3://汉宁窗
        for (int i = 0; i <= M; i++) {
            x[i] = x[i] * 0.5 * (1-cos(2*pi*i/M));
        }
        break;
    case 4://哈明窗
        for (int i = 0; i <= M; i++) {
            x[i] = x[i] * (0.54 - 0.46*cos(2 * pi * i / M));
        }
        break;
    case 5://布莱克曼窗
        for (int i = 0; i <= M; i++) {
            x[i] = x[i] * (0.42 - 0.5 * cos(2 * pi * i / M) + 0.08 * cos(4 * pi * i / M));
        }
        break;
    }
    for (int i = M + 1; i < N; i++) { //补零
        x[i] = 0;
    }
}
void Rader(double *x, std::complex<double>* X,int N)
{
    int cur_rev = 0;
    int k = N / 2;  //权系数初始值
    X[0].real(x[cur_rev]); 
    X[0].imag(0);
    for (int j = 1; j <= N - 1; j++) {
        if (cur_rev < k) { //j-1的倒序数的最高位为0，说明从j-1到j只有最低位从0变1，没有进位
            cur_rev = cur_rev + k;//把j-1对应的倒序数的最高位从0变成1得到当前j的倒序数

        }
        else {
            while (cur_rev >= k) {    //j-1的倒序数的最高位为1，则从j-1到j有进位
                cur_rev = cur_rev - k;//把j-1的倒序数最高位从1变为0
                k = k / 2;            //次高位为1，则更新权系数，这个循环实际上在求从j-1到j进了几位
            }
            cur_rev = cur_rev + k;    //最高位为0时跳出循环，把最高位置1得到当前倒序数
            k = N / 2;//还原权系数值
        }
        X[j].real(x[cur_rev]);
        X[j].imag(0);
    }
 }
std::complex<double>* FFT(std::complex<double>* X, std::complex<double>* X_temp,int N)
{
   
    std::complex<double>* W_N = (std::complex<double> *)malloc((int)(sizeof(std::complex<double>) * N/2));
    for (int i = 0; i < (int)(N / 2); i++) {
        W_N[i] = std::complex<double>(cos(2 * pi * i / N), -sin(2 * pi *i / N));
    }
    int m = (int)round(log(N) / log(2));//蝶形层数
    for (int i=1; i <= m; i++) { //m层
        for (int j=0; j < (int)pow(2,m - i); j++) { //每一层做2^i点DFT的次数
            for (int k = 0; k < (int)pow(2,i - 1); k++) { //DFT的前半部分
               X[k + j * (int)pow(2, i) + (int)pow(2, i - 1)] = W_N[k*(int)pow(2,m-i)]* X[k + j * (int)pow(2, i) + (int)pow(2, i - 1)];
               X_temp[k] = X[k + j * (int)pow(2, i)];
               X[k + j * (int)pow(2, i)] = X[k + j * (int)pow(2, i)]+X[k + j * (int)pow(2, i) + (int)pow(2, i - 1)];
            }
            for (int k = (int)pow(2, i - 1); k <(int)pow(2,i); k++) { //DFT的后半部分
                X[k + j * (int)pow(2, i)] = -X[k + j * (int)pow(2, i)] + X_temp[k - (int)pow(2, i - 1)];
            }
        }

    }
    return X;

}
std::complex<double>* DFT(complex<double>* X_DFT,double* x,int N)
{
    
    for (int i = 0; i < N; i++) { //N点DFT
        X_DFT[i] = std::complex<double>{ 0,0 };
        for (int j = 0; j < N; j++) { //每点进行N次复乘
            X_DFT[i] = X_DFT[i] + std::complex<double>{cos(2 * pi * i * j / N),-sin(2 * pi * i * j / N)}*std::complex<double>{ x[j], 0 };
        }
        //cout << X_DFT[i] << endl;

    }
    return X_DFT;

}
int main(void)
{   
    int M;//窗长
    int N;//变换点数
    double fs;//采样率
    int window_type;//窗函数类型
    clock_t start, end;//用于计时
    printf("请输入窗长M：");
    scanf_s("%d", &M);
    printf("请输入变换点数N：");
    scanf_s("%d", &N);
    printf("请输入采样率fs：");
    scanf_s("%lf", &fs);
    printf("请输入窗函数类型：");
    scanf_s("%d", &window_type);
    double* x = (double*)malloc(sizeof(double) * N);
    std::complex<double>* X = (std::complex<double> *)malloc(sizeof(std::complex<double>) * N);
    std::complex<double>* X_temp = (std::complex<double> *)malloc((int)(sizeof(std::complex<double>) * N / 2));
    std::complex<double>* X_DFT = (std::complex<double> *)malloc(sizeof(std::complex<double>) * N);
    Sample(x, M, fs);
    Windowing(x,M,N,window_type);//加窗
    start = clock();
    X_DFT=DFT(X_DFT,x,N);
    end = clock();
    cout << "DFT用时：" << 1000 * (float)(end - start) / CLOCKS_PER_SEC<< "ms" << endl;
    start = clock();
    Rader(x,X,N);//得到倒位序的x[n]
    X=FFT(X,X_temp,N);
    end = clock();
    cout << "FFT用时：" << 1000 * (float)(end - start) / CLOCKS_PER_SEC<< "ms" << endl;
    ofstream p;
    string filenamep = "FFT_M" + (string)to_string(M) + "_N"+ (string)to_string(N) +"_Window"+ (string)to_string(window_type) +"_fs"+ (string)to_string((int)fs) +".csv";
    p.open(filenamep,ios::out);
    p << "Real" << "," << "Imag" << endl;
    for (int i = 0; i < N; i++) {
        p << (string)to_string(X[i].real()) << "," << (string)to_string(X[i].imag()) << endl;
    }
    p.close();
    ofstream q;
    string filenameq = "DFT_M" + (string)to_string(M) + "_N" + (string)to_string(N) + "_Window" + (string)to_string(window_type) + "_fs" + (string)to_string((int)fs) + ".csv";;
    q.open(filenameq, ios::out);
    q << "Real" << "," << "Imag" << endl;
    for (int i = 0; i < N; i++) {
        q << (string)to_string(X_DFT[i].real()) << "," << (string)to_string(X_DFT[i].imag()) << endl;
    }
    q.close();

    
}