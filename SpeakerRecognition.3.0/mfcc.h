//
//  mfcc.h
//  SpeakerRecognition.2.0
//
//  Created by Nancy Fan on 16/5/17.
//  Copyright © 2016年 Nancy Fan. All rights reserved.
//

#ifndef mfcc_h
#define mfcc_h

#define _CRT_SECURE_NO_WARNINGS

#include<iostream>
#include<fstream>
#include<string>
#include<cstdio>
#include<cmath>
#include<vector>
#include <complex>
#include <bitset>
//#include <conio.h>
//#include<fstream.h>

//#define DEBUG
#define RUN

using namespace std;



const int FS = 16;
const int FrmLen = 20 * FS;  //窗口长度
const unsigned long FFTLen = 512; //参与FFT运算的512个数据
const double PI = 3.1415926536;
const int FiltNum = 25;   //滤波器组数，一共25组
const int PCEP = 13;     //MFCC的阶数，最后得到的关于的13个MFCC的系数
const int GmmNum = 8;
const double infinite = 1e+32;



typedef struct _TWavHeader
{
    int rId;    //标志符（RIFF）
    int rLen;   //数据大小,包括数据头的大小和音频文件的大小
    int wId;    //格式类型（"WAVE"）
    int fId;    //"fmt"
    
    int fLen;   //Sizeof(WAVEFORMATEX)
    
    short wFormatTag;       //编码格式，包括WAVE_FORMAT_PCM，WAVEFORMAT_ADPCM等
    short nChannels;        //声道数，单声道为1，双声道为2
    int nSamplesPerSec;   //采样频率
    int nAvgBytesPerSec;  //每秒的数据量
    short nBlockAlign;      //块对齐
    short wBitsPerSample;   //WAVE文件的采样大小
    int fact; //0x6661 6374
    int wSampleLength; //4
    int factdata;
    int dId; //0x64617461
    int dbyte;
    //int dId;              //"data"
    //int wSampleLength;    //音频数据的大小
} TWavHeader;


typedef struct _GmmPara
{
    double weight;  //存放权值
    double mean[PCEP];  //存放均值
    double variance[PCEP];  //存放方差
}GmmPara;


int getmfcc(char *path, string Name);
void InitHamming();
void HammingWindow(short* buf, float* data);
float GetSTE(short* data);
int GetZcr(short *data);
void compute_fft(float *buffer, vector<complex<float> >& vecList);
void FFT(const unsigned long & ulN, vector<complex<float> >& vecList); //FFT的实际程序
void display(const unsigned long & ulN, vector<complex<float> >& vecList);
void InitFilt(float *FiltCoe1, float *FiltCoe2, int *Num); //初始化滤波器
void CFilt(float *spdata, float *FiltCoe1, float *FiltCoe2, int *Num, float *En, vector<complex<float> >& vecList);
void MFCC(float *En); //计算MFCC的13个系数
void cluster();
void EM();
float CalculateLikelihood();
void WritePara();
int get_mfcc(int mode);
size_t myfread(void *buffer, size_t size, size_t count, FILE *stream, long long line);
int vad(FILE * INFILE, FILE* OUTFILE);

#endif /* mfcc_h */
