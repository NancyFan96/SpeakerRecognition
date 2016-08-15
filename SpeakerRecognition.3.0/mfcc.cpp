#include "mfcc.h"
#include <iostream>

double Hamming[FrmLen];

vector<float>xishu;

GmmPara  para[GmmNum];
extern ofstream logfile;

size_t myfread(void *buffer, size_t size, size_t count, FILE *stream, long long line)
{
    size_t cnt = 0;
    short mybuffer;
    while(cnt < count && fread(&mybuffer, size, 1, stream)!=0){
        if(mybuffer* mybuffer >line){
            ((short*)buffer)[cnt]=mybuffer;
            cnt++;
            //cout <<cnt<<":  "<< mybuffer <<endl;
        }
    }
    
    return cnt;
}
int vad(FILE * INFILE, FILE* OUTFILE)
{
    return 0;
}

int getmfcc(char *path, string Name)
{
    TWavHeader waveheader;
    FILE *sourcefile;
    short buffer[FrmLen];
    float data[FrmLen];  //加窗后得到的数据
    int index = 0, count = 0;
    float energy = 0.0, sum = 0.0;
    float FiltCoe1[FFTLen / 2 + 1];  //左系数
    float FiltCoe2[FFTLen / 2 + 1];  //右系数
    int Num[FFTLen / 2 + 1];     //决定每个点属于哪一个滤波器，一般而言，每个点会包含在相邻的两个滤波器中，这里是与该点相关的第二个滤波器
    float En[FiltNum + 1];         //频带能量
    
    vector<complex<float> > vecList;
    ofstream outfile1(Name);          //建立音频文件所对应的.txt文件，用于存储提取MFCC的内容
    sourcefile = fopen(path, "rb");     //打开要提取的音频文件
    
    if(sourcefile == NULL){
        cout << "找不到音频文件:"<<path<<endl;
        return 1;
    }
    
    
    
    
    //get vad.wav
    
    long long aveEn = 0;
    long long line = 0;
    int cnt = 0;
    int rlen = 0;
    char vadfilename[128]="./vad/";

    rlen = strlen(Name.c_str());
    while(Name.c_str()[rlen-1]!='/')
        rlen--;
    strcat(vadfilename, Name.c_str()+rlen);
    
    rlen = strlen(vadfilename);
    vadfilename[rlen-3]='w';
    vadfilename[rlen-2]='a';
    vadfilename[rlen-1]='v';
    
    fread(&waveheader, sizeof(struct _TWavHeader), 1, sourcefile);
    while (rlen = fread(buffer, sizeof(short), FrmLen, sourcefile), rlen )
    {
        long add = 0;
        for(int i = 1; i <= rlen;i++)
            add = add*(i-1)/i+(buffer[i]*buffer[i])/i;
        cnt+=rlen;
        aveEn = aveEn*(cnt-1)/cnt+add/cnt;
    }
    cnt<<=1;
    //cout << "readin cnt = "<<cnt<<endl;
    cnt=0;
    FILE * outfile = fopen(vadfilename, "wb");
    
    fseek(sourcefile, sizeof(struct _TWavHeader),0);
    fwrite(&waveheader, sizeof(struct _TWavHeader), 1, outfile);
    line = 0.02*aveEn;
    
    while ((rlen = myfread(buffer, sizeof(short), FrmLen, sourcefile, line)),rlen)
    {
        rlen = fwrite(buffer, sizeof(short), FrmLen, outfile);
        cnt+=rlen;
    }
    
    cnt <<=1;
    //cout << "write cnt = "<<cnt<<endl;
    cnt = cnt + 58;
    //cout << "new flen = "<<cnt<<endl;
    int * intbuff = new int;
    *intbuff = cnt - 8;//rlen
    fseek(outfile, sizeof(int), 0);
    fwrite(intbuff, sizeof(int), 1, outfile);
    *intbuff = (cnt - 56)/4;//factdata
    fseek(outfile, sizeof(int)*11, 0);
    fwrite(intbuff, sizeof(int), 1, outfile);
    *intbuff = cnt - 56;//dbyte
    fseek(outfile, sizeof(int)*13, 0);
    fwrite(intbuff, sizeof(int), 1, outfile);
    fclose(outfile);
    fclose(sourcefile);
    
    sourcefile = fopen(vadfilename, "rb");
    //get ave energy
    /*long long aveEn = 0;
    long long line = 0;
    int cnt = 0;
    int rlen = 0;
    
    while (rlen = fread(buffer, sizeof(short), FrmLen, sourcefile), rlen )
    {
        long add = 0;
        for(int i = 1; i <= rlen;i++)
            add = add*(i-1)/i+(buffer[i]*buffer[i])/i;
        cnt+=rlen;
        aveEn = aveEn*(cnt-1)/cnt+add/cnt;
    }
     
    fseek(sourcefile, sizeof(struct _TWavHeader),0);
    
    line = 0.01*aveEn;
     */
    fread(&waveheader, sizeof(struct _TWavHeader), 1, sourcefile);
    //从文件sourcefile中读入1个_TWavHeader字节的数据项到waveheader所确定的地址中
    InitHamming();//初始化汉明窗
    InitFilt(FiltCoe1, FiltCoe2, Num); //初始化MEL滤波系数
    
    /*int fread( void *buffer, size_t size, size_t num, FILE *stream );
     函数fread()读取[num]个对象(每个对象大小为size(大小)指定的字节数),
     并把它们替换到由buffer(缓冲区)指定的数组. 数据来自给出的输入流.
     函数的返回值是读取的内容数量...
     */
    
    while (fread(buffer, sizeof(short), FrmLen, sourcefile) == FrmLen)
    {
 
        HammingWindow(buffer, data);
        compute_fft(data, vecList);
        CFilt(data, FiltCoe1, FiltCoe2, Num, En, vecList);
        MFCC(En);
        vecList.clear();
        index++;
    }

    
    outfile1 << index << endl;
    logfile << index << endl;
    
    int length = xishu.size();
    for (int i = 0; i<length; ++i)
    {
        outfile1 << xishu[i] << ' ';
        if ((i + 1) % 13 == 0)
            outfile1 << endl;
    }
    //::MessageBox(NULL, "已经提取好了MFCC特征", "恭喜!", NULL);
    return 0;
    
}

void InitHamming()   //初始化汉明窗,目前常用的窗函数是Hamming窗(即升余弦窗)
{                    //W(n)=0.54-0.46*cos(2*PI*n)/(N-1),   0<=n<=N-1
    float twopi;
    int i;
    twopi = 8.0F*atan(1.0F);  //求2*Pi
    for (i = 0; i<FrmLen; i++)
    {
        Hamming[i] = (float)(0.54 - 0.46*cos((float)i*twopi / (float)(FrmLen - 1)));
    }
}

void HammingWindow(short* buf, float* data) //加窗后的信号y(n)=x(n)*w(n),x(n)为帧信号,w(n)为窗函数
{
    int i;
    for (i = 0; i<FrmLen; i++)
    {
        data[i] = buf[i] * Hamming[i];
    }
}

float GetSTE(short* data)
{
    int i;
    float STE = 0.0;
    for (i = 0; i<FrmLen; i++)
    {
        STE += (data[i] * data[i]);
    }
    STE /= (float)FrmLen;
    //	cout<<STE<<' ';
    return STE;
}

int GetZcr(short *data)
{
    int pre = data[0];
    int zero = 0;
    for (int i = 1; i<FrmLen; ++i)
    {
        if ((pre>0 && data[i]<0) || pre<0 && data[i]>0)
            zero++;
        pre = data[i];
    }
    return zero;
}

void compute_fft(float *data, vector<complex<float> >& vecList)
{
    for (int i = 0; i<FFTLen; ++i)
    {
        if (i<FrmLen)
        {
            complex<float> temp(data[i]);
            vecList.push_back(temp);
        }
        else
        {
            complex<float> temp(0);
            vecList.push_back(temp);
        }
    }
    FFT(FFTLen, vecList);
}
/*void compute_fft(double *data)
 {
 for(int i=0;i<FFTLen;++i)
 {
 if(i<FrmLen)
 vecList.push_back(data[i]);
 else
 vecList.push_back(0);
 }
 FFT(FFTLen,vecList);
 }*/
void FFT(const unsigned long & ulN, vector<complex<float> >& vecList)
{
    //得到幂数
    
    unsigned long ulPower = 0; //幂数
    unsigned long ulN1 = ulN - 1;
    while (ulN1 > 0)
    {
        ulPower++;
        ulN1 /= 2;
    }
    //反序
    
    bitset<sizeof(unsigned long) * 8> bsIndex; //二进制容器
    unsigned long ulIndex; //反转后的序号
    unsigned long ulK;
    for (unsigned long p = 0; p < ulN; p++)
    {
        ulIndex = 0;
        ulK = 1;
        bsIndex = bitset<sizeof(unsigned long) * 8>(p);
        for (unsigned long j = 0; j < ulPower; j++)
        {
            ulIndex += bsIndex.test(ulPower - j - 1) ? ulK : 0;
            ulK *= 2;
        }
        
        if (ulIndex > p)
        {
            complex<float> c = vecList[p];
            vecList[p] = vecList[ulIndex];
            vecList[ulIndex] = c;
        }
    }
    
    //计算旋转因子
    
    vector<complex<float> > vecW;
    for (unsigned long i = 0; i < ulN / 2; i++)
    {
        vecW.push_back(complex<float>(cos(2 * i * PI / ulN), -1 * sin(2 * i * PI / ulN)));
    }
    
    /*for(unsigned long m = 0; m < ulN / 2; m++)
     {
     cout<< "\nvW[" << m << "]=" << vecW[m];
     } */
    
    //计算FFT
    
    unsigned long ulGroupLength = 1; //段的长度
    unsigned long ulHalfLength = 0; //段长度的一半
    unsigned long ulGroupCount = 0; //段的数量
    complex<float> cw; //WH(x)
    complex<float> c1; //G(x) + WH(x)
    complex<float> c2; //G(x) - WH(x)
    for (unsigned long b = 0; b < ulPower; b++)
    {
        ulHalfLength = ulGroupLength;
        ulGroupLength *= 2;
        for (unsigned long j = 0; j < ulN; j += ulGroupLength)
        {
            for (unsigned long k = 0; k < ulHalfLength; k++)
            {
                cw = vecW[k * ulN / ulGroupLength] * vecList[j + k + ulHalfLength];
                c1 = vecList[j + k] + cw;
                c2 = vecList[j + k] - cw;
                vecList[j + k] = c1;
                vecList[j + k + ulHalfLength] = c2;
            }
        }
    }
}

void display(const unsigned long & ulN, vector<complex<float> >& vecList)
{
    ofstream outfile("out.txt");
    outfile << "\n\n===========================Display The Result=========================" << endl;
    for (unsigned long d = 0; d < ulN; d++)
    {
        outfile << "X(" << d << ")\t\t\t = " << vecList[d] << endl;
    }
}

/*
 设置滤波器参数
 输入参数：无
 输出参数：*FiltCoe1---三角形滤波器左边的系数
 *FiltCoe2---三角形滤波器右边的系数
 *Num     ---决定每个点属于哪一个滤波器
 */
void InitFilt(float *FiltCoe1, float *FiltCoe2, int *Num)
{
    int i, j;
    float Freq;
    int FiltFreq[FiltNum + 1] = { 0,100,200,300,400,500,600,700,800,900,1000,
        1149,1320,1516,1741,2000,2297,2639,3031,3482,4000,
        4595,5278,6063,6964,8001 };//滤波器的中心频率
    int BW[FiltNum + 1] = { 100,100,100,100,100,100,100,100,100,100,124,
        160,184,211,242,278,320,367,422,484,556,
        639,734,843,969,1112 };//滤波器的带宽
    for (i = 0; i <= FFTLen / 2; i++)
    {
        Num[i] = 0;
    }
    
    for (i = 0; i <= FFTLen / 2; i++)
    {
        Freq = FS * 1000.0F * (float)(i) / (float)(FFTLen);
        for (j = 0; j <FiltNum; j++)
        {
            if (Freq >= (float)FiltFreq[j] && Freq <= (float)FiltFreq[j + 1])
            {
                Num[i] = j;
                if (j == 0)
                {
                    FiltCoe1[i] = 0.0F;
                }
                else
                {
                    FiltCoe1[i] = ((float)(FiltFreq[j] + BW[j]) - Freq) / (float)(BW[j]);
                }
                FiltCoe2[i] = (Freq - (float)(FiltFreq[j + 1] - BW[j + 1])) / (float)(BW[j + 1]);
                FiltCoe1[i] = FiltCoe1[i] * FiltCoe1[i];
                FiltCoe2[i] = FiltCoe2[i] * FiltCoe2[i];
                break;
            }
        }
    }
    
}


/*
 根据滤波器参数计算频带能量
 输入参数：*spdata  ---预处理之后的一帧语音信号
 *FiltCoe1---三角形滤波器左边的系数
 *FiltCoe2---三角形滤波器右边的系数
 *Num     ---决定每个点属于哪一个滤波器
 
 输出参数：*En      ---输出对数频带能量
 */
void CFilt(float *spdata, float *FiltCoe1, float *FiltCoe2, int *Num, float *En, vector<complex<float> >& vecList)
{
    
    float temp = 0;
    int id, id1, id2;
    
    for (id = 0; id <= FiltNum; id++)
    {
        En[id] = 0.0F;
    }
    for (id = 0; id < FFTLen / 2; id++)
    {
        temp = vecList[id].real()*vecList[id].real() + vecList[id].imag()*vecList[id].imag();
        id1 = Num[id];
        id2 = id1 + 1;
        En[id1] = En[id1] + FiltCoe1[id] * temp;
        En[id2] = En[id2] + FiltCoe2[id] * temp;
    }
    for (id = 1; id <= FiltNum; id++)
    {
        if (En[id] != 0)
            En[id] = (float)log(En[id]);
    }
}

/*
 计算MFCC系数
 输入参数：*En ---对数频带能量
 */

void MFCC(float *En)
{
    int idcep, iden;
    float Cep[13];
    
    for (idcep = 0; idcep < PCEP; idcep++)
    {
        Cep[idcep] = 0.0;
        
        for (iden = 1; iden <= FiltNum; iden++)
        {
            Cep[idcep] = Cep[idcep] + En[iden] * (float)cos((idcep + 1) * (iden - 0.5F) * PI / (FiltNum));
        }
        Cep[idcep] = Cep[idcep] / 10.0F;
        xishu.push_back(Cep[idcep]);
    }
}

void cluster()
{
    int i, j, k, iter = 0, index;
    int TotalFrame = xishu.size() / PCEP;
    cout << TotalFrame << endl;
    float OldMean[GmmNum][PCEP], NewMean[GmmNum][PCEP], variance[GmmNum][PCEP];
    for (i = 0; i<GmmNum; ++i)
        for (j = 0; j<PCEP; ++j)
        {
            index = i*PCEP + j;
            OldMean[i][j] = xishu[index];
        }
    
    for (i = 0; i<GmmNum; ++i)
        for (j = 0; j<PCEP; ++j)
        {
            NewMean[i][j] = 0;
        }
    float sum, temp, min, sum1;
    int key;
    float TempMean[GmmNum][PCEP];
    int temp1[GmmNum];
    
    /*	for(i=0;i<GmmNum;++i)
     {
     for(j=0;j<PCEP;++j)
     cout<<OldMean[i][j]<<' ';
     cout<<endl;
     }*/
    while (iter<1000)
    {
        iter++;
        for (i = 0; i<GmmNum; ++i)
        {
            temp1[i] = 0;
            for (j = 0; j<PCEP; ++j)
                TempMean[i][j] = 0.0;
        }
        for (i = 0; i<TotalFrame; ++i)
        {
            min = infinite;
            for (j = 0; j<GmmNum; ++j)
            {
                sum = 0.0;
                for (k = 0; k<PCEP; ++k)
                {
                    index = i*PCEP + k;
                    temp = xishu[index] - OldMean[j][k];
                    sum += temp*temp;
                }
                if (sum<min)
                {
                    min = sum;
                    key = j;
                }
            }
            temp1[key]++;
            for (j = 0; j<PCEP; ++j)
            {
                index = i*PCEP + j;
                TempMean[key][j] += xishu[index];
            }
        }
        for (i = 0; i<GmmNum; ++i)
        {
            if (temp1[i] != 0)
            {
                for (j = 0; j<PCEP; ++j)
                    NewMean[i][j] = TempMean[i][j] / temp1[i];
            }
            if (temp1[i] == 0)
            {
                for (j = 0; j<PCEP; ++j)
                    NewMean[i][j] = OldMean[i][j];
            }
        }
        sum = 0.0;
        for (i = 0; i<GmmNum; ++i)
        {
            sum1 = 0.0;
            for (j = 0; j<PCEP; ++j)
            {
                temp = NewMean[i][j] - OldMean[i][j];
                sum1 += temp*temp;
            }
            sum1 = sqrt(sum1);
            sum += sum1;
        }
        if (sum<0.0001)
            break;
        for (i = 0; i<GmmNum; ++i)
            for (j = 0; j<PCEP; ++j)
            {
                OldMean[i][j] = NewMean[i][j];
                NewMean[i][j] = 0;
            }
    }
    cout << iter << endl;
    //	cout<<iter<<endl;
    /*	ofstream outfile2("out2.txt");
     for(i=0;i<GmmNum;++i)
     {
     outfile2<<temp1[i]<<endl;
     //	for(j=0;j<PCEP;++j)
     //		outfile2<<OldMean[i][j]<<' ';
     //	outfile2<<endl;
     }*/
    for (i = 0; i<GmmNum; ++i)
        temp1[i] = 0;
    
    for (i = 0; i<GmmNum; ++i)
        for (j = 0; j<PCEP; ++j)
            variance[i][j] = 0.0;
    
    for (i = 0; i<TotalFrame; ++i)
    {
        min = infinite;
        for (j = 0; j<GmmNum; ++j)
        {
            sum = 0.0;
            for (k = 0; k<PCEP; ++k)
            {
                index = i*PCEP + k;
                temp = xishu[index] - OldMean[j][k];
                sum += temp*temp;
            }
            if (sum<min)
            {
                min = sum;
                key = j;
            }
        }
        temp1[key]++;
        for (j = 0; j<PCEP; ++j)
        {
            index = i*PCEP + j;
            variance[key][j] += (xishu[index] - OldMean[key][j])*(xishu[index] - OldMean[key][j]);
        }
    }
    for (i = 0; i<GmmNum; ++i)
    {
        if (temp1[i]>0)
        {
            para[i].weight = ((float)temp1[i]) / TotalFrame;
            for (j = 0; j<PCEP; ++j)
                para[i].variance[j] = variance[i][j] / temp1[i];
        }
        else
        {
            para[i].weight = 1e-10;
            for (j = 0; j<PCEP; ++j)
                para[i].variance[j] = 0.001;
        }
        for (j = 0; j<PCEP; ++j)
            para[i].mean[j] = OldMean[i][j];
    }
}

void EM()
{
    float result1 = 0.0, result2;
    int TotalFrame = xishu.size() / PCEP;
    int i, j, k, index, count = 0;
    float temp, temp1, temp2, sum, sum1;
    double poslimit = 1e-100;
    float twopi = 8.0F*atan(1.0F);
    float mid[GmmNum], mid1[GmmNum], mid2[GmmNum][PCEP], mid3[GmmNum][PCEP];
    result2 = CalculateLikelihood();
    cout << result2 << endl;
    while (fabs(result1 - result2) / TotalFrame>0.0005)
    {
        count++;
        for (i = 0; i<GmmNum; ++i)
        {
            mid1[i] = 0.0;
            for (j = 0; j<PCEP; ++j)
            {
                mid2[i][j] = 0.0;
                mid3[i][j] = 0.0;
            }
        }
        for (i = 0; i<TotalFrame; ++i)
        {
            sum1 = 0.0;
            for (j = 0; j<GmmNum; ++j)
            {
                sum = 0.0;
                temp2 = 1.0;
                for (k = 0; k<PCEP; ++k)
                {
                    index = i*PCEP + k;
                    temp = xishu[index] - para[j].mean[k];
                    temp1 = temp*temp;
                    temp1 = temp1 / (2 * para[j].variance[k]);
                    sum -= temp1;
                    temp2 /= sqrt(para[j].variance[k]);
                }
                sum = exp(sum)*temp2 / sqrt(twopi);
                if (sum == 0)
                    sum = poslimit;
                mid[j] = sum*para[j].weight;
                sum1 += mid[j];
            }
            for (j = 0; j<GmmNum; ++j)
            {
                mid1[j] += mid[j] / sum1;
                for (k = 0; k<PCEP; ++k)
                {
                    index = i*PCEP + k;
                    mid2[j][k] += mid[j] * xishu[index] / sum1;
                    mid3[j][k] += mid[j] * (xishu[index] - para[j].mean[k])*(xishu[index] - para[j].mean[k]) / sum1;
                }
            }
        }
        for (i = 0; i<GmmNum; ++i)
        {
            para[i].weight = mid1[i] / TotalFrame;
            //		if(para[i].weight<0.0001)
            //			para[i].weight=0.0001;
            for (j = 0; j<PCEP; ++j)
            {
                para[i].mean[j] = mid2[i][j] / mid1[i];
                para[i].variance[j] = mid3[i][j] / mid1[i];
            }
        }
        result1 = result2;
        result2 = CalculateLikelihood();
        cout << result2 << endl;
    }
    cout << count << endl;
}

float CalculateLikelihood()
{
    int i, j, k, index;
    float temp, temp1, temp2, sum, sum1, sum2 = 0.0;
    double poslimit = 1e-100;
    float twopi = 8.0F*atan(1.0F);
    int TotalFrame = xishu.size() / PCEP;
    for (i = 0; i<TotalFrame; ++i)
    {
        sum1 = 0.0;
        for (j = 0; j<GmmNum; ++j)
        {
            sum = 0.0;
            temp2 = 1.0;
            for (k = 0; k<PCEP; ++k)
            {
                index = i*PCEP + k;
                temp = xishu[index] - para[j].mean[k];
                temp1 = temp*temp;
                temp1 = temp1 / (2 * para[j].variance[k]);
                sum -= temp1;
                temp2 /= sqrt(para[j].variance[k]);
            }
            sum = exp(sum)*temp2 / sqrt(twopi);
            if (sum == 0)
                sum = poslimit;
            sum1 += sum*para[j].weight;
        }
        //		cout<<sum1<<endl;
        sum2 += log(sum1);
        //		cout<<sum2<<endl;
    }
    return sum2;
}

void WritePara()
{
    
    int i, j;
    ofstream outfile("guzhen.txt");
    for (i = 0; i<GmmNum; ++i)
    {
        outfile << para[i].weight << ' ';
        outfile << endl;
        for (j = 0; j<PCEP; ++j)
        {
            outfile << para[i].mean[j] << ' ';
        }
        outfile << endl;
        for (j = 0; j<PCEP; ++j)
        {
            outfile << para[i].variance[j] << ' ';
        }
        outfile << endl;
        
    }
}


int get_mfcc(int mode)              //0 success
{
    extern ofstream tofile;
    extern char filename[128];
    extern bool flag;
    char path[128], name[128];
    char trainfilename[128]="./train/";
    char testfilename[128]="./test/";
    int end = 127;
    int start = 127;
    int namelen = 0;
    char ch;
    
    
    cout << "请输入要提取特征值的输入wav路径" << endl;
    cin >> path;
          
    if((path[0]==0)||
       ((path[0]=='Q'||path[0]=='K')&&path[1]==0)){
        flag = false;
        return 0;
    }
    
    end = strlen(path);
    start = end;
    while(path[start]!='/' && start!=0)start--;
    if(path[start]=='/') start++;
    namelen = end - start;
    strncpy(name, path + start , namelen+1);
    name[namelen-3]='t';
    name[namelen-2]='x';
    name[namelen-1]='t';
    
    xishu.clear();
    
    strcpy(filename, name);
    if(mode == 0||mode == 1){
        strcat(trainfilename, name);
        logfile << path <<"   "<< name<<"   ";
        return getmfcc(path, trainfilename);
    }
    else{
        strcat(testfilename, name);
        logfile << path <<"   "<< name<<"   ";
        return getmfcc(path,testfilename);
    }
}




