#define _CRT_SECURE_NO_WARNINGS
#include"hmm.h"
#include<iostream>
#include<cmath>
#include<math.h>
#include<fstream>
using namespace std;
#define MAXN 4000 //◊Ó∂‡ ˝æ›∏ˆ ˝
#define N 5//“˛≤ÿ
#define M 64//ø…π€≤Ï

//»´æ÷±‰¡ø ∂‘trainŒƒº˛kmeansæ€¿‡∫ÛµƒΩ·π˚
double kmeans[M][13];
void init_km()
{
	ifstream inkm;
    inkm.open("KMresult.txt", ios::in); //wenjianming
	for (int i = 0; i < M; i++)
		for (int j = 0; j < 13; j++)
			inkm >> kmeans[i][j];
	inkm.close();
}
//º∆À„¡Ω∏ˆœÚ¡ø÷Æº‰µƒ≈∑ œæ‡¿Î
double caldis(double x1[], double x2[])
{
	double sum = 0;
	for (int i = 0; i < 13; i++)
		sum += (x1[i] - x2[i])*(x1[i] - x2[i]);
	return sqrt(sum);
}
//∑µªÿœÚ¡øx1 Ù”⁄ƒƒ∏ˆ¿‡
int belongto(double x1[])
{
	double mind = caldis(x1, kmeans[0]);
	int minI = 0;
	for (int i = 1; i < 64; i++)
	{
		double tmpd = caldis(x1, kmeans[i]);
		if (tmpd < mind)
		{
			mind = tmpd;
			minI = i;
		}
	}
	return (minI + 1);
}
int hmm_test(int t, double xt[][13])
{
    extern ifstream tifile;
    extern int NFILES;
    
	init_km();
    int O[MAXN], id = -1;
	HMM * hmm[500];//最多500个wav
    char file_to_check[128]="./hmm/h";
    char namebuf[128];
    extern char pname[128];
    

    tifile.open("table.txt");
    tifile.getline(namebuf, 120);
    
	for (int f = 0; f < NFILES; ++f)
	{
        tifile >> id;
        tifile >> namebuf;
        tifile >> pname;
        
        
        strcat(file_to_check, namebuf);
		FILE* fp = fopen(file_to_check, "r");
		hmm[f] = new HMM;
		ReadHMM(fp, hmm[f]);
        file_to_check[7]=0;
        fclose(fp);
	}
    tifile.close();

	int i = 1;
	while (i <= t)
	{
		O[i] = belongto(xt[i]);
		i++;
	}
	int T_test = i;
	int q[MAXN] = { 0 };
	int ** psi = new int*[T_test + 1];
	for (int j = 0; j < T_test + 1; j++)
		psi[j] = new int[T_test + 1];
	double pf;
	double ** delta = new double*[T_test + 1];
	for (int j = 0; j < T_test + 1; j++)
		delta[j] = new double[T_test + 1];
	
	//’“µΩ◊Ó¥Ûµƒ∏≈¬ “‘º∞∂‘”¶µƒ◊Óø…ƒ‹πÏº£
	double pmax = 0;
	int answer = -1;
    cout << "呐咪呐咪哄……"<< endl;
    cout << "id probability:"<< endl;
    	for (int f = 0; f < NFILES; f++)
	{
		Viterbi(hmm[f], T_test - 1, O, delta, psi, q, &pf);
		cout << f << "   " << pf << endl;
		if (pf > pmax)
		{
			pmax = pf;
			answer = f;
			//£∫√ø∏ˆƒ£–Õœ¬∂º”√viterbi«Û◊Ó¥Û∏≈¬ £¨‘⁄≤ªÕ¨ƒ£–Õ÷–—°‘Ò∏≈¬ ◊Ó¥Ûµƒ
		}
	}

    tifile.open("table.txt");
    tifile.getline(namebuf, 120);
    
    for (int f = 0; f <= answer; ++f)
    {
        tifile >> id;
        tifile >> namebuf;
        tifile >> pname;
    }
    tifile.close();
    
    
    return answer;

}
void hmm_train(char * infile, char* hmmfile)
{
	init_km();
	int T = 0;//ø…π€≤Ï◊¥Ã¨–Ú¡–÷–◊¥Ã¨ ˝
	int O[MAXN];
	double x[MAXN][13];
	int sample_num = 0;

	HMM * hmm[2];
	hmm[0] = new HMM;
	hmm[1] = new HMM;

    FILE * fp = fopen(infile, "r");//单个训练

	//int weigh[13] = { 1,2,1,1,1,1,2,1,2,1,1,1,1 };
	while (fscanf(fp,"%d",&T)>0)//
	{
		if (T == 0)
			break;
		T++;
		sample_num++;
		//for (int k = 0; k < 13; k++) cin >> x[0][k];//∆˙÷√≤ª”√£¨»›“◊”–‘Î…˘
		for (int i = 1; i< T; ++i)
		{
			for (int k = 0; k < 13; k++)
			{
				 fscanf(fp,"%lf",&(x[i][k]));
			}
			O[i] = belongto(x[i]);
		}
        
		double ** alpha = new double*[T + 1];
		for (int j = 0; j < T + 1; j++)
			alpha[j] = new double[T + 1];
		double ** beta = new double*[T + 1];
		for (int j = 0; j < T + 1; j++)
			beta[j] = new double[T + 1];
		double ** gamma = new double*[T + 1];
		for (int j = 0; j < T + 1; j++)
			gamma[j] = new double[T + 1];

		double pi, pf;
		int pniter;

		int f = sample_num == 1 ? 0 : 1;
		InitHMM(hmm[f], N, M, 1);
		BaumWelch(hmm[f], T - 1, O, alpha, beta, gamma, &pniter, &pi, &pf);
		delete alpha;
		delete beta;
		delete gamma;

		//Ω´–¬hmm÷–µƒA,B,pi÷µº”µΩchmmµƒ∂‘”¶Œª÷√…œ
		if (f == 0)
			continue;

		for (int n = 1; n <= N; ++n)
		{
			//add A
			for (int n1 = 1; n1 <= N; ++n1)
			{
				hmm[0]->A[n][n1] += hmm[1]->A[n][n1];
			}
			//add B
			for (int m = 1; m <= M; ++m)
			{
				hmm[0]->B[n][m] += hmm[1]->B[n][m];
			}
			//add pi
			hmm[0]->pi[n] += hmm[1]->pi[n];
		}
        
	}
    fclose(fp);
	

	//«Û≥ˆ∆Ωæ˘µƒ◊™“∆æÿ’Û∫Õ∏≈¬ ∑÷≤ºœÚ¡ø
	for (int n = 1; n <= N; ++n)
	{
		//avg A
		for (int n1 = 1; n1 <= N; ++n1)
		{
			hmm[0]->A[n][n1] /= sample_num;
		}
		//avg B
		for (int m = 1; m <= M; ++m)
		{
			hmm[0]->B[n][m] /= sample_num;
		}
		//avg pi
		hmm[0]->pi[n] /= sample_num;
	}

	fp = fopen(hmmfile, "w");
	PrintHMM(fp, hmm[0]);
    fclose(fp);

}