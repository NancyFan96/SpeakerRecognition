/*
**      File:   backward.cpp
**      功能：给定观察值序列和HMM模型，利用前向后向算法
**            求取其概率
*/

#include <stdio.h>
#include "hmm.h"

/***************************************************************************
** 函数名称：ForwardWithScale
** 功能：后前算法估计参数（带比例因子修正）
** 参数：phmm：指向HMM的指针
**       T：观察值序列的长度
**       O：观察值序列
**       beta：运算中用到的临时数组
**       scale：比例因子数组
**       pprob：返回值，所要求的概率
*****************************************************************************/
void ForwardWithScale(HMM *phmm, int T, int *O, double **alpha, 
	double *scale, double *pprob)
{
	int     i, j;   /* 状态指示 */
	int     t;      /* 时间下标 */
	double sum;
 
 
	/* 1. 初始化 */
	scale[1] = 0.0;
	for(i = 1; i <= phmm->N; i ++)
	{
		alpha[1][i] = phmm->pi[i]*(phmm->B[i][O[1]]);
		scale[i] += alpha[1][i];
	}
	for (i = 1; i <= phmm->N; i++)
		alpha[1][i] = 1.0/scale[1]; 
 
	/* 2. 递归 */
	for (t = 1; t < T; t++) 
	{
		scale[t+1] = 0.0;
		for (j = 1; j <= phmm->N; j++) 
		{
			sum = 0.0;
			for (i = 1; i <= phmm->N; i++)
				sum += phmm->A[i][j] * alpha[t][i];
			alpha[t+1][j] = sum*(phmm->B[j][O[t+1]]);
			scale[t+1]+=alpha[t+1][j];
		}
		for(j = 1; j <= phmm->N; j ++)
			alpha[t+1][j]/=scale[t+1];
	}
		/* 3. 终止 */
	*pprob = 0.0;
	for (t = 1; t<T; t++)
		*pprob += log(scale[t]);

}
