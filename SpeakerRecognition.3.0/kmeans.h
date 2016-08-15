//
//  kmeans.h
//  SpeakerRecognition.2.0
//
//  Created by Nancy Fan on 16/5/17.
//  Copyright © 2016年 Nancy Fan. All rights reserved.
//

#ifndef kmeans_h
#define kmeans_h

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include<time.h>
#define k 64//簇的数目
using namespace std;
//存放元组的属性信息
typedef vector<double> Tuple;//存储每条数据记录

double getDistXY(const Tuple& t1, const Tuple& t2);
int clusterOfTuple(Tuple means[], const Tuple& tuple);
double getVar(vector<Tuple> clusters[], Tuple means[]);
Tuple getMeans(const vector<Tuple>& cluster);
void print(const vector<Tuple> clusters[]);
void KMeans(vector<Tuple>& tuples);
int do_kmeans();



#endif /* kmeans_h */
