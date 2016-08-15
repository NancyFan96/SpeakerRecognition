#include <iostream>
#include "kmeans.h"

int dataNum;//数据集中数据记录数目
int dimNum;//每条记录的维数

//计算两个元组间的欧几里距离
double getDistXY(const Tuple& t1, const Tuple& t2)
{
    double sum = 0;
    for (int i = 1; i <= dimNum; ++i)
    {
        sum += (t1[i] - t2[i]) * (t1[i] - t2[i]);
    }
    return sqrt(sum);
}

//根据质心，决定当前元组属于哪个簇
int clusterOfTuple(Tuple means[], const Tuple& tuple) {
    double dist = getDistXY(means[0], tuple);
    double tmp;
    int label = 0;//标示属于哪一个簇
    for (int i = 1; i<k; i++) {
        tmp = getDistXY(means[i], tuple);
        if (tmp<dist) { dist = tmp; label = i; }
    }
    return label;
}
//获得给定簇集的平方误差
double getVar(vector<Tuple> clusters[], Tuple means[]) {
    double var = 0;
    for (int i = 0; i < k; i++)
    {
        vector<Tuple> t = clusters[i];
        for (int j = 0; j< t.size(); j++)
        {
            var += getDistXY(t[j], means[i]);
        }
    }
    //cout<<"sum:"<<sum<<endl;
    return var;
    
}
//获得当前簇的均值（质心）
Tuple getMeans(const vector<Tuple>& cluster) {
    
    int num = cluster.size();
    Tuple t(dimNum + 1, 0);
    for (int i = 0; i < num; i++)
    {
        for (int j = 1; j <= dimNum; ++j)
        {
            t[j] += cluster[i][j];
        }
    }
    for (int j = 1; j <= dimNum; ++j)
        t[j] /= num;
    return t;
    //cout<<"sum:"<<sum<<endl;
}

void print(const vector<Tuple> clusters[])
{
    for (int lable = 0; lable<k; lable++)
    {
        cout << "第" << lable + 1 << "个簇：" << endl;
        vector<Tuple> t = clusters[lable];
        for (int i = 0; i<t.size(); i++)
        {
            cout << i + 1 << ".(";
            for (int j = 0; j <= dimNum; ++j)
            {
                cout << t[i][j] << ", ";
            }
            cout << ")\n";
        }
    }
}

void KMeans(vector<Tuple>& tuples) {
    vector<Tuple> clusters[k];//k个簇
    Tuple means[k];//k个中心点
    int i = 0;
    //cout << tuples.size() << endl;
    //一开始随机选取k条记录的值作为k个簇的质心（均值）
    srand(time(NULL));
    for (i = 0; i<k;) {
        int iToSelect = rand() % tuples.size();
        if (means[i].size() == 0)
        {
            for (int j = 0; j <= dimNum; ++j)
            {
                means[i].push_back(tuples[iToSelect][j]);
            }
            ++i;
        }
    }
    
    int lable = 0;
    //根据默认的质心给簇赋值
    for (i = 0; i != tuples.size(); ++i) {
        lable = clusterOfTuple(means, tuples[i]);
        clusters[lable].push_back(tuples[i]);
    }
    double oldVar = -1;
    double newVar = getVar(clusters, means);
    
#ifdef DEBUG_K
    cout << "初始的的整体误差平方和为：" << newVar << endl;
#endif
    
    int t = 0;
    while (abs(newVar - oldVar) >= 1) //当新旧函数值相差不到1即准则函数值不发生明显变化时，算法终止
    {
#ifdef DEBUG_K
        cout << "第 " << ++t << " 次迭代开始：" << endl;
#endif
        for (i = 0; i < k; i++) //更新每个簇的中心点
        {
            means[i] = getMeans(clusters[i]);
        }
        oldVar = newVar;
        newVar = getVar(clusters, means); //计算新的准则函数值
        for (i = 0; i < k; i++) //清空每个簇
        {
            clusters[i].clear();
        }
        //根据新的质心获得新的簇
        for (i = 0; i != tuples.size(); ++i) {
            lable = clusterOfTuple(means, tuples[i]);
            clusters[lable].push_back(tuples[i]);
        }
#ifdef DEBUG_K
        cout << "此次迭代之后的整体误差平方和为：" << newVar << endl;
#endif
    }
#ifdef DEBUG_K
    cout << "最终迭代之后的整体误差平方和为：" << newVar << endl;
    cout << "The result is:\n";
    system("pause");
#endif
    ofstream outfile;
    outfile.open("KMresult.txt", ios::out);
    for (int zz = 0; zz < 64; zz++)
    {
        for (int kk = 1; kk <= 13; kk++)
        {
            outfile << means[zz][kk] << " ";
        }
        outfile << endl;
    }
#ifdef DEBUG_K
    print(clusters);
#endif
}
int do_kmeans() {
    
    char fname[256]="kmeans_src.txt";
    
    //cout << "请输入存放数据的文件名： ";
    //cin >> fname;
    cout << "请依次输入: 维数dimNum 样本数目dataNum" << endl;
    cin >> dimNum >> dataNum;
    ifstream infile(fname);
    if (!infile) {
        cout << "不能打开输入的文件" << fname << endl;
        return 0;
    }
    vector<Tuple> tuples;
    //从文件流中读入数据  
    for (int i = 0; i<dataNum && !infile.eof(); ++i)
    {
        string str;
        getline(infile, str);
        istringstream istr(str);
        Tuple tuple(dimNum + 1, 0);//第一个位置存放记录编号，第2到dimNum+1个位置存放实际元素  
        tuple[0] = i + 1;
        for (int j = 1; j <= dimNum; ++j)
        {
            istr >> tuple[j];
        }
#ifdef DEBUG_K
        cout << i << endl;
#endif
        tuples.push_back(tuple);
    }
#ifdef DEBUG_K
    cout << endl << "开始聚类" << endl;
#endif
    KMeans(tuples);
    //system("pause");
    return 0;
}