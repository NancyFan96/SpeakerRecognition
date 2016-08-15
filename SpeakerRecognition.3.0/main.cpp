#define _CRT_SECURE_NO_WARNINGS
//#include <Windows.h>
#include <iostream>
#include <fstream>
//#include <conio.h>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <string>
#include <vector>
#include <list>
#include "hmm.h"
#include "kmeans.h"
#include "mfcc.h"

using namespace std;


//#define DEBUG
//#define DEBUG_K
//#define INOUT
//#define LONGLOG
//#define VAD

#define min(a,b) ((a<b)? a:b)
#define MAXN 4000 //最多数据个数
char cmode;
int mode = 0; //全局变量，标记当前是训练模式还是测试模式：true-测试，false——训练
int id = 1;

ofstream logfile;           //log文档专用
ofstream tofile, OFILE;     //t- 对照列表汇总专用
ifstream tifile, IFILE;
char trainfilename[128]="./train/";
char hmmfilename[128]="./hmm/h";
char testfilename[128]="./test/";
char filename[128];
char pname[128];
int NFILES = 0;
bool flag = true; //执行标识
char ch;

//主函数
int main(int argc, char *argv[])
{
#ifdef INOUT
    cout<<"请输入模式： 0-keams, 1-train, 2-test"<<endl;
    cin >> ch;
    if(ch == '1'){
        freopen("./io/intrain.txt", "r", stdin);
        freopen("./io/outtrain.txt", "w", stdout);

    }
    if(ch =='2'){
        freopen("./io/intest.txt", "r", stdin);
        freopen("./io/outtest.txt", "w", stdout);
    }
    if(ch =='0'){
        freopen("./io/inkmeans.txt", "r", stdin);
        freopen("./io/outkmeans.txt", "w", stdout);
    }
    
#endif
    
#ifdef LONGLOG
    logfile.open("mylog.txt",ios::app|ios::out);
#endif
#ifndef LONGLOG
    logfile.open("mylog.txt",ios::out);
#endif
    if (!logfile) {
        cout << "不能打开日志文件 mylog.txt！" << endl;
        return 9;
    }
    time_t t = time(0);
    char tmptime[64];
    strftime( tmptime, sizeof(tmptime), "%Y/%m/%d %X %A",localtime(&t) );
    logfile << tmptime << endl;
    logfile << "样本wav名  txt文件名  行数" << endl;

    
    
    
    /********************************* HMM train and test START *********************************/
    
    srand(time(NULL));

    while(1){
#ifdef VAD
        while(get_mfcc(0)!=0){
            cout <<"MFCC提取错误!"<< endl;
            cin.clear();
            cin.sync();
        }
        if(flag == false) {
            cout <<"MFCC exit！"<<endl;
            return 0;
        }
#endif
        
        cout << endl << endl;
        cout << "0:    训练准备，Kmeans"<<endl;
        cout<<  "1:    训练模式"<<endl;
        cout<<  "2:    测试模式"<<endl;
        cout<<  "Q:    退出"<< endl;
        cout << endl;
        cin >> cmode;
        //训练文件，内有13维特征向量
        //HMM模型
        
    
        flag = true;
        /********************************* 合并MFCC文件，由K-means聚类 ********************************/
        if (cmode == '0')
        {
            mode = cmode - '0';
            long long totalLineN = 0;
            cout <<"训练准备！ 输入K，将生成K-means聚类"<< endl;
            cout <<"是否重新写入Kmeans_src.txt?Y/N"<<endl;
            cin >> ch;
            if(ch=='Y'){
                OFILE.open("kmeans_src.txt");
                OFILE.close();
            }
            //else
              //  OFILE.open("kmeans_src.txt",ios::app|ios::out);
            
            
            while(1){
                
                // filename 为样例对应的mfcc文件 同名txt
                while(get_mfcc(0)!=0){
                    cout <<"MFCC提取错误!"<< endl;
                    cin.clear();
                    cin.sync();
                }
                
                if(flag == false) {
                    OFILE.close();
                    do_kmeans();
                    cout <<"K-means聚类完成！"<<endl;
                    cout <<"退出准备模式"<<endl;
                    break;
                }
                
                
                char linebuf[1024];
                int lineN = 0;
               
                OFILE.open("kmeans_src.txt",ios::app|ios::out);
                strcat(trainfilename, filename);
                IFILE.open(trainfilename);
                if(!IFILE){
                    cout<<"不能添入"<<trainfilename<<endl;
                    return 8;
                }
                IFILE >> lineN;
                IFILE.get();
                /*cout<<"是否自定义汇入规模（推荐值1000）：Y/N?"<<endl;
                cin >> ch;
                if(ch=='Y'){
                    cout << "输入自定义汇入规模（推荐值1000）："<<endl;
                    cin >>lineN;
                }
                 */
                
                lineN = min(1000,lineN);
                while((IFILE.getline(linebuf, 1023)))
                    OFILE<<linebuf<<endl;
                    
                    
                IFILE.close();
                OFILE.close();
                trainfilename[8]=0;
                totalLineN += lineN;
                
                cout << "kmeans_src.txt已更新, 总行数："<< totalLineN << endl;                
            }
            
            
        }
        else if(cmode == '1')//训练模式
        {
            mode = cmode - '0';
            cout <<"训练ing！ 输入Q退出训练模式"<< endl;
            cout <<"是否重置训练库 table.txt?Y/N"<<endl;
            cin >> ch;
            if(ch=='Y'){
                tofile.open("table.txt");
                tofile << "id  filename  name" << endl;
                tofile.close();
            }
            else{
                char bbuf[128];
                char bbuf_cpy[128];
                tifile.open("table.txt");
                while(tifile.getline(bbuf, 128)){
                    strcpy(bbuf_cpy, bbuf);
                }
                sscanf(bbuf_cpy, "%d",&NFILES);
                NFILES++;
            }
           
            while(1){                
                // filename 为样例对应的mfcc文件 同名txt
                while(get_mfcc(1)!=0){
                    cout <<"MFCC提取错误!"<< endl;
                    cin.clear();
                    cin.sync();
                }
                
                if(flag == false) {
                    cout << "当前训练库规模:" << NFILES << endl;
                    cout <<"退出训练模式"<<endl;
                    break;
                }
               
                
                strcat(hmmfilename, filename);
                strcat(trainfilename, filename);
                hmm_train(trainfilename, hmmfilename);
                
                cout << "请输入人名"<<endl;
                cin >> pname;
            
                tofile.open("table.txt",ios::app|ios::out);
                tofile << NFILES << "   " << filename<< "   " << pname<<endl;
                tofile.close();
                cout<<"学习了ID为 "<<NFILES<<"的模型（Y^-^Y）"<< endl;
                
                NFILES++;
                //trainfilename="./train/"
                //hmmfilename = "./hmm/h";
                trainfilename[8]=0;
                hmmfilename[7]=0;
            }
        }
        else if(cmode == '2')//测试模式
        {
            mode = cmode - '0';
            cout <<"测试ing~; 输入Q退出测试模式"<< endl;
            cout << "请输入训练库规模："<< endl;
           
                cin>>NFILES;
          
            

            while(1){
                double xm[MAXN][13];
                int lineN = 0;
                int id = 0;
                
                
                
                // filename 为样例对应的mfcc文件 同名txt
                while(get_mfcc(mode)!=0){
                    cout <<"MFCC提取错误!"<< endl;
                    cin.clear();
                    cin.sync();
                }
                if(flag == false) {
                    cout <<"退出测试模式【猜对的话再来找我哦"<<endl;
                    break;
                }
               
                
                strcat(testfilename, filename);
                IFILE.open(testfilename);
                if(!IFILE){
                    cout<<"不能打开"<<testfilename<<endl;
                    return 8;
                }
                IFILE >> lineN;
                cout << "当前测试音频的特征规模为: "<< lineN <<endl;
#ifndef INOUT
                cout << "是否自定义测试规模（推荐值100）：Y/N?"<<endl;
                cin >> ch;
                if(ch=='Y'){
                    cout << "输入自定义测试规模（推荐值100）："<<endl;
                    cin >>lineN;
                }
#endif
              
                
#ifdef INOUT
                //出于结果精度考虑
                lineN = min(120, lineN);
#endif
                for(int i = 0; i < lineN; i++){
                    for(int j = 0; j < 13; j++){
                        IFILE >>xm[i][j];
                    }
                }
                IFILE.close();
                testfilename[7]=0;//testfilename="./test/"
                id = hmm_test(lineN, xm);
                
                
                char name[128];
                int tmp = strlen(filename);
                strcpy(name, filename);
                while(name[tmp-1]!='.') tmp--;
                name[tmp-1]=0;
                
                cout << "我猜你 "<<name<<" 是 ┑(￣Д ￣)┍"<< endl;
                cout << "id = "<<id<<"   name :"<< pname << endl;
               
            }
        }
        /********************************* HMM train and test END *********************************/
       else if(cmode=='Q')
           break;
    }
    logfile<<endl<<endl;
    logfile.close();
    
#ifdef INOUT
    fclose(stdin);
    fclose(stdout);
#endif
    
    return 0;
}

