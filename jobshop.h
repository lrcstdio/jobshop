//
// Created by huangzhen on 18-5-14.
//
#ifndef JOBSHOP_H
#define JOBSHOP_H
#ifndef JSP_EVOLUTION_H
#define JSP_EVOLUTION_H

const double r1=5;
const double r2=-0.42;
//交叉比例，种群数量，迭代次数,最大工件数，最大机器数，最大工序数
const double cr = 0.5;
const int ps = 150;
const int limit=500000;
const int maxn=20+1,maxm=10+1;
const int maxp=maxn*maxm;
//分别表示对应机器，所需时间,对应工件，在对应工件中的相对工序
struct process{
    int m,t,i,j;
};
class jobshop
{
public:
    int n,m,len,ans;                 //工件数，机器数，工序数，最小完成时间
    double ttime;                           //程序运行时间
    process pc[maxp];                //工序数组
    int genes[ps+2][maxp];                  //基因组
    int tgenes[ps+2][maxp];                 //交叉过程中存储父代
    int theBestGene[maxp];           //最优基因
    int  startpro[maxn];                    //每个工件开始的工序号
    bool book1[maxp],book2[maxp];           //记录交叉过程中是否重复
    int TMJ[maxm][maxp];
    int TT[maxm][maxp];              //解码矩阵
    int b[maxn],fa[maxp],d[maxm],p[maxn];   //辅助数组
    double q[ps];                           //适值数组
    double mr1,mr2;                         //两种变异概率
    int cur,st,et;                  //检修开始时间和结束时间
public:
    jobshop();                 //构造函数
    int decode(int *);         //解码
    void lsort(int *);         //局部排序
    void input();              //输入
    void init();               //初始化种群基因组
    void crossover();          //交叉
    void mutate();             //变异
    void rebuild();            //重建
    void output(double);       //输出
    void iter();               //迭代
};

#endif // JOBSHOP_H
#endif // JSP_EVOLUTION_H
