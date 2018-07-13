#include "jobshop.h"
#include <cstdio>
#include <iostream>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <QDebug>
using namespace std;
jobshop::jobshop()
{
    len=0;
    mr1=0.5;
    mr2=0.9;
    ans=0x7fffffff;
    st=0;
    et=0;
    memset(pc,0,sizeof(pc));
    memset(TMJ,0,sizeof(TMJ));
    memset(genes,0,sizeof(genes));
    memset(theBestGene,0,sizeof(theBestGene));
}
void jobshop::input()
{
    int j,k,a,b,cnt=0;//按照格式读入
    FILE *f;
    f=fopen("input.txt","r+");
    fscanf(f,"%d %d",&n,&m);
    for(k=0;k<n;k++)
    {
        j = 0;
        startpro[k] = cnt + 1;
        while (fscanf(f, "%d %d", &a, &b)) {
            //qDebug()<<a<<b;
            pc[++len].m = a;
            pc[len].t = b;
            pc[len].i = k;
            pc[len].j = j++;
            TMJ[a][len] = b;
            ++cnt;
            char ch;
            if ((ch = fgetc(f)) == EOF || ch == '\n')
                break;
        }
    }
    fclose(f);
}
int jobshop::decode(int *g)
{
    memset(TT,0,sizeof(TT));
    memset(d,0,sizeof(d));   //对应机器上的最大完成时间
    memset(p,0,sizeof(p));   //对应工件上的最大完成时间

    for(int j=0;j<len;j++)//枚举染色体中的每个基因
    {
        int k=g[j];
        TT[pc[k].m][k] = max(d[pc[k].m], p[pc[k].i]) + pc[k].t;//由机器约束和工件约束得到当前工序的结束时间
        d[pc[k].m] = p[pc[k].i] = TT[pc[k].m][k];//更新机器约束和工件约束
    }

    int maxt=0;
    for(int i=0;i<m;i++)
        maxt = max(maxt,d[i]);//得到当前解
    if(maxt<ans)//更新全局最优解
    {
        ans=maxt;
        memcpy(theBestGene,g,len*sizeof(int));
    }
    return maxt;
}
#define TIME 1
void jobshop::rebuild()
{
    srand((unsigned)time(NULL));
    for(int k=0;k<ps;k++)
    {
        memcpy(genes[k],theBestGene,len*sizeof(int));
        int t=rand()%(len*len);//以随机程度改变当前较优解
        while(t--) {
            int i = rand() % len;//随机选取两个基因互换
            int j = rand() % len;
            int tmp = genes[k][i];
            genes[k][i] = genes[k][j];
            genes[k][j]=tmp;
        }
        lsort(genes[k]);//部分排序保证解合法
    }
}
void jobshop::iter() {
    clock_t startTime;
    startTime=clock();//开始计时
    init();//初始化种群
    int index=0,preans=0,cnt=0;
    int rhd=5000;
    while((double)(clock()-startTime)/CLOCKS_PER_SEC<TIME && ++index<limit) {//在迭代次数内和时间范围内
        crossover();//交叉
        mutate();//变异

        if(preans==ans)
            cnt++;
        else {
            cnt=0;preans=ans;
        }
        if(cnt==rhd)//如果迭代rhd次没有更新全局最优解,则重建
        {
            rebuild();
            mr1=0.5;
            cnt=0;
        }
//        if(index%10000==0) output((double(clock()-startTime)/CLOCKS_PER_SEC));
    }
    output((double(clock()-startTime)/CLOCKS_PER_SEC));//输出结果
}
void jobshop::crossover()
{
    int nn=0,n1,n2,m1,m2,cnt=1;
    srand((unsigned)time(0));
    memcpy(tgenes,genes,sizeof(genes));

    double maxq=-100;
    for(int i=0;i<ps;i++)
    {
        q[i]=r1*exp(r2*decode(genes[i]));//计算种群每个染色体的适应度
        if(q[i]>maxq) {maxq=q[i];nn=i;}//寻找当代种群的最优解
    }
    memcpy(genes[0],tgenes[nn],len*sizeof(int));//最优解直接进入下一代
    int lhd=(int)ps*(1-cr);
    while(cnt<lhd)//选取部分直接进入下一代
    {
        do{
            n1=rand()%ps;
            n2=rand()%ps;
         } while(n1==n2);
    nn=(q[n1]>q[n2])?n1:n2;
    memcpy(genes[cnt],tgenes[nn],len*sizeof(int));
    lsort(genes[cnt++]);
    }
    while(cnt<ps)//种群互相交叉得到
    {
        do{
            n1=rand()%lhd;
            n2=rand()%lhd;
        } while(n1==n2);//选取两个个体

        memset(book1,0,sizeof(book1));
        memset(book2,0,sizeof(book2));

        do{
          m1=rand()%len;
          m2=rand()%len;
        }while(m1==m2);//顺序交叉选取交叉部分
        if(m1>m2) {int t;t=m1;m1=m2;m2=t;}//保证m1<m2

        for(int i=m1;i<=m2;i++)
        {
          genes[cnt][i]=genes[n1][i];//交叉部分直接复制
          for(int j=0;j<len;j++)
            if(genes[n2][j]==genes[n1][i]) {book2[j]=true;break;}//把已经复制的基因标记
          genes[cnt+1][i]=genes[n2][i];//同上
          for(int j=0;j<len;j++)
            if(genes[n1][j]==genes[n2][i]) {book1[j]=true;break;}
        }

        int j1=0,j2=0;
        for(int i=0;i<m1;i++)
        {
          while(book2[j2]) j2++;//依次复制胜余部分，保证相对顺序
          genes[cnt][i]=genes[n2][j2];j2++;
          while(book1[j1]) j1++;
          genes[cnt+1][i]=genes[n1][j1];j1++;
        }
        for(int i=m2+1;i<len;i++)
        {
          while(book2[j2]) j2++;//同上
          genes[cnt][i]=genes[n2][j2];j2++;
          while(book1[j1]) j1++;
          genes[cnt+1][i]=genes[n1][j1];j1++;
        }
        lsort(genes[cnt]);//部分排序保证解合法
        lsort(genes[cnt+1]);
        cnt+=2;
  }
}
void jobshop::mutate() {
    srand((unsigned)time(NULL));
    bool flag;

    double sum=0;
    for(int i=0;i<ps;i++) {sum+=q[i];q[i]*=ps;}
    if(mr1-0.3>0) mr1-=0.0001;//变异概率随着迭代次数减少
    for(int k=0;k<ps;k++){
        flag=false;
        if(q[k]>sum) flag=true;

        double p=rand()/double(RAND_MAX);
        if((flag&&p<mr1) || (!flag&&p<mr2))
        {
            int i=rand()%len;//随机选择两个位置交换
            int j=rand()%len;
            int tmp=genes[k][i];
            genes[k][i]=genes[k][j];
            genes[k][j]=tmp;
        }
        lsort(genes[k]);//重排保证解合法
    }
}
void jobshop::output(double t)
{
    ttime=t;
    decode(theBestGene);//对最有解解码
    FILE *f;
    f=fopen("output1.txt","w");

    for(int i=0;i<m;i++)
    {
        fprintf(f,"M%d", i);
        for(int j=0;j<len;j++)
        {
            int k=theBestGene[j];
            if(pc[k].m == i)
                fprintf(f," (%d,%d-%d,%d)",TT[i][k]-pc[k].t, pc[k].i, pc[k].j, TT[i][k]);//分别对应该工序的开始时间,所属工件
                                                                                        // 对应工件的工序和结束时间
        }
        fputc('\n',f);
    }
    fprintf(f,"Time Used:%0.3lfs\n",t);
    fprintf(f,"End Time:%d\n", ans);
    fclose(f);
}
void jobshop::init()
{
    srand((unsigned)time(NULL));
    for (int i = 0; i < len; ++i) {
        fa[i]=i+1;
    }//生成染色体模板
    for(int k=0;k<ps;k++)
    {
        memcpy(genes[k],fa,sizeof(fa));
        int t=len*len;
        while(t--) {//乱排保证初始种群分布均匀，即足够乱
            int i = rand() % len;//随机选择两个基因足够乱
            int j = rand() % len;
            int tmp = genes[k][i];
            genes[k][i] = genes[k][j];
            genes[k][j]=tmp;
        }
        lsort(genes[k]);//重新部分排序保证染色体合法
    }
}
void jobshop::lsort(int *a)
{
    memset(b,0,sizeof(b));//清0
    for(int i=0;i<len;i++)
    {
        int t=a[i];
        a[i]=(b[pc[t].i]++)+startpro[pc[t].i];//保证每个工件给个工序相对有序
    }
}



