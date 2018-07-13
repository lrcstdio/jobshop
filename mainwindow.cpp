#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFile>
#include <QTextStream>
#include <QPainter>
#include <QKeyEvent>
#include <QDebug>
#include <ctime>
#include <QImage>
#include <QPixmap>
#include <QTime>
#include <QTimer>
#include <cstdlib>
#include <iostream>
#include <algorithm>
using namespace std;
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    valid=false;
    cnt=0;
    ma=0;
    ta=0;
    t=0;
    count=0;
    isEnd=false;
    currentime=0;
    ot=new QTimer;
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::updateTT()
{
   /*
    TT为开始时间
    updateTT函数只需要在需要检修时调用一次即可
    updateTT函数需要传参
    检修前需要先对最优基因解码
    */
        memset(js.d,0,sizeof(js.d));
        memset(js.p,0,sizeof(js.p));
        for(int j=0;j<js.len;j++)
        {
            int k=js.theBestGene[j];
            js.TT[js.pc[k].m][k]-=js.pc[k].t;
        }
        bool flag =false;
        for(int j=0;j<js.len;j++)//枚举染色体中的每个基因-
        {
            int k=js.theBestGene[j];
            if(js.pc[k].m == js.cur && js.TT[js.pc[k].m][k]>=js.st && js.TT[js.pc[k].m][k]<js.et)
            {
                js.TT[js.pc[k].m][k]=js.et;
                flag=true;
            }
            else if(js.pc[k].m == js.cur && js.TT[js.pc[k].m][k] < js.st &&
                    js.TT[js.pc[k].m][k] + js.pc[k].t >js.st)
            {
                js.et+=js.TT[js.pc[k].m][k] + js.pc[k].t - js.st;
                js.st=js.TT[js.pc[k].m][k] + js.pc[k].t;
                check[count].st=js.st;
                check[count].et=js.et;
            }
            if(flag && js.TT[js.pc[k].m][k] < max(js.d[js.pc[k].m], js.p[js.pc[k].i]))
                js.TT[js.pc[k].m][k] = max(js.d[js.pc[k].m], js.p[js.pc[k].i]);
            js.d[js.pc[k].m]=js.p[js.pc[k].i]=js.TT[js.pc[k].m][k]+js.pc[k].t;
        }
        int maxt=0;
        for(int i=0;i<js.m;i++)
            maxt = max(maxt,js.d[i]);//得到当前解
        js.ans = maxt;
        for(int j=0;j<js.len;j++)
        {
            int k=js.theBestGene[j];
            js.TT[js.pc[k].m][k]+=js.pc[k].t;
        }
}

void MainWindow::output2()
{
    FILE *f;
    f=fopen("output2.txt","w");

    for(int i=0;i<js.m;i++)
    {
        int fixcur=-1,last=0;
        for(int k=1;k<=count;k++)
            if(check[k].m==i)
            {
                fixcur=k;
                break;
            }
        fprintf(f,"M%d", i);
        for(int j=0;j<js.len;j++)
        {
            int k=js.theBestGene[j];
            if(js.pc[k].m == i)
            {
                if(fixcur!=-1&&check[fixcur].et<=(js.TT[i][k]-js.pc[k].t))
                {
                    fprintf(f," (%d,检修,%d)",check[fixcur].st,check[fixcur].et);
                    fixcur=-1;
                }
                fprintf(f," (%d,%d-%d,%d)",js.TT[i][k]-js.pc[k].t, js.pc[k].i, js.pc[k].j,js.TT[i][k]);//分别对应该工序的开始时间,所属工件
                last=js.TT[i][k];
            }                                                                           // 对应工件的工序和结束时间
        }
        if(fixcur!=-1&&check[fixcur].st>=last)
            fprintf(f," (%d,检修,%d)",check[fixcur].st,check[fixcur].et);
        fputc('\n',f);
    }
    QTime tm;
    tm.start();
    fprintf(f,"Time Used:%.3lfs\n",(double)timestart.secsTo(tm));
    fprintf(f,"End Time:%d\n", js.ans);
    fclose(f);
}
void MainWindow::on_loadorder_clicked()
{
    js.input();
    QFile file("input.txt");
    ui->order_display->setReadOnly(true);
    QTextStream out(&file);
    if (file.open(QIODevice::ReadOnly)){
       ui->order_display->setText(out.readAll());
    }
    file.close();
    ui->loadorder->setEnabled(false);
}

void MainWindow::printGantt()
{
//    qDebug()<<"yes"<<currentime;
    currentime+=6;
    QPixmap pix(currentime,340);
    pix.fill(Qt::lightGray);
    QPainter painter;
    painter.begin(&pix);
    int y=-10;
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.setPen(QColor(Qt::black));
    double sc=1;
    painter.drawLine(QPoint(0,10),QPoint(0,340));
    painter.drawLine(QPoint(0,340),QPoint(1000,340));
    srand((unsigned)time(NULL));
    int last=0;
    for(int i=0;i<js.m;i++)
    {
        y+=30;
        for(int j=0;j<js.len;j++)
        {
            int k=js.theBestGene[j];
            if(js.pc[k].m == i)
             {
                  int cntp=rand()%3;
                  while(cntp==last)
                      cntp=rand()%3;
                  last=cntp;
                  if(cntp==1)
                        painter.setBrush(QColor(255,255,0));
                  else if(cntp==2)
                        painter.setBrush(QColor(0,255,0));
                  else
                        painter.setBrush(QColor(0,255,255));
                  int start=js.TT[i][k]-js.pc[k].t;
                  int end=js.TT[i][k];
                  QRectF rect((int)start*sc,y,(int)(end-start)*sc,20);
                  painter.drawRect(rect);
                  QRectF rec(-10+(int)start*sc,y,(int)(end-start)*sc+20,20);
                  painter.drawText(rec, Qt::AlignHCenter|Qt::AlignVCenter,
                  QString::number(js.pc[k].i)+"-"+QString::number(js.pc[k].j));
             }
        }
    }
    painter.setBrush(QColor(255,0,0));
    for(int i=0;i<=count;i++)
    {
        if(check[i].et<=currentime)
        {
            QRectF rect(check[i].st,(check[i].m+1)*30-10,check[i].et-check[i].st,20);
            painter.drawRect(rect);
        }
        else
        {
            if(check[i].st<currentime)
            {
                QRectF rect(check[i].st,(check[i].m+1)*30-10,currentime-check[i].st,20);
                painter.drawRect(rect);
            }
        }
    }
    painter.end();
    Pix2=pix;
    update();
    if(currentime>=js.ans)
    {
        ot->stop();
        ui->showfixans->setText(QString::number(js.ans));
        output2();
    }
}
void MainWindow::on_start_clicked()
{
    js.iter();
    QPixmap pix(1300,340);
    pix.fill(Qt::lightGray);
    QPainter painter;
    painter.begin(&pix);
    int y=-10;
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.setPen(QColor(Qt::black));
    js.decode(js.theBestGene);
    double sc=1;
    painter.drawLine(QPoint(0,10),QPoint(0,340));
    painter.drawLine(QPoint(0,340),QPoint(1300,340));
    srand((unsigned)time(NULL));
    int last=0;
    for(int i=0;i<js.m;i++)
    {
        y+=30;
        for(int j=0;j<js.len;j++)
        {
            int k=js.theBestGene[j];
            if(js.pc[k].m == i)
             {
                  int cntp=rand()%3;
                  while(cntp==last)
                      cntp=rand()%3;
                  last=cntp;
                  if(cntp==1)
                        painter.setBrush(QColor(255,255,0));
                  else if(cntp==2)
                        painter.setBrush(QColor(0,255,0));
                  else
                        painter.setBrush(QColor(0,255,255));
                  int start=js.TT[i][k]-js.pc[k].t;
                  int end=js.TT[i][k];
                  QRectF rect((int)start*sc,y,(int)(end-start)*sc,20);
                  painter.drawRect(rect);
                  QRectF rec(-10+(int)start*sc,y,(int)(end-start)*sc+20,20);
                  painter.drawText(rec, Qt::AlignHCenter|Qt::AlignVCenter,
                  QString::number(js.pc[k].i)+"-"+QString::number(js.pc[k].j));
             }
        }
    }
    painter.end();
    Pix1=pix;
    update();
    ui->time_used->setText(QString::number(js.ttime));
    ui->end_time->setText(QString::number(js.ans));
}

void MainWindow::paintEvent(QPaintEvent *event)
{
    Q_UNUSED(event);
    QPainter painter(this);
    painter.drawPixmap(60,10,Pix1);
    painter.drawPixmap(60,420,Pix2);
}
void MainWindow::keyPressEvent(QKeyEvent *event)
{
   qDebug()<<event->text();
   if(event->key()==Qt::Key_C)
   {
       s+="C";
       valid=true;
   }
   else if(valid&&(event->key()==Qt::Key_Enter||event->key()==Qt::Key_Return))
   {
       s+="\n";
       valid=false;
       t=0;
       cnt=0;
       //qDebug()<<"C"<<ma<<ta;
       QTime ttime;
       ttime.start();
       js.st=timestart.secsTo(ttime)*3;
       js.cur=ma;
       js.et=js.st+ta;
       check[++count].m=ma;
       check[count].st=js.st;
       check[count].et=js.et;
       ma=0;
       ta=0;
//       qDebug()<<"start_fix"<<js.cur<<js.st<<js.et;
       updateTT();
   }
   else if(event->key()==Qt::Key_Shift)
   {
       t=0;
       cnt++;
       s+=" ";
   }
   else
   {
       t=event->text().toInt();
       if(cnt==1)
          ma=t;
       else
       {
           ta=ta*10+t;
       }
       s+=event->text();
   }
   ui->fix_display->setText(s);
}

void MainWindow::on_model_clicked()
{
    timestart.start();
    connect(ot,SIGNAL(timeout()),this,SLOT(printGantt()));
    ot->start(2000);
}
