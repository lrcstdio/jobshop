#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "jobshop.h"
#include <QImage>
#include <QTime>
struct Check{           //检修结构体
  int m;
  int st;
  int et;
};
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    jobshop js;                     //js实现各种jobshop功能
    int currentime;                 //加工当前时间
    QPixmap Pix1;                   //绘图容器
    QPixmap Pix2;
    QTime  timestart;               //模拟开始时间
    Check check[maxm];              //记录检修命令
    bool isEnd;                     //是否结束模拟过程
    QTimer *ot;                     //计时器
public:
    void updateTT();                //更新TT解码数组
    void output2();                 //输出output2
public slots:
    void on_loadorder_clicked();    //订单载入
    void printGantt();              //动态绘制甘特图
    void on_start_clicked();        //jobshop计算函数
    void on_model_clicked();        //模拟加工启动函数
private:
    Ui::MainWindow *ui;
    int t,cnt,count;    //辅助变量
    int ma,ta;          //检修机器和时间
    bool valid;         //是否为检修指令
    QString s;          //记录检修指令
protected:
    void paintEvent(QPaintEvent *event);            //重绘事件
    void keyPressEvent(QKeyEvent *event);           //键盘事件
};

#endif // MAINWINDOW_H
