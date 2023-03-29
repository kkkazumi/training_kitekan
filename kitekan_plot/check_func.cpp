#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

double func(double factor, double mental, int func_num);
double inv_norm(int num, double val);

#define FACT_TYPE 10
#define SIGNAL_TYPE 4
#define N 100

int main(int argc, const char* argv[]){//(argv[0]=./estimate),argv[1]=usrname

  std::string output_csv_file_path="func_check.csv";
  std::ofstream out_file(output_csv_file_path);

  int func_num[10][1] = {{0},{1},{2},{3},{4},{5},{6},{7},{8}};
  for(int sit_num=0;sit_num<8;sit_num++){
    for(int emo_num=0;emo_num<1;emo_num++){
      for(int m=0;m<3;m++){
        double m_conv = (double)m*0.5+0.01;
        printf("s-%d,s-%d,m-%f\n",sit_num,emo_num,m_conv);
        for(int f=0;f<100;f++){
          //f_val[f] = (double)f/100.0;
          double value=func(inv_norm(func_num[sit_num][emo_num],(double)f/100.0),m_conv,func_num[sit_num][emo_num]);//(1e+1);//因子から表情推定値を出す。mentalは基底関数のパラメータ
          //s_val[f]=value;
          if(f<N-1){
            out_file<<value<<",";
          }else{
            out_file<<value<<"\n";
          }
        }
      }
    }
  }
}
