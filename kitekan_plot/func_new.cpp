#include <math.h>
#include <iostream>
#include <fstream>

//encourage, symp, teasing, unrelated, no action,

int func_num[10][1] = {{0},{1},{2},{3},{4},{5},{6},{7},{8}};

double sig(double factor,double a,double b,double c){
  return 1.0/(1.0+exp(a*(b*factor-c)));
}

double gauss(double factor, double a, double b, double c){
  return exp(-1.0*pow((a*factor-b),2.0)/c);
}

double inv_down(double factor, double a, double b){
  double A = pow(2.0,a);
  double B = A - 1.0;
  double C = b * factor + 1;
  double D = pow(C,a);
  return A/(B*D) - 1.0/B;
}

double inv_up(double factor,double a,double b){
  double A = pow(2.0,a);
  double B = A - 1.0;
  double C = b * factor + 1;
  double D = pow(C,a);
  return (D-1.0)/B;
}

double quadra_func(double factor, double a, double b, double c){
  return a*factor*factor + b*factor + c;
}

//1〜10種類ある予測関数番号をもらうと計算をする関数を作ります
double func(double factor, double mental, int func_num){
  double a,b,c;
  double ret;
  //time of trial ##
  switch(func_num){
    case 0:
      //time of trial happy
      a = 2.03E-06 * mental;
      b = -5.95E-05 * mental;
      c = 2.56E-03 * mental;
      ret = quadra_func(factor,a,b,c);
      break;
    case 1:
      //win-lose		
      a = 0.00013507* mental;
      b = 0.00013532* mental;
      c = 0.00205905* mental;
      ret = quadra_func(factor,a,b,c);
      break;
    case 2:
      //num_enc			
      a = -4.20E-04 * mental;
      b = -3.36E-05 * mental;
      c = 1.23E-03 * mental;
      ret = quadra_func(factor,a,b,c);
      break;
    case 3:
      //num_symp			
      a = 0.00469428 * mental;
      b = -0.00357062 * mental;
      c = 0.00137434 * mental;
      ret = quadra_func(factor,a,b,c);
      break;
    case 4:
      //num_teas			
      a = -0.03822286 * mental;
      b = 0.04032191* mental;
      c = 0.00111175* mental;
      ret = quadra_func(factor,a,b,c);
      break;
    case 5:
      //num_un-related			
      a = -0.04513776 * mental;
      b = 0.04547926 * mental;
      c = 0.00079288 * mental;
      ret = quadra_func(factor,a,b,c);
      break;
    case 6:
      //num_nobehave			
      a = 0.00221025 * mental;
      b = -0.0013246* mental;
      c = 0.00052344* mental;
      ret = quadra_func(factor,a,b,c);
      break;
    case 7:
      //motion_intensity			
      a = 0.01275412 * mental;
      b = -0.00737356 * mental;
      c = 0.00203659 * mental;
      ret = quadra_func(factor,a,b,c);
      break;
  }
  return ret;
}

double inv_norm(int num, double val){
  double norm_val[40] = {80,1,1,1,1,1,1,1};
  val = val*norm_val[num];
  return val;
}
