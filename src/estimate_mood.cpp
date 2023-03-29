#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <cmath>

#include <Eigen/Core>
#include <Eigen/LU>

#include <sys/stat.h>
#include <sys/types.h>

#define N 2
#define MAX_TRAIN_DATA_LEN 145

using namespace std;
using namespace Eigen;

int minElement(double* array, size_t size);

void pseudo_inverse(MatrixXd *O, MatrixXd P);

void read_file(char *filename, int size2, double **data);

double func(double factor, double mental, int func_num);
double func2(double factor, double mental, int func_num);
double inv_norm(int num, double val);

typedef Matrix<double,N,8,RowMajor> Mat_phi;

void setPhi(double (&factor)[N][8],double mental[N], int emo_num, Mat_phi *Phi);
void get_para(int id_num, int set_num, int t, MatrixXd *W1_para,char *dir_path);
void f_mm(MatrixXd W1_para, double (&factor)[N][8], double mental[N], double *output);
//mental model
double fE(MatrixXd W1_para, double (&factor)[N][8], double mental[N], double sig_total[4]);

void out_tmppara(char *filename, MatrixXd W1_para);

void grad_descent(double M_org[N], double* M,MatrixXd W1_para,double (&factor)[N][8], double sig_total[1],int id_num,int test_num,int moto_count,char *dir_path,int set_num);

int main(int argc, const char* argv[]){//(argv[0]=./estimate),argv[1]=usrname

  double *_factor[N],*_signal[N];
  double factor_size[N][8],signal_size[N][1];//データ数は１５個、要素数は７こ、
  double factor[N][8], signal[N][1];
  for(int i=0;i<N;i++) _factor[i] = factor_size[i];
  for(int i=0;i<N;i++) _signal[i] = signal_size[i];

  char factor_name[100];
  char signal_name[100];
  char out_name[100];
  char face_name[100];
  char err_log[100];
  char dir_path[100];

  //for(int id=0;id<9;id++){ //commented out for debug 20230306
    int id = 0;//for debug 20230306
    int id_num = id+1;
    std::cout<<"id_num: "<<id_num<<std::endl;
    for(int _set_num =0; _set_num<4; _set_num++){
      int set_num = 5*(_set_num+1);
      for(int t=0; t<MAX_TRAIN_DATA_LEN; t++){
        sprintf(dir_path,"../sep_data/%d/",id_num);

        sprintf(factor_name,"%sfactor_test%d-%d.csv",dir_path,set_num,t);
        sprintf(signal_name,"%sface_test%d-%d.csv",dir_path,set_num,t);
        sprintf(out_name,"%soutput_test/TESTestimated_%d-%d.csv",dir_path,set_num,t);//コマンドライン引数のmental.csvのタイトルをファイル名にする
          //sprintf(out_name,"%soutput/TESTestimated_%d-%d.csv",dir_path,set_num,t);//コマンドライン引数のmental.csvのタイトルをファイル名にする//commented out for debug 20230306

        printf("debug:before readfile;%s\n",factor_name);

        read_file(factor_name, 8,_factor);
        read_file(signal_name, 4,_signal);

        memcpy(factor,*_factor,sizeof(factor));
        memcpy(signal,*_signal,sizeof(signal));

        double m_org[N];//Mデータは1次元１５このでーた
        double mental[N];//Mデータは1次元１５このでーた
        //srand((unsigned) time(NULL));
        for(int ment_num=0;ment_num<N;ment_num++){
          m_org[ment_num]={0};//初期値random
          //m_org[ment_num]=(rand()%100-50.0)*2;//初期値random
        }
        memcpy(mental,m_org,sizeof(m_org));

        double dEdM[N] = {0};
        int errflg=0;
        double err=10;

        std::cout << std::fixed;
        int count=0;

        char tmp_w_filename[50];
        sprintf(tmp_w_filename,"%soutput_test/weight_zero_trial4-2_%d-%d.csv",dir_path,set_num,t);
        char mental_out_file[50];
        sprintf(mental_out_file,"%soutput_test/mental_log_trial4-2_%d-%d.csv",dir_path,set_num,t);

        ofstream init_out(tmp_w_filename,ios::app);
        ofstream mental_out(mental_out_file,ios::app);
        ofstream write_mood(out_name,ios::app);

        init_out<<id_num<<","<<count<<"==================="<<std::endl;

        MatrixXd W1_para(8,1);
        //int t =0;
        printf("debug:set_num%d\n",set_num);
        get_para(id_num,set_num,t,&W1_para,dir_path);

        double sig_total[1]={0};
        for(int emo_num=0;emo_num<1;emo_num++){
          for(int col_num1 = 0; col_num1<N; col_num1++){
            sig_total[emo_num] += signal[col_num1][emo_num];
          }
        }

        int count_mgrad=0;
        double E_pre =fE(W1_para,factor,m_org,sig_total);

        grad_descent(m_org,mental,W1_para,factor,sig_total,id_num,t,count,dir_path,set_num);

        //S誤を求める
        err = fE(W1_para, factor, mental,sig_total) - E_pre;

        init_out<<id_num<<","<<","<<count<<"==================="<<std::endl;
        //out_tmppara(tmp_w_filename, W1_para);
        //err_ofs<<err<<std::endl;

        memcpy(m_org,mental,sizeof(mental));
        if(err<=0.00000001){
          errflg=1;
        }else if(count>=5000000){//KENSYU: loop count
          errflg=1;
        }
        for(int i=0;i<N;i++){
          write_mood<<mental[i]<<std::endl;
        }
      }
    }
  //} //commented out for debug 20230306 (to for loop in user id)
}

void read_file(char *filename, int size2, double **data){
    int i,j;
    double hoga;
    double hage;

    std::string data_f[size2];
    std::stringstream ss;

    std::string line;

    std::ifstream ifs;

    ifs.open(filename);
    j=0;
    while(getline(ifs,line)){
      std::stringstream ss(line);
      for(i=0;i<size2;i++){
        getline(ss,data_f[i],',');
        hoga=std::stod(data_f[i]);
        hage=std::stod(data_f[i]);
        data[j][i]=hage;
      }
      j++;
    }

}

void grad_descent(double M_org[N], double* M,MatrixXd W1_para,double (&factor)[N][8], double sig_total[1],int id_num,int test_num,int moto_count,char *dir_path,int set_num){
  double dEdM[N]={0};
  double err=100;
  double Msize[N];
  int count=0;
  char mental_out_file[50];
  sprintf(mental_out_file,"%soutput_test/mental_log_trial5_%d-%d.csv",dir_path,set_num,test_num);//for debug 20230306
  //sprintf(mental_out_file,"%soutput/mental_log_trial5_%d-%d.csv",dir_path,set_num,t);//commented out for debug 20230306
  ofstream mental_out(mental_out_file,ios::app);

  while(err>0.0000001){
    double E_pre  = fE(W1_para, factor, M_org,sig_total);
    //std::cout<<"E_pre: "<<E_pre<<", ";
    for(int i=0;i<N;i++){
      double tmp_M = M_org[i];
      double E_cand[3] = {0};
      double delta[3] = {0};
      for(int dm = -1; dm<2;dm++){
        M[i] = tmp_M + dm*10;
        //E_cand[dm+1] = fE(M);
        E_cand[dm+1] = fE(W1_para, factor, M,sig_total);
        //std::cout<<"E_cand["<<dm+1<<"]: "<<E_cand[dm+1]<<std::endl;
        delta[dm+1] = (E_cand[dm+1] - E_pre);
      }
      int index = minElement(delta,3);
      dEdM[i] = -delta[index];
      M[i]=tmp_M;
    }
    err=0;
    for(int j=0;j<N;j++){
      M[j] = M_org[j]+dEdM[j];
      err+=(dEdM[j])*(dEdM[j])/2.0;
    }
    //std::cout<<std::setprecision(4)<<"err:"<<err<<std::endl;
    memcpy(M_org,M,sizeof(Msize));

    cout<<"gd_descent("<<id_num<<"-"<<test_num<<"-"<<moto_count<<"),"<<std::setprecision(10)<<err<<","<<count<<",";
    mental_out<<"gd_descent("<<id_num<<"-"<<test_num<<"-"<<moto_count<<"),"<<std::setprecision(10)<<err<<","<<count<<",";
    for(int m_check=0;m_check<N;m_check++){
      cout<<std::setprecision(2)<<M[m_check]<<", ";
      mental_out<<std::setprecision(2)<<M[m_check]<<", ";
    }
    cout<<std::endl;
    mental_out<<std::endl;
    count++;
  }
  cout<<std::endl;
}



int minElement(double* array, size_t size){
  int min = array[0];
  int index = 0;
  //std::cout<<"element check"<<std::endl;
  //std::cout<<array[0]<<",";
  for (size_t i = 1; i < size; ++i) {
  //  std::cout<<array[i]<<",";
    if (min < array[i]) {
      min = array[i];
      index = i;
    }
  }
  if((abs(array[0] - array[1])<0.0001) && (abs(array[2] - array[1])<0.0001))  index = 1;
  //std::cout<<"index: "<<index<<std::endl;
  return index;
}


void setPhi(double (&factor)[N][8],double mental[N], int emo_num, Mat_phi *Phi){

  Mat_phi Phi1;//計画行列のこと
  double m_conv = 0.0;

  int func_num[8][1] = {{0},{1},{2},{3},{4},{5},{6},{7}};
  //int func_num[10][4]={{0,1,2,3},{4,5,6,7},{8,9,10,11},//もとのやつパート１
  //    {12,13,14,15},{16,17,18,19},{20,21,22,23},//7つの状況に対する予測関数、各情動4つ分ずつある
  //    {24,25,26,27},{28,29,30,31},{32,33,34,35},{36,37,38,39}};//もとのやつパート２

  for(int data_num=0;data_num<N;data_num++){
    m_conv = 9.98/(1.0 + exp(-0.5 * mental[data_num]))+0.01;
    for(int sit_num=0;sit_num<8;sit_num++){
      if(emo_num==0){
        Phi1(data_num,sit_num)=func(inv_norm(func_num[sit_num][emo_num],factor[data_num][sit_num]),m_conv,func_num[sit_num][emo_num]);//(1e+1);//因子から表情推定値を出す。mentalは基底関数のパラメータ
      }else if(emo_num==1){
        Phi1(data_num,sit_num)=func(inv_norm(func_num[sit_num][emo_num],factor[data_num][sit_num]),m_conv,func_num[sit_num][emo_num]);//(1e+1);
      }else if(emo_num==2){
        Phi1(data_num,sit_num)=func(inv_norm(func_num[sit_num][emo_num],factor[data_num][sit_num]),m_conv,func_num[sit_num][emo_num]);//(1e+1);
      }else if(emo_num==3){
        Phi1(data_num,sit_num)=func(inv_norm(func_num[sit_num][emo_num],factor[data_num][sit_num]),m_conv,func_num[sit_num][emo_num]);//(1e+1);
      }
      //多分ここは大丈夫。
    }
  }
  *Phi = Phi1;
}

void get_para(int id_num, int set_num, int t , MatrixXd *W1_para,char *dir_path){

  char hap_in[100];//,sup_in[100],ang_in[100],sad_in[100];
  double ment_est_para[8][1] = { 0 };//現在のM推定するときに使う(mode=1)

	sprintf(hap_in,"%soutput/hap_weight%d-%d.csv",dir_path,set_num,t);//for debug 20230306
	//sprintf(hap_in,"%soutput/hap_weight%d-%d.csv",dir_path,set_num,t);//commented out for debug 20230306

  ifstream hap_model_in(hap_in);//happyのおもみを読み込み

  MatrixXd W1(8,1);

  for (int i = 0; i< 8; i++){
    hap_model_in >> ment_est_para[i][0];
    std::cout<<ment_est_para[i][0];
  }

  for(int emo_num=0;emo_num<1;emo_num++){
    for(int sit_num=0;sit_num<8;sit_num++){
      W1(sit_num,emo_num)=ment_est_para[sit_num][emo_num];
    }
  }
  *W1_para = W1;

}

void f_mm(MatrixXd W1_para, double (&factor)[N][8], double mental[N], double *output){

  Mat_phi Phi1;//,Phi2,Phi3,Phi4;//計画行列のこと
  setPhi(factor, mental,0,&Phi1);
  //setPhi(factor, mental,1,&Phi2);
  //setPhi(factor, mental,2,&Phi3);
  //setPhi(factor, mental,3,&Phi4);

  for(int col_num1=0;col_num1<N;col_num1++){
    for(int emo_num=0;emo_num<1;emo_num++){//まだ閉じてない
      for(int sit_num=0;sit_num<8;sit_num++){
        //switch(emo_num){
        //  case 0:
        output[0] += W1_para(sit_num,emo_num)*Phi1(col_num1,sit_num);
        //    break;
        //  case 1:
        //    output[1] += W1_para(sit_num,emo_num)*Phi2(col_num1,sit_num);
        //    break;
        //  case 2:
        //    output[2] += W1_para(sit_num,emo_num)*Phi3(col_num1,sit_num);
        //    break;
        //  case 3:
        //    output[3] += W1_para(sit_num,emo_num)*Phi4(col_num1,sit_num);
            //std::cout<<"col_num1: "<<col_num1<<",";
            //std::cout<<"emo_num: "<<emo_num<<",";
            //std::cout<<"sit_num: "<<sit_num<<",";
            //std::cout<<"output check: "<<output[emo_num]<<std::endl;
            //std::cout<<"para check: "<<W1_para(sit_num,emo_num)<<std::endl;
        //    break;
        }
      }
    }
}
//mental model

double fE(MatrixXd W1_para, double (&factor)[N][8], double mental[N], double sig_total[1]){
  double output=0;
  double total[1]={0};

  f_mm(W1_para, factor, mental, total);
  for(int emo_num=0;emo_num<1;emo_num++){
    output+= 0.1*total[emo_num]-sig_total[emo_num];
  }
  return output;
}
