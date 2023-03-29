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

#define N 30
#define SET_NUM_LEN 145

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
void get_para(double (&factor)[N][8], double (&signal)[N][1], double mental[N], MatrixXd *W1_para);
void f_mm(MatrixXd W1_para, double (&factor)[N][8], double mental[N], double *output, char *filename);
//mental model
double fE(MatrixXd W1_para, double (&factor)[N][8], double mental[N], double sig_total[1], char *filename);

void out_tmppara(char *filename, MatrixXd W1_para);

//void grad_descent(double M_org[N], double* M,MatrixXd W1_para,double (&factor)[N][8], double sig_total[1],int id_num,int moto_count,char* filename, char* dir_path, int set_num);
void grad_descent(double M_org[N], double* M,MatrixXd W1_para,double (&factor)[N][8], double sig_total[1],int id_num,int moto_count,char* filename, char* dir_path);

int main(int argc, const char* argv[]){//(argv[0]=./estimate),argv[1]=usrname

  double *_factor[N],*_signal[N];
  double factor_size[N][8],signal_size[N][1];//データ数は１５個、要素数は７こ、
  double factor[N][8], signal[N][1];
  for(int i=0;i<N;i++) _factor[i] = factor_size[i];
  for(int i=0;i<N;i++) _signal[i] = signal_size[i];

  char factor_name[100];
  char signal_name[100];
  char face_name[200];
  char err_log[200];
  char dir_path[50];
  MatrixXd W1_para(8,1);

  int id_num=0;
  if(argc > 1){
    std::string id_input_str = argv[1];
    id_num=stoi(id_input_str);
    printf("id is %d\n",id_num);
  }

	sprintf(dir_path,"../data/%d/",id_num);

	sprintf(factor_name,"%sfactor_train.csv",dir_path);
	//sprintf(factor_name,"%sfactor_train%d-%d.csv",dir_path,N,set_num);
	std::cout<<factor_name<<std::endl;

	sprintf(signal_name,"%sface_train.csv",dir_path);
	//sprintf(signal_name,"%sface_train%d-%d.csv",dir_path,N,set_num);
	sprintf(err_log,"%soutput/err_phi.csv",dir_path);
	//sprintf(err_log,"%soutput/err_phi%d-%d.csv",dir_path,N,set_num);

	std::ofstream err_ofs(err_log);

	read_file(factor_name, 8,_factor);
	read_file(signal_name, 4,_signal);

	memcpy(factor,*_factor,sizeof(factor));
	memcpy(signal,*_signal,sizeof(signal));

	double m_org[N];
	double mental[N];
	srand((unsigned) time(NULL));
	for(int ment_num=0;ment_num<N;ment_num++){
	m_org[ment_num]={0};
	}
	memcpy(mental,m_org,sizeof(m_org));

	double dEdM[N] = {0};

	int errflg=0;
	double err=10;

	std::cout << std::fixed;
	int count=0;

	char tmp_w_filename[200];
	char mental_filename[200];
	sprintf(tmp_w_filename,"%soutput/weight_zero_trial00.csv",dir_path);
	//sprintf(tmp_w_filename,"%soutput/weight_zero_trial00%d-%d.csv",dir_path,N,set_num);
	sprintf(mental_filename,"%soutput/mental_log_trial00.csv",dir_path);
	//sprintf(mental_filename,"%soutput/mental_log_trial00%d-%d.csv",dir_path,N,set_num);
	ofstream init_out(tmp_w_filename,ios::app);
	ofstream mental_out(mental_filename,ios::app);

	init_out<<id_num<<","<<count<<"==================="<<std::endl;

	while(errflg==0){
		std::cout<<"err"<<std::setprecision(10)<<err<<" count:"<<count<<", M:";
		mental_out<<"err"<<std::setprecision(10)<<err<<" count:"<<count<<", M:";

		for(int m_check=0;m_check<N;m_check++){
			std::cout<<std::setprecision(2)<<mental[m_check]<<", ";
			mental_out<<std::setprecision(2)<<mental[m_check]<<", ";
		}
		std::cout<<std::endl;
		mental_out<<std::endl;

		//MatrixXd W1_para(10,4);
		get_para(factor, signal, mental, &W1_para);

		double sig_total[1]={0};

		for(int emo_num=0;emo_num<1;emo_num++){
			for(int col_num1 = 0; col_num1<N; col_num1++){
				sig_total[emo_num] += signal[col_num1][emo_num];
			}
		}

		int count_mgrad=0;
		//TODO
		char some_filename[200];
		sprintf(some_filename,"%soutput/signal_out.csv",dir_path);
		double E_pre =fE(W1_para,factor,m_org,sig_total,some_filename);

		grad_descent(m_org,mental,W1_para,factor,sig_total,id_num,count,some_filename,dir_path);
		//grad_descent(m_org,mental,W1_para,factor,sig_total,id_num,count,some_filename,dir_path,set_num);

		//S誤を求める
		//TODO
		err = abs(fE(W1_para, factor, mental,sig_total,some_filename) - E_pre);

		init_out<<id_num<<","<<count<<"==================="<<std::endl;
		out_tmppara(tmp_w_filename, W1_para);
		err_ofs<<err<<std::endl;
		std::cout<<id_num<<","<<count<<std::endl;

		count++;
		memcpy(m_org,mental,sizeof(mental));
		if(err<=0.0001){
			errflg=1;
		}else if(count>=50000){//KENSYU: loop count
			errflg=1;
		}
			printf("checkpoint err%f\n",err);
		}

	char hap[200],sup[100],ang[100],sad[100],face[100];
	sprintf(hap,"%soutput/hap_weight.csv",dir_path);
	//sprintf(hap,"%soutput/hap_weight%d-%d.csv",dir_path,N,set_num);

	ofstream log1(hap,ios::app);
	for(int factor_num=0;factor_num<8;factor_num++){
		log1<<W1_para(factor_num,0)<<endl;
	}
	out_tmppara(tmp_w_filename, W1_para);

	char sigtar_filename[200];
	sprintf(sigtar_filename,"%soutput/filename.csv",dir_path);
	//sprintf(sigtar_filename,"%soutput/filename%d-%d.csv",dir_path,N,set_num);
	double total[1]={0};
	f_mm(W1_para, factor, mental, total,sigtar_filename);
}


void out_tmppara(char *filename, MatrixXd W1_para){
  ofstream out_para(filename,ios::app);
  for(int e_num=0;e_num<1;e_num++){
    for(int f_num=0;f_num<8;f_num++){
      out_para<<W1_para(f_num,e_num)<<",";
    }
    out_para<<std::endl;
  }
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
        getline(ss,data_f[i],'\t');
        hoga=std::stod(data_f[i]);
        hage=std::stod(data_f[i]);
        data[j][i]=hage;
      }
      j++;
    }

}

//void grad_descent(double M_org[N], double* M,MatrixXd W1_para,double (&factor)[N][8], double sig_total[1],int id_num,int moto_count, char* filename, char *dir_path,int set_num){
void grad_descent(double M_org[N], double* M,MatrixXd W1_para,double (&factor)[N][8], double sig_total[1],int id_num,int moto_count, char* filename, char *dir_path){
  double dEdM[N]={0};
  double err=100;
  double Msize[N];
  int count=0;
  char mental_log_second[200];
  sprintf(mental_log_second,"%soutput/mental_log_trial2.csv",dir_path);
  ofstream mental_out(mental_log_second,ios::app);

  while(err>0.001){
    //TODO	  
    double E_pre  = fE(W1_para, factor, M_org,sig_total,filename);
    //std::cout<<"E_pre: "<<E_pre<<", ";
    for(int i=0;i<N;i++){
      double tmp_M = M_org[i];
      double E_cand[3] = {0};
      double delta[3] = {0};
      for(int dm = -1; dm<2;dm++){
        M[i] = tmp_M + pow(10.0,dm);
        E_cand[dm+1] = fE(W1_para, factor, M,sig_total,filename);
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

    mental_out<<"gd_descent("<<id_num<<"-"<<"-"<<moto_count<<"),"<<std::setprecision(10)<<err<<","<<count<<",";
    for(int m_check=0;m_check<N;m_check++){
      //cout<<std::setprecision(2)<<M[m_check]<<", ";
      mental_out<<std::setprecision(2)<<M[m_check]<<", ";
    }

    //cout<<std::endl;
    mental_out<<std::endl;
    count++;
    if(count>10000){
      break;
    }

  }

  //cout<<std::endl;
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
  //int func_num[8][1] = {{0},{4},{8},{12},{16},{20},{24},{28},{32},{36}};
  
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

void get_para(double (&factor)[N][8], double (&signal)[N][1], double mental[N], MatrixXd *W1_para){
  Mat_phi Phi1;//,Phi2,Phi3,Phi4;//計画行列のこと  //行優先の行列

  setPhi(factor, mental,0,&Phi1);

  MatrixXd Sig_tar(N,1);//目標値tに相当します//最後の１列は全部１

  for(int d_num=0;d_num<N;d_num++){
    for(int e_num=0;e_num<1;e_num++){
      Sig_tar(d_num,e_num)=signal[d_num][e_num];
    }
  }

  MatrixXd Phi1_t(8,N);
  MatrixXd PP1(8,8);
  MatrixXd PP1_inv(8,8);

  PP1=Phi1.transpose()*Phi1;

  pseudo_inverse(&PP1_inv,PP1);

  W1_para->col(0)=PP1_inv*Phi1.transpose()*Sig_tar.col(0);//どの因子が表情表出に大きな影響を及ぼすかの重み

}

void f_mm(MatrixXd W1_para, double (&factor)[N][8], double mental[N], double *output, char *filename){

  Mat_phi Phi1;//,Phi2,Phi3,Phi4;//計画行列のこと
  setPhi(factor, mental,0,&Phi1);

  //make new file
  //char sigtar_name[50];
  //sprintf(sigtar_name, "signal_out_%d.csv",N);
  std::ofstream sigtar_out(filename,std::ios::out|std::ios::trunc);
  //std::ofstream sigtar_out(sigtar_name,std::ios::out|std::ios::trunc);
  for(int col_num1=0;col_num1<N;col_num1++){
    for(int emo_num=0;emo_num<1;emo_num++){//まだ閉じてない
      for(int sit_num=0;sit_num<8;sit_num++){
        output[0] += W1_para(sit_num,emo_num)*Phi1(col_num1,sit_num);
      }
      //print output[emo_num]
      sigtar_out<<output[emo_num];
      if(emo_num<3){
        sigtar_out<<",";
      }else if(emo_num==3){
        sigtar_out<<std::endl;
      }

    }
    //std::cout<<std::endl;
  }

}
//mental model

double fE(MatrixXd W1_para, double (&factor)[N][8], double mental[N], double sig_total[1],char *filename){
  double output=0;
  double total[1]={0};

  f_mm(W1_para, factor, mental, total, filename);
  for(int emo_num=0;emo_num<1;emo_num++){
    output+= 0.1*total[emo_num]-sig_total[emo_num];
  }
  return output;
}
