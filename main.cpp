#include <iostream>
#include <complex>
#include <math.h>
#include <fstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <ctime>

using namespace std;

double Q_erf_D ( double parm ) {///////误差函数的导数

    double const PI = acos(double(-1));
    double result_Q_erf_D = 2 * ( 1 / (sqrt(PI)) ) * exp(-1*(parm*parm)) ;
    return result_Q_erf_D ;

}

double fun_vert_x ( double vec_a[] , double vec_b[] ) {/////点乘

    double lin_c = 0 ;
    for ( int i = 0 ; i < 3 ; i ++ ) {
        //
        lin_c = lin_c + ( vec_a[i] * vec_b[i] ) ;
    }
    return lin_c;
}

void fun_vert_j ( double vec_c[] , double vec_a[], double vec_b[] ) {////向量减法

    for ( int i = 0 ; i < 3 ; i ++ ) {
        //
        vec_c[i] = vec_a[i] - vec_b[i] ;
    }
}

void fun_matrix_ab ( double matrix_n[][3], double vec_a[] , double vec_b[] ) {/////向量转置乘=矩阵

    for ( int i = 0 ; i < 3 ; i ++ ) {
        for ( int j = 0; j < 3 ; j ++ ) {
            matrix_n[i][j] = vec_a[i] * vec_b[j] ;
        }
    }
}

void fun_matrix_mAB ( double matrix_n[][3] , double matrix_A[][3] , double matrix_B[][3] )  {//////矩阵相乘

    for ( int i = 0 ; i < 3 ; i ++ ) {
        for ( int j = 0; j < 3 ; j ++ ) {
            //
            for ( int k = 0 ; k < 3 ; k++ ) {
                matrix_n[i][j] = matrix_n[i][j] + ( matrix_A[i][k] * matrix_B[k][j] ) ;
            }
        }
    }
}

void fun_matrix_vector ( double vec_c[], double matrix_n[][3], double vec_a[] ) {///////矩阵乘向量

    for ( int i = 0 ; i < 3 ; i ++ ) {
        for ( int j = 0 ; j < 3 ; j ++ ) {
            //
            vec_c[i] = vec_c[i] + (matrix_n[i][j] * vec_a[j] ) ;
        }
    }
}

void fun_matrix_mAJIAB ( double matrix_n[][3] , double matrix_A[][3] , double matrix_B[][3] )  {//////矩阵相加

    for ( int i = 0 ; i < 3 ; i ++ ) {
        for ( int j = 0; j < 3 ; j ++ ) {
            //
            matrix_n[i][j] =  ( matrix_A[i][j] + matrix_B[i][j] ) ;
        }
    }
}

//////////////FElnFX
double fun_FEln ( double a1 , double a2 , double a0 , double N_Success , double x_Success[] , double y_Success[] , double N_Failure, double x_Failure[] , double y_Failure[] ) {

    double FEln = 0 ;
    for ( int i = 0 ; i < N_Success ; i ++ ) {
        FEln = FEln - log((1+erf((a1*x_Success[i])+(a2*y_Success[i])+(a0)))/2) ;
//        cout<<i<<" "<<x_Success[i]<<" "<<y_Success[i]<<" "<<(a1*x_Success[i])+(a2*y_Success[i])+(a0)<<" "<<(1+erf((a1*x_Success[i])+(a2*y_Success[i])+(a0)))/2<<" "<<FEln<<endl;
    }
    for ( int i = 0 ; i < N_Failure ; i ++ ) {
        FEln = FEln - log( 1 - ( ((1+erf((a1*x_Failure[i])+(a2*y_Failure[i])+(a0)))/2) ) );
    }
    return FEln ;
}
/////////////////梯度FElnFX
void fun_DFElnD ( double d_vect[] , double a1 , double a2 , double a0 , double N_Success , double x_Success[] , double y_Success[] , double N_Failure, double x_Failure[] , double y_Failure[] ) {

    //cout<<d_vect[0]<<" "<<d_vect[1]<<" "<<d_vect[2]<<" "<<a1<<a2<<a0<<" "<<x_Success[0]<<" "<<y_Success[0]<<" "<<x_Failure[0]<<" "<<y_Failure[0]<<" "<<N_Success<<" "<<N_Failure<<endl;
    double d_1 = 0;
    double d_2 = 0;
    double d_0 = 0;
    for ( int i = 0 ; i < N_Success ; i ++ ) {
        d_1 = d_1 + ((0.5*(Q_erf_D((a1*x_Success[i])+(a2*y_Success[i])+(a0)))*x_Success[i])/( (1+erf((a1*x_Success[i])+(a2*y_Success[i])+(a0)))/2.0 ));
        //cout<<i<<" "<<(a1*x_Success[i])+(a2*y_Success[i])+(a0)<<" "<<(0.5*Q_erf_D((a1*x_Success[i])+(a2*y_Success[i])+(a0))*x_Success[i])<<endl;
        d_2 = d_2 + ((0.5*(Q_erf_D((a1*x_Success[i])+(a2*y_Success[i])+(a0)))*y_Success[i])/( (1+erf((a1*x_Success[i])+(a2*y_Success[i])+(a0)))/2 ));
        d_0 = d_0 + ((0.5*(Q_erf_D((a1*x_Success[i])+(a2*y_Success[i])+(a0))))/( (1+erf((a1*x_Success[i])+(a2*y_Success[i])+(a0)))/2 ));
    //    cout<<d_1<<" "<<d_2<<" "<<d_0<<" "<<x_Success[i]<<" "<<y_Success[i]<<" "<<((0.5*(Q_erf_D((a1*x_Success[i])+(a2*y_Success[i])+(a0)))*x_Success[i])/( (1+erf((a1*x_Success[i])+(a2*y_Success[i])+(a0)))/2.0 ))<<endl;
    }
    for ( int i = 0 ; i < N_Failure ; i ++ ) {
        d_1 = d_1 + ((-0.5*(Q_erf_D((a1*x_Failure[i])+(a2*y_Failure[i])+(a0)))*x_Failure[i])/( ( 1 - (1+erf((a1*x_Failure[i])+(a2*y_Failure[i])+(a0)))/2 ) ));
        d_2 = d_2 + ((-0.5*(Q_erf_D((a1*x_Failure[i])+(a2*y_Failure[i])+(a0)))*y_Failure[i])/( ( 1 - (1+erf((a1*x_Failure[i])+(a2*y_Failure[i])+(a0)))/2 ) ));
        d_0 = d_0 + ((-0.5*(Q_erf_D((a1*x_Failure[i])+(a2*y_Failure[i])+(a0))))/( ( 1 - (1+erf((a1*x_Failure[i])+(a2*y_Failure[i])+(a0)))/2 ) ));
   //     cout<<d_1<<" "<<d_2<<" "<<d_0<<" "<<x_Failure[i]<<" "<<y_Failure[i]<<" "<<((-0.5*(Q_erf_D((a1*x_Failure[i])+(a2*y_Failure[i])+(a0)))*x_Failure[i])/( ( 1 - (1+erf((a1*x_Failure[i])+(a2*y_Failure[i])+(a0)))/2 ) ))<<endl;
    }
    d_vect[0] = d_1;
    d_vect[1] = d_2;
    d_vect[2] = d_0;
}
double fun_DEFln_landa ( double d1 , double d2, double d0, double a1 , double a2 , double a0 , double N_Success , double x_Success[] , double y_Success[] , double N_Failure, double x_Failure[] , double y_Failure[] ) {
 //   cout<<d1<<" "<<d2<<" "<<d0<<" "<<a1<<" "<<a2<<" "<<a0<<endl;
    double d_k_a1 = 0;
    double d_k_a2 = 0;
    double d_k_a0 = 0;
    for ( int i = 0 ; i < N_Success ; i ++ ) {
        d_k_a1 = d_k_a1 + ((0.5*(Q_erf_D((a1*x_Success[i])+(a2*y_Success[i])+(a0)))*x_Success[i]*d1)/( (1+erf((a1*x_Success[i])+(a2*y_Success[i])+(a0)))/2.0 ));
    //    cout<<d_k_a1<<" "<<((0.5*(Q_erf_D((a1*x_Success[i])+(a2*y_Success[i])+(a0)))*x_Success[i]*d1)/( (1+erf((a1*x_Success[i])+(a2*y_Success[i])+(a0)))/2.0 ))<<endl;
        d_k_a2 = d_k_a2 + ((0.5*(Q_erf_D((a1*x_Success[i])+(a2*y_Success[i])+(a0)))*y_Success[i]*d2)/( (1+erf((a1*x_Success[i])+(a2*y_Success[i])+(a0)))/2 ));
        d_k_a0 = d_k_a0 + ((0.5*(Q_erf_D((a1*x_Success[i])+(a2*y_Success[i])+(a0)))*d0)/( (1+erf((a1*x_Success[i])+(a2*y_Success[i])+(a0)))/2 ));
    }
    for ( int i = 0 ; i < N_Failure ; i ++ ) {
        d_k_a1 = d_k_a1 + ((-0.5*(Q_erf_D((a1*x_Failure[i])+(a2*y_Failure[i])+(a0)))*x_Failure[i]*d1)/( ( 1 - (1+erf((a1*x_Failure[i])+(a2*y_Failure[i])+(a0)))/2 ) ));
    //    cout<<d_k_a1<<" "<<((-0.5*(Q_erf_D((a1*x_Failure[i])+(a2*y_Failure[i])+(a0)))*x_Failure[i]*d1)/( ( 1 - (1+erf((a1*x_Failure[i])+(a2*y_Failure[i])+(a0)))/2 ) ))<<" "<<( ( 1 - (1+erf((a1*x_Failure[i])+(a2*y_Failure[i])+(a0)))/2 ) )<<endl;
        d_k_a2 = d_k_a2 + ((-0.5*(Q_erf_D((a1*x_Failure[i])+(a2*y_Failure[i])+(a0)))*y_Failure[i]*d2)/( ( 1 - (1+erf((a1*x_Failure[i])+(a2*y_Failure[i])+(a0)))/2 ) ));
        d_k_a0 = d_k_a0 + ((-0.5*(Q_erf_D((a1*x_Failure[i])+(a2*y_Failure[i])+(a0)))*d0)/( ( 1 - (1+erf((a1*x_Failure[i])+(a2*y_Failure[i])+(a0)))/2 ) ));
    }
 //   cout<<d_k_a1 + d_k_a2 + d_k_a0<<" "<<endl;
    return (d_k_a1 + d_k_a2 + d_k_a0) ;

}

int main()
{
    ofstream OUT_print;
    OUT_print.open("out.txt",ios::trunc);
    cout << "Hello world!" << endl;
    double aZ = 0.0 ;
    double resul_erf ;
    resul_erf = erf (aZ) ;///////////////////////////////误差函数
    cout<<"erf("<<aZ<<")="<<resul_erf<<endl;
    //erf
    double res2 ;
    res2 = Q_erf_D(2) ;//////////////////////////////////////误差函数导数
    cout<<res2<<endl;
    //erf'
    //matrix
    double vect_a[3] = {1,2,3} ;
    double vect_b[3] = {2,3,4} ;
    double p;
    p = fun_vert_x(vect_a,vect_b);
    cout<<"a*b= "<<p<<endl;
    //////////////////////////////////////////////////向量点乘
    double vect_c[3] ;
    fun_vert_j(vect_c,vect_a,vect_b);
    for (int i = 0 ; i < 3 ; i ++) {
        cout<<vect_c[i]<<" ";
    }
    cout<<endl;
    //////////////////////////////////////////////////向量减法
    double matrix_double[3][3];
    fun_matrix_ab(matrix_double,vect_a,vect_b);
    for ( int i = 0 ; i < 3 ; i ++ ) {
        for ( int j = 0; j < 3 ; j ++ ) {
            cout<<*(*(matrix_double+i)+j)<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
    /////////////////////////////////////////////////向量转置乘向量=矩阵
    double matrix_A[3][3] = {{1,2,3},{2,3,4},{3,4,5}};
    double matrix_B[3][3] = {{1,2,3},{2,3,4},{3,4,5}};
    double matrix_n[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    fun_matrix_mAB( matrix_n,matrix_A,matrix_B );
    for ( int i = 0 ; i < 3 ; i ++ ) {
        for ( int j = 0; j < 3 ; j ++ ) {
            cout<<*(*(matrix_n+i)+j)<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
    ///////////////////////////////////////////////////矩阵相乘
    double v_c[3] = {0,0,0};
    double v_a[3] = {1,1,1};
    fun_matrix_vector( v_c,matrix_n,v_a );
    for ( int i = 0 ; i < 3 ; i++ ) {
        cout<<v_c[i]<<" ";
    }
    cout<<endl;
    //////////////////////////////////////////////////矩阵乘向量
    double cs_ln = log(2.7182818284);
    cout<<cs_ln<<endl;
    //main/////////////////////////////////////////////////////////////////////main/////////////////
    cout.precision(16);
    ifstream OpenFile("A_data.txt", ios::in);
    if (OpenFile.fail()){
        cout<<"Can not open target file"<<endl;
        return 0;
    }
    std::string lineStr;
    std::vector<double> x_A;
    std::vector<double> y_A;
    int k1 = 0 ;
    if (OpenFile){
        //

        while (getline(OpenFile,lineStr)){
            //
            k1 ++ ;
            string x_k1,y_k1;
            x_k1=lineStr.substr(0,11);
            y_k1=lineStr.substr(13,22);
            x_A.push_back(atof(x_k1.c_str()));
            y_A.push_back(atof(y_k1.c_str()));
            cout<<x_A[k1-1]<<" "<<y_A[k1-1]<<endl;
        }

    }
    cout<<k1<<endl;
    OpenFile.close();//////////////////////////////////////AAAAAAAAAAA
    double x_AA[k1];
    double y_AA[k1];
    for ( int i = 0 ; i < k1 ; i ++ ) {
        x_AA[i] = x_A[i];
        y_AA[i] = y_A[i];
    }
    ifstream OpenFileB("B_data.txt", ios::in);
    if (OpenFileB.fail()){
        cout<<"Can not open target file"<<endl;
        return 0;
    }
    std::string lineStrB;
    std::vector<double> x_B;
    std::vector<double> y_B;
    int k2 = 0 ;
    if (OpenFileB){
        //

        while (getline(OpenFileB,lineStrB)){
            //
            k2 ++ ;
            string x_k2,y_k2;
            x_k2=lineStrB.substr(0,11);
            y_k2=lineStrB.substr(13,22);
            x_B.push_back(atof(x_k2.c_str()));
            y_B.push_back(atof(y_k2.c_str()));
            cout<<x_B[k2-1]<<" "<<y_B[k2-1]<<endl;
        }

    }
    cout<<k2<<endl;
    OpenFileB.close();//////////////////////////////////////BBBBBBBBBBB
    double x_BB[k1];
    double y_BB[k1];
    for ( int i = 0 ; i < k2 ; i ++ ) {
        x_BB[i] = x_B[i];
        y_BB[i] = y_B[i];
    }
    ////构建-ELN 以及导函数
    //////////////////读取数据/////////////////////////////////////////////////////////////////////
    double a1 = -0.00072;
    double a2 = 0.2;
    double a0 = 0.0015;
    double H_matrix[3][3] = {{1,0,0},{0,1,0},{0,0,1}} ;
    double epslong = 0.01;//收敛极限
    double d_IN[3] = {0,0,0};
    cout<<d_IN[0]<<d_IN[1]<<d_IN[2]<<" "<<a1<<a2<<a0<<" "<<k1<<endl;
    cout<<x_AA[96]<<" "<<y_AA[96]<<" "<<k2<<" "<<x_BB[0]<<" "<<y_BB[0]<<endl;
    fun_DFElnD( d_IN, a1 , a2 , a0 , k2 , x_BB , y_BB , k1 , x_AA , y_AA  );
    ///////////////////////测试fun_-FX
    double cs_fun_FX = 0 ;
    cs_fun_FX = fun_FEln(a1,a2,a0,k2,x_BB,y_BB,k1,x_AA,y_AA);
    cout<<"-fx= "<<cs_fun_FX<<endl;
    //////////////////////////////////epsion求解
    cout<<d_IN[0]<<" "<<d_IN[1]<<" "<<d_IN[2]<<endl;
    double fan2_epslong = (d_IN[0]*d_IN[0]) + (d_IN[1]*d_IN[1]) + (d_IN[2]*d_IN[2]) ;
    cout<<fan2_epslong<<endl;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////+while
    while ( fan2_epslong > epslong ) {
    ////////////////////////////////步长迭代
        double epslong2 = 0.01;
        double landa = 0.00001 ;//步长
        double next_a1,next_a2,next_a0;
        next_a1 = a1 + (landa*d_IN[0]);
        next_a2 = a2 + (landa*d_IN[1]);
        next_a0 = a0 + (landa*d_IN[2]);
        cout<<next_a1<<" "<<next_a2<<" "<<next_a0<<endl;
        double PDPDPD = 0;
        PDPDPD = fun_DEFln_landa(d_IN[0],d_IN[1],d_IN[2],next_a1,next_a2,next_a0,k2,x_BB,y_BB,k1,x_AA,y_AA);
        cs_fun_FX = fun_FEln(next_a1,next_a2,next_a0,k2,x_BB,y_BB,k1,x_AA,y_AA);
        cout<<PDPDPD<<"fx="<<cs_fun_FX<<endl;
        while ( (PDPDPD*PDPDPD) > epslong2 ) {
            landa = landa + PDPDPD*0.00000001 ;
            next_a1 = a1 + (landa*d_IN[0]);
            next_a2 = a2 + (landa*d_IN[1]);
            next_a0 = a0 + (landa*d_IN[2]);
            PDPDPD = fun_DEFln_landa(d_IN[0],d_IN[1],d_IN[2],next_a1,next_a2,next_a0,k2,x_BB,y_BB,k1,x_AA,y_AA);
            cs_fun_FX = fun_FEln(next_a1,next_a2,next_a0,k2,x_BB,y_BB,k1,x_AA,y_AA);
      //      cout<<PDPDPD<<" fx= "<<cs_fun_FX<<endl;
        }
    //    cout<<"步长："<<landa<<endl;
        next_a1 = a1 + (landa*d_IN[0]);
        next_a2 = a2 + (landa*d_IN[1]);
        next_a0 = a0 + (landa*d_IN[2]);
        double deta[3] = {0,0,0};
        deta[0] = next_a1 - a1;
        deta[1] = next_a2 - a2;
        deta[2] = next_a0 - a0;
        double d_IN_next[3] = {0,0,0};
        double gama[3] = {0,0,0};
        fun_DFElnD( d_IN_next, next_a1 , next_a2 , next_a0 , k2 , x_BB , y_BB , k1 , x_AA , y_AA  );
        gama[0] = ( -1 * d_IN_next[0] ) + d_IN[0];
        gama[1] = ( -1 * d_IN_next[1] ) + d_IN[1];
        gama[2] = ( -1 * d_IN_next[2] ) + d_IN[2];
   //     cout<<"deta: "<<deta[0]<<" "<<deta[1]<<" "<<deta[2]<<" gama: "<<gama[0]<<" "<<gama[1]<<" "<<gama[2]<<endl;
        double H_matrix2[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        double csH_matrix1[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        double csH_matrix2[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        fun_matrix_ab(csH_matrix1,deta,deta);
        double csH_M2_xishu = 0;
        csH_M2_xishu = fun_vert_x(deta,gama);
   //     cout<<csH_M2_xishu<<endl;
        for ( int i = 0 ; i < 3 ; i ++ ) {
            for ( int j = 0 ; j < 3 ; j ++ ) {
                csH_matrix1[i][j] = csH_matrix1[i][j] / csH_M2_xishu ;
            }
        }
  //      cout<<"csH_matrix1:"<<endl;
  //      cout<<csH_matrix1[0][0]<<" "<<csH_matrix1[0][1]<<" "<<csH_matrix1[0][2]<<endl;
  //      cout<<csH_matrix1[1][0]<<" "<<csH_matrix1[1][1]<<" "<<csH_matrix1[1][2]<<endl;
  //      cout<<csH_matrix1[2][0]<<" "<<csH_matrix1[2][1]<<" "<<csH_matrix1[2][2]<<endl;
        double csH_matrix2_A[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        double csH_matrix2_B[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        double csH_matrix2_C[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        fun_matrix_ab(csH_matrix2_A,gama,gama);
  //      cout<<csH_matrix2_A[0][0]<<" "<<csH_matrix2_A[0][1]<<" "<<csH_matrix2_A[0][2]<<endl;
  //      cout<<csH_matrix2_A[1][0]<<" "<<csH_matrix2_A[1][1]<<" "<<csH_matrix2_A[1][2]<<endl;
  //      cout<<csH_matrix2_A[2][0]<<" "<<csH_matrix2_A[2][1]<<" "<<csH_matrix2_A[2][2]<<endl;
        fun_matrix_mAB(csH_matrix2_B,H_matrix,csH_matrix2_A);
   //     cout<<csH_matrix2_B[0][0]<<" "<<csH_matrix2_B[0][1]<<" "<<csH_matrix2_B[0][2]<<endl;
   //     cout<<csH_matrix2_B[1][0]<<" "<<csH_matrix2_B[1][1]<<" "<<csH_matrix2_B[1][2]<<endl;
   //     cout<<csH_matrix2_B[2][0]<<" "<<csH_matrix2_B[2][1]<<" "<<csH_matrix2_B[2][2]<<endl;
        fun_matrix_mAB(csH_matrix2_C,csH_matrix2_B,H_matrix);
   //     cout<<csH_matrix2_C[0][0]<<" "<<csH_matrix2_C[0][1]<<" "<<csH_matrix2_C[0][2]<<endl;
   //     cout<<csH_matrix2_C[1][0]<<" "<<csH_matrix2_C[1][1]<<" "<<csH_matrix2_C[1][2]<<endl;
   //     cout<<csH_matrix2_C[2][0]<<" "<<csH_matrix2_C[2][1]<<" "<<csH_matrix2_C[2][2]<<endl;
        double csV_xishu[3] = {0,0,0};
        fun_matrix_vector(csV_xishu,H_matrix,gama);
   //     cout<<csV_xishu[0]<<" "<<csV_xishu[1]<<" "<<csV_xishu[2]<<endl;
        double csH_M2_xishu2 = 0;
        csH_M2_xishu2 = fun_vert_x(gama,csV_xishu);
   //     cout<<csH_M2_xishu2<<endl;
        for ( int i = 0 ; i < 3 ; i ++ ) {
            for ( int j = 0 ; j < 3 ; j ++ ) {
                csH_matrix2[i][j] = ( -1 * csH_matrix2_C[i][j] ) / csH_M2_xishu2 ;
            }
        }
   //     cout<<"csH_matrix2: "<<endl;
   //     cout<<csH_matrix2[0][0]<<" "<<csH_matrix2[0][1]<<" "<<csH_matrix2[0][2]<<endl;
   //     cout<<csH_matrix2[1][0]<<" "<<csH_matrix2[1][1]<<" "<<csH_matrix2[1][2]<<endl;
   //     cout<<csH_matrix2[2][0]<<" "<<csH_matrix2[2][1]<<" "<<csH_matrix2[2][2]<<endl;
        fun_matrix_mAJIAB(H_matrix2,csH_matrix1,csH_matrix2);
   //     cout<<H_matrix2[0][0]<<" "<<H_matrix2[0][1]<<" "<<H_matrix2[0][2]<<endl;
   //     cout<<H_matrix2[1][0]<<" "<<H_matrix2[1][1]<<" "<<H_matrix2[1][2]<<endl;
   //     cout<<H_matrix2[2][0]<<" "<<H_matrix2[2][1]<<" "<<H_matrix2[2][2]<<endl;
        double H_matrix_Next[3][3] = {{0,0,0},{0,0,0},{0,0,0}} ;
        fun_matrix_mAJIAB(H_matrix_Next,H_matrix,csH_matrix2);
    //    cout<<"H_matrix_Next: "<<endl;
    //    cout<<H_matrix_Next[0][0]<<" "<<H_matrix_Next[0][1]<<" "<<H_matrix_Next[0][2]<<endl;
    //    cout<<H_matrix_Next[1][0]<<" "<<H_matrix_Next[1][1]<<" "<<H_matrix_Next[1][2]<<endl;
    //    cout<<H_matrix_Next[2][0]<<" "<<H_matrix_Next[2][1]<<" "<<H_matrix_Next[2][2]<<endl;
        double D_HD_lin[3] = {0,0,0};
        fun_matrix_vector ( D_HD_lin, H_matrix_Next , d_IN_next );
    //    cout<<D_HD_lin[0]<<" "<<D_HD_lin[1]<<" "<<D_HD_lin[2]<<endl;
        d_IN[0] = D_HD_lin[0];
        d_IN[1] = D_HD_lin[1];
        d_IN[2] = D_HD_lin[2];
        for ( int i = 0 ; i < 3 ; i ++ ) {
            for ( int j = 0 ; j < 3 ; j ++ ) {
                H_matrix[i][j] = H_matrix_Next[i][j];
            }
        }
        a1 = next_a1;
        a2 = next_a2;
        a0 = next_a0;
        fan2_epslong = (d_IN[0]*d_IN[0]) + (d_IN[1]*d_IN[1]) + (d_IN[2]*d_IN[2]);
        OUT_print<<a1<<" "<<a2<<" "<<a0<<" fan2_epslong :"<<fan2_epslong<<endl;
        OUT_print<<"FX= "<<fun_FEln(a1,a2,a0,k2 , x_BB , y_BB , k1 , x_AA , y_AA)<<endl;
    }///////////////////////////////////////////////////////////////////////////////////////////////////////////////////while
    //
    OUT_print.close();
    return 0;
}
