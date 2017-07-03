//inverse discrete fourier transform iDFT 1D algorithm source code implementation
//odwrotna dyskretna transformacja fouriera iDFT 1D c++ kod źródłowy

//haven't try it with other function that cos(x)+jsin(x)=sin(x+pi/2)+jsin(x)

//author marcin matysek (r)ewertyn.PL

 #include <iostream>
#include "conio.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <complex>

using namespace std;


//complex number method:
 void fun_fourier_transform_DFT_method5_full_complex(int N,double *table2);
//fun_fourier_transform_DFT_method6_full_complex(int N,std::complex<double> tab[])

void fun_inverse_fourier_transform_DFT_method5_full_complex_inverse(int N,double *table2);

//these others work only  if audio samples are not complex numbers but result of DFT  it shows complex numbers Hzthat but numbers are combinet to normal number
void fun_fourier_transform_DFT_method1(int N,double *table2);
void fun_inverse_fourier_transform_DFT_method1(int N,double *table2);

void fun_fourier_transform_DFT_method2(int N,double *table2);
void fun_inverse_fourier_transform_DFT_method2(int N,double *table2);

void fun_fourier_transform_DFT_method3(int N,double *table2);
void fun_inverse_fourier_transform_DFT_method3(int N,double *table2);

void fun_fourier_transform_DFT_method4(int N,double *table2);
void fun_inverse_fourier_transform_DFT_method4(int N,double *table2);

//inverse_method4 works only witch discrete _method4
//inverse_method2 works only witch discrete _method2
//inverse_method3 works only witch discrete _method3

static double diffclock(clock_t clock2,clock_t clock1)
{
    double diffticks=clock1-clock2;
    double diffms=(diffticks)/(CLOCKS_PER_SEC/1000);
    return diffms;
}
int main()
{
 int N;
    //if N==period of signal in table tab[] then resolution = 1 Hz

     //tab[A] A=2*N;
//tab[0 to N-1]  -re numbers
//tab[N-1 to 2N-1]  -im numbers

    // period= 16  samples=16
N=16;
    double tab[32]={-0.923879533,0.382683432,1.03153013,0.923879533,0.923879533,1.465075634,1.796896994,0.923879533,-0.923879533,
    -2.230442498,-1.796896994,-0.158512669,0.923879533,0.382683432,-1.03153013,-1.689246397};

   // period=12  samples=12
   //N=12;
   // double tab[32]={-0.923879533,0.70664666,0.996551596,0.923879533,1.659378744,1.369473808,-0.923879533,
   // -2.29335334,-0.735499212,0.923879533,-0.072672064,-1.630526192};

//tab[20]=3.86;//im number


double time2;

    for(int i=0;i<2;i++)
    {
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab[j+i*N]*1000)/1000<<"  ";
    }
    cout<<endl;
    }
    cout<<endl;

    clock_t start = clock();
    //fun_fourier_transform_DFT_method1(N,tab);
    fun_fourier_transform_DFT_method5_full_complex(N,tab);
    time2=diffclock( start, clock() );

    for(int i=0;i<2;i++)
    {
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab[j+i*N]*1000)/1000<<"  ";
    }
    cout<<endl;
    }
    cout<<endl;

    //fun_inverse_fourier_transform_DFT_method1(N,tab);
    fun_inverse_fourier_transform_DFT_method5_full_complex_inverse(N,tab);

    for(int i=0;i<2;i++)
    {
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab[j+i*N]*1000)/1000<<"  ";
    }
    cout<<endl;
    }
    cout<<endl;

    system("pause");
    return 0;
}


void fun_fourier_transform_DFT_method5_full_complex(int N,double *table2)
{
    int i=N,j=N;
    double (*(table1)[2]);
    table1[0]=&table2[0];
    table1[1]=&table2[0+N];
    const double pi=3.141592653589793238462;

    double** table3 = new double*[2];
    for( int i = 0; i < 2; ++i)
    {
      table3[i] = new double[N];
    }
for (int i=0;i<2;i++)
{
    for(int j=0;j<N;j++)
    {
       table3[i][j]=0;
    }
}
for (int i=0;i<N;i++)
{
    for(int j=0;j<N;j++)
    {
//complex number method:
           table3[0][i]=table3[0][i]+table1[0][j]*cos(i*j*2*pi/(float)N);
          table3[1][i]=table3[1][i]-table1[0][j]*sin(i*j*2*pi/(float)N);
          table3[0][i]=table3[0][i]-table1[1][j]*sin(i*j*2*pi/(float)N)*-1;//im*im
          table3[1][i]=table3[1][i]+table1[1][j]*cos(i*j*2*pi/(float)N);

//other that methot is for B from inverse transform
       //   table3[0][i]=table3[0][i]+table1[0][j]*cos(i*j*2*pi/(float)N);
       // table3[0][i]=table3[0][i]-table1[1][j]*sin(i*j*2*pi/(float)N);
      //    table3[1][i]=table3[1][i]-table1[1][j]*cos(i*j*2*pi/(float)N);
      //  table3[1][i]=table3[1][i]-table1[0][j]*sin(i*j*2*pi/(float)N);

//that methot is for C from inverse transform
        //   table3[0][i]=table3[0][i]+table1[0][j]*cos(i*j*2*pi/(float)N);
       // table3[0][i]=table3[0][i]+table1[1][j]*sin(i*j*2*pi/(float)N);
      //    table3[1][i]=table3[1][i]+table1[1][j]*cos(i*j*2*pi/(float)N);
      //  table3[1][i]=table3[1][i]+table1[0][j]*sin(i*j*2*pi/(float)N);

//that methot is for D from inverse transform
      //     table3[0][i]=table3[0][i]+table1[0][j]*cos(i*j*2*pi/(float)N);
       //   table3[0][i]=table3[0][i]-table1[1][j]*sin(i*j*2*pi/(float)N);
       //   table3[1][i]=table3[1][i]+table1[1][j]*cos(i*j*2*pi/(float)N);
       //   table3[1][i]=table3[1][i]-table1[0][j]*sin(i*j*2*pi/(float)N);

    }
}
    for(int j=0;j<N;j++)
    {
      table1[0][j] =table3[0][j];
      table1[1][j] =table3[1][j];
    }

   for( int i = 0; i < 2; ++i)
   {
       delete[] table3[i];
    }
    delete[] table3;

}

void fun_inverse_fourier_transform_DFT_method5_full_complex_inverse(int N,double *table2)
{
    int i=N,j=N;
    double (*(table1)[2]);
    table1[0]=&table2[0];
    table1[1]=&table2[0+N];
    const double pi=3.141592653589793238462;

    double** table3 = new double*[2];
    for( int i = 0; i < 4; ++i)
    {
      table3[i] = new double[N];

    }
for (int i=0;i<4;i++)
{
    for(int j=0;j<N;j++)
    {
       table3[i][j]=0;
    }
}
for (int i=0;i<N;i++)
{
    for(int j=0;j<N;j++)
    {
//complex number method:
          table3[0][i]=table3[0][i]+table1[0][j]*cos(i*j*2*pi/(float)N);
          table3[1][i]=table3[1][i]+table1[0][j]*sin(i*j*2*pi/(float)N);
          table3[0][i]=table3[0][i]+table1[1][j]*sin(i*j*2*pi/(float)N)*-1;//im*im
          table3[1][i]=table3[1][i]+table1[1][j]*cos(i*j*2*pi/(float)N);

//other that methot is for B from discrete transform
          //table3[0][i]=table3[0][i]-table1[0][j]*cos(i*j*2*pi/(float)N);
       // table3[0][i]=table3[0][i]+table1[1][j]*sin(i*j*2*pi/(float)N);
       //   table3[1][i]=table3[1][i]-table1[1][j]*cos(i*j*2*pi/(float)N);
      //  table3[1][i]=table3[1][i]-table1[0][j]*sin(i*j*2*pi/(float)N);

//that methot is for C from discrete transform
      //    table3[0][i]=table3[0][i]+table1[0][j]*cos(i*j*2*pi/(float)N);
      //  table3[0][i]=table3[0][i]+table1[1][j]*sin(i*j*2*pi/(float)N);
      //    table3[1][i]=table3[1][i]+table1[1][j]*cos(i*j*2*pi/(float)N);
     //   table3[1][i]=table3[1][i]+table1[0][j]*sin(i*j*2*pi/(float)N);

 //that methot is for D from discrete transform
         // table3[0][i]=table3[0][i]-table1[0][j]*cos(i*j*2*pi/(float)N);
         // table3[0][i]=table3[0][i]+table1[1][j]*sin(i*j*2*pi/(float)N);
         //table3[1][i]=table3[1][i]-table1[1][j]*cos(i*j*2*pi/(float)N);
         // table3[1][i]=table3[1][i]+table1[0][j]*sin(i*j*2*pi/(float)N);

    }
}
    for(int j=0;j<N;j++)
    {
      table1[0][j] =table3[0][j]/N;
      table1[1][j] =table3[1][j]/N;
    }

   for( int i = 0; i < 4; ++i)
   {
       delete[] table3[i];
    }
    delete[] table3;

}

void fun_fourier_transform_DFT_method1(int N,double *table2)
{
    int i=N,j=N;
    double (*(table1)[2]);
    table1[0]=&table2[0];
    table1[1]=&table2[0+N];
    const double pi=3.141592653589793238462;

    double** table3 = new double*[2];
    for( int i = 0; i < 2; ++i)
    {
      table3[i] = new double[N];
    }
for (int i=0;i<2;i++)
{
    for(int j=0;j<N;j++)
    {
       table3[i][j]=0;
    }
}
for (int i=0;i<N;i++)
{
    for(int j=0;j<N;j++)
    {

          table3[0][i]=table3[0][i]+table1[0][j]*cos(i*j*2*pi/(float)N);
        table3[1][i]=table3[1][i]+table1[0][j]*sin(i*j*2*pi/(float)N);

    }
}
    for(int j=0;j<N;j++)
    {
      table1[0][j] =table3[0][j]+table3[1][j];
      table1[1][j] =0;
    }

   for( int i = 0; i < 2; ++i)
   {
       delete[] table3[i];
    }
    delete[] table3;

}

void fun_inverse_fourier_transform_DFT_method1(int N,double *table2)
{
    int i=N,j=N;
    double (*(table1)[2]);
    table1[0]=&table2[0];
    table1[1]=&table2[0+N];
    const double pi=3.141592653589793238462;

    double** table3 = new double*[2];
    for( int i = 0; i < 2; ++i)
    {
      table3[i] = new double[N];

    }
for (int i=0;i<2;i++)
{
    for(int j=0;j<N;j++)
    {
       table3[i][j]=0;
    }
}
for (int i=0;i<N;i++)
{
    for(int j=0;j<N;j++)
    {

          table3[0][i]=table3[0][i]+table1[0][j]*cos(i*j*2*pi/(float)N);
        table3[1][i]=table3[1][i]+table1[0][j]*sin(i*j*2*pi/(float)N);

    }
}
    for(int j=0;j<N;j++)
    {
      table1[0][j] =(table3[0][j]+table3[1][j])/N;//tylko tym się różni od normalnej transformnacji
      table1[1][j] =0;
    }

   for( int i = 0; i < 2; ++i)
   {
       delete[] table3[i];
    }
    delete[] table3;

}

void fun_fourier_transform_DFT_method2(int N,double *table2)
{
    int i=N,j=N;
    double (*(table1)[2]);
    table1[0]=&table2[0];
    table1[1]=&table2[0+N];
    const double pi=3.141592653589793238462;

    double** table3 = new double*[2];
    for( int i = 0; i < 2; ++i)
    {
      table3[i] = new double[N];

    }
for (int i=0;i<2;i++)
{
    for(int j=0;j<N;j++)
    {
       table3[i][j]=0;
    }
}
for (int i=0;i<N;i++)
{
    for(int j=0;j<N;j++)
    {

          table3[0][i]=table3[0][i]+table1[0][j]*cos(i*j*2*pi/(float)N);
        table3[1][i]=table3[1][i]+table1[0][j]*sin(i*j*2*pi/(float)N);

    }
}
    for(int j=0;j<N;j++)
    {
      table1[0][j] =table3[0][j];
      table1[1][j] =table3[1][j];
    }

   for( int i = 0; i < 2; ++i)
   {
       delete[] table3[i];
    }
    delete[] table3;

}

void fun_inverse_fourier_transform_DFT_method2(int N,double *table2)
{
    int i=N,j=N;
    double (*(table1)[2]);
    table1[0]=&table2[0];
    table1[1]=&table2[0+N];
    const double pi=3.141592653589793238462;

    double** table3 = new double*[4];
    for( int i = 0; i < 4; ++i)
    {
      table3[i] = new double[N];
    }
for (int i=0;i<4;i++)
{
    for(int j=0;j<N;j++)
    {
       table3[i][j]=0;
    }
}
for (int i=0;i<N;i++)
{
    for(int j=0;j<N;j++)
    {

          table3[0][i]=table3[0][i]+table1[0][j]*cos(i*j*2*pi/(float)N);
        table3[1][i]=table3[1][i]+table1[0][j]*sin(i*j*2*pi/(float)N);
          table3[2][i]=table3[2][i]+table1[1][j]*cos(i*j*2*pi/(float)N);
        table3[3][i]=table3[3][i]+table1[1][j]*sin(i*j*2*pi/(float)N);

    }
}
    for(int j=0;j<N;j++)
    {
      table1[0][j] =(table3[0][j]+table3[1][j]+table3[2][j]+table3[3][j])/N;//tylko tym się różni od normalnej transformnacji
      table1[1][j] = 0;
    }

   for( int i = 0; i < 4; ++i)
   {
       delete[] table3[i];
    }
    delete[] table3;

}

void fun_fourier_transform_DFT_method3(int N,double *table2)
{
    int i=N,j=N;
    double (*(table1)[2]);
    table1[0]=&table2[0];
    table1[1]=&table2[0+N];
    const double pi=3.141592653589793238462;

    double** table3 = new double*[2];
    for( int i = 0; i < 2; ++i)
    {
      table3[i] = new double[N];

    }
for (int i=0;i<2;i++)
{
    for(int j=0;j<N;j++)
    {
       table3[i][j]=0;
    }
}
for (int i=0;i<N;i++)
{
    for(int j=0;j<N;j++)
    {

          table3[0][i]=table3[0][i]+table1[0][j]*cos(i*j*2*pi/(float)N);
        table3[1][i]=table3[1][i]+table1[0][j]*sin(i*j*2*pi/(float)N);

    }
}
    for(int j=0;j<N;j++)
    {
      table1[0][j] =table3[0][j];
      table1[1][j] =table3[1][j];
    }

   for( int i = 0; i < 2; ++i)
   {
       delete[] table3[i];
    }
    delete[] table3;

}

void fun_inverse_fourier_transform_DFT_method3(int N,double *table2)
{
    int i=N,j=N;
    double (*(table1)[2]);
    table1[0]=&table2[0];
    table1[1]=&table2[0+N];
    const double pi=3.141592653589793238462;

    double** table3 = new double*[4];
    for( int i = 0; i < 4; ++i)
    {
      table3[i] = new double[N];

    }
for (int i=0;i<4;i++)
{
    for(int j=0;j<N;j++)
    {
       table3[i][j]=0;
    }
}
for (int i=0;i<N;i++)
{
    for(int j=0;j<N;j++)
    {

          table3[0][i]=table3[0][i]+table1[0][j]*cos(i*j*2*pi/(float)N);
        table3[1][i]=table3[1][i]+table1[0][j]*sin(i*j*2*pi/(float)N);
          table3[2][i]=table3[2][i]+table1[1][j]*cos(i*j*2*pi/(float)N);
        table3[3][i]=table3[3][i]+table1[1][j]*sin(i*j*2*pi/(float)N);

    }
}
    for(int j=0;j<N;j++)
    {
      table1[0][j] =(table3[0][j]+table3[3][j])/N;
      table1[1][j] =(table3[1][j]+table3[2][j])/N;

      if(table1[0][j]>=0)
      {
      table1[0][j]=table1[0][j]*table1[0][j];
      }
      else
      {
      table1[0][j]=table1[0][j]*table1[0][j]*-1;
      }

      if(table1[1][j]>=0)
      {
      table1[1][j]=table1[1][j]*table1[1][j];
      }
      else
      {
      table1[1][j]=table1[1][j]*table1[1][j]*-1;
      }

      if(table1[0][j]+table1[1][j]>=0)
      {
      table1[0][j]=sqrt(fabs(table1[0][j]+table1[1][j]));
      }
      else
      {
      table1[0][j]=sqrt(fabs(table1[0][j]+table1[1][j]))*-1;
      }

      table1[1][j]=0;
    }

   for( int i = 0; i < 4; ++i)
   {
       delete[] table3[i];
    }
    delete[] table3;

}
void fun_fourier_transform_DFT_method4(int N,double *table2)
{
    int i=N,j=N;
    double (*(table1)[2]);
    table1[0]=&table2[0];
    table1[1]=&table2[0+N];
    const double pi=3.141592653589793238462;

    double** table3 = new double*[2];
    for( int i = 0; i < 2; ++i)
    {
      table3[i] = new double[N];

    }
for (int i=0;i<2;i++)
{
    for(int j=0;j<N;j++)
    {
       table3[i][j]=0;
    }
}
for (int i=0;i<N;i++)
{
    for(int j=0;j<N;j++)
    {

          table3[0][i]=table3[0][i]+table1[0][j]*cos(i*j*2*pi/(float)N);
        table3[1][i]=table3[1][i]+table1[0][j]*sin(i*j*2*pi/(float)N);

    }
}
    for(int j=0;j<N;j++)
    {
      table1[0][j] =table3[0][j]-table3[1][j];
      table1[1][j] =0;
    }

   for( int i = 0; i < 2; ++i)
   {
       delete[] table3[i];
    }
    delete[] table3;

}

void fun_inverse_fourier_transform_DFT_method4(int N,double *table2)
{
    int i=N,j=N;
    double (*(table1)[2]);
    table1[0]=&table2[0];
    table1[1]=&table2[0+N];
    const double pi=3.141592653589793238462;

    double** table3 = new double*[2];
    for( int i = 0; i < 2; ++i)
    {
      table3[i] = new double[N];

    }
for (int i=0;i<2;i++)
{
    for(int j=0;j<N;j++)
    {
       table3[i][j]=0;
    }
}
for (int i=0;i<N;i++)
{
    for(int j=0;j<N;j++)
    {

          table3[0][i]=table3[0][i]+table1[0][j]*cos(i*j*2*pi/(float)N);
        table3[1][i]=table3[1][i]+table1[0][j]*-1*sin(i*j*2*pi/(float)N);

    }
}
    for(int j=0;j<N;j++)
    {
      table1[0][j] =(table3[0][j]+table3[1][j])/N;
      table1[1][j] =0;
    }

   for( int i = 0; i < 2; ++i)
   {
       delete[] table3[i];
    }
    delete[] table3;

}

void fun_fourier_transform_DFT_method6_full_complex(int N,std::complex<double> tab[])
{
    const double pi=3.141592653589793238462;
    std::complex<double> tab2[16]={};    // tab2[]==N
    std::complex<double>  w[1]={{1,1}};
    std::complex<double>  w2[1]={{1,1}};



for (int i=0;i<N;i++)
{
    for(int j=0;j<N;j++)
    {
          //complex number method:
          w[0].real()=cos(i*j*2*pi/(float)N);
          w[0].imag()=(-sin(i*j*2*pi/(float)N));
          tab2[i]=tab2[i]+tab[j]*w[0];

    }
}
/*
    for(int j=0;j<N;j++)
    {
      tab[j].real() =tab2[j].real()*2/N;
      tab[j].imag() =tab2[j].imag()*2/N;
    }
*/
}

//http://inverse-fast-fourier-transform-fft.blogspot.com/
