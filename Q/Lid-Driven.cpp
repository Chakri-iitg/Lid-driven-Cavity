#include <iostream>
#include<fstream>
#include<cmath>
#include<bits/stdc++.h>

using namespace std;

int main()
{
  cout<<fixed<<setprecision(6);
   int i,j,m,n,p;
   double dx,dy,x,y;

  cout<<"Enter the values of Grid points MxN: "<<endl;
  cin>>m>>n;


    p=(m-2)*(n-2);

    dx=1.0/(m-1);
    dy=1.0/(n-1);

    double beta =(dx/dy);
    double d=pow(beta,2);
    double k=(1.0/(2.0*(1+d)));
    double Re ;

    cout<<"Enter the value of Reynolds number: "<<endl;
    cin>>Re;

    double psi_n[m][n],omega_n[m][n],psi_o[m][n],omega_o[m][n];
    double u[m][n],v[m][n];
    double error_psi=0.0,error_omega=0.0;
    int iteration=0;


    for (j=0;j<n;j++){
        for(i=0;i<m;i++){

            if(j==0){

                u[i][j]=0.0;
                v[i][j]=0.0;
                psi_n[i][j]=0.0;

            }
            else if (j==(n-1)){

                 u[i][j]=1.0;
                 v[i][j]=0.0;
                 psi_n[i][j]=0.0;
                }
            else if(i==0){

                 u[i][j]=0.0;
                 v[i][j]=0.0;
                 psi_n[i][j]=0.0;
            }
            else if(i==(m-1)){

                 u[i][j]=0.0;
                 v[i][j]=0.0;
                 psi_n[i][j]=0.0;
                }
            else {
                 u[i][j]=0.0;
                 v[i][j]=0.0;
                 psi_n[i][j]=0.0;
                }
            }
        }

       for(j=0;j<n;j++){
        for(i=0;i<m;i++){


            if(j==0){

            omega_n[i][j]=(2.0/pow(dy,2))*(psi_n[i][j]-psi_n[i][j+1]);

            }
            else if (j==(n-1)){

              omega_n[i][j]=((2.0/pow(dy,2))*(psi_n[i][j]-psi_n[i][j-1]))-((2.0/dy)*u[i][j]);

                }
            else if(i==0){

               omega_n[i][j]=(2.0/pow(dx,2))*(psi_n[i][j]-psi_n[i+1][j]);

            }
            else if(i==(m-1)){

                omega_n[i][j]=(2.0/pow(dx,2))*(psi_n[i][j]-psi_n[i-1][j]);


                }
            else {

                 omega_n[i][j]=0.0;

                }
            }
        }


       do{

        for(j=0;j<n;j++){
            for(i=0;i<m;i++){

                psi_o[i][j]=psi_n[i][j];
                omega_o[i][j]=omega_n[i][j];

                }
            }


         for(j=1;j<(n-1);j++){
            for(i=1;i<(m-1);i++){


                psi_n[i][j]=(k*(psi_n[i+1][j]+psi_n[i-1][j]+(d*(psi_n[i][j+1]+psi_n[i][j-1]))+(pow(dx,2)*omega_n[i][j])));

                }
            }



        for(j=1;j<(n-1);j++){
            for(i=1;i<(m-1);i++){


              omega_n[i][j]=(k*(((1.0-((psi_n[i][j+1]-psi_n[i][j-1])*((beta*Re)/4.0)))*omega_n[i+1][j])
                             +((1.0+((psi_n[i][j+1]-psi_n[i][j-1])*((beta*Re)/4.0)))*omega_n[i-1][j])
                             +((1.0+((psi_n[i+1][j]-psi_n[i-1][j])*(Re/(4.0*beta))))*(d*omega_n[i][j+1]))
                             +((1.0-((psi_n[i+1][j]-psi_n[i-1][j])*(Re/(4.0*beta))))*(d*omega_n[i][j-1]))));

                }
            }



             for(j=0;j<n;j++){
          for(i=0;i<m;i++){


            if(j==0){

            omega_n[i][j]=(2.0/pow(dy,2))*(psi_n[i][j]-psi_n[i][j+1]);

            }
            else if (j==(n-1)){

              omega_n[i][j]=((2.0/pow(dy,2))*(psi_n[i][j]-psi_n[i][j-1]))-((2.0/dy)*u[i][j]);

                }
            else if(i==0){

               omega_n[i][j]=(2.0/pow(dx,2))*(psi_n[i][j]-psi_n[i+1][j]);

            }
            else if(i==(m-1)){

                omega_n[i][j]=(2.0/pow(dx,2))*(psi_n[i][j]-psi_n[i-1][j]);


                }

            }
        }


         error_psi=0.0;
         error_omega=0.0;


         for(j=1;j<(n-1);j++) {
            for(i=1;i<(m-1);i++){

                error_psi=error_psi+pow((psi_n[i][j]-psi_o[i][j]),2.0);
                error_omega=error_omega+pow((omega_n[i][j]-omega_o[i][j]),2.0);

                }
            }


          error_psi=sqrt(error_psi/p);
          error_omega=sqrt(error_omega/p);

            cout<<"Iteration-"<<iteration<<"\t\t";
            cout<<error_psi<<"\t\t"<<error_omega<<endl;

           iteration++;
        }while(error_psi>pow(10,-6)||error_omega>pow(10,-6));


        for(j=1;j<(n-1);j++){
           for(i=1;i<(m-1);i++){

            u[i][j]=(0.5/dy)*(psi_n[i][j+1]-psi_n[i][j-1]);

            v[i][j]=(-0.5/dx)*(psi_n[i+1][j]-psi_n[i-1][j]);



          }
         }

    ofstream stream("224103307_Assignment3-StreamFunction.dat");
     stream<<fixed<<setprecision(3);
    stream<<"TITLE = STREAM FUNCTION"<<endl<<"VARIABLES = \"x\", \"y\", \"psi\""<<endl;
    stream<<"Zone T = \"BLOCK1\", i=100, j=100, F=POINT\n\n"<<endl;
    for(int i=0; i<=m; i++)
    {
        for(int j=0; j<=n; j++)
        {
             x= i*dx;
             y= j*dy;
            stream<<x<<"   "<<y<<"   "<<psi_n[i][j]<<endl;
        }
    }
    ofstream vorticity("224103307_Assignment3-Vorticity.dat");
    vorticity<<fixed<<setprecision(3);
    vorticity<<"TITLE = Vorticity"<<endl<<"VARIABLES = \"x\", \"y\", \"omega\""<<endl;
    vorticity<<"Zone T = \"BLOCK1\", i=100, j=100, F=POINT\n\n"<<endl;
    for(int i=0; i<=m; i++)
    {
        for(int j=0; j<=n; j++)
        {
             x= i*dx;
             y= j*dy;
            vorticity<<x<<"   "<<y<<"   "<<omega_n[i][j]<<endl;
        }
    }
    ofstream x_velocity("224103307_Assignment3-U.dat");
    ofstream y_velocity("224103307_Assignment3-V.dat");

    x_velocity<<fixed<<setprecision(3);
    x_velocity<<"TITLE = U velocity"<<endl<<"VARIABLES = \"x\", \"y\", \"U\""<<endl;
    x_velocity<<"Zone T = \"BLOCK1\", i=100, j=100, F=POINT\n\n"<<endl;

     y_velocity<<fixed<<setprecision(3);
    y_velocity<<"TITLE = V velocity"<<endl<<"VARIABLES = \"x\", \"y\", \"V\""<<endl;
    y_velocity<<"Zone T = \"BLOCK1\", i=100, j=100, F=POINT\n\n"<<endl;

    for(int i=0; i<=m; i++)
    {
        for(int j=0; j<=n; j++)
        {
             x= i*dx;
             y= j*dy;
            x_velocity<<x<<"   "<<y<<"   "<<u[i][j]<<endl;
            y_velocity<<x<<"   "<<y<<"   "<<v[i][j]<<endl;

        }
    }

    ofstream plot1("3307-100u.dat");
    ofstream plot2("3307-100v.dat");
    plot1<<fixed<<setprecision(4);
    plot2<<fixed<<setprecision(4);

    for(int i=0; i<=m; i++)
    {    x= i*dx;
        plot1<<x<<" "<<v[i][50]<<endl;

    }
        for(int j=0; j<=n; j++)
        {

             y= j*dy;
            plot2<<u[50][j]<<"   "<<y<<endl;

        }












    return 0;
    }
