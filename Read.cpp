#include <iostream>
#include <cmath>
#include <complex>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <sstream>
#include <cmath>
#include "basic.h"
using namespace std;
double *h_x,*h_z;
double *J_x,*J_z;
double *K_x;
int n,p,groundstatenum,*flag;
dcomp iota;

int readfunc(string filename)
{
    string line,line1,word;
    int i,j,counter, counter1=0;
    ifstream infile;
    infile.open(filename);

    //infile.open("/home/home3/student/bo561338/Desktop/vrinda/18spin/hamiltonian_" + std::to_string(filenumber) + ".dat");

    if(!infile){
        cout<<"Cannot open input file"<<endl;
    }
    else{

        while(getline(infile,line)){
            if (line[1]=='!'){
                continue;
            }
            else{

                //cout<<line<<endl;
                //string check=;
                word=" NUMBER OF SITES";
                if(line.find(word)!=std::string::npos){
                    getline(infile,line);

                    line1=line[1];
                    line1+=line[2];
                    //cout<<line1<<endl;
                    n=atoi(line1.c_str());
                    p=1<<n;

                    h_x=new double [n];
                    h_z=new double [n];
                    J_x=new double [n*n];
                    J_z=new double [n*n];
                    K_x=new double [n*n*n];
                    flag=new int [n*n];
                    //cout<<n<<endl;
                    for(i=0;i<n;i++){
                        for(j=0;j<n;j++){
                            flag[j+n*i]=0;
                        }}
                    for (i=0;i<n;i++)
                    {
                        h_z[i]=0;
                        for (j=0;j<n;j++)
                        {
                            J_z[j+n*i]=0;

                        }
                    }
                }


                word=" HEISENBERG INTERACTION :";
                if (line.find(word)!=std::string::npos){
                    //cout<<"YES"<<endl;
                  do{
                    getline(infile,line);
                    if (line.find("END")==std::string::npos){
                        counter=0;
                        //cout<<line<<endl;
                        line1=line[1];
                        line1+=line[2];
                        i=atoi(line1.c_str())-1;
                        line1=line[6];
                        line1+=line[7];
                        j=atoi(line1.c_str())-1;
                        stringstream ss(line);
                        while(getline(ss,word,',')){
                            counter++;
                            
                            if (counter==5){
                                J_z[j+n*i]=atof(word.c_str());
                                flag[j+n*i]=1;
                                //flag[j+n*i]=1;
                            }
                        }
                        //cout<<"#J_z["<<i<<"]["<<j<<"]= "<<J_z[i][j]<<"   ";
                        //cout<<"i= "<<i<<" j= "<<j<<endl;
                    }

                } while(line.find("END")==std::string::npos);

                }
                 word=" SPIN :";
                if(line.find(word)!=std::string::npos){
                    counter1++;

                    //cout<<line<<endl;
                     getline(infile,line);
                     //cout<<line<<endl;
                        for (int k=0;k<line.length();k++){

                            if (line[k]=='D'){
                                line[k]='E';
                            }
                        }
                        //cout<<line<<endl;
                        counter=0;
                        stringstream ss(line);
                        while(getline(ss,word,',')){
                            counter++;
                            //cout<<word<<"       ";
                            //cout<<"counter= "<<counter<<"counter1= "<<counter1<<endl;
                            if (counter==7){
                                //cout<<atof(word.c_str());
                                h_z[counter1-1]=atof(word.c_str());
                                //cout<<"h_z["<<counter1-1<<"]= "<<h_z[counter1-1]<<endl;
                            }
                        }
                        //cout<<"#h_z["<<counter1-1<<"]= "<<h_z[counter1-1]<<"     ";

                }
                word=" GROUND STATE NUMBER";
                counter=0;
                if(line.find(word)!=std::string::npos){
                    stringstream ss(line);
                    while(getline(ss,word,',')){
                        counter++;
                        if(counter==1){
                        counter1=0;
                        stringstream ss2(word);
                        while(getline(ss2,word,'=')){
                            counter1++;
                            if (counter1==2){
                                groundstatenum=pow(2,n)-atof(word.c_str())-1;

                            }

                        }
                    }
                }

                }
            }
        }
    }
    return 0;

}
