
#include <iostream>
#include "TF1.h"
#include "TApplication.h"
#include "TCanvas.h"
//#include "TSystem.h"
#include <cmath>

using namespace std;

double gauss(double x, double mu, double sigma){
    return exp( -0.5*(x-mu)*(x-mu)/(sigma*sigma))/(sqrt(2*M_PI)*sigma);
}

Double_t myfunction(Double_t* x, Double_t* par){
    double delta_b = par[0];
    double sigma = par[1];
    double A = par[2];
    int n_clock = (int) par[3]; // nclocks * 25 ns = recovery time
    double sum = 0 ;
    for(unsigned int i=0;i<n_clock;i++){
        sum+=A*gauss(*x,delta_b+delta_b/n_clock*i, sigma); // supose a linear evolution of the baseline
        //try to directly put the formula of gaussian instead of calling a function
    }

    sum+=(1-A*n_clock)*gauss(*x,0,sigma);

    //add a term for all events outside the HIP recovery time
    // gaussian centered one 0 with a width = sigma
    // amplitude is defined with respect to A
    // amplitude  = 1-f_HIP*n_clocks
    return sum;
}

int main(int argc, char** argv){

    TApplication* theApp = new TApplication("App", &argc, argv);
    double freq = 1E-3;
    double rtime = 7e-4; // recovery time
    double lhctime = 25e-6; // time btw 2 bunch crossing
    double phip = 1e-3;  // probability of HIP
    double delta_b = -15; // variation of the baseline due to HIP
    double sigma = 5;     // rms of the baseline

    int n_clock = rtime/lhctime;
    TF1* f = new TF1("tot",myfunction,-100,20,4);
    /*
    f->SetParameter(0,delta_b);
    f->SetParameter(1, sigma);
    f->SetParameter(2, phip);
    f->SetParameter(3, (double)n_clock);
    */
    f->SetParameters(delta_b,sigma,phip,(double)n_clock);
    f->SetParNames("delta_b","sigma","prob","nclock");
    Double_t par[] = {delta_b, sigma, phip, (double) n_clock};
    double x = 20;
    cout<<gauss(1,0,1)<<endl;
    cout<<myfunction(&x,par)<<endl;
    //cout<<f->Eval(-5)<<endl; 

    f->Draw();

    TF1* g = new TF1("g","gaus",0,10);
    g->SetParameter(0,10);
    g->SetParameter(1,5);
    g->SetParameter(2,1);
    //g->Draw();
    //c2->Draw();

    theApp->Run();
    delete theApp;
}



