#pragma once
#include "complex.hpp"
#include "FFT.hpp"

void cutoff(const int N, double T, comp G[N], double magnification, double frequency);

void Fourier(Sound* sound,double startTime, double endTime){
    const int N=2048;//maxfrequency = N/2T
    comp G[N];
    
    FFT(sound,startTime,endTime,N,G);
    
    //cutoff(N,endTime-startTime,G,0,300);
    
    IFFT(sound,startTime,endTime,N,G);
    //FFT(sound,startTime,endTime,N,G);
}

void cutoff(const int N, double T, comp G[N], double magnification, double frequency){
    if(frequency>N/(2*T)) return;
    int nearest_index=frequency*T+0.5;//round
    
    int width=10;
    for(int i=-width;i<=width;i++){
        G[nearest_index+i].re*=magnification  ;G[nearest_index-i].im*=magnification;
        G[N-nearest_index+i].re*=magnification;G[N-nearest_index-i].im*=magnification;
    }
}
