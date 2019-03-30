#pragma once
#include "complex.hpp"
#include "FFT.hpp"
void shift_wave(const int N, comp G[N], int time);

void Fourier(Sound* sound,double startTime, double endTime){
    if(startTime>=endTime) return;
    if(startTime<0) return;
    if(endTime>sound->length) return;
    
    const int N=2048;//maxfrequency = N/2T
    comp G[N];
    
    printf("max frequency = %f\n",N/(2*(endTime-startTime)));
    printf("interval frequency = %f\n",1/(2*(endTime-startTime)));
    
    static bool draw=true;
    
    
    FFT(sound,startTime,endTime,N,G,draw);
    
    
    IFFT(sound,startTime,endTime,N,G);
    FFT(sound,startTime,endTime,N,G,draw);
    draw=false;
}
