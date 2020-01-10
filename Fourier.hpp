#pragma once
#include "complex.hpp"
#include "FFT.hpp"

void Fourier(Sound* sound,double startTime, double endTime,bool draw,int drawTime){
    if(startTime>=endTime) return;
    if(startTime<0) return;
    if(endTime>sound->length) return;
    
    const int N=2048;//maxfrequency = N/2T
    comp G[N];
    
    printf("max frequency = %f\n",N/(2*(endTime-startTime)));
    printf("interval frequency = %f\n",1/(2*(endTime-startTime)));
    
    
    FFT(sound,startTime,endTime,N,G,draw,drawTime);
    
    
    IFFT(sound,startTime,endTime,N,G);
    //FFT(sound,startTime,endTime,N,G,draw,drawTime);
}
