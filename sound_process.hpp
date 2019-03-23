#pragma once

void synthesize(Sound* src1,Sound* src2,Sound* dst,double startTime,double endTime){
    for (int i = startTime*src1->samplingFrequency; i < endTime*src1->samplingFrequency; i++) {
        dst->waveData[i]=(src1->waveData[i]+src2->waveData[i])/2;
    }
}

void copy(Sound* src, Sound* dst, double startTime_src, double startTime_dst, double duration){
    for(int i=0;i<duration*src->samplingFrequency;i++){
        int t_src=startTime_src*src->samplingFrequency+i;
        int t_dst=startTime_dst*dst->samplingFrequency+i;
        dst->waveData[t_dst]=src->waveData[t_src];
    }
}
