#pragma once

void synthesize(Sound* src1,Sound* src2,Sound* dst,double startTime,double endTime){
    for (int i = startTime*src1->samplingFrequency; i < endTime*src1->samplingFrequency; i++) {
        if(i<0 || i>src1->samplingFrequency*src1->length || i>src2->samplingFrequency*src2->length){
            break;
        }
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

void fast_forward(Sound* src, Sound* dst, double startTime, double endTime, double magnification){
    for(int i=startTime*src->samplingFrequency/magnification;i<endTime*src->samplingFrequency/magnification;i++){
        dst->waveData[i]=src->waveData[(int)(i*magnification)];
    }
}

void eliminate_vocal(Sound* soundL, Sound* soundR,Sound* dst,double balance) {
	if (soundR->length != soundL->length) {
		printf("load deffernt size file!\n");
		return;
	}
	if (soundR->samplingFrequency != soundL->samplingFrequency) {
		printf("load deffernt sampling frequency file!\n");
		return;
	}
	if (balance < -1 || balance > 1) {
		printf("balance should be in range from -1 to 1!\n");
		return;
	}

	dst->New(soundR->samplingFrequency, soundR->length);

	for (int i = 0; i < soundR->samplingFrequency*soundR->length; i++) {
		dst->waveData[i] = soundL->waveData[i] * ((-balance + 1) / 2) - soundR->waveData[i] * ((balance + 1) / 2);
	}
}
