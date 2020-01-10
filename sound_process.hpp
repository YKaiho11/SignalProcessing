#pragma once

void synthesize(Sound* src1,Sound* src2,Sound* dst,double startTime,double endTime){
    for (int i = startTime*src1->samplingFrequency; i < endTime*src1->samplingFrequency; i++) {
        if(i<0 || i>src1->samplingFrequency*src1->length || i>src2->samplingFrequency*src2->length){
            break;
        }
        dst->waveData[i]=src1->waveData[i]+src2->waveData[i];
		if (dst->waveData[i] > 32767) dst->waveData[i] = 32767;
		if (dst->waveData[i] <- 32767) dst->waveData[i] = -32767;
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
        signed short db=soundL->waveData[i] * ((-balance + 1) / 2.0) - soundR->waveData[i] * ((balance + 1) / 2.0);
        if(soundR->amplitude<db) db=soundR->amplitude;
        if(soundR->amplitude<-db) db=-soundR->amplitude;
        dst->waveData[i] = db;
        for(int j=0;j<100;j++){
            if(i==(int)(soundR->samplingFrequency*soundR->length*j/100.0))
                printf("%d%% done!\n",j);
        }
	}
}

void stereo2monoral(Sound* stereo,Sound* R,Sound* L){
    L->New(stereo->samplingFrequency,stereo->length/2);
    R->New(stereo->samplingFrequency,stereo->length/2);
    
    for(int i=0;i<stereo->samplingFrequency*stereo->length;i++){
        if(i%2==0){
            R->waveData[i/2]=stereo->waveData[i];
        }
        else{
            L->waveData[(i-1)/2]=stereo->waveData[i];
        }
    }
}

Sound* make_eliminatedSound(char* path_to_stereo){
    Sound* L=new Sound();
    Sound* R=new Sound();
    Sound* stereo=new Sound();
    stereo->input(path_to_stereo);
    
    stereo2monoral(stereo,R,L);
    
    Sound* dst=new Sound();
    
    eliminate_vocal(L,R,dst,0);
    
    return dst;
}

Sound* cut_first(char* path) {
	Sound* src = new Sound();
	src->input(path);

	int start;
	for (int i = 0; i < src->length * src->samplingFrequency; i++) {
		if (src->waveData[i] != 0) {
			start = i;
			break;
		}
	}

	Sound* dst = new Sound();
	dst->New(src->samplingFrequency, src->length - src->samplingFrequency*start);
	for (int i = 0; i < dst->length*dst->samplingFrequency; i++) {
		dst->waveData[i] = src->waveData[(int)(i + start * src->samplingFrequency)];
	}

	return dst;
}