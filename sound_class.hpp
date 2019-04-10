#pragma once

class Sound{
public:
    int samplingFrequency;//Hz
    double length;//second
    signed short* waveData;
    int amplitude=32767;//量子化ビット数は16
    
    Sound(){}
    ~Sound(){
        delete waveData;
    }
    
    void New(int f, double l){
        samplingFrequency=f;
        length=l;
        
        waveData = new signed short[samplingFrequency*length];
        for(int i=0;i<samplingFrequency*length;i++)
            waveData[i]=0;
    }
    
    void play(){
        //OpenALの下準備　おまじない的な
        ALCdevice *device = alcOpenDevice(NULL);
        ALCcontext *context = alcCreateContext(device, NULL);
        alcMakeContextCurrent(context);
        
        //バッファ(保存領域)とソース(音源)を宣言
		ALuint buffer[2];
		ALuint source;
        //それを生成
        alGenBuffers(2, buffer);
        alGenSources(1, &source);
        
        //fail safe
		for (int i = 0; i < samplingFrequency*length; i++) {
			if (waveData[i] > amplitude)
				waveData[i] = amplitude;
			if (waveData[i] < -amplitude)
				waveData[i] = -amplitude;
		}
        
        //バッファに音源データを入れる
		alBufferData(buffer[0], AL_FORMAT_MONO16, &waveData[0], length*samplingFrequency * sizeof(signed short), samplingFrequency);
        //ソースにバッファを適用
        alSourcei(source, AL_BUFFER, buffer[0]);
        //ループ再生をON
        //alSourcei(source, AL_LOOPING, AL_TRUE);
        //ソースを再生！
        alSourcePlay(source);
        
        //ここで一時停止
        //system("PAUSE");//available only by WindowsOS
        cv::waitKey(length*1000);//fuction of openCV
        
        // バッファの破棄
        alDeleteBuffers(2, buffer);
        // ソースの破棄
        alDeleteSources(1, &source);
        
        //OpenALの後始末
        alcMakeContextCurrent(nullptr);
        alcDestroyContext(context);
        alcCloseDevice(device);
    }
    
    /********some types of wave********/
    void sine_curve(int frequency,double startTime,double endTime){
        for (int i = startTime*samplingFrequency; i < endTime*samplingFrequency; i++) {
            if(i>length*samplingFrequency) break;
            waveData[i] = amplitude * sin(2 * M_PI*i * frequency / samplingFrequency);
        }
    }
    
    void saw_curve(int frequency, double startTime, double endTime){
        bool up=true;
        double val=-amplitude;
        int count=0;
        for(int i=startTime*samplingFrequency;i<endTime*samplingFrequency;i++){
            if(i>length*samplingFrequency) break;
            if(val>amplitude){
                val=amplitude;
            }
            if(val<-amplitude){
                val=-amplitude;
            }
            if(count==samplingFrequency/(2*frequency)){
                if(up) up=false;
                else up=true;
                count=0;
            }
            waveData[i]=val;
            if(up) val+=4*frequency*amplitude/samplingFrequency;
            else val-=4*frequency*amplitude/samplingFrequency;
            count++;
        }
    }
    
    void square_curve(int frequency,double startTime,double endTime){
        int count=0;
        bool up=true;
        for(int i=startTime*samplingFrequency;i<endTime*samplingFrequency;i++){
            if(i>length*samplingFrequency) break;
            if(count==samplingFrequency/(2*frequency)){
                if(up) up=false;
                else up=true;
                count=0;
            }
            waveData[i]=amplitude/2*(up?1:-1);
            count++;
        }
    }
    /****************/
    
    
    /********some effects********/
    void amplify(double startTime,double endTime,double gain){
        for(int i=startTime*samplingFrequency;i<endTime*samplingFrequency;i++){
            if(i>length*samplingFrequency) break;
            waveData[i]*=gain;
            if(waveData[i]>amplitude) waveData[i]=amplitude;
            if(waveData[i]<-amplitude) waveData[i]=-amplitude;
        }
    }
    
    void attack(double startTime, double endTime){
        for(int i=startTime*samplingFrequency;i<endTime*samplingFrequency;i++){
            if(i>length*samplingFrequency) break;
            waveData[i]*=(i-startTime*samplingFrequency)/((endTime-startTime)*samplingFrequency);
        }
    }
    
    void decay(double startTime, double endTime){
        for(int i=startTime*samplingFrequency;i<endTime*samplingFrequency;i++){
            if(i>length*samplingFrequency) break;
            waveData[i]*=(endTime*samplingFrequency-i)/((endTime-startTime)*samplingFrequency);
        }
    }
    
    void smoothing(double startTime,double endTime,int n){
        signed short* tmp;
        tmp = new signed short[samplingFrequency*length];
        
        for(int i=startTime*samplingFrequency;i<endTime*samplingFrequency;i++){
            if(i>length*samplingFrequency) break;
            double sum=0;
            int count=0;
            for(int j=-n;j<=n;j++){
                if(i+j>length*samplingFrequency) break;
                if(i+j<0) continue;
                sum+=waveData[i+j];
                count++;
            }
            double average=sum/count;
            tmp[i]=average;
        }
        
        for(int i=startTime*samplingFrequency;i<endTime*samplingFrequency;i++){
            if(i>length*samplingFrequency) break;
            waveData[i]=tmp[i];
        }
        delete tmp;
    }
    /****************/
    
    
    void draw_wave(double startTime,double endTime){
        int width=800;
        int height=400;
        Mat win(Size(width,height),CV_8U,Scalar::all(255));
        for(int i=0;i<width;i++){
            rectangle(win,Point(i,height/2+(int)(1.0*height/2/amplitude*waveData[(int)((i+startTime*samplingFrequency)*(endTime*samplingFrequency-startTime*samplingFrequency)/width)])),Point(i,height/2),Scalar(0));
        }
        imshow("win",win);
        waitKey();
    }
    
    void input(char filename[64]){
        WAV_PRM prm_in;
        double *data_in;
        
        //wavファイルの読み込み
        data_in = audio_read(&prm_in, filename);
        
        samplingFrequency=prm_in.fs;
        if(prm_in.bits==16) amplitude=32767;
        length=1.0*prm_in.L/samplingFrequency;
        
        waveData = new signed short[(int)(samplingFrequency*length)];
        for(int i=0;i<samplingFrequency*length;i++)
            waveData[i]=0;
        
        for(int i=0;i<prm_in.L;i++){
            waveData[i]=data_in[i]*amplitude;
        }
        
        //メモリ解放
        free(data_in);
    }
    
    void output(char filename[64]){
        //変数宣言
        WAV_PRM prm_out;
        double *data_out;
        int i;
        
        //パラメータコピー
        prm_out.fs =  samplingFrequency;
        if(amplitude==32767) prm_out.bits = 16;
        prm_out.L =  (int)(length*samplingFrequency);
        
        //データコピー(実際にはこの代わりにエフェクト処理をかける)
        data_out = (double*)calloc(prm_out.L, sizeof(double)); //メモリの確保

        for (i = 0;i < prm_out.L;i++){
            data_out[i] = waveData[i]*1.0/amplitude;
        }
        
        //書き込み
        audio_write(data_out, &prm_out, filename);
        
        //メモリ解放
        free(data_out);
    }
};
