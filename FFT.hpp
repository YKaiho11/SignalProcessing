#ifndef FFT_HPP
#define FFT_HPP

#include <complex>
#include <cmath>
#include <utility>

typedef std::complex<double> complex_t;

void F(int N, int s, int q, int d, complex_t* x) {
	// N : 系列長
	// s : ストライド(偶数成分や奇数成分の歩幅)
	// q : ブロックの開始位置
	// d : ブロックの開始位置のビット反転位置計算用変数
	// x : FFT する系列(入出力)
	const int m = N / 2;
	const double theta0 = 2 * M_PI / N;

	if (N > 1) {
		for (int p = 0; p < m; p++) {
			const complex_t wp = complex_t(cos(p*theta0), -sin(p*theta0));
			const complex_t a = x[q + p + 0];
			const complex_t b = x[q + p + m];
			x[q + p + 0] = a + b;
			x[q + p + m] = (a - b) * wp;
		}
		F(N / 2, 2 * s, q + 0, d + 0, x); // 偶数位置成分
		F(N / 2, 2 * s, q + m, d + s, x); // 奇数位置成分
	}
	else if (q > d) std::swap(x[q], x[d]); // ビット反転並べ替え

}


void fft(int N, complex_t* x) {// フーリエ変換
// N : 系列長
// x : フーリエ変換する系列(入出力)
	F(N, 1, 0, 0, x);
}

void ifft(int N, complex_t* x) {// 逆フーリエ変換
// N : 系列長
// x : 逆フーリエ変換する系列(入出力)
	for (int p = 0; p < N; p++) x[p] = conj(x[p]);
	F(N, 1, 0, 0, x);
	for (int k = 0; k < N; k++) x[k] = conj(x[k]) / double(N);
}

double hamming_window(const int N, int i) {
	return 0.54 - 0.46*cos(2 * M_PI*i / N);
}


void FFT(Sound* sound, double startTime, double endTime, const int N, comp* G, bool drawResult,int drawTime){
    int i;
    complex_t* x;
    x=(complex_t*)malloc(N*sizeof(complex_t));
    
    for(i=0;i<N;i++){
        x[i].real(sound->waveData[(int)(sound->samplingFrequency*((endTime-startTime)/N*i+startTime))]*hamming_window(N,i));
        x[i].imag(0);
    }
    
    fft(N,x);
    
    for(i=0;i<N;i++){
        G[i].re=x[i].real();
        G[i].im=x[i].imag();
    }
    
    //max frequency
    double max=0;
    for(i=0;i<N/2;i++){
        if(abs(G[i])>max) max=abs(G[i]);
    }
    
    if(drawResult){
        /****drawing****/
        int width=800;
        int height=400;
        int height_offset=50;
        Mat win(Size(width,height+height_offset),CV_8U,Scalar::all(252));
        for(i=0;i<width;i++){
            rectangle(win,Point(i,height-abs(G[i*N/2/width])*height/max),Point(i,height),Scalar(0));
        }

        for(i=0;i<5;i++){
            char c[16];
            sprintf(c,"|%d",(int)(1.0*i*N/2/(endTime-startTime)/5));
            putText(win, c, Point(i*width/5,height+height_offset/2), FONT_HERSHEY_SIMPLEX, 1.2, Scalar(0,0,200), 1, 0);
        }
        
        imshow("FFT result",win);
        //imwrite("output_FFT.png",win);
        waitKey(drawTime);
    }
    
    free(x);
}

void IFFT(Sound* sound, double startTime, double endTime, const int N, comp* G){
    complex_t* x;
    x=(complex_t*)malloc(N*sizeof(complex_t));
    
    for(int i=0;i<N;i++){
        x[i].real(G[i].re);
        x[i].imag(G[i].im);
    }
    
    ifft(N,x);
    
    bool isFirst=true;
    double val_pre;
    for(int i=0;i<N;i++){
        double val_tmp=x[i].real()/hamming_window(N,i);
		sound->waveData[(int)(sound->samplingFrequency*((endTime - startTime) / N * i + startTime))] = val_tmp;
        
        
        if(isFirst) {
            isFirst=false;
            val_pre=val_tmp;
        }
        else{
            int a=(int)(sound->samplingFrequency*((endTime-startTime)/N*(i-1)+startTime))+1;
            int b=(int)(sound->samplingFrequency*((endTime-startTime)/N*i+startTime));
            for(int j=a;j<=b;j++){
                if(a!=b){
                    sound->waveData[j]=1.0*(j-a)*(val_tmp-val_pre)/(b-a)+val_pre;
                }
            }
            val_pre=val_tmp;
        }
    }
    free(x);
}


#endif