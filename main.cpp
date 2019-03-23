#include <iostream>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#pragma (lib,"opencv_world342.lib")
#include <OpenAL/al.h>
#include <OpenAL/alc.h>
#pragma comment(lib, "OpenAL32.lib")
using namespace std;
using namespace cv;

#define SAMPLINGFREQUENCY 44100
#define _USE_MATH_DEFINES

#include "wav_inoutput.hpp"
#include "sound_class.hpp"
#include "sound_process.hpp"
#include "Fourier.hpp"

int main() {
    srand((unsigned int)time(NULL));
    int i,j;
    
    Sound* sound1;
    sound1=new Sound();
    sound1->input("output.wav");
    sound1->play();
    delete sound1;
    
    Sound* sound2;
    sound2=new Sound();
    sound2->New(SAMPLINGFREQUENCY,1);
    sound2->sine_curve(3000,0,1);
    sound2->play();
    //delete sound2;
    
    Sound* sound3;
    sound3=new Sound();
    sound3->New(SAMPLINGFREQUENCY,1);
    sound3->sine_curve(1090,0,1);
    sound3->amplify(0,1,0.8);
    sound3->play();
    //delete sound3;
    
    Sound* sound4;
    sound4=new Sound();
    sound4->New(SAMPLINGFREQUENCY,1);
    sound4->sine_curve(440,0,1);
    sound4->play();
    //delete sound5;
    
    Sound* sound5;
    sound5=new Sound();
    sound5->New(SAMPLINGFREQUENCY,1);
    synthesize(sound2,sound3,sound5,0,1);
    synthesize(sound5,sound4,sound5,0,1);
    sound5->play();
    Fourier(sound5,0,0.1);
    for(i=1;i<10;i++)
        copy(sound5,sound5,0,0.1*i,0.1);
    sound5->play();
    delete sound2;delete sound3;delete sound4;delete sound5;
    
    return 0;
}
