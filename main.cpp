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
#include "tone_list.hpp"

int main() {
    Tone t;
    srand((unsigned int)time(NULL));
    int i,j;
    
    Sound* sound1;
    sound1=new Sound();
    sound1->input("piano/C#5.wav");
    int count=50;
    for(i=5;i<count;i++){
        Fourier(sound1,sound1->length/count*i,sound1->length/count*(i+1));
    }
    sound1->play();
    delete sound1;
    return 0;
    
    destroyAllWindows();
    return 0;
}
