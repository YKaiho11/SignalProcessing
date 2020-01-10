#include "pch.h"
#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS 
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#pragma comment(lib,"opencv_world342.lib")
#include <al.h>
#include <alc.h>
#pragma comment(lib, "OpenAL32.lib")
using namespace std;
using namespace cv;

#define SAMPLINGFREQUENCY 44100

#include "wav_inoutput.hpp"
#include "sound_class.hpp"
#include "sound_process.hpp"
#include "Fourier.hpp"
#include "tone_list.hpp"
#include "NMF.hpp"
#include "cNMF.hpp"
//todo test cut_first in "sound_process.hpp"

int main() {
	srand((unsigned int)time(NULL));

	/**Sound* sound1;
	sound1 = new Sound();
	sound1->New(SAMPLINGFREQUENCY, 2);
	sound1->sine_curve(440, 0, 1);

	Sound* sound2;
	sound2 = new Sound();
	sound2->New(SAMPLINGFREQUENCY, 2);
	sound2->sine_curve(1000, 0.5, 1.5);

	Sound* sound3;
	sound3 = new Sound();
	sound3->New(SAMPLINGFREQUENCY, 2);
	sound3->sine_curve(675, 1, 2);

	Sound* sound4;
	sound4 = new Sound();
	sound4->New(SAMPLINGFREQUENCY, 2);
	synthesize(sound1, sound4, sound4, 0, 1);
	sound4->amplify(0, 1, 2);
	synthesize(sound2, sound4, sound4, 0.5, 1.5);
	sound4->amplify(1, 1.5, 2);
	synthesize(sound3, sound4, sound4, 1, 2);
	sound4->amplify(1.5, 2, 2);
	sound4->smoothing(0, 2, 5);

	cNMF(sound1, 0, 2, 100, 10);

	delete sound1; delete sound2; delete sound3; delete sound4;**/

	/**/Sound* sound5;
	sound5 = new Sound();
	sound5->input((char*)"../piano/C3.wav");
	Sound* sound6;
	sound6 = new Sound();
	sound6->input((char*)"../piano/E3.wav");
	Sound* sound7;
	sound7 = new Sound();
	sound7->input((char*)"../piano/G3.wav");
	synthesize(sound6, sound7, sound7, 0, 10);
	synthesize(sound5, sound7, sound7, 0, 10);
	cNMF(sound7, 0, sound7->length, 100, 10);

	delete sound5; delete sound6; delete sound7;/**/
	destroyAllWindows();
	return 0;
}
