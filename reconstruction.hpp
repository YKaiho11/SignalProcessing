#pragma once
void add_wave(int selectNo, int N, int K, int I, double** H, double **U, comp*** V, Sound* sound);

void reconstruct(int N,int K,int I,double** H,double** U,comp*** V,Sound* base_sound) {
	int n, k, i;
	Sound* sound;
	sound = new Sound();
	sound->New(base_sound->samplingFrequency, base_sound->length);

	char c;
	int selectNo;
	bool finish = false;

	while (!finish) {
		while (true) {
			printf("select your number, q to quit\n");
			c = waitKey();
			selectNo = c - '0';//convert char to int

			if (c == 'q') {
				finish = true;
				break;
			}

			if (!(0 <= selectNo && selectNo < K)) {
				printf("you can select %d to %d\n", 0, K);
				continue;
			}
			else break;
		}

		if (finish) {
			printf("finish synthesizing sound\n");
			break;
		}

		printf("you select %d\n", selectNo);
		printf("calculating...\n");

		add_wave(selectNo, N, K, I, H, U, V, sound);
		printf("IFFT done!\n");
	}

	for (i = 0; i < I; i++) {
		sound->draw_wave(sound->length*i/I, sound->length*(i+1) / I);
	}
	waitKey();
}

void add_wave(int selectNo, int N, int K, int I, double** H, double **U, comp*** V, Sound* sound) {
	int n, i;
	Sound* tmp_sound;
	tmp_sound = new Sound();
	tmp_sound->New(sound->samplingFrequency, sound->length);

	for (i = 0; i < I; i++) {
		comp* G;
		G = new comp[N];
		for (n = 0; n < N; n++) {
			comp c = H[n][selectNo] * U[selectNo][i] * V[n][selectNo][i];
			G[n] = c;
		}

		IFFT(tmp_sound, tmp_sound->length*i / I, tmp_sound->length*(i + 1) / I, N, G);
		delete[] G;
	}

	for (i = 0; i < tmp_sound->length*tmp_sound->samplingFrequency; i++) {
		sound->waveData[i] += tmp_sound->waveData[i];
	}
	
	delete tmp_sound;
}