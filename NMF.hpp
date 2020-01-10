#pragma once
const double interval = 0.1;

double Euclid(double x, double x_);
double KL(double x, double x_);
double IS(double x, double x_);

int decrease_base(double **V, int N, int K, int *index, double threshold);
void merge_base(double **V, double **T, int N, int K, int I, int *index, int &decreaseNo);
void addList(int val, int *index, int decreaseNo);
void NMF_calc(double **X, double **V, double **T, int N, int K, int I);

void show_base(const int width, const int height, int K, int N, int decreaseNo, int decreaseIndex[], double **V,String win_name);
void show_base_log(const int width, const int height, int K, int N, int decreaseNo, int decreaseIndex[], double **V);
void show_weight(const int width, const int height, int K, int I, int decreaseNo, int decreaseIndex[], double **T, String win_name);

void NMF(Sound* sound, double startTime, double endTime, int I, int K) {
	int i, n, k;
	//X=V*T
	//N*I=N*K K*I
	//k:base number,i:length
	int N = 2048;
	const int option = 2;//0,1,2
	comp* G;
	G = new comp[N];

	double** X; //N*I
	double** X_;//N*I
	double** V; //N*K
	double** T; //K*I

	/**** allocate memory****/
	X = new double*[N];
	X_ = new double*[N];
	V = new double*[N];
	T = new double*[K];
	if (X == nullptr || X_ == nullptr || V == nullptr || T == nullptr) {
		printf("failed to allocate memory\n");
		return;
	}
	for (i = 0; i < N; i++) {
		X[i] = new double[I];
		X_[i] = new double[I];
		V[i] = new double[K];
		if (X[i] == nullptr || X_[i] == nullptr || V[i] == nullptr) {
			printf("failed to allocate memory\n");
			return;
		}
	}
	for (i = 0; i < K; i++) {
		T[i] = new double[I];
		if (T[i] == nullptr) {
			printf("failed to allocate memory\n");
			return;
		}
	}

	//generate X
	for (i = 0; i < I; i++) {
		FFT(sound, (endTime - startTime)*i / I + startTime, (endTime - startTime)*(i + 1) / I + startTime, N, G, false, 0);
		double max_s = 0;
		for (n = 0; n < N; n++) {
			if (max_s < pow(abs(G[n]), 1))max_s = pow(abs(G[n]), 1);
		}

		if (max_s == 0) {
			for (n = 0; n < N; n++) {
				X[n][i] = 0.00000000000000001;
			}
		}
		else {
			for (n = 0; n < N; n++) {
				X[n][i] = pow(abs(G[n]), 1) / max_s;
			}
		}
	}

	/***********initialize*********/
	for (n = 0; n < N; n++) {
		for (k = 0; k < K; k++) {
			V[n][k] = 1.0*rand() / RAND_MAX;
		}
	}

	for (k = 0; k < K; k++) {
		for (i = 0; i < I; i++) {
			T[k][i] = 1.0*rand() / RAND_MAX;
		}
	}

	//main calculaiton
	NMF_calc(X, V, T, N, K, I);

	int *decreaseIndex;
	decreaseIndex = new int[K];
	for (k = 1; k < K; k++) decreaseIndex[k] = 0;
	decreaseIndex[0] = -1;
	int decreaseNo = decrease_base(V, N, K, decreaseIndex, 0.07);
	merge_base(V, T, N, K, I, decreaseIndex, decreaseNo);

	const int width = 600;
	const int height = 600;

	//show result
	show_base(width, height, K, N, decreaseNo, decreaseIndex, V,"NMF base");
	show_weight(width, height, K, I, decreaseNo, decreaseIndex, T,"NMF weight");

	/*****free memory ****/
	for (i = 0; i < N; i++) {
		delete[] X[i];
		delete[] X_[i];
		delete[] V[i];
	}
	for (i = 0; i < K; i++) {
		delete[] T[i];
	}
	delete[] X;
	delete[] X_;
	delete[] V;
	delete[] T;
	delete[] G;
	delete[] decreaseIndex;
}

void NMF_calc(double** X, double** V, double** T, int N, int K, int I) {
	int n, k, i;
	double **X_ = new double*[N];//N*I
	for (n = 0; n < N; n++) X_[n] = new double[I];
	int option = 2;

	double distance = 0, distance_pre;
	static bool isFirst = true;
	while (true) {
		//multiplication of V and T
		for (n = 0; n < N; n++) {
			for (i = 0; i < I; i++) {
				X_[n][i] = 0;
				for (k = 0; k < K; k++) {
					X_[n][i] += V[n][k] * T[k][i];
				}
			}
		}

		//calculate distance
		distance_pre = distance;
		distance = 0;
		for (n = 0; n < N; n++) {
			for (i = 0; i < I; i++) {
				switch (option) {
				case 0:
					distance += Euclid(X[n][i], X_[n][i]);
					break;
				case 1:
					distance += KL(X[n][i], X_[n][i]);
					break;
				case 2:
					distance += IS(X[n][i], X_[n][i]);
					break;
				}
			}
		}

		//change of distance
		if (isFirst) isFirst = false;
		else {
			//printf("%f,%f\n", distance, distance_pre);
			printf("calculating NMF... d=%f    distance=%f\n", abs(distance_pre - distance),distance);
			if (abs(distance_pre - distance) < 0.03) break;/************************************************************SET THRESHOLD*****/
		}

		//update V
		double sum1;
		double sum2;
		for (n = 0; n < N; n++) {
			for (k = 0; k < K; k++) {
				sum1 = 0;
				sum2 = 0;
				for (int ii = 0; ii < I; ii++) {
					switch (option) {
					case 0:
						sum1 += X[n][ii] * T[k][ii];
						sum2 += X_[n][ii] * T[k][ii];
						break;
					case 1:
						sum1 += X[n][ii] / X_[n][ii] * T[k][ii];
						sum2 += T[k][ii];
						break;
					case 2:
						sum1 += X[n][ii] * T[k][ii] / (X_[n][ii] * X_[n][ii]);
						sum2 += T[k][ii] / X_[n][ii];
						break;
					}
				}
				switch (option) {
				case 0:
				case 1:
					V[n][k] *= sum1 / sum2;
					break;
				case 2:
					V[n][k] *= sqrt(sum1 / sum2);
					break;
				}
			}
		}

		//update T
		for (k = 0; k < K; k++) {
			for (i = 0; i < I; i++) {
				sum1 = 0;
				sum2 = 0;
				for (int j = 0; j < N; j++) {
					switch (option) {
					case 0:
						sum1 += X[j][i] * V[j][k];
						sum2 += X_[j][i] * V[j][k];
						break;
					case 1:
						sum1 += X[j][i] / X_[j][i] * V[j][k];
						sum2 += V[j][k];
						break;
					case 2:
						sum1 += X[j][i] * V[j][k] / (X_[j][i] * X_[j][i]);
						sum2 += V[j][k] / X_[j][i];
						break;
					}
				}
				switch (option) {
				case 0:
				case 1:
					T[k][i] *= sum1 / sum2;
					break;
				case 2:
					T[k][i] *= sqrt(sum1 / sum2);
					break;
				}
			}
		}
	}

	for (n = 0; n < N; n++) delete[] X_[n];
	delete[] X_;
}

void addList(int val,int *index,int decreaseNo){
	if (decreaseNo == 0) {
		index[0] = val;
		return;
	}

	int i, j;
	if (val < index[0]) {
		i = 0;
	}
	else {
		for (i = 1; i < decreaseNo; i++) {
			if (index[i-1] < val && val < index[i]) break;
		}
	}

	for (j = decreaseNo; j > i; j--) {
		index[j] = index[j - 1];
	}
	index[i] = val;
}

void merge_base(double **V, double **T, int N, int K, int I, int *index, int &decreaseNo) {
	int *peak;
	peak = new int[K];
	int n, k, k_, i;
	double max;
	//find each peak
	for (k = 0; k < K; k++) {
		bool judge = true;
		for (i = 0; i < decreaseNo; i++) {
			if (index[i] == k) {
				judge = false;
				break;
			}
		}

		if (!judge) {
			peak[k] = -1;
		}
		else {
			peak[k] = 0;
			max = V[0][k];
			for (n = 1; n < N / 2; n++) {
				if (V[n][k] > max) {
					max = V[n][k];
					peak[k] = n;
				}
			}
		}
	}

	for (k = 0; k < K; k++) {
		for (k_ = k + 1; k_ < K; k_++) {
			bool judge = true;
			for (i = 0; i < decreaseNo; i++) {
				if (index[i] == k || index[i] == k_) {
					judge = false;
					break;
				}
			}

			if (judge) {
				if (peak[k] == peak[k_]) {
					addList(k_, index, decreaseNo);
					decreaseNo++;

					for (i = 0; i < I; i++) {
						T[k][i] += T[k_][i];
					}
				}
			}
		}
	}


	delete[] peak;
}


int decrease_base(double **V,int N,int K,int *index, double threshold) {
	double sum_effective;
	double max_v;
	int n, k, decreaseNo;

	decreaseNo = 0;
	for (k = 0; k < K; k++) {
		sum_effective = 0;
		max_v = 0;
		for (n = 0; n < N; n++) {
			if (max_v < V[n][k]) max_v = V[n][k];
			sum_effective += V[n][k];
		}

		if (sum_effective / (max_v*N) > threshold || max_v==0) {
			index[decreaseNo] = k;
			decreaseNo++;
		}
	}
	return decreaseNo;
}


double Euclid(double x, double x_) {
	return pow(abs(x - x_), 2);
}
double KL(double x, double x_) {
	return x*log(x / x_) - x + x_;
}
double IS(double x, double x_) {
    if(x==0 || x_==0) return 0;
	return x / x_ - log(x / x_) - 1;
}

void show_base(const int width,const int height,int K,int N,int decreaseNo,int decreaseIndex[],double **V,String win_name) {
	//todo: draw histogram as log
	int k, n, i;
	int base_width = width;
	int base_height = height / K;
	Mat base(Size(base_width, (K - decreaseNo) * base_height), CV_8UC3, Scalar::all(0));

	int index = 0;
	for (k = 0; k < K; k++) {//V[n][k]
		if (decreaseIndex[index] == k) {
			index++;
			continue;
		}

		double max_v = 0;
		for (n = 0; n < N / 2; n++) {
			if (max_v < V[n][k]) max_v = V[n][k];
		}

		n = 0;
		for (i = 0; i < base_width; i++) {
			int n0 = n;
			while (n < 1.0*i / base_width * N / 2) {
				n++;
			}
			int n1 = n;

			if (n1 - n0 != 0) {
				double average = 0;
				for (int j = n0; j < n1; j++) {
					average += V[j][k];
				}
				average /= n1 - n0;
				int color[3] = { (int)(255 - average / max_v * 255),0,(int)(average / max_v * 255) };
				for (int ii = i; ii < i + 1; ii++) {
					for (int jj = base_height * (k - index); jj < base_height * (k - index + 1); jj++) {
						base.data[3 * (jj*base.cols + ii)] = color[0];
						base.data[3 * (jj*base.cols + ii) + 1] = color[1];
						base.data[3 * (jj*base.cols + ii) + 2] = color[2];
					}
				}
			}
			else {

			}
		}
	}
	imshow(win_name, base);
	waitKey(1);
}


void show_base_log(const int width, const int height, int K, int N, int decreaseNo, int decreaseIndex[], double **V) {
	//todo: complete this function
	int k, n, i;
	int base_width = width;
	int base_height = height / K;
	Mat base(Size(base_width, (K - decreaseNo) * base_height), CV_8UC3, Scalar::all(0));

	int index = 0;
	for (k = 0; k < K; k++) {//V[n][k]
		if (decreaseIndex[index] == k) {
			index++;
			continue;
		}

		double max_v = 0;
		for (n = 0; n < N / 2; n++) {
			if (max_v < V[n][k]) max_v = V[n][k];
		}

		n = 0;
		for (i = 0; i < base_width; i++) {
			int n0 = n;
			while (n < 1.0*i / base_width * N / 2) {
				n++;
			}
			int n1 = n;

			if (n1 - n0 != 0) {
				double average = 0;
				for (int j = n0; j < n1; j++) {
					average += V[j][k];
				}
				average /= n1 - n0;
				int color[3] = { (int)(255 - average / max_v * 255),0,(int)(average / max_v * 255) };
				for (int ii = i; ii < i + 1; ii++) {
					for (int jj = base_height * (k - index); jj < base_height * (k - index + 1); jj++) {
						base.data[3 * (jj*base.cols + ii)] = color[0];
						base.data[3 * (jj*base.cols + ii) + 1] = color[1];
						base.data[3 * (jj*base.cols + ii) + 2] = color[2];
					}
				}
			}
			else {

			}
		}
	}
	imshow("base", base);
	waitKey(1);
}

void show_weight(const int width,const int height,int K,int I,int decreaseNo,int decreaseIndex[],double **T,String win_name) {
	int k, i;
	int index = 0;
	int base_height = height / K;
	int width_feature = width / I;
	Mat feature(Size(I * width_feature, (K - decreaseNo) * base_height), CV_8UC3, Scalar::all(0));
	for (k = 0; k < K; k++) {//T[K][I]
		if (decreaseIndex[index] == k) {
			index++;
			continue;
		}
		double max_t = 0;
		for (i = 0; i < I; i++) {
			if (max_t < T[k][i]) max_t = T[k][i];
		}

		for (i = 0; i < I; i++) {
			int color[3] = { (int)(255 - T[k][i] / max_t * 255),0,(int)(T[k][i] / max_t * 255) };
			for (int ii = i * width_feature; ii < (i + 1)*width_feature; ii++) {
				for (int jj = (k - index) * base_height; jj < (k - index + 1)*base_height; jj++) {
					feature.data[3 * (jj*feature.cols + ii)] = color[0];
					feature.data[3 * (jj*feature.cols + ii) + 1] = color[1];
					feature.data[3 * (jj*feature.cols + ii) + 2] = color[2];
				}
			}
		}
	}
	imshow(win_name, feature);
	waitKey(1);
}