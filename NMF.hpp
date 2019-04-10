#pragma once
const double interval = 0.1;
#define EPSILON 0.00001

double Euclid(double x, double x_);
double KL(double x, double x_);
double IS(double x, double x_);

void NMF(Sound* sound, double startTime, double endTime, const int I) {
	int i, n, k;
	//X=TV
	//N*I=N*K K*I
	//k:base number,i:length
	const int N = 2048;
	const int K = 10;
	const int option = 2;//0,1,2
	comp G[N];

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
		FFT(sound, (endTime - startTime)*i / I + startTime, (endTime - startTime)*(i + 1) / I + startTime, N, G, false);
		double max_s = 0;
		for (n = 0; n < N; n++) {
			if (max_s < pow(abs(G[n]), 2))max_s = pow(abs(G[n]), 2);
		}

		for (n = 0; n < N; n++) {
			X[n][i] = pow(abs(G[n]), 2) / max_s;
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
			//printf("d=%f\n", abs(distance_pre - distance));
			if (abs(distance_pre - distance) < EPSILON) break;
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

	//show result
	int base_width = 800;
	int base_height = 10;
	Mat base(Size(base_width, K * base_height), CV_8UC3, Scalar::all(0));
	Mat feature(Size(I * 10, K * 10), CV_8UC3, Scalar::all(0));

	for (k = 0; k < K; k++) {//V[n][k]
		double max_v = 0;
		for (n = 0; n < N; n++) {
			if (max_v < V[n][k]) max_v = V[n][k];
		}

		n = 0;
		for (i = 0; i < base_width; i++) {
			int n0 = n;
			while (n < 1.0*i / base_width * N) {
				n++;
			}
			int n1 = n;

			if (n1 - n0 != 0) {
				double average = 0;
				for (int j = n0; j < n1; j++) {
					average += V[j][k];
				}
				average /= n1 - n0;
				int color[3] = { 255 - average / max_v * 255,0,average / max_v * 255 };
				for (int ii = i; ii < i + 1; ii++) {
					for (int jj = base_height * k; jj < base_height * (k + 1); jj++) {
						base.data[3 * (jj*base.cols + ii)] = color[0];
						base.data[3 * (jj*base.cols + ii)+1] = color[1];
						base.data[3 * (jj*base.cols + ii)+2] = color[2];
					}
				}
			}
			else {

			}
		}
	}
	imshow("base", base);
	waitKey();

	/*****free memory ****/
	for (i = 0; i < N; i++) {
		delete X[i];
		delete X_[i];
		delete V[i];
	}
	for (i = 0; i < K; i++) {
		delete T[i];
	}
	delete X;
	delete X_;
	delete V;
	delete T;
}


double Euclid(double x, double x_) {
	return pow(abs(x - x_), 2);
}
double KL(double x, double x_) {
	return x*log(x / x_) - x + x_;
}
double IS(double x, double x_) {
	return x / x_ - log(x / x_) - 1;
}