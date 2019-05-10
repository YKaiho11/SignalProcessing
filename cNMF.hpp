#pragma once
#include "NMF.hpp"
//References

//(1)COMPLEX NMF:A NEW SPARSE REPRESENTATION FOR ACOUSTIC SIGNALS
//Hirokazu Kameoka et al.

//(2)New Methods of Complex Matrix Factorization for Single-Channel Source Separation and Analysis
//Brian John King

//Fx,t=sum Hk,x*Uk,t*Vk,x,t
//this x,k,t <-> my n,k,i

void cNMF(Sound* sound, double startTime, double endTime, int I, int K) {
	//Yni~Fni=sum(l)  Hnk*Uki*exp(j*Vnki)
	//n,k,i->k,l,m
	int n, k, i;

	int N = 2048;
	double lambda = 1;//0<ƒÉ
	double p = 1;//0<p<=2

	comp* G;
	comp **Y;//N*I
	comp **F;//N*I
	double **H;//N*K
	double **U;//K*I
	comp ***V;//N*K*I
	double **V_;//K*I
	double ***B;//N*K*I
	comp ***X;//N*K*I

	//allocate memory
	G = new comp[N];

	Y = new comp*[N];
	for (n = 0; n < N; n++) Y[n] = new comp[I];

	F = new comp*[N];
	for (n = 0; n < N; n++) F[n] = new comp[I];

	H = new double*[N];
	for (n = 0; n < N; n++) H[n] = new double[K];

	U = new double*[K];
	for (k = 0; k < K; k++) U[k] = new double[I];

	V = new comp**[N];
	for (n = 0; n < N; n++) {
		V[n] = new comp*[K];
		for (k = 0; k < K; k++) V[n][k] = new comp[I];
	}

	V_ = new double*[K];
	for (k = 0; k < K; k++) {
		V_[k] = new double[I];
	}

	B = new double**[N];
	for (n = 0; n < N; n++) {
		B[n] = new double*[K];
		for (k = 0; k < K; k++) B[n][k] = new double[I];
	}

	X = new comp**[N];
	for (n = 0; n < N; n++) {
		X[n] = new comp*[K];
		for (k = 0; k < K; k++) X[n][k] = new comp[I];
	}
    
    //if failed to allocate memory
    if(Y==nullptr || F==nullptr || H==nullptr || U==nullptr || V==nullptr || V_==nullptr || B==nullptr || X==nullptr){
        printf("failed to allocate memory\n");
        return;
    }
    
    for(n=0;n<N;n++){
        if(Y[n]==nullptr || F[n]==nullptr || H[n]==nullptr || V[n]==nullptr || X[n]==nullptr){
            printf("failed to allocate memory\n");
            return;
        }
        for(k=0;k<K;k++){
            if(V[n][k]==nullptr || B[n][k]==nullptr || X[n][k]==nullptr){
                printf("failed to allocate memory\n");
                return;
            }
        }
    }
    
    for(k=0;k<K;k++){
        if(U[k]==nullptr || V_[k]==nullptr){
            printf("failed to allocate memory\n");
            return;
        }
    }

	//generate Y
	for (i = 0; i < I; i++) {
		FFT(sound, (endTime - startTime)*i / I + startTime, (endTime - startTime)*(i + 1) / I + startTime, N, G, false);
		double max_s = 0;
		for (n = 0; n < N; n++) {
			if (max_s < abs(G[n])) max_s = pow(abs(G[n]), 1);
		}

		if (max_s == 0) {
			for (n = 0; n < N; n++) {
				Y[n][i].re = 0.00000000000000001;
                Y[n][i].im = 0.00000000000000001;
			}
		}
		else {
			for (n = 0; n < N; n++) {
				Y[n][i].re = G[n].re / max_s;
				Y[n][i].im = G[n].im / max_s;
			}
		}
	}

	//initialize H,U,V,V_,B
	for (n = 0; n < N; n++) {
		for (k = 0; k < K; k++) {
            H[n][k] = 1.0*rand() / RAND_MAX + 0.5;
		}
	}
	for (k = 0; k < K; k++) {
		for (i = 0; i < I; i++) {
            U[k][i] = 1.0*rand() / RAND_MAX + 0.5;
            V_[k][i] = 1.0*rand() / RAND_MAX + 0.5;
		}
	}
	for (n = 0; n < N; n++) {
		for (k = 0; k < K; k++) {
			for (i = 0; i < I; i++) {
				//if successful, below should be deleted
				/**V[n][k][i].re = 2 * (double)rand() / RAND_MAX - 1;
				V[n][k][i].im = 2 * (double)rand() / RAND_MAX - 1;
				double normalizing = abs(V[n][k][i]);
				V[n][k][i].re /= normalizing;
				V[n][k][i].im /= normalizing;**/

				comp c = Y[n][i] / abs(Y[n][i]);
				V[n][k][i] = c;

				do {
					B[n][k][i] = (double)rand() / RAND_MAX;
				} while (B[n][k][i] == 0);
			}
		}
	}

	//set sum of B as 1
	for (n = 0; n < N; n++) {
		for (i = 0; i < I; i++) {
			double sum = 0;
			for (k = 0; k < K; k++) {
				sum += B[n][k][i];
			}
			for (k = 0; k < K; k++) {
				B[n][k][i] /= sum;
			}
		}
	}

	//pre calculation using nonnegatice matrices factorization (NMF)
	double **Y_ = new double*[N];
	for (n = 0; n < N; n++) Y_[n] = new double[I];
	for (n = 0; n < N; n++)
		for (i = 0; i < I; i++)
			Y_[n][i] = abs(Y[n][i]);
	NMF_calc(Y_, H, U, N, K, I);
	for (n = 0; n < N; n++) delete[] Y_[n];
	delete[] Y_;


	/********** main culculation **********/
	double distance = 0, distance_pre;
	static bool isFirst = true;
	while (true) {

		//make F
		for (n = 0; n < N; n++) {
			for (i = 0; i < I; i++) {
				F[n][i].re = 0; F[n][i].im = 0;
				for (k = 0; k < K; k++) {
					comp c = F[n][i] + H[n][k] * U[k][i] * V[n][k][i];
					F[n][i] = c;
				}
			}
		}


		//calculate distance
		distance_pre = distance;
		distance = 0;
		for (n = 0; n < N; n++) {
			for (i = 0; i < I; i++) {
				distance += pow(abs(Y[n][i] - F[n][i]), 2);
			}
		}
		for (k = 0; k < K; k++) {
			for (i = 0; i < I; i++) {
				distance += 2 * lambda*pow(U[k][i], p);
			}
		}

		//change of distance
		if (isFirst) isFirst = false;
		else {
			//printf("%f,%f\n", distance, distance_pre);
			printf("calculating CMF... d=%f    distance=%f\n", abs(distance_pre - distance),distance);
			if (abs(distance_pre - distance) < 0.001) break;/*********************************************  SET THRESHOLD **************************/
		}

		//update B
		for (n = 0; n < N; n++) {
			for (i = 0; i < I; i++) {
				double sum = 0;
				for (k = 0; k < K; k++) {
					sum += H[n][k] * U[k][i];
				}
				if (sum == 0) {
					for (k = 0; k < K; k++) {
						B[n][k][i] = 1.0 / K;
					}
				}
				else {
					for (k = 0; k < K; k++) {
						B[n][k][i] = H[n][k] * U[k][i] / sum;
					}
				}
			}
		}

		//update X
		for (n = 0; n < N; n++) {
			for (k = 0; k < K; k++) {
				for (i = 0; i < I; i++) {
					comp c2 = H[n][k] * U[k][i] * V[n][k][i] + B[n][k][i] * (Y[n][i] - F[n][i]);
					X[n][k][i] = c2;
				}
			}
		}

		//update V
		for (n = 0; n < N; n++) {
			for (k = 0; k < K; k++) {
				for (i = 0; i < I; i++) {
					if (abs(X[n][k][i]) == 0) {
						V[n][k][i].re = 0;
						V[n][k][i].im = 0;
					}
					else {
						comp c = X[n][k][i] / abs(X[n][k][i]);
						V[n][k][i] = c;
					}
				}
			}
		}

		//update H
		/**/
		for (n = 0; n < N; n++) {
			double sum = 0;
			for (k = 0; k < K; k++) {
				double sum1 = 0;
				double sum2 = 0;
				for (i = 0; i < I; i++) {
					if (B[n][k][i] != 0) {
						sum1 += U[k][i] * abs(X[n][k][i]) / B[n][k][i];
						sum2 += U[k][i] * U[k][i] / B[n][k][i];
					}
				}
				if (sum2 == 0) {
					H[n][k] = 0;
				}
				else {
					H[n][k] = sum1 / sum2;
				}
				sum += H[n][k];
			}

			if (sum != 0) {
				for (k = 0; k < K; k++) {
					H[n][k] /= sum;
				}
			}
		}
		/**/

		//update U
		/**/
		for (k = 0; k < K; k++) {
			for (i = 0; i < I; i++) {
				double sum1 = 0;
				double sum2 = 0;
				for (n = 0; n < N; n++) {
					if (B[n][k][i] != 0) {
						sum1 += H[n][k] * abs(X[n][k][i]) / B[n][k][i];
						sum2 += H[n][k] * H[n][k] / B[n][k][i];
					}
				}
				U[k][i] = sum1 / (sum2 + lambda * p*pow(U[k][i], p - 2));//U<->V_
			}
		}
		/**/
	}
	/********************/



	int *decreaseIndex;
	decreaseIndex = new int[K];
	for (k = 1; k < K; k++) decreaseIndex[k] = 0;
	decreaseIndex[0] = -1;
	int decreaseNo = decrease_base(H, N, K, decreaseIndex, 0.11);
	printf("decrease = %d\n", decreaseNo);
	merge_base(H, U, N, K, I, decreaseIndex, decreaseNo);
	printf("decrease + merge = %d\n", decreaseNo);
    printf("remain = %d\n",K-decreaseNo);
	if (decreaseNo == K) {
		printf("No Base!\n");
		waitKey();
		return;
	}

	const int width = 600;
	const int height = 600;
	//show result
	int base_width = width;
	int base_height = height / K;
	Mat base(Size(base_width, (K - decreaseNo) * base_height), CV_8UC3, Scalar::all(0));

	int index = 0;
	for (k = 0; k < K; k++) {//V[n][k]
		if (decreaseIndex[index] == k) {
			index++;
			continue;
		}
		printf("k=%d\n", k);
		printf("index=%d\n", index);

		double max_v = 0;
		for (n = 0; n < N / 2; n++) {
			if (max_v < H[n][k]) max_v = H[n][k];
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
					average += H[j][k];
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

	index = 0;
	int width_feature = width / I;
	Mat feature(Size(I * width_feature, (K - decreaseNo) * base_height), CV_8UC3, Scalar::all(0));
	for (k = 0; k < K; k++) {//T[K][I]
		if (decreaseIndex[index] == k) {
			index++;
			continue;
		}
		double max_t = 0;
		for (i = 0; i < I; i++) {
			if (max_t < U[k][i]) max_t = U[k][i];
		}

		for (i = 0; i < I; i++) {
			int color[3] = { (int)(255 - U[k][i] / max_t * 255),0,(int)(U[k][i] / max_t * 255) };
			for (int ii = i * width_feature; ii < (i + 1)*width_feature; ii++) {
				for (int jj = (k - index) * base_height; jj < (k - index + 1)*base_height; jj++) {
					feature.data[3 * (jj*feature.cols + ii)] = color[0];
					feature.data[3 * (jj*feature.cols + ii) + 1] = color[1];
					feature.data[3 * (jj*feature.cols + ii) + 2] = color[2];
				}
			}
		}
	}
	imshow("feature", feature);
	waitKey();


	//free memory
	for (n = 0; n < N; n++) {
		delete[] Y[n];
		delete[] F[n];
		delete[] H[n];
		for (k = 0; k < K; k++) {
			delete[] V[n][k];
			delete[] B[n][k];
			delete[] X[n][k];
		}
		delete[] V[n];
		delete[] B[n];
		delete[] X[n];
	}
	for (k = 0; k < K; k++) {
		delete[] U[k];
		delete[] V_[k];
	}
	delete[] G;
	delete[] Y;
	delete[] F;
	delete[] H;
	delete[] U;
	delete[] V;
	delete[] V_;
	delete[] B;
	delete[] X;
}
