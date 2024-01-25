#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <fstream>
#include <Eigen>
#include <Eigen/FFT>


//int main() {
//    //Eigenを用いたFFT
//    double dt = 0.01;           //サンプリング周期
//    size_t n = 256;             //データ数
//    double df = 1 / (dt * n);   //周波数周期
//    double f0 = 10, f1 = 20;    //テスト用周波数設定
//
//    Eigen::VectorXd t(n), y(n), y_ifft(n), amp(n), f(n);//EigenClassベクトルの作成
//    Eigen::VectorXcd y_fft(n);//EigenClass複素ベクトルの作成
//
//    //仮データ作成
//    for (int i = 0; i < n; i++) {
//        f(i) = df * i;
//        t(i) = dt * i;
//        y(i) = 2 * std::sin(2 * M_PI * f0 * t(i)) + std::sin(2 * M_PI * f1 * t(i));
//    }
//
//    //FFT
//    Eigen::FFT<double> fft;
//    fft.fwd(y_fft, y);
//    amp = (y_fft / (n / 2)).cwiseAbs();
//
//    //CSV書き込み処理
//    std::ofstream dataFile("../data.csv");
//    for (int i = 0; i < n/2; i++) {
//        //横軸周波数、縦軸振幅
//        dataFile << f[i] << "," << amp[i] << std::endl;
//        std::cout << f[i] << "," << amp[i] << std::endl;
//    }
//
//    return 0;
//}


std::vector<std::complex<double>> FFT(std::vector<std::complex<double>> X);

//一度に2^N個のデータを読み込んで表示
int main() {
    int data_N = 256;
    double dt = 0.01;
    double df = 1 / (dt * data_N);   //周波数周期
    std::vector<double> time(data_N), freq(data_N);
    std::vector<std::complex<double>> val(data_N);
    std::vector<std::complex<double>> val_fft(data_N);
    std::vector<double> val_fft_real(data_N);

    double f1 = 10;
    double f2 = 20;


    for (int i = 0; i < data_N; i++) {
        time[i] = dt * i;
        freq[i] = df * i;
        val[i] = std::sin(2 * M_PI * f1 * time[i]) + std::sin(2 * M_PI * f2 * time[i]);
    }

    val_fft = FFT(val);

    for (int i = 0; i < data_N; i++) {
        val_fft_real[i] = std::abs(val_fft[i]);
    }

    std::ofstream dataFile("../data.csv");
    for (int i = 0; i < data_N / 2; i++) {
        dataFile << freq[i] << "," << val_fft_real[i] << std::endl;
        std::cout << freq[i] << "," << val_fft_real[i] << std::endl;
    }

    return 0;
}

std::vector<std::complex<double>> FFT(std::vector<std::complex<double>> X) {
    const int N = X.size();
    if (N == 1) return X;
    std::vector<std::complex<double>> even(N / 2), odd(N / 2);
    for (int i = 0; i < N / 2; i++) {
        even[i] = X[i * 2];
        odd[i] = X[i * 2 + 1];
    }
    even = FFT(even);
    odd = FFT(odd);
    for (int i = 0; i < N / 2; i++) {
        odd[i] *= std::polar(1.0, 2 * M_PI * i / N);
        X[i] = even[i] + odd[i];
        X[N / 2 + i] = even[i] - odd[i];
    }
    return X;
}
