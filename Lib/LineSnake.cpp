#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "ImgVec.hpp"
#define PI 3.14159265
using namespace std;

vector<double> calculate_Eint1(vector<int> &y)
{
	int length = y.size();

	double d = 0;
	for (int i = 0; i < length; i++)
	{
		d += i == 0 ? abs(y.back() - y[i]) : abs(y[i - 1] - y[i]);
	}
	d /= length;

	vector<double> Eint1(length);
	double _max;
	double denominator, numerator;
	for (int i = 0; i < length; i++)
	{
		_max = 0;
		for (int j = -1; j < 2; j++)
		{
			denominator = (i == length - 1) ? pow(d - abs(y[i] + j - y[0]), 2) : pow(d - abs(y[i] + j - y[i + 1]), 2);
			_max = max(_max, denominator);
		}
		Eint1[i] = (i == length - 1) ? pow(d - abs(y[i] - y[0]), 2) / _max : pow(d - abs(y[i] - y[i + 1]), 2) / _max;
	}
	return Eint1;
}

vector<double> calculate_Eint2(vector<int> &y)
{
	int length = y.size();
	double _max, temp;
	vector<double> Eint2(length);

	for (int i = 0; i < length; i++)
	{
		_max = 0;
		if (i == 0)
		{
			for (int j = -1; j < 2; j++)
			{
				temp = pow(y.back() - 2 * (y[i] + j) + y[i + 1], 2);
				_max = max(temp, _max);
			}
			Eint2[i] = pow(y.back() - 2 * y[i] + y[i + 1], 2) / _max;
		}
		else if (i == length - 1)
		{
			for (int j = -1; j < 2; j++)
			{
				temp = pow(y[i - 1] - 2 * (y[i] + j) + y[0], 2);
				_max = max(temp, _max);
			}
			Eint2[i] = pow(y[i - 1] - 2 * y[i] + y[0], 2) / _max;
		}
		else
		{
			for (int j = -1; j < 2; j++)
			{
				temp = pow(y[i - 1] - 2 * (y[i] + j) + y[i + 1], 2);
				_max = max(temp, _max);
			}
			Eint2[i] = pow(y[i - 1] - 2 * y[i] + y[i + 1], 2) / _max;
		}
	}
	return Eint2;
}

vector<double> calculate_Eimg(MyImage<double> &I, MyImage<double> &Ibg, vector<int> &x, vector<int> &y)
{
	int length = y.size();

	vector<double> Eimg(length);

	vector<double> G(length);
	for (int i = 0; i < length; i++)
		G[i] = (I.getPixel(y[i] - 1, x[i]) - I.getPixel(y[i] + 1, x[i])) / 2;

	vector<double> C(length);
	for (int i = 0; i < length; i++)
		C[i] = G[i] / (Ibg.getPixel(y[i], x[i]) + 1e-8);

	// double _minC = *min_element(C.begin(), C.end());
	// double _maxC = *max_element(C.begin(), C.end());
	// for (int i = 0; i < length; i++)
	// 	Eimg[i] = (_minC - C[i]) / (_maxC - _minC);
	return C;
}

vector<double> calculate_dist(vector<int> &y)
{
	int length = y.size();
	vector<double> dist(length);
	for (int i = 1; i < length; i++)
		dist[i] = abs(y[i] - y[i-1]);
	dist[0] = dist.back();
	return dist;
}

vector<double> calculate_Econ(vector<int> &y)
{
	int length = y.size();
	double y_mean, D;
	y_mean = 0;
	for (auto num : y)
	{
		y_mean += num;
	}
	y_mean /= length;
	D = 0;
	for (auto num : y)
	{
		D += abs(num - y_mean);
	}
	D /= length;

	vector<double> Econ(length);
	double _max, temp;
	for (int i = 0; i < length; i++)
	{
		_max = 0;
		for (int j = -1; j < 2; j++)
		{
			temp = pow(D - abs(y[i] + j - y_mean), 2);
			_max = max(temp, _max);
		}
		Econ[i] = -pow(D - abs(y[i] - y_mean), 2) / _max;
	}
	return Econ;
}

double calculate_E(MyImage<double> &I, MyImage<double> &Ibg, vector<int> &x, vector<int> &y, double alpha, double beta, double gamma, double sigma, double dist)
{
	vector<double> Eint1 = calculate_Eint1(y);
	vector<double> Eint2 = calculate_Eint2(y);
	vector<double> Eimg = calculate_Eimg(I, Ibg, x, y);
	vector<double> Econ = calculate_Econ(y);
	vector<double> E_dist = calculate_dist(y);
	double E = 0;
	int length = y.size();
	for (int i = 0; i < length; i++)
		E += alpha * Eint1[i] + beta * Eint2[i] + gamma * Eimg[i] + sigma * Econ[i] + dist * E_dist[i];
	return E / length;
}

vector<int> iterate(MyImage<double> &I, MyImage<double> &Ibg,vector<int>&y, double alpha, double beta, double gamma, double sigma, double dist)
{
	vector<int> x(36);
	for (int i = 0; i < 36; i++)
		x[i] = i * 10;
	int length = x.size();
	double init_E = calculate_E(I, Ibg, x, y, alpha, beta, gamma, sigma, dist);
	double E_old, E_new;
	double minE = 0;
	double tempE;
	int minEj;
	int it = 0;
	E_old = init_E;
	E_new = E_old + 10;
	cout << "init_E:" << init_E << endl;
	while (abs(E_new - E_old) > 0.01)
	{
		it++;
		E_old = E_new;
		for (int i = 0; i < length; i++)
		{
			if (i == 0)
				minE = E_old;
			else
				minE = minE;
			minEj = 0;
			for (int j = -5; j < 21; j++)
			{
				vector<int> y_new = y;
				y_new[i] += j;
				tempE = calculate_E(I, Ibg, x, y_new, alpha, beta, gamma, sigma, dist);
				if (tempE < minE)
				{
					minE = tempE;
					minEj = j;
				}
			}
			y[i] += minEj;
			y[i] = y[i] < 23 ? 23 : (y[i] > 124 ? 124 : y[i]);
			cout << "It" << it << "\tNo." << i << "\tminE:" << minE << "\tminEj:" << minEj << endl;
		}
		E_new = minE;
	}
	return y;
}

extern "C" __declspec(dllexport) void active_contour(double *polar_img, double *gaussian_img, double *input_y, double *output_y, int rows, int cols, double alpha, double beta, double gamma, double sigma,double dist)
{
	MyImage<double> I(rows, cols);
	I.data = polar_img;
	MyImage<double> Ibg(rows, cols);
	Ibg.data = gaussian_img;
	vector<int> y_ori(36);
	for(int i=0;i<36;i++)
		y_ori[i] = input_y[i];
	vector<int> y = iterate(I, Ibg,y_ori, alpha, beta, gamma, sigma,dist);

	for (int i = 0; i < 36; i++)
		output_y[i] = y[i];
	I.data = nullptr;
	Ibg.data = nullptr;
}