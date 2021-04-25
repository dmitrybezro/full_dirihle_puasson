#pragma once
#include "AllMethods.h"
#include "MethodOutside.h"

#ifndef MSG_H
#define MSG_H

using namespace std;

class MSG :
	public MethodOutside
{
private:
	TMatrix<double> H;
	double a, bet, Ahh;
public:
	MSG() : MethodOutside() {}
	MSG(int N, int M, TVector<double> XBorder, TVector<double> YBorder, int TASKNumber) : MethodOutside(
		N, M, XBorder, YBorder, TASKNumber) {
		a = 0;
		Ahh = 1;
		H = TMatrix<double>(m + 1);
		for (int j = 0; j <= m; j++)
			H[j] = TVector<double>(n + 1);
	}

	virtual void SetParametrs(); //Параметр метода
	virtual double Runner(); //Итерационный метод
	
	virtual TVector<double> MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations, char Name);
	virtual TVector<double> MethodError(double eps, int MaxIterations);
};

#endif
