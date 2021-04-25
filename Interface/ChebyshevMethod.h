#pragma once
#include "AllMethods.h"
#include "MethodInside.h"

#ifndef CHEBYSHEVMETHOD_H
#define CHEBYSHEVMETHOD_H

using namespace std;

class ChebyshevMethod :
	public MethodInside
{
private:
	int K; //Количество параметров 
	TVector<double> Tau;
	// optimal здесь включает выключает сортировку параметров
public:
	ChebyshevMethod() : MethodInside() {}
	ChebyshevMethod(int N, int M, TVector<double> XBorder, TVector<double> YBorder, int TASKNumber) : MethodInside(
		N, M, XBorder, YBorder, TASKNumber) {
		Tau = TVector<double>(1);
		optimal = false;
		K = 8;
	}

	virtual double GetParametr(); //Отправить значение параметра
	virtual void SetUserParametr(double Parametr); //Уставить количество параметров

	virtual void SetParametrs(); //Параметр метода

	virtual double Runner(int index); //Итерационный метод

	virtual TVector<double> MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations, char Name); //Точность решения, для основных задач
	virtual TVector<double> MethodError(double eps, int MaxIterations); //Погрешность решения, для тестовых задач

};

#endif