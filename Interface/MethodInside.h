#pragma once
#include "AllMethods.h"

#ifndef METHODINSIDE_H
#define METHODINSIDE_H

class MethodInside :
	public AllMethods
{
protected:
	bool optimal;
public:
	MethodInside() : AllMethods() {}
	MethodInside(int N, int M, TVector<double> XBorder, TVector<double> YBorder, int TASKNumber) : AllMethods(
		N, M, XBorder, YBorder, TASKNumber) {
	}

	virtual void Optimal(bool wish); //Разрешить или рапретить использоване оптимального параметра

	virtual TVector<double> MethodError(double eps, int MaxIterations); //Погрешность решения, для тестовых задач
};

#endif
