#pragma once

#ifndef ALLMETHODS_H
#define ALLMETHODS_H

#include <string>
#include <fstream>
#include <iostream>
#include <math.h>
#include "TMatrix.h"

using namespace std;

class AllMethods
{
protected:
	TMatrix<double> V, R, F; //Решение задачи, невязка и значение функции в правой части соответственно
	TVector<double> xBorder, yBorder; //Границы по оси x и y соответственно
	double h, k; //Шаг по оси x и y соответственно
	double hE, kE, A; //Вспомогательные данные
	int TaskNumber, n, m; //Номер задачи, размерность сетки
	double pi;

public:
	AllMethods();
	AllMethods(int N, int M, TVector<double> XBorder, TVector<double> YBorder, int TASKNumber); //Конструктор инициализации

	double F_Function(double x, double y); //Функция в правой части (-f(x,y))
	void FunctionInicialisation();
	double ExactSolution(double x, double y);//Точные решения
	double XInicialConditions(double x, int Num); //Начальные условия в виде функций параллельно оси X
	double YInicialConditions(double y, int Num); //Начальные условия в виде функций параллельно оси Y
	void Inicialisation(); //Заполнение начальных условий в наше решение

	virtual void SetParametrs(); //Параметр метода

	virtual double Runner(); //Итерационный метод

	void GetRes(double **mas); //Записать решение в динамический массив

	virtual TVector<double> MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations, char Name); //Точность решения, для основных задач
	virtual TVector<double> MethodError(double eps, int MaxIterations); //Погрешность решения, для тестовых задач

	void VectorNevyazki();
	double NevyazkaInf(); //Невязка по норме бесконечность
	double NevyazkaEvkl(); //Евклидова норма невязки

	void ChangeGrid(int N, int M); //Изменить сетку
	void ChangeBorders(TVector<double> XBorder, TVector<double> YBorder); //Изменить границу
	void ChangeTask(int TASKNumber); //Изменить задачу
	void Reset();//Обнуление результатов

	void XInterpolation();
	void YInterpolation();
	void Average();

	void SaveGrid(string s);
	void SaveData(string s);

	virtual void SetUserParametr(double Parametr) { }
	virtual double GetParametr() { return 0; }
	virtual void Optimal(bool wish) {}

};

#endif
