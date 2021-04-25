#include "MethodInside.h"


void  MethodInside::Optimal(bool wish)
{
	optimal = wish;
}


TVector<double> MethodInside::MethodError(double eps, int MaxIterations)
{
	//Считаем точное решение
	TMatrix<double> u(m + 1);
	double x, y = yBorder[0];
	for (int j = 0; j <= m; j++)
	{
		u[j] = TVector<double>(n + 1);
		x = xBorder[0];
		for (int i = 0; i <= n; i++)
		{
			u[j][i] = ExactSolution(x, y);
			x += h;
		}
		y += k;
	}

	double error = 0, accurancy = eps + 1; //Погрешность решения и точость метода соответственно
	double sup; //Переменная помощник
	int IterationsCount = 0; //Количество итераций

	if (optimal)
		SetParametrs();

	while ((accurancy > eps) && (IterationsCount < MaxIterations))
	{
		accurancy = Runner();
		IterationsCount++;
	}

	SaveData("MainSolutE.txt"); //Заполнили файл численным решением
	SaveGrid("SupSolutE.txt");  //Заполнили файл начальной информацией
	SaveGrid("DifferenceE.txt"); //Заполняем файл начальной информацией

	ofstream SupSolut("SupSolutE.txt", ios::app), Difference("DifferenceE.txt", ios::app);//Запись в файл будет продолжена с последнего элемента

	int ix, jy;
	for (int j = 0; j <= m; j++)
		for (int i = 0; i <= n; i++)
		{
			sup = fabs(V[j][i] - u[j][i]);
			SupSolut << u[j][i] << endl;
			Difference << sup << endl;
			if (sup > error)
			{
				error = sup;
				ix = i;
				jy = j;
			}
		}

	VectorNevyazki();
	sup = NevyazkaEvkl();
	TVector<double> result(6);
	result[0] = accurancy; //Точность метода
	result[1] = IterationsCount; //Количество итераций
	result[2] = error; //Погрешность решения
	result[3] = sup; //Евклидова норма невязки
	result[4] = xBorder[0] + h * ix; //Значение x в самой плохой точке
	result[5] = yBorder[0] + k * jy; //Значение y в самой плохой точке

	return result;
}