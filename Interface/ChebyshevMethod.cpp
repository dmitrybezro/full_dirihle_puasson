#include "ChebyshevMethod.h"


void ChebyshevMethod::SetUserParametr(double Parametr)
{
	K = (int)Parametr;
}


double ChebyshevMethod::GetParametr()
{
	return K;
}


void ChebyshevMethod::SetParametrs()
{
	double Max, Min;
	Min = -4 * hE*pow(sin(pi / (2.0*n)), 2) - 4 * kE*pow(sin(pi / (2.0*m)), 2);
	Max = -4 * hE*pow(cos(pi / (2.0*n)), 2) - 4 * kE*pow(cos(pi / (2.0*m)), 2);

	Tau = TVector<double>(K);

	int step = K / 2;
	TVector<int> sup(K);
	//
	if (optimal)
	{
		for (int i = 2; i <= K; i *= 2)
		{
			for (int j = step; j < K; j += 2 * step)
			{
				sup[j] = (i - 1) - sup[j - step];
			}
			step /= 2;
		}
	}
	else {
		for (int i = 0; i < K; i++)
			sup[i] = i;
	}

	//
	for (int i = 0; i < K; i++)
	{
		Tau[i] = 2 / (Max + Min + (Max - Min)*cos(pi*(2 * sup[i] + 1) / (2.0*K)));
	}
}


double ChebyshevMethod::Runner(int index)
{
	double temp;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
			R[j][i] = A * V[j][i] + hE * (V[j][i + 1] + V[j][i - 1]) + kE * (V[j + 1][i] + V[j - 1][i]) + F[j][i];


	double accurancy = 0;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			temp = V[j][i];
			//
			V[j][i] = V[j][i] - Tau[index] * R[j][i];
			//
			temp = fabs(V[j][i] - temp);
			if (temp > accurancy)
				accurancy = temp;
		}
	return accurancy;
}


TVector<double> ChebyshevMethod::MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations, char Name)
{
	ChebyshevMethod Solution(2 * n, 2 * m, xBorder, yBorder, TaskNumber);
	switch (Name)
	{
	case 'A':
		Solution.Average();
		break;
	case 'X':
		Solution.XInterpolation();
		break;
	case 'Y':
		Solution.YInterpolation();
		break;
	}
	SetParametrs();
	Solution.SetParametrs();

	TVector<double> accurancy(2); //Точность метода
	//Чтобы зайти в цикл
	for (int i = 0; i < 2; i++)
		accurancy[i] = 1 + eps[i];

	TVector<int> IterationsCount(2);
	//
	int TauIndex = 0;
	while ((accurancy[0] > eps[0]) && (IterationsCount[0] < MaxIterations[0]))
	{
		if (TauIndex < K - 1)
		{
			Runner(TauIndex);
			TauIndex++;
		}
		else
		{
			accurancy[0] = Runner(TauIndex);
			TauIndex = 0;
			IterationsCount[0]++;
		}
	}
	TauIndex = 0;
	while ((accurancy[1] > eps[1]) && (IterationsCount[1] < MaxIterations[1]))
	{
		if (TauIndex < K - 1)
		{
			Runner(TauIndex);
			TauIndex++;
		}
		else
		{
			accurancy[1] = Solution.Runner(TauIndex);
			TauIndex = 0;
			IterationsCount[1]++;
		}
	}
	//

	SaveData("MainSolutA.txt"); //Заполнили файл численным решением
	Solution.SaveData("SupSolutA.txt");  //Заполнили файл с численным решением на вспомогательной сетке
	SaveGrid("DifferenceA.txt"); //Заполняем файл начальной информацией

	ofstream Difference("DifferenceA.txt", ios::app); //Запись в файл будет продолжена с последнего элемента

	double error = 0; //Точность решения
	double sup; //Вспомогательная переменная
	int ix, jy;
	for (int j = 0; j <= m; j++)
		for (int i = 0; i <= n; i++)
		{
			sup = fabs(V[j][i] - Solution.V[2 * j][2 * i]);
			Difference << sup << endl;
			if (sup > error)
			{
				error = sup;
				ix = i;
				jy = j;
			}
		}

	double temp;
	VectorNevyazki();
	Solution.VectorNevyazki();

	sup = NevyazkaEvkl();
	temp = Solution.NevyazkaEvkl();
	TVector<double> result(9);
	result[0] = accurancy[0]; //Точность метода на сетке (n+1,m+1)
	result[1] = IterationsCount[0]; //Количество итераций на сетке (n+1,m+1)
	result[2] = accurancy[1]; //Точность метода на сетке (2n+1,2m+1)
	result[3] = IterationsCount[1]; //Количество итераций на сетке (2n+1,2m+1)
	result[4] = error; //Точность решения
	result[5] = sup; //Евклидова норма невязки на основной сетке
	result[6] = temp; //Евклидова норма невязки на вспомогательной сетке
	result[7] = xBorder[0] + ix * h; //Значение x в самой плохой точке
	result[8] = yBorder[0] + jy * k; //Значение y в самой плохой точке

	return result;
}


TVector<double> ChebyshevMethod::MethodError(double eps, int MaxIterations)
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

	SetParametrs();

	int TauIndex = 0;
	while ((accurancy > eps) && (IterationsCount < MaxIterations))
	{
		if (TauIndex < K - 1)
		{
			Runner(TauIndex);
			TauIndex++;
		}
		else
		{
			accurancy = Runner(TauIndex);
			TauIndex = 0;
			IterationsCount++;
		}
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
