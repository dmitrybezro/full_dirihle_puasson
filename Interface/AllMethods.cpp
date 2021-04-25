#include "AllMethods.h"


AllMethods::AllMethods()
{
	pi = 3.14;
}


AllMethods::AllMethods(int N, int M, TVector<double> XBorder, TVector<double> YBorder, int TASKNumber)
{
	pi = 2 * asin(1);
	TaskNumber = TASKNumber;
	n = N; m = M;

	xBorder = XBorder; yBorder = YBorder;

	h = (xBorder[xBorder.Size() - 1] - xBorder[0]) / (1.0*n);
	k = (yBorder[yBorder.Size() - 1] - yBorder[0]) / (1.0*m);
	hE = -1 / pow(h, 2);
	kE = -1 / pow(k, 2);
	A = -2 * (hE + kE);

	V = TMatrix<double>(m + 1);
	for (int j = 0; j <= m; j++)
		V[j] = TVector<double>(n + 1);
	R = V;
	F = V;
	FunctionInicialisation();
	Inicialisation();
}


double AllMethods::F_Function(double x, double y)
{
	double res;
	switch (TaskNumber)
	{
	case 1:
		res = -pow(pi, 2)*sin(pi*x*y)*(pow(x, 2) + pow(y, 2));
		break;
	case 2:
		res = -4;
		break;
	case 3:
		res = exp(x*y)*(pow(x, 2) + pow(y, 2));
		break;
	case 4:
		res = exp(-x * pow(y, 2));
		break;
	default:
		res = 0;
	}
	return res;
}


void AllMethods::FunctionInicialisation()
{
	double x, y = yBorder[0] + k;
	for (int j = 1; j < m; j++)
	{
		x = xBorder[0] + h;
		for (int i = 1; i < n; i++)
		{
			F[j][i] = F_Function(x, y);
			x += h;
		}
		y += k;
	}
}


double AllMethods::ExactSolution(double x, double y)
{
	double res;
	switch (TaskNumber)
	{
	case 1:
		res = sin(pi*x*y);
		break;
	case 2:
		res = 1 - pow(x - 1, 2) - pow(y - 0.5, 2);
		break;
	case 3:
		res = exp(x*y);
		break;
	default:
		res = 0;
	}
	return res;
}


double AllMethods::XInicialConditions(double x, int Num)
{
	double res;
	switch (TaskNumber)
	{
	case 1:
		switch (Num)
		{
		case 1:
			res = sin(yBorder[0] * pi*x);
			break;
		case 2:
			res = sin(yBorder[yBorder.Size() - 1] * pi*x);
			break;
		}
		break;
	case 2:
		switch (Num)
		{
		case 1:
			res = 0.75 - pow(x - 1, 2);
			break;
		case 2:
			res = 0.75 - pow(x - 1, 2);
			break;
		}
		break;
	case 3:
		switch (Num)
		{
		case 1:
			res = 1;
			break;
		case 2:
			res = exp(x);
			break;
		}
		break;
	case 4:
		switch (Num)
		{
		case 1:
			res = (x - 1)*(x - 2);
			break;
		case 2:
			res = x * (x - 1)*(x - 2);
			break;
		}
		break;
	default:
		res = 0;
	}
	return res;
}


double AllMethods::YInicialConditions(double y, int Num)
{
	double res;
	switch (TaskNumber)
	{
	case 1:
		switch (Num)
		{
		case 1:
			res = sin(xBorder[0] * pi*y);
			break;
		case 2:
			res = sin(xBorder[xBorder.Size() - 1] * pi*y);
			break;
		}
		break;
	case 2:
		switch (Num)
		{
		case 1:
			res = -pow(y - 0.5, 2);
			break;
		case 2:
			res = -pow(y - 0.5, 2);
			break;
		}
		break;
	case 3:
		switch (Num)
		{
		case 1:
			res = exp(-y);
			break;
		case 2:
			res = 1;
			break;
		}
		break;
	case 4:
		switch (Num)
		{
		case 1:
			res = (y - 2)*(y - 3);
			break;
		case 2:
			res = y * (y - 2)*(y - 3);
			break;
		}
		break;
	default:
		res = 0;
	}
	return res;
}


void AllMethods::Inicialisation()
{
	double x, y = yBorder[0];
	for (int j = 0; j <= m; j++)
	{
		V[j][0] = YInicialConditions(y, 1);
		V[j][n] = YInicialConditions(y, 2);
		y += k;
	}
	x = xBorder[0] + h;
	for (int i = 1; i < n; i++)
	{
		V[0][i] = XInicialConditions(x, 1);
		V[m][i] = XInicialConditions(x, 2);
		x += h;
	}
}


void AllMethods::SetParametrs()
{

}


double AllMethods::Runner()
{
	return 0;
}


void AllMethods::GetRes(double **mas)
{
	for (int j = 0; j <= m; j++)
		for (int i = 0; i <= n; i++)
			mas[j][i] = V[j][i];
}


TVector<double> AllMethods::MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations, char Name) //Точность решения
{
	return 0;
}


TVector<double> AllMethods::MethodError(double eps, int MaxIterations) //Погрешность решения
{
	return 0;
}


void AllMethods::VectorNevyazki()
{
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
			R[j][i] = A * V[j][i] + hE * (V[j][i - 1] + V[j][i + 1]) + kE * (V[j - 1][i] + V[j + 1][i]) + F[j][i];
}


double AllMethods::NevyazkaInf()
{
	double res = 0, temp;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			temp = fabs(R[j][i]);
			if (temp > res)
				res = temp;
		}

	return res;
}


double AllMethods::NevyazkaEvkl()
{
	double res = 0;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
			res += pow(R[j][i], 2);

	return sqrt(res);
}


void AllMethods::ChangeGrid(int N, int M)
{
	n = N;
	m = M;
	h = (xBorder[xBorder.Size() - 1] - xBorder[0]) / (1.0*n);
	k = (yBorder[yBorder.Size() - 1] - yBorder[0]) / (1.0*m);
	hE = -1 / pow(h, 2);
	kE = -1 / pow(k, 2);
	A = -2 * (hE + kE);
	V = TMatrix<double>(m + 1);
	for (int j = 0; j <= m; j++)
		V[j] = TVector<double>(n + 1);
	R = V;
	F = V;
	FunctionInicialisation();
	Inicialisation();
}


void AllMethods::ChangeBorders(TVector<double> XBorder, TVector<double> YBorder)
{
	xBorder = XBorder;
	yBorder = YBorder;
	h = (xBorder[xBorder.Size() - 1] - xBorder[0]) / (1.0*n);
	k = (yBorder[yBorder.Size() - 1] - yBorder[0]) / (1.0*m);
	hE = -1 / pow(h, 2);
	kE = -1 / pow(k, 2);
	A = -2 * (hE + kE);
	V = TMatrix<double>(m + 1);
	for (int j = 0; j <= m; j++)
		V[j] = TVector<double>(n + 1);
	R = V;
	F = V;
	FunctionInicialisation();
	Inicialisation();
}


void AllMethods::ChangeTask(int TASKNumber)
{
	TaskNumber = TASKNumber;
	V = TMatrix<double>(m + 1);
	for (int j = 0; j <= m; j++)
		V[j] = TVector<double>(n + 1);
	R = V;
	F = V;
	FunctionInicialisation();
	Inicialisation();
}


void AllMethods::Reset()
{
	V = TMatrix<double>(m + 1);
	for (int j = 0; j <= m; j++)
		V[j] = TVector<double>(n + 1);
	R = V;
	F = V;
	FunctionInicialisation();
	Inicialisation();
}


void AllMethods::XInterpolation()
{
	double x, a, b;
	a = xBorder[0]; b = xBorder[xBorder.Size() - 1];
	for (int j = 1; j < m; j++)
	{
		x = xBorder[0] + h;
		for (int i = 1; i < n; i++)
		{
			V[j][i] = ((x - a)*V[j][n] - (x - b)*V[j][0]) / (b - a);
			x += h;
		}
	}
}


void AllMethods::YInterpolation()
{
	double y, c, d;
	c = yBorder[0]; d = yBorder[yBorder.Size() - 1];
	for (int i = 1; i < n; i++)
	{
		y = yBorder[0] + k;
		for (int j = 1; j < m; j++)
		{
			V[j][i] = ((y - c)*V[m][i] - (y - d)*V[0][i]) / (d - c);
			y += k;
		}
	}
}


void AllMethods::Average()
{
	double temp;
	XInterpolation();
	double y, c, d;
	c = yBorder[0]; d = yBorder[yBorder.Size() - 1];
	for (int i = 1; i < n; i++)
	{
		y = yBorder[0] + k;
		for (int j = 1; j < m; j++)
		{
			temp = ((y - c)*V[m][i] - (y - d)*V[0][i]) / (d - c);
			V[j][i] = (V[j][i] + temp) / 2;
			y += k;
		}
	}
}


void AllMethods::SaveGrid(string s)
{
	ofstream file(s);

	file << n << endl << m << endl;

	for (int i = 0; i < xBorder.Size(); i++)
		file << xBorder[i] << endl;
	file << h << endl;

	for (int i = 0; i < yBorder.Size(); i++)
		file << yBorder[i] << endl;
	file << k << endl;

	file.close();

}


void AllMethods::SaveData(string s)
{
	ofstream file(s);

	file << n << endl << m << endl;

	for (int i = 0; i < xBorder.Size(); i++)
		file << xBorder[i] << endl;
	file << h << endl;

	for (int i = 0; i < yBorder.Size(); i++)
		file << yBorder[i] << endl;
	file << k << endl;

	for (int j = 0; j <= m; j++)
		for (int i = 0; i <= n; i++)
			file << V[j][i] << endl;
	
	file.close();
	
}

