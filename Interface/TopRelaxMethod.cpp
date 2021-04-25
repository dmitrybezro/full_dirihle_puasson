#include "TopRelaxMethod.h"


void TopRelaxMethod::SetUserParametr(double Parametr)
{
	w = Parametr;
}


double TopRelaxMethod::GetParametr()
{
	return w;
}


void TopRelaxMethod::SetParametrs()
{
	double temp1, temp2;
	temp1 = pow(sin(pi / (2.0*n)), 2);
	temp2 = pow(sin(pi / (2.0*m)), 2);

	temp2 = 2 * pow(k, 2)*temp1 + 2 * pow(h, 2)*temp2;
	temp1 = pow(h, 2) + pow(k, 2);
	temp2 = temp2 / temp1;
	w = 2 / (1 + sqrt(temp2*(2 - temp2))); //����������� ��������
	
}


double TopRelaxMethod::Runner()
{
	double temp, prev;
	double errorM = 0; //�������� ������
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			prev = V[j][i];
			//
			temp = A * prev + hE * (V[j][i - 1] + V[j][i + 1]) + kE * (V[j - 1][i] + V[j + 1][i]);
			V[j][i] = prev - w * (temp + F[j][i]) / A;
			//
			temp = fabs(V[j][i] - prev);
			if (temp > errorM)
				errorM = temp;
		}
	return errorM;
}


TVector<double> TopRelaxMethod::MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations, char Name)
{
	TopRelaxMethod Solution(2 * n, 2 * m, xBorder, yBorder, TaskNumber);
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
	if (optimal)
	{
		SetParametrs();
		Solution.SetParametrs();
	}
	else {
		Solution.Optimal(false);
		Solution.SetUserParametr(w);
	}

	TVector<double> accurancy(2); //�������� ������
	//����� ����� � ����
	for (int i = 0; i < 2; i++)
		accurancy[i] = 1 + eps[i];

	TVector<int> IterationsCount(2);
	while ((accurancy[0] > eps[0]) && (IterationsCount[0] < MaxIterations[0]))
	{
		accurancy[0] = Runner();
		IterationsCount[0]++;
	}
	while ((accurancy[1] > eps[1]) && (IterationsCount[1] < MaxIterations[1]))
	{
		accurancy[1] = Solution.Runner();
		IterationsCount[1]++;
	}

	SaveData("MainSolutA.txt"); //��������� ���� ��������� ��������
	Solution.SaveData("SupSolutA.txt");  //��������� ���� � ��������� �������� �� ��������������� �����
	SaveGrid("DifferenceA.txt"); //��������� ���� ��������� �����������

	ofstream Difference("DifferenceA.txt", ios::app); //������ � ���� ����� ���������� � ���������� ��������

	double error = 0; //�������� �������
	double sup; //��������������� ����������
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
	temp = Solution.NevyazkaEvkl();
	sup = NevyazkaEvkl();
	TVector<double> result(9);
	result[0] = accurancy[0]; //�������� ������ �� ����� (n+1,m+1)
	result[1] = IterationsCount[0]; //���������� �������� �� ����� (n+1,m+1)
	result[2] = accurancy[1]; //�������� ������ �� ����� (2n+1,2m+1)
	result[3] = IterationsCount[1]; //���������� �������� �� ����� (2n+1,2m+1)
	result[4] = error; //�������� �������
	result[5] = sup; //��������� ����� ������� �� �������� �����
	result[6] = temp; //��������� ����� ������� �� ��������������� �����
	result[7] = xBorder[0] + ix * h; //�������� x � ����� ������ �����
	result[8] = yBorder[0] + jy * k; //�������� y � ����� ������ �����

	return result;
}

