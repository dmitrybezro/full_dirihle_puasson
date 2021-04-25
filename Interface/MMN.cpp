#include "MMN.h"


void MMN::SetParametrs()
{
	double Arr = 0, ArAr = 0;
	double temp;
	for (int j = 1; j < m; j++)
	{
		for (int i = 1; i < n; i++)
		{
			temp = A * R[j][i] + hE * (R[j][i + 1] + R[j][i - 1]) + kE * (R[j + 1][i] + R[j - 1][i]);
			Arr += temp * R[j][i];
			ArAr += pow(temp, 2);
		}
	}
	Tau = Arr / ArAr;
}


double MMN::Runner()
{
	double temp;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
			R[j][i] = A * V[j][i] + hE * (V[j][i + 1] + V[j][i - 1]) + kE * (V[j + 1][i] + V[j - 1][i]) + F[j][i];

	SetParametrs();
	double accurancy = 0;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			temp = V[j][i];
			//
			V[j][i] = V[j][i] - Tau * R[j][i];
			//
			temp = fabs(V[j][i] - temp);
			if (temp > accurancy)
				accurancy = temp;
		}
	return accurancy;
}


TVector<double> MMN::MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations, char Name)
{
	MMN Solution(2 * n, 2 * m, xBorder, yBorder, TaskNumber);
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

	TVector<double> accurancy(2); //�������� ������
	//����� ����� � ����
	for (int i = 0; i < 2; i++)
		accurancy[i] = 1 + eps[i];

	TVector<int> IterationsCount(2);
	//
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
	//

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

	sup = NevyazkaEvkl();
	temp = Solution.NevyazkaEvkl();
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

