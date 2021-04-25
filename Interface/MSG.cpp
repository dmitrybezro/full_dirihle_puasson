#include "MSG.h"


void MSG::SetParametrs()
{
	double b = 0;
	double temp = 0, res = 0;

	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			res = A * H[j][i] + hE * (H[j][i + 1] + H[j][i - 1]) + kE * (H[j + 1][i] + H[j - 1][i]);
			temp += res * R[j][i];
		}

	b = temp / Ahh;

	H = R * (-1) + H * b;

	Ahh = 0; temp = 0;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			res = A * H[j][i] + hE * (H[j][i + 1] + H[j][i - 1]) + kE * (H[j + 1][i] + H[j - 1][i]);
			temp += R[j][i] * H[j][i];
			Ahh += res * H[j][i];
		}
	a = -temp / Ahh;
	
}


double MSG::Runner()
{
	double temp;
	bet = 0;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			R[j][i] = A * V[j][i] + hE * (V[j][i + 1] + V[j][i - 1]) + kE * (V[j + 1][i] + V[j - 1][i]) + F[j][i];
			bet += fabs(R[j][i]);
		}

	double accurancy = 0;
	SetParametrs();
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			temp = V[j][i];
			//
			V[j][i] = V[j][i] + a * H[j][i];
			//
			temp = fabs(V[j][i] - temp);
			if (temp > accurancy)
				accurancy = temp;
		}
	return accurancy;
}


TVector<double> MSG::MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations, char Name)
{
	MSG Solution(2 * n, 2 * m, xBorder, yBorder, TaskNumber);
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

	double sup; //��������������� ����������
	TVector<double> accurancy(2); //�������� ������
	//����� ����� � ����
	for (int i = 0; i < 2; i++)
		accurancy[i] = 1 + eps[i];

	TVector<int> IterationsCount(2);

	//
	while ((accurancy[0] > eps[0]) && (IterationsCount[0] < MaxIterations[0]))
	{
		sup = Runner();
		if (bet <= pow(10, -12))
			break;
		accurancy[0] = sup;
		IterationsCount[0]++;
	}

	while ((accurancy[1] > eps[1]) && (IterationsCount[1] < MaxIterations[1]))
	{
		sup = Solution.Runner();
		if (Solution.bet <= pow(10, -12))
			break;
		accurancy[1] = sup;
		IterationsCount[1]++;
	}
	//

	SaveData("MainSolutA.txt"); //��������� ���� ��������� ��������
	Solution.SaveData("SupSolutA.txt");  //��������� ���� � ��������� �������� �� ��������������� �����
	SaveGrid("DifferenceA.txt"); //��������� ���� ��������� �����������

	ofstream Difference("DifferenceA.txt", ios::app); //������ � ���� ����� ���������� � ���������� ��������

	double error = 0; //�������� �������
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


TVector<double> MSG::MethodError(double eps, int MaxIterations)
{
	//������� ������ �������
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

	double error = 0, accurancy = eps + 1; //����������� ������� � ������� ������ ��������������
	double sup; //���������� ��������
	int IterationsCount = 0; //���������� ��������

	while ((accurancy > eps) && (IterationsCount < MaxIterations))
	{
		sup = Runner();
		if (bet <= pow(10, -12))
			break;
		accurancy = sup;
		IterationsCount++;
	}

	SaveData("MainSolutE.txt"); //��������� ���� ��������� ��������
	SaveGrid("SupSolutE.txt");  //��������� ���� ��������� �����������
	SaveGrid("DifferenceE.txt"); //��������� ���� ��������� �����������

	ofstream SupSolut("SupSolutE.txt", ios::app), Difference("DifferenceE.txt", ios::app);//������ � ���� ����� ���������� � ���������� ��������

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
	result[0] = accurancy; //�������� ������
	result[1] = IterationsCount; //���������� ��������
	result[2] = error; //����������� �������
	result[3] = sup; //��������� ����� �������
	result[4] = xBorder[0] + h * ix; //�������� x � ����� ������ �����
	result[5] = yBorder[0] + k * jy; //�������� y � ����� ������ �����

	return result;
}
