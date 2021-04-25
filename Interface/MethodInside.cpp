#include "MethodInside.h"


void  MethodInside::Optimal(bool wish)
{
	optimal = wish;
}


TVector<double> MethodInside::MethodError(double eps, int MaxIterations)
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

	if (optimal)
		SetParametrs();

	while ((accurancy > eps) && (IterationsCount < MaxIterations))
	{
		accurancy = Runner();
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