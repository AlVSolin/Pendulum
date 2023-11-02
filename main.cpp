#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <fstream>

#include <Eigen/Dense> // ���������� Eigen ��� ������� ���.



int main() {
	// ���������� ����������� ����������.
	const double m = 1.0, L = 5.0; // m -- ����� �������, L -- ����� �������.
	const double x_start = 3, y_start = -4; // ��������� ���������� �������.
	const double vx_start = - 0.8, vy_start = - 0.6; // vx_start -- �������� ��������� �������� �� OX, vy_start -- �� OY.
	const double h = 0.0001; // h -- ��� �� ������� (������ - 0.0001, ��������� - 0.001).
	const double eps = 0.0001; // ��������, ���������� �� �������� ������� ������� ��������� (������ - 0.0001, ��������� - 0.001).
	const double y_step = 0.000001; // ��� �� 'y' ��� ��������� ������� ��������� x^2 + y^2 = L^2 (������ - 0.000001, ��������� - 0.0001).
	const double h2 = h * h;
	std::vector<double> x(3); // ���������� ��� �������� x_0, x_1 � x_2.
	std::vector<double> y(3); // ���������� ��� �������� y_0, y_1 � y_2.
	double t_end, t; // t_end -- �������� ������ �������, t -- ������� ������ �������.
	double g; // ������� �������� ������� g(t).
	std::ofstream output_file("output.txt"); // ���� ��� ������ ����������.

	// ��������� ����������� ������ �� ������������ (��� ��������� ��������� ������������� ����� ������ �� �����������).
	std::cout << "Please enter the end time:" << std::endl;
	std::cout << "End time is... "; std::cin >> t_end;

	// ������������� ��������� �������.
	t = 0;
	x[1] = x_start;
	y[1] = y_start;
	x[0] = x[1] - vx_start * h;
	y[0] = y[1] - vy_start * h;

	// ���������� ������ ��������� A * (x T) = b.
	Eigen::Matrix2d A(2, 2);
	Eigen::Vector2d b(2);
	Eigen::Vector2d r(2);

	int k = -1; // ��������������� �������.

	clock_t time_start = clock();

	// �������� ����, ����������� �������� x(t) � y(t).
	while (t < t_end) {
		g = 9.81 + 0.05 * sin(2*M_PI*t);

		// ����������� ���������� �������� y_2 � ����������� �� ����, ������ ��� �������� ��� ���.
		if (k == -1) {
			y[2] = y[1] - 0.01;
		}
		else {
			y[2] = y[1] - 10 * fabs(y[1] - y[0]);
		}

		// ��������� ���������� ������������� ������ A � b ���������� ��������� A * (x T) = b.
		A(0, 0) = m/h2; A(0, 1) = x[0] / L; b(0) = m * (2 * x[1] - x[0]) / h2;
		A(1, 0) = 0; A(1, 1) = y[0] / L; b(1) = m * (2 * y[1] - y[0]) / h2 - m * g - m * y[2] / h2;

		// ��������� ���������� x_2 ��� �������� ���� �������� y_2.
		x[2] = (b(0) - A(0, 1) * b(1) / A(1, 1)) / A(0, 0);

		//r = A.colPivHouseholderQr().solve(b); x[2] = r(0);

		// ���� ���������� ���������� y_2 ������, ��� ��� ������� ������� A * (x T) = b ���� x_2 � y_2 ������������� ������� x^2 + y^2 = L^2.
		while (fabs(x[2] * x[2] + y[2] * y[2] - L * L) > eps) {
			y[2] += y_step;

			// ���������� ������������� ������ A � b ���������� ��������� A * (x T) = b.
			A(0, 0) = m / h2; A(0, 1) = x[0] / L; b(0) = m * (2 * x[1] - x[0]) / h2;
			A(1, 0) = 0; A(1, 1) = y[0] / L; b(1) = m * (2 * y[1] - y[0]) / h2 - m * g - m * y[2] / h2;

			// ���������� x_2 ��� �������� ���� �������� y_2.
			x[2] = (b(0) - A(0, 1) * b(1) / A(1, 1)) / A(0, 0);

			// r = A.colPivHouseholderQr().solve(b); x[2] = r(0);
		}

		// ����� ������ �� �����:
		if (k % 100 == 0) { std::cout << "x = " << x[2] << "\t" << "y = " << y[2] << "\t" << "t = " << t << std::endl; }

		// ����� ������ � ����:
		if (k % 100 == 0) { output_file << x[2] << "\t" <<  y[2] << "\t" << t << std::endl; }

		x[0] = x[1]; x[1] = x[2]; y[0] = y[1]; y[1] = y[2];
		t += h;
		k++;
	}

	clock_t time_end = clock();

	std::cout << "Work time is " << (double)(time_end - time_start) / CLOCKS_PER_SEC << std::endl;
	
	output_file.close();

	return(0);
}