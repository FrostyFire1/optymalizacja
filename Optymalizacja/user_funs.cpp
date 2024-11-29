﻿#define _USE_MATH_DEFINES
#include"user_funs.h"
#include <cmath>
matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
	int n = get_len(Y[0]);
	double teta_max = Y[1](0, 0);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	double m = 1, l = 0.5, b = 0.5, g = 9.81;
	double I = m*pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1))*ud2(0) - m*g*l*sin(Y(0)) - b*Y(1)) / I;
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2) {
	return -cos(0.1 * x(0)) * exp(-pow((0.1 * x(0) - 2 * 3.14), 2)) + 0.002 * pow((0.1 * x(0)), 2);
}

matrix df1(double t, matrix Y, matrix ud1, matrix ud2) {
    matrix dY(Y);  // dY to macierz przechowująca pochodne VA', VB' i TB'

    // Parametry stałe
    double a = 0.98;            // wsp. lepkości
    double b = 0.63;            // wsp. zwężenia strumienia
    double g = 9.81;            // przysp. ziemskie
    double PA = 0.5;            // pole pow. zbiornika A
    double PB = 1;              // pole pow. zbiornika B
    double DB = 0.00365665;     // wielkość otworu w zbiorniku B
    double Fin = 0.01;          // ilość wody wpływającej do B (z rury)
    double Tin = 20;            // temp. wody wpływającej
    double TA = 90.0;           // temp. wody w zbiorniku A
    double DA = ud1(0);         // wielkość otworu w zbiorniku A (parametr z ud2)

    // Stan w danym momencie, podany przez wektor Y
    double VA = Y(0);           // objętość w zbiorniku A
    double VB = Y(1);           // objętość w zbiorniku B
    double TB = Y(2);           // temperatura w zbiorniku B

    // Obliczenie przepływu wody z A i B
    double FAout = VA > 0 ? a * b * DA * sqrt(2 * g * VA / PA) : 0;
    double FBout = VB > 0 ? a * b * DB * sqrt(2 * g * VB / PB) : 0;

    // Obliczenie pochodnych VA', VB' i TB'
    dY(0) = -FAout;                                 // VA' (zmiana objętości w A)
    dY(1) = FAout + Fin - FBout;                    // VB' (zmiana objętości w B)
    dY(2) = FAout / VB * (TA - TB) + Fin / VB * (Tin - TB); // TB' (zmiana temperatury w B)

    return dY;
}

matrix ff2T(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	matrix Y0 = matrix(3, new double[3] {5, 1, 20});
	matrix* Y = solve_ode(df1, 0, 1, 2000, Y0, x);
	int n = get_len(Y[0]);
	double max = Y[1](0, 2);
	for (int i = 0; i < n; i++) {
		if (max < Y[1](i, 2)) {
			max = Y[1](i, 2);
		}
	}
	y = abs(max - 50);
	return y[0];
}

matrix df2(double t, matrix Y, matrix ud1, matrix ud2) {
	matrix dY(Y);
	double l = 1.0; // długość ramienia robota
	double mr = 1.0; // masa ramienia robota
	double mc = 5.0; // masa ciężarka
	double b = 0.5; // Współczynnik tarcia
	double I = (1.0 / 3.0) * (mr * pow(l, 2)) + mc * pow(l, 2);//moment bezwładności

	double Mt = ud2(0) * (ud1(0) - Y(0)) + ud2(1) * (ud1(1) - Y(1));
	
	dY(0) = Y(1);
	dY(1) = (Mt - b * Y(1)) / I;
	return dY;
}

matrix ff2ActuallyT(matrix x, matrix ud1, matrix ud2) {
	double x1 = x(0);
	double x2 = x(1);

	double result = pow(x1, 2) + pow(x2, 2) - cos(2.5 * 3.14 * x1) - cos(2.5 * 3.14 * x2) + 2;

	return matrix(1, 1, result);
}
matrix ff2ActuallyR(matrix x, matrix ud1, matrix ud2) {
	double y = 0;
	matrix Y0(2, 1), Yref(2, new double[2] {3.14, 0});
	matrix* Y = solve_ode(df2, 0, 0.1, 100, Y0, Yref, x);
	int n = get_len(Y[0]);
	for (int i = 0; i < n; i++) {
		y += 10 * pow(Yref(0) - Y[1](i, 0), 2);
		y += pow(Yref(1) - Y[1](i, 1), 2);
		y += pow(x(0)*(Yref(0) - Y[1](i,0)) + x(1)*(Yref(1) - Y[1](i,1)),2);
	}
	return y * 0.1;
}

matrix ff3T_outside(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = (sin(M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2)))) / (M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2)));

	double a = m2d(ud1);
	double c = m2d(ud2);

	//g1
	if (-x(0) + 1 > 0)
		y = y + c * pow(-x(0) + 1, 2);
	//g2
	if (-x(1) + 1 > 0)
		y = y + c * pow(-x(1) + 1, 2);
	//g3
	if (norm(x) - a > 0)
		y = y + c * pow(norm(x) - a, 2);

	return y;
}

matrix ff3T_inside(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = (sin(M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2)))) / (M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2)));

	double a = m2d(ud1);
	double c = m2d(ud2);

	//g1
	if (-x(0) + 1 > 0)
		y = 1E10;
	else
		y = y - c / (-x(0) + 1);
	//g2
	if (-x(1) + 1 > 0)
		y = 1E10;
	else
		y = y - c / (-x(1) + 1);
	//g3
	if (norm(x) - a > 0)
		y = 1E10;
	else
		y = y - c / (norm(x) - a);

	return y;
}

matrix ff3R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0(4, new double[4] {0.0, x(0), 100, 0});
	matrix* Y = solve_ode(df3, 0.0, 0.01, 7.0, Y0, ud1, x(1));
	int n = get_len(Y[0]);
	int i50 = 0;
	int i_0 = 0;
	for (int i = 0; i < n; ++i)
	{
		if (abs(Y[1](i, 2) - 50.0) < abs(Y[1](i50, 2) - 50.0))
			i50 = i;
		if (abs(Y[1](i, 2)) < abs(Y[1](i_0, 2)))
			i_0 = i;
	}

	y = -Y[1](i_0, 0);

	if (abs(x(0)) - 10 > 0)
		y = y + ud2 * pow(abs(x(0)) - 10, 2);
	if (abs(x(1)) - 15 > 0)
		y = y + ud2 * pow(abs(x(1)) - 15, 2);
	if (abs(Y[1](i50, 0) - 5.0) - 0.5 > 0)
		y = y + ud2 * pow(abs(Y[1](i50, 0) - 5.0) - 0.5, 2);

	return y;
}

matrix df3(double t, matrix Y, matrix ud1, matrix ud2)
{
	//Wektor zmian po czasie
	matrix dY(4, 1);

	//Zmienne dane
	double x = Y(0);
	double v_x = Y(1);
	double y = Y(2);
	double v_y = Y(3);

	//Dane zadania
	double C = ud1(0);
	double rho = ud1(1);
	double r = ud1(2);
	double m = ud1(3);
	double g = ud1(4);

	double s = M_PI * pow(r, 2);

	double D_x = (1.0 / 2.0) * C * rho * s * v_x * abs(v_x);
	double D_y = (1.0 / 2.0) * C * rho * s * v_y * abs(v_y);
	double FM_x = rho * v_y * m2d(ud2) * M_PI * pow(r, 3);
	double FM_y = rho * v_x * m2d(ud2) * M_PI * pow(r, 3);

	dY(0) = v_x;
	dY(1) = (-D_x - FM_x) / m;
	dY(2) = v_y;
	dY(3) = ((-m * g) - D_y - FM_y) / m;

	return dY;
}