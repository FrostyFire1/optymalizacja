/*********************************************
Kod stanowi uzupełnienie materiałów do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udostępniony na licencji CC BY-SA 3.0
Autor: dr inż. Łukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include"opt_alg.h"

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab2();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	//system("pause");
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;
	int Nmax = 10000;
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
	solution opt;
	a(0) = -1;
	a(1) = 2;

	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Wahadlo
	Nmax = 1000;
	epsilon = 1e-2;
	lb = 0;
	ub = 5;
	double teta_opt = 1;

	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	ofstream Sout("symulacja_lab0.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}

void lab1()
{
	/*
	std::ofstream theory("teoria.csv");
	theory << "x(0);a;b;Liczba wywołań funkcji celu;x*;y*;Liczba wywołań funkcji celu;Minimum lokalne/globalne;x*;y*;Liczba wywołań funkcji celu;Minimum lokalne/globalne\n";
	srand(time(NULL));
	//zadanie teoretyczne
	
	double* res = new double[2] { 0, 0 };
	double x0 = 50, d = 5, alpha = 1.75;
	int Nmax = 10000;
	double min = -100, max = 100, maxRand = max - min;

	//double a = 50, b = 70;
	double epsilon = 0.00001;
	double gamma = 0.000001;
	solution wynik;
	for (int i = 0; i < 100; i++)
	{
		x0 = (static_cast<double>(rand()) / RAND_MAX) * maxRand + min;
		cout << x0 << endl;
		res = expansion(ff1T, x0, d, alpha, Nmax);
		//cout << res[0] << endl << res[1] << endl << solution::f_calls << endl << endl;
		//Sout << "x" << res[0] << ";" << "x" << res[1] << ";" << "x" << solution::f_calls << "\n";
		printf("%f, %f, %f, %d\n", x0, res[0], res[1], solution::f_calls);
		theory << format("{};{};{};{};",x0, res[0], res[1], solution::f_calls);
		//printf("============================================\n");
		wynik = fib(ff1T, res[0], res[1], epsilon);
		cout << wynik << endl;
		theory << format("{};{};{};-;", wynik.x(0), wynik.y(0), wynik.f_calls);

		//printf("============================================\n");
		wynik = lag(ff1T, res[0], res[1], epsilon, gamma, Nmax);
		//Sout << "x" << wynik.x << "x" << wynik.y << "x" << wynik.f_calls << "\n";
		cout << wynik << endl;
		theory << format("{};{};{};-\n", wynik.x(0), wynik.y(0), wynik.f_calls);
		//printf("============================================\n");

	}
	theory.close();
	
	*/
	//zadanie praktyczne
	/*
	double* res = new double[2] { 0, 0 };
	double da = 0.005, delta_da = 0.002;
	double alpha = 1.5, epsilon = 0.0001, gamma = 0.000001;
	int nmax = 1000;

	res = expansion(ff2T, da, delta_da, alpha, nmax);
	solution wynik;
	//wynik = fib(ff2T, res[0], res[1], epsilon); 
	//cout << "Metoda Fibonacciego: " << endl;
	//cout << "Optymalna wielosc otworu D_A: " << wynik.x << "\nMaksymalna temperatura wody w zbiorniku: " << wynik.y + 50 << "|" << wynik.y << "\nLiczna wywolan fukcji: " << wynik.f_calls << "\nExit flag: " << wynik.flag << endl;;
	//cout << wynik;
	wynik = lag(ff2T, res[0], res[1], epsilon, gamma, nmax);
	cout << "Metoda Lagrangea: " << endl;
	cout << "Optymalna wielosc otworu D_A: " << wynik.x << "\nMaksymalna temperatura wody w zbiorniku: " << wynik.y + 50 << "|" << wynik.y << "\nLiczna wywolan fukcji: " << wynik.f_calls << "\nExit flag: " << wynik.flag << endl;
	//cout << wynik;

	//symlacja 
	
	// Warunki początkowe
	matrix Y0(3, 1);
	Y0(0) = 5.0;   // Początkowa objętość w zbiorniku A (VA)
	Y0(1) = 1.0;   // Początkowa objętość w zbiorniku B (VB)
	Y0(2) = 20.0;  // Początkowa temperatura w zbiorniku B (TB)

	// Czas symulacji
	double t0 = 0.0;            // Początkowy czas
	double tend = 2000.0;       // Końcowy czas symulacji
	double dt = 1.0;            // Krok czasowy

	matrix ud1(1, 1), ud2;
	ud1(0) = m2d(wynik.x);
	
	matrix* S = solve_ode(df1, t0, dt, tend, Y0, ud1, ud2);

	int N = static_cast<int>(floor((tend - t0) / dt) + 1);
	std::vector<double> t_values(N);
	std::vector<double> VA_values(N);
	std::vector<double> VB_values(N);
	std::vector<double> TB_values(N);

	for (int i = 0; i < N; ++i) {
		t_values[i] = S[0](i);
		VA_values[i] = S[1](i, 0); // Pierwsza kolumna to VA
		VB_values[i] = S[1](i, 1); // Druga kolumna to VB
		TB_values[i] = S[1](i, 2); // Trzecia kolumna to TB
	}

	std::ofstream file("symulacja_lab1.csv");
	file << "Czas;VA;VB;TB\n";
	for (int i = 0; i < N; ++i) {
		file << t_values[i] << "x;" << VA_values[i] << "x;" << VB_values[i] << "x;" << TB_values[i] << "x\n";
	}
	file.close();
	*/
}

void lab2()
{
	int num_optimalizatoins = 100;
	int max_iterations = 10000;
	double epsilon = 1e-6;
	double alphaHJ = 0.5;
	double alphaRos = 1.5;
	double beta = 0.5;
	std::vector<double> step_sizes = { 1.0, 0.5, 0.1 };
	std::vector<std::vector<double>> results(num_optimalizatoins, std::vector<double>(step_sizes.size()));
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(-1.0, 1.0);

	std::ofstream file("results_table1.txt");
	//make header each other column in excel
	file << "step_size,optim_number,x1,x2,x_HJ1,x_HJ2,y,f_calls_HJ,is_global_min_HJ,x_Rosen.1,x_ROsen2,y_Rosen,f_calls_Rosen,is_global_min_Rosen\n";



	for (double step_size : step_sizes) {
		for(int i=0;i<num_optimalizatoins;i++){
			//losowanie punktu startowego
			matrix x0(2, 1);
			x0(0) = dis(gen);
			x0(1) = dis(gen);
			
			// Parametry dla Rosenbrocka
			matrix s0(2, 1);
			s0(0) = step_size;
			s0(1) = step_size;

			// Parametry dla Hooke’a-Jeevesa
			solution Xopt_HJ = HJ(ff3T, x0, step_size,alphaHJ, epsilon, max_iterations);
			bool is_global_min_HJ = (Xopt_HJ.y(0) < epsilon); // Sprawdzamy, czy osiągnięto minimum globalne

			// Parametry dla Rosenbrocka
			solution Xopt_Rosen = Rosen(ff3T, x0, s0, alphaRos, beta, epsilon, max_iterations);
			bool is_global_min_Rosen = (Xopt_Rosen.y(0) < epsilon); // Sprawdzamy, czy osiągnięto minimum globalne

			// niech wypisze mi cout dane ktore mialy by byc zapisane do pliku
			std::cout << "Podejscie dla step_size: " << step_size << " i numeru optymalizacji: " << i+1 << std::endl;
			std::cout << "step_size: " << step_size << std::endl;
			std::cout << "optim_number: " << i << std::endl;
			std::cout << "x0: " << x0(0) << std::endl;
			std::cout << "y0: " << x0(1) << std::endl;

			std::cout << "x1_HJ: " << Xopt_HJ.x(0) << std::endl;
			std::cout << "x2_HJ: " << Xopt_HJ.x(1) << std::endl;
			std::cout << "y_HJ: " << Xopt_HJ.y(0) << std::endl;
			std::cout << "f_calls_HJ: " << Xopt_HJ.f_calls << std::endl;
			std::cout << "is_global_min_HJ: " << is_global_min_HJ << std::endl;

			std::cout << "x1_Rosen: " << Xopt_Rosen.x(0) << std::endl;
			std::cout << "x2_Rosen: " << Xopt_Rosen.x(1) << std::endl;
			std::cout << "y_Rosen: " << Xopt_Rosen.y(0) << std::endl;

			std::cout << "f_calls_Rosen: " << Xopt_Rosen.f_calls << std::endl;
			std::cout << "is_global_min_Rosen: " << is_global_min_Rosen << std::endl;

			//zapis do pliku
			file << step_size << "," << i << "," << x0(0) << "," << x0(1) << "," << Xopt_HJ.x(0) << "," << Xopt_HJ.x(1) << "," << Xopt_HJ.y(0) << "," << Xopt_HJ.f_calls << "," << is_global_min_HJ << "," << Xopt_Rosen.x(0) << "," << Xopt_Rosen.x(1) << "," << Xopt_Rosen.y(0) << "," << Xopt_Rosen.f_calls << "," << is_global_min_Rosen << "\n";


		}
	}

	file.close();
	std::cout << "Zakończono obliczenia\n";



}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
