#include <iostream>
#include <vector>
#include <valarray>
#include <complex>
#include <string>
#include <stdlib.h> 
#include <algorithm> 
#include <chrono> 
#include <fstream>


using namespace std::chrono;
using namespace std;

// Multiplication of large integers using Schoenhage-Strassen algorithm
// http://www.cs.rug.nl/~ando/pdfs/Ando_Emerencia_multiplying_huge_integers_using_fourier_transforms_paper.pdf

typedef complex<double> Complex;
typedef valarray<Complex> PolynomialCoefficients;

const double PI = 3.141592653589793238460;

void FFT(PolynomialCoefficients& a) {
	// Cooley-Tukey algorithm: radix-2 decimation-in-time Fast Fourier Transform
	const size_t N = a.size();
	if (N <= 1) return;

	PolynomialCoefficients even = a[slice(0, N / 2, 2)];
	PolynomialCoefficients odd = a[slice(1, N / 2, 2)];

	FFT(even);
	FFT(odd);

	for (size_t k = 0; k < N / 2; ++k)
	{
		Complex t = polar(1.0, -2 * PI*k / N) * odd[k];
		// polar(rho, theta) = rho(cos(theta) + isin(theta)) = rho*e^(i*theta)
		// polar(1.0, -2*PI*k/N) = e^i(-2*PI*k/N)
		a[k] = even[k] + t;
		a[k + N / 2] = even[k] - t;
	}
}


void iFFT(PolynomialCoefficients& a) {
	//inverse Fourier transform
	a = a.apply(conj); // conjugate
	FFT(a);
	a = a.apply(conj);
	a /= a.size();
}


valarray<int> multiply(PolynomialCoefficients a, PolynomialCoefficients b, bool show_details) {
	PolynomialCoefficients c;
	FFT(a);
	FFT(b);
	c = a * b;
	iFFT(c);

	if (show_details) {
		cout << endl << "FFT(a):" << endl;
		for (long i = 0; i < a.size(); i++) {
			if (i != 0) cout << ", ";
			cout << a[i];
		}
		cout << endl;

		cout << endl << "FFT(b):" << endl;
		for (long i = 0; i < b.size(); i++) {
			if (i != 0) cout << ", ";
			cout << b[i];
		}
		cout << endl;

		cout << endl << "iFFT( FFT(a).FFT(b) ):" << endl;
		for (long i = 0; i < c.size(); i++) {
			if (i != 0) cout << ", ";
			cout << c[i];
		}
		cout << endl;
	}

	valarray<int> res(a.size());
	for (long i = 0; i < a.size(); i++) {
		res[i] = round(real(c[i]));
	}

	if (show_details) {
		cout << endl << res[a.size() - 1] << "(10)^" << a.size() - 1 << "+ ";
		for (long i = 0; i < a.size() - 1; i++) {
			if (i != 0) cout << "+ ";
			cout << res[i] << "(10)^" << a.size() - 2 - i;
		}
		cout << endl << endl;
	}

	for (long i = a.size() - 1; i > 0; i--) {
		res[i] += res[i - 1] % 10;
		res[i - 1] = int(res[i - 1] / 10);
		if (res[i] > 9) res[i - 1] += int(res[i] / 10);
		res[i] %= 10;
	}

	return res;
}


void executeFastMultiplication(string& a, string& b, bool show_details) {

	long n = a.size();
	long m = b.size();
	long d = n - m;

	if (d > 0) {
		for (long i = 0; i < d; i++) b.insert(0, "0");
	}
	if (d < 0) {
		d = -d;
		for (long i = 0; i < d; i++) a.insert(0, "0");
	}

	long k = max(n, m);
	long p = 2;
	while (p < 2 * k) p *= 2;
	PolynomialCoefficients a_coeff(p);
	PolynomialCoefficients b_coeff(p);

	for (long i = 0; i < p - k; i++) a_coeff[i] = 0.0;
	for (long i = 0; i < p - k; i++) b_coeff[i] = 0.0;
	long t = p - k;
	for (int i = 0; i < k; i++) a_coeff[t + i] = a[i] - '0';
	for (int i = 0; i < k; i++) b_coeff[t + i] = b[i] - '0';

	valarray<int> c(p);

	auto start = high_resolution_clock::now();
	c = multiply(a_coeff, b_coeff, show_details);
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);

	cout << "a*b = ";

	bool started_printing = false;
	for (int j = 0; j < p; j++) {
		if (c[j] == 0 && !started_printing) {
			continue;
		}
		else {
			started_printing = true;
			cout << c[j];
		}

	}
	cout << endl;
	cout << "This was computed using Schoenhage-Strassen algorithm in " << duration.count() << " milliseconds" << endl << endl;
}


struct bigInt {
	vector<int> a;

	bigInt() {
	}

	bigInt(vector<int> b) {
		a = b;
	}

	bigInt(int b) {
		while (b) {
			a.insert(a.begin(), b % 10);
			b = int(b / 10);
		}
	}

	bigInt operator + (const bigInt &v) const {
		// make copies of the two vectors, s and t, then manipulate those copies
		vector<int> s = a;
		vector<int> t = v.a;
		int num_of_digits = max(s.size(), t.size());
		while (s.size() < num_of_digits)
			s.insert(s.begin(), 0);
		while (t.size() < num_of_digits)
			t.insert(t.begin(), 0);
		bool carry_one = false;
		for (int i = 0; i < num_of_digits; i++) {
			int q = t[num_of_digits - 1 - i] + s[num_of_digits - 1 - i];
			if (carry_one)
				q++;
			t[num_of_digits - 1 - i] = q % 10;
			if (q > 9) {
				carry_one = true;
			}
			else carry_one = false;
		}
		if (carry_one)
			t.insert(t.begin(), 1);
		return bigInt(t);
	}

	bigInt operator * (const bigInt &v) const {
		bigInt product;
		for (long i = 0; i < a.size(); ++i) {
			for (long j = 0; j < v.a.size(); ++j) {
				vector<int> p(i + j, 0);
				int h = a[i] * v.a[j];
				while (h) {
					p.insert(p.begin(), h % 10);
					h = int(h / 10);
				}
				/*
				for (auto x : p) cout << x;
				cout << endl;
				*/
				product = product + bigInt(p);
			}
		}
		return product;
	}

};


ostream& operator<<(ostream& os, const bigInt& b)
{
	for (auto i : b.a) os << i;
	return os;
}


void executeMultiplication(string& a, string& b) {

	vector<int> vector_a;
	vector<int> vector_b;
	for (auto x : a) vector_a.insert(vector_a.begin(), x - '0');
	for (auto y : b) vector_b.insert(vector_b.begin(), y - '0');
	bigInt A = bigInt(vector_a);
	bigInt B = bigInt(vector_b);


	auto start = high_resolution_clock::now();
	bigInt C = A * B;
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);

	cout << "a*b = ";
	cout << C << endl;
	cout << endl;
	cout << "This was computed using a naive algorithm in " << duration.count() << " milliseconds" << endl << endl;
}


string readFile(string& fileName) {
	string a;
	string x;

	ifstream inFile;

	inFile.open(fileName);
	if (!inFile) {
		cout << "Unable to open file " << fileName << endl;
	}
	while (inFile >> x) {
		a.append(x);
	}
	return a;
}


string getUserInputa() {
	string a;
	cout << "input a: " << endl;
	cin >> a;
	return a;
}


string getUserInputb() {
	string b;
	cout << "input b: " << endl;
	cin >> b;
	return b;
}


void showAbout() {
	cout << "Author: Michael Angel" << endl;
	cout << "email: mike00632@gmail.com" << endl;
	cout << "version: 1.0" << endl;
}


void showHelp() {
	cout << "Commands:" << endl;
	cout << "- \"multiply\": slowly multiply two input numbers" << endl;
	cout << "- \"fastmultiply\": quickly multiply two input numbers" << endl;
	cout << "- \"fastmultiplydetails\": show detailed calculations of fast multiplication" << endl;
	cout << "- \"a*b\": open files a.txt and b.txt and slow multiply" << endl;
	cout << "- \"a*bfast\": open files a.txt and b.txt and fast multiply" << endl;
	cout << "- \"about\": author and version" << endl;
	cout << "- \"help\"" << endl;
	cout << "- \"exit\"" << endl;
}


void giveIntro() {
	cout << "***********************************************************************" << endl;
	cout << "Multiplication of large intergers using Schoenhage-Strassen algorithm" << endl;
	cout << "***********************************************************************" << endl;
	cout << "Commands:" << endl;
	cout << "- \"multiply\": slowly multiply two input numbers" << endl;
	cout << "- \"fastmultiply\": quickly multiply two input numbers" << endl;
	cout << "- \"fastmultiplydetails\": show detailed calculations of fast multiplication" << endl;
	cout << "- \"a*b\": open files a.txt and b.txt and slow multiply" << endl;
	cout << "- \"a*bfast\": open files a.txt and b.txt and fast multiply" << endl;
	cout << "- \"about\": author and version" << endl;
	cout << "- \"help\"" << endl;
	cout << "- \"exit\"" << endl;
}


int main() {
	giveIntro();
	string user_input;
	while (user_input != "exit") {
		cin >> user_input;
		if (user_input == "fastmultiply" || user_input == "fm") {
			string a = getUserInputa();
			string b = getUserInputb();
			executeFastMultiplication(a, b, false);
		}
		if (user_input == "fastmultiplydetails" || user_input == "fmd") {
			string a = getUserInputa();
			string b = getUserInputb();
			executeFastMultiplication(a, b, true);
		}
		if (user_input == "multiply" || user_input == "m") {
			string a = getUserInputa();
			string b = getUserInputb();
			executeMultiplication(a, b);
		}
		if (user_input == "a*b") {
			string file1 = "a.txt";
			string file2 = "b.txt";
			string a = readFile(file1);
			string b = readFile(file2);
			executeMultiplication(a, b);
		}
		if (user_input == "a*bfast") {
			string file1 = "a.txt";
			string file2 = "b.txt";
			string a = readFile(file1);
			string b = readFile(file2);
			executeFastMultiplication(a, b, false);
		}
		if (user_input == "about") showAbout();
		if (user_input == "help") showHelp();
	}
}