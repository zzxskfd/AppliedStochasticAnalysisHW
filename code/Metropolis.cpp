#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <random>
#include <ctime>
#include <fstream>
#include <cmath>
#include <list>
#include <algorithm>
using namespace std;

inline double sg_rand() {
	// generate a single-precision random number uniformly distributed in [0, 1).
	return double(rand()) / (RAND_MAX + 1);
}

inline double db_rand() {
	// generate a double-precision random number uniformly distributed in [0, 1).
	return (double(rand()) * (RAND_MAX + 1) + double(rand())) / (RAND_MAX + 1)
			/ (RAND_MAX + 1);
}

enum Methods {
	combined, original, kmc, wolff, sw
};
const char* METHOD_NAMES[5] = {"combined", "original", "kmc", "wolff", "sw"};

typedef pair<int, int> Cd;
class IsingBoard {
private:
	int N;
	bool *s; //states, or sigma, or spins, with one filling row&column, size N*N
	double beta, h, J, Kb;
	Methods method;
	int* nbsum;	// neighbor sum, size N*N
	int kmc_n[2][5];
	// the index of 10 classes, class (i, j): state i with j nearby spins up
	double c2p[2][5]; // p_j = min(1, exp(-beta*dH_j))
	double last_time_passed, total_time_passed;
	double H, M;
	bool keep_H, keep_M;
	vector<double> spcf;
	const int spcf_len;
	bool keep_spcf;
	inline double bool2sign(bool s) {
		return s ? 1 : -1;
	}
	inline double boolMbool(bool s1, bool s2) {
		return (s1 ^ s2) ? -1 : 1;
	}
	// Coordinates and indexes
	inline Cd u(Cd cd) {
		return Cd((cd.first - 1 + N) % N, cd.second);
	}
	inline Cd d(Cd cd) {
		return Cd((cd.first + 1) % N, cd.second);
	}
	inline Cd l(Cd cd) {
		return Cd(cd.first, (cd.second - 1 + N) % N);
	}
	inline Cd r(Cd cd) {
		return Cd(cd.first, (cd.second + 1) % N);
	}
	vector<Cd> nbs_in_cd(Cd cd) {
		vector<Cd> nb = { u(cd), d(cd), l(cd), r(cd) };
		return nb;
	}
	vector<Cd> nbs_in_cd(int ind) {
		Cd cd = ind2cd(ind);
		vector<Cd> nb = { u(cd), d(cd), l(cd), r(cd) };
		return nb;
	}
	vector<int> nbs_in_ind(Cd cd) {
		vector<int> nb = { cd2ind(u(cd)), cd2ind(d(cd)), cd2ind(l(cd)), cd2ind(
				r(cd)) };
		return nb;
	}
	vector<int> nbs_in_ind(int ind) {
		Cd cd = ind2cd(ind);
		vector<int> nb = { cd2ind(u(cd)), cd2ind(d(cd)), cd2ind(l(cd)), cd2ind(
				r(cd)) };
		return nb;
	}
	inline int cd2ind(Cd cd) {
		return cd.first * N + cd.second;
	}
	inline Cd ind2cd(int ind) {
		return Cd(ind / N, ind % N);
	}
public:
	// H = - J * sum(s[i]s[j], i and j are adjacent) - h * sum(s[i])
	// M = sum(s[i])
	IsingBoard(int _N, double _beta, double _h, double _J, double _Kb,
			Methods _method, bool _keep_H = true, bool _keep_M = true,
			int _spcf_len = 0) :
			N(_N), beta(_beta), h(_h), J(_J), Kb(_Kb), method(_method), keep_H(
					_keep_H), keep_M(_keep_M), spcf_len(_spcf_len) {
		last_time_passed = 0;
		total_time_passed = 0;
		s = new bool[N * N];
		nbsum = new int[N * N];
		keep_spcf = (spcf_len > 0);
		initS();
		update_p();
		calcVars();
	}
	~IsingBoard() {
		delete s;
		delete nbsum;
	}

	void initS() {
//		for (int i = 0; i < N; ++i)
//			for (int j = 0; j < N; ++j)
//				*(s + i * (N+1) + j) = rand() % 2;		// Uniformly Random
////				*(s+i*(N+1)+j) = (i+j)%2;					// Alternating Directions
////				*(s+i*(N+1)+j) = 1;							// All Positive
//		for (int i = 0; i <= N; ++i) {
//			*(s + i * (N + 1) + N) = *(s + i * (N + 1));
//			*(s + N * (N + 1) + i) = *(s + i);
//		}
		for (int i = 0; i < N * N; ++i)
			s[i] = rand() % 2;
	}

	double get_beta() {
		return beta;
	}
	double get_H() {
		if (!keep_H) {
			cout << "Warning: not keeping H!" << endl;
		}
		return H;
	}
	double get_M() {
		if (!keep_M) {
			cout << "Warning: not keeping H!" << endl;
		}
		return M;
	}
	vector<double> get_SPCF() {
		return spcf;
	}
	double get_total_time_passed() {
		return total_time_passed;
	}
	double get_last_time_passed() {
		return last_time_passed;
	}
	void reset_time_passed() {
		total_time_passed = 0;
	}

	void set_beta(double _beta) {
		beta = _beta;
		update_p();
	}

	void calcVars() {
		calc_nbsum();
		if (keep_H || keep_M)
			calcHM();
		if (keep_spcf) {
			if (spcf_len > N) {
				cout << "Warning: spcf_len > N";
			}
			spcf = vector<double>(spcf_len, 0);
			calcSPCF();
		}
		if (method == Methods::kmc)
			calcKMC();
	}

	void calcHM() {
		// This function must be used under updated nbsum.
		double H1 = 0, H2 = 0;
		for (int i = 0; i < N * N; ++i) {
			H1 += bool2sign(s[i]) * (nbsum[i] - 2);
			H2 += bool2sign(s[i]);
		}
		H = -J * H1 - h * H2;
		M = H2;
	}

	void update_p() {
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 5; ++j) {
				double dH = 2 * J * bool2sign(i) * (2 * j - 4)
						+ 2 * h * bool2sign(i);
				c2p[i][j] = min(1.0, exp(-beta * dH));
			}
//		for (int i = 0; i < 2; ++i)
//			for (int j = 0; j < 5; ++j)
//				cout << c2p[i][j] << " ";
//		cout << endl;
//		system("pause");
	}

	void calc_nbsum() {
		for (int i = 0; i < N * N; ++i)
			nbsum[i] = 0;
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < N; ++j) {
				int tmp = s[i * N + j];
				nbsum[i * N + (j + 1) % N] += tmp;
				nbsum[i * N + (j - 1 + N) % N] += tmp;
				nbsum[((i + 1) % N) * N + j] += tmp;
				nbsum[((i - 1 + N) % N) * N + j] += tmp;
			}
	}

	void calcKMC() {
		for (int ci = 0; ci < 2; ++ci)
			for (int cj = 0; cj < 5; ++cj)
				kmc_n[ci][cj] = 0;
		for (int i = 0; i < N * N; ++i)
			++kmc_n[s[i]][nbsum[i]];
//		for (int ci = 0; ci < 2; ++ci)
//			for (int cj = 0; cj < 5; ++cj)
//				cout << kmc_n[ci][cj] << " ";
//		cout << endl;
//		system("pause");
	}

	void calcSPCF() {
		// spatial correlation function
		for (int i = 0; i < spcf_len; ++i) {
			spcf[i] = 0;
			for (int j = 0; j < N; ++j)
				for (int k = 0; k < N; ++k)
					spcf[i] += boolMbool(s[N * j + k], s[N * j + (k + i) % N]);
//			cout << spcf[i] << " ";
		}
//		cout << endl;
//		system("pause");
	}

	void flip(int ind) {
		Cd cd = ind2cd(ind);
		vector<int> nb = nbs_in_ind(cd);
		bool ci = s[ind];
		int cj = nbsum[ind];
		int dnb = ci ? -1 : 1;
		// update H, M, spcf
		if (keep_H) {
			H = H + 2 * J * bool2sign(ci) * (2 * cj - 4)
					+ 2 * h * bool2sign(ci);
		}
		if (keep_M) {
			M = M - 2 * bool2sign(ci);
		}
		if (keep_spcf) {
			for (int i = 1; i < spcf_len; ++i)
				spcf[i] = spcf[i]
						- 2 * bool2sign(ci)
								* (bool2sign(
										s[cd.first * N + (cd.second + i) % N])
										+ bool2sign(
												s[cd.first * N
														+ (cd.second - i + N)
																% N]));
		}
		// update kmc_bin
		if (method == Methods::kmc) {
			--kmc_n[ci][cj];
			++kmc_n[!ci][cj];
			for (int i = 0; i < 4; ++i) {
				int tmp = nb[i];
				--kmc_n[s[tmp]][nbsum[tmp]];
				++kmc_n[s[tmp]][nbsum[tmp] + dnb];
			}
		}
		// update nbsum
		for (int i = 0; i < 4; ++i)
			nbsum[nb[i]] += dnb;
		// flip
		s[ind] = !s[ind];
	}

//	inline void flip(Cd cd) {
//		flip(cd2ind(cd));
//	}

	void flip(list<int> inds) {
		// for small size, update one-by-one, else update as a whole
		if (int(inds.size()) < N * 4) {
			for (list<int>::iterator it = inds.begin(); it != inds.end(); ++it)
				flip(*it);
		} else {
			for (list<int>::iterator it = inds.begin(); it != inds.end(); ++it)
				s[*it] = !s[*it];
			calcVars();
		}
	}

	void transform() {
		if (method == Methods::original) {
			int i = rand() % N;
			int j = rand() % N;
			Cd cd = Cd(i, j);
			int ind = cd2ind(cd);
			bool ci = s[ind];
			int cj = nbsum[ind];
			double A = c2p[ci][cj];

			double judge = db_rand();
			if (judge < A)
				// accept
				flip(ind);
			last_time_passed = 1.0;
		} else if (method == Methods::kmc) {
			double Q[11];
			Q[0] = 0;
			int ci, cj, cind;
			for (cind = 0; cind < 10; ++cind) {
				ci = cind / 5;
				cj = cind % 5;
				Q[cind + 1] = Q[cind] + kmc_n[ci][cj] * c2p[ci][cj];
			}
			double judge = db_rand();
			for (cind = 0; cind < 10; ++cind)
				if (Q[cind + 1] / Q[10] > judge)
					break;
			ci = cind / 5;
			cj = cind % 5;
//			if (kmc_n[ci][cj] == 0 || cind < 0 || cind > 9) {
//				cout << cind << endl;
//				for (int i = 0; i < 11; ++i)
//					cout << Q[i] << " ";
//				cout << endl;
//				cout << judge << endl;
//				system("pause");
//			}
			int ind_in_c = db_rand() * kmc_n[ci][cj];
			int ind;
			for (ind = 0; ind < N * N; ++ind) {
				if (s[ind] == ci && nbsum[ind] == cj)
					--ind_in_c;
				if (ind_in_c < 0)
					break;
			}
//			if (ind_in_c != -1) {
//				cout << cind << " " << ind_in_c << " " << ind << endl;
//				system("pause");
//			}
			flip(ind);
			last_time_passed = -1.0 / Q[10] * log(1 - db_rand()) * N * N;
		} else if (method == Methods::wolff) {
			Cd cd = Cd(rand() % N, rand() % N);
			bool checked[N * N];
			for (int i = 0; i < N * N; ++i)
				checked[i] = false;
			list<int> inds;
			inds.push_back(cd2ind(cd));
			checked[cd2ind(cd)] = true;
			const bool init_s = s[inds.front()];
			const double p = 1.0 - exp(-2.0 * beta * J);
			for (list<int>::iterator inds_it = inds.begin();
					inds_it != inds.end(); ++inds_it) {
				int tmp = *inds_it;
				vector<int> tmp_nb_ind = nbs_in_ind(tmp);
				for (int j = 0; j < 4; ++j)
					if (!checked[tmp_nb_ind[j]]) {
						// bond with probability 1 - exp(-2 * beta * J)
						if (s[tmp_nb_ind[j]] == init_s) {
							if (sg_rand() < p) {
								// bonded with 1
								inds.push_back(tmp_nb_ind[j]);
								checked[tmp_nb_ind[j]] = true;
							}
						} else {
							// different color, no need to re-search
							checked[tmp_nb_ind[j]] = true;
						}
					}
			}
			flip(inds);
			last_time_passed = N * N * 4.0 / (1.0 + exp(-30.0 * beta + 13.0));
		} else if (method == Methods::sw) {
			const double p = 1.0 - exp(-2.0 * beta * J);
			bool bonds[N * N][4];
			for (int i = 0; i < N * N; ++i) {
				vector<int> nb = nbs_in_ind(i);
				bonds[i][0] = (s[i] == s[nb[0]]) && (sg_rand() < p);
				bonds[nb[0]][1] = bonds[i][0];
				bonds[i][2] = (s[i] == s[nb[2]]) && (sg_rand() < p);
				bonds[nb[2]][3] = bonds[i][2];
			}
			list<int> inds_to_flip;
			bool checked[N * N];
			for (int i = 0; i < N * N; ++i)
				checked[i] = false;
			for (int ind = 0; ind < N * N; ++ind) {
				if (checked[ind])
					continue;
				list<int> inds;
				inds.push_back(ind);
				checked[ind] = true;
				for (list<int>::iterator inds_it = inds.begin();
						inds_it != inds.end(); ++inds_it) {
					int tmp = *inds_it;
					vector<int> tmp_nb_ind = nbs_in_ind(tmp);
					for (int j = 0; j < 4; ++j)
						if (!checked[tmp_nb_ind[j]] && bonds[tmp][j]) {
							inds.push_back(tmp_nb_ind[j]);
							checked[tmp_nb_ind[j]] = true;
						}
				}
				if (rand() % 2)
					inds_to_flip.splice(inds_to_flip.end(), inds);
//				for (list<int>::iterator it = inds_to_flip.begin();it!=inds_to_flip.end();++it)
//					cout << *it << " ";
//				cout << endl;
//				system("pause");
			}
			flip(inds_to_flip);
			last_time_passed = 4 * N * N;
		} else {
			cout << "method not supported: " << method << endl;
			exit(-1);
		}
		total_time_passed += last_time_passed;
	}

	void transform_for_time(double _total_time) {
		reset_time_passed();
		bool _keep_H = keep_H;
		bool _keep_M = keep_M;
		bool _keep_spcf = keep_spcf;
		keep_H = false;
		keep_M = false;
		keep_spcf = false;
		while (total_time_passed < _total_time)
			transform();
		reset_time_passed();
		keep_H = _keep_H;
		keep_M = _keep_M;
		keep_spcf = _keep_spcf;
		calcVars();
	}

	void printGrid(const int max_i = 11, const int max_j = 11) {
		int mi = min(max_i, N), mj = min(max_j, N);
		for (int i = 0; i < mi; ++i) {
			for (int j = 0; j < mj; ++j) {
				bool tmp = *(s + i * N + j);
				cout << (tmp ? " 1" : "-1") << " ";
			}
			cout << endl;
		}
	}

};

const Methods METHOD = Methods::kmc;
const int EXP_N = 32;
// (a) find the critical T*
void work_a() {
	const double h = 0, J = 1, Kb = 1;
	const double T_star = 2.0 * J / Kb / log(1.0 + sqrt(2.0));
	int CASENUM = 41;
	int N = EXP_N;
	int NN = N * N;
	const double warmup_time = 1e7;
	const double short_warmup_time = 5e6;
	const double test_time = 1e7;
	clock_t time_stt;

	// Set betas for testing
	double betas[CASENUM];
	for (int i = 0; i < CASENUM; ++i)
		betas[i] = 1.0 / (T_star + 0.8 * (double(i) / CASENUM - 0.5));

	double Us[CASENUM];
	double Cs[CASENUM];
	cout << "Initializing..." << endl;
	IsingBoard board(N, betas[0], h, J, Kb, METHOD, true, false, 0);
	double beta, H, ttp, ltp;

	cout << endl << "Warming up......" << endl;
//	Warm-up only once: make sure betas are somewhat continuous
	board.transform_for_time(warmup_time);
	for (int casenum = 0; casenum < CASENUM; ++casenum) {
		time_stt = clock();
		beta = betas[casenum];
		board.set_beta(beta);
		cout << "beta: " << beta << " T: " << 1.0 / beta << endl;
		// Short warm-up
		board.transform_for_time(short_warmup_time);
		double EH = 0, EH2 = 0;
		for (ttp = 0; ttp < test_time;) {
			board.transform();
			H = board.get_H();
			ttp = board.get_total_time_passed();
			ltp = board.get_last_time_passed();
			EH = EH / ttp * (ttp - ltp) + H / ttp * ltp;
			EH2 = EH2 / ttp * (ttp - ltp) + H * H / ttp * ltp;
		}

		Us[casenum] = EH / NN;
		Cs[casenum] = Kb * beta * beta * (EH2 - EH * EH) / NN;

		board.printGrid(16, 16);
		cout << endl;
		cout << "U: " << Us[casenum] << endl;
		cout << "C: " << Cs[casenum] << endl;
		cout << "time elapsed: " << clock() - time_stt << "ms" << endl;
		cout << endl;
	}

	cout << "Wrting into result_a.txt..." << endl;
	ofstream fout;
	fout.open("result_a.txt");
	for (int casenum = 0; casenum < CASENUM; ++casenum)
		fout << betas[casenum] << ", " << Us[casenum] << ", " << Cs[casenum]
				<< endl;
	fout.close();

}
// (b) plot M with regard to T and h
void work_b() {
	const double J = 1, Kb = 1;
	const double T_star = 2.0 * J / Kb / log(1.0 + sqrt(2.0));
	int H_CASENUM;
	if (METHOD == Methods::wolff || METHOD == Methods::sw)
		H_CASENUM = 1;
	else
		H_CASENUM = 20;
	int CASENUM = 41;
	int N = EXP_N;
	int NN = N * N;
	const double warmup_time = 1e7;
	const double short_warmup_time = 5e6;
	const double test_time = 1e7;
	clock_t time_stt;

	// Set betas and hs for testing
	double betas[CASENUM];
	for (int i = 0; i < CASENUM; ++i)
		betas[i] = 1.0 / (T_star + 1.0 * (double(i) / CASENUM - 0.5));

	double hs[H_CASENUM];
	if (METHOD == Methods::wolff || METHOD == Methods::sw)
		hs[0] = 0;
	else
		for (int i = 0; i < H_CASENUM; ++i)
			hs[i] = 0.4 * (double(i) / H_CASENUM - 0.5);

	double Ms[H_CASENUM][CASENUM];
	double h, beta, M, ttp, ltp;
	for (int h_casenum = 0; h_casenum < H_CASENUM; ++h_casenum) {
		h = hs[h_casenum];
		cout << "Initializing..." << endl;
		IsingBoard board(N, betas[0], h, J, Kb, METHOD, false, true, 0);
		cout << endl << "Warming up......" << endl;
		//	Warm-up only once: make sure betas are somewhat continuous
		board.transform_for_time(warmup_time);
		for (int casenum = 0; casenum < CASENUM; ++casenum) {
			time_stt = clock();
			beta = betas[casenum];
			board.set_beta(beta);
			cout << "h: " << h << " beta: " << beta << " T: " << 1.0 / beta
					<< endl;
			// Short warm-up
			board.transform_for_time(short_warmup_time);
			double EM = 0;
			for (ttp = 0; ttp < test_time;) {
				board.transform();
				M = board.get_M();
				ttp = board.get_total_time_passed();
				ltp = board.get_last_time_passed();
				EM = EM / ttp * (ttp - ltp) + M / ttp * ltp;
			}
			Ms[h_casenum][casenum] = EM / NN;

			board.printGrid(16, 16);
			cout << endl;
			cout << "m: " << Ms[h_casenum][casenum] << endl;
			cout << "time elapsed: " << clock() - time_stt << "ms" << endl;
			cout << endl;
		}
	}

	cout << "Wrting into result_b.txt..." << endl;
	ofstream fout;
	fout.open("result_b.txt");
	fout << "NaN";
	for (int casenum = 0; casenum < CASENUM; ++casenum)
		fout << ", " << betas[casenum];
	fout << endl;
	for (int h_casenum = 0; h_casenum < H_CASENUM; ++h_casenum) {
		fout << hs[h_casenum];
		for (int casenum = 0; casenum < CASENUM; ++casenum)
			fout << ", " << Ms[h_casenum][casenum];
		fout << endl;
	}
	fout.close();

}
// (c) spatial correlation function
void work_c() {
	const double h = 0, J = 1, Kb = 1;
	const double T_star = 2.0 / log(1.0 + sqrt(2.0));
	const int CASENUM = 41;
	const int N = EXP_N;
	const int NN = N * N;
	const int spcf_len = N / 3 + 1;

	const double warmup_time = 1e7;
	const double short_warmup_time = 5e6;
	const double test_time = 1e7;
	clock_t time_stt;
//	int time_records[CASENUM];

	// Set betas for testing
	double betas[CASENUM];
	for (int i = 0; i < CASENUM; ++i)
		betas[i] = 1.0 / (T_star + 0.05 + 0.5 * (double(i) / CASENUM));
	cout << "Initializing..." << endl;
	IsingBoard board(N, betas[0], h, J, Kb, METHOD, false, false, spcf_len);
	double beta, ltp, ttp;

	cout << endl << "Warming up......" << endl;
//	Warm-up only once: be sure betas are somewhat continuous
	board.transform_for_time(warmup_time);
	cout << endl;

//	double ksis[CASENUM];
	double spcf[CASENUM][spcf_len];
	for (int casenum = 0; casenum < CASENUM; ++casenum) {
		time_stt = clock();
		beta = betas[casenum];
		board.set_beta(beta);
		cout << "beta: " << beta << " T: " << 1.0 / beta << endl;
		for (int j = 0; j < spcf_len; ++j)
			spcf[casenum][j] = 0;
		// Short warm-up
		board.transform_for_time(short_warmup_time);
		for (ttp = 0; ttp < test_time;) {
			board.transform();
			ttp = board.get_total_time_passed();
			ltp = board.get_last_time_passed();
			vector<double> tmp_spcf = board.get_SPCF();
			for (int j = 0; j < spcf_len; ++j)
				spcf[casenum][j] = spcf[casenum][j] / ttp * (ttp - ltp)
						+ tmp_spcf[j] / ttp * ltp;
		}
		for (int j = 0; j < spcf_len; ++j)
			spcf[casenum][j] /= NN;

//		// linear regression
//		double t1 = 0, t2 = 0, t3 = 0, t4 = 0, tmp;
//		for (int i = 1; i < spcf_len; ++i) {
//			t1 += i * i;
//			t2 += i;
//			tmp = log(abs(spcf[casenum][i]) + 1e-10);
//			t3 += i * tmp;
//			t4 += tmp;
//		}
////		cout << t1 << " " << t2 << " " << t3 << " " << t4 << endl;
////		system("pause");
//		ksis[casenum] = (t3 * spcf_len - t2 * t4) / (t1 * spcf_len - t2 * t2);
//		ksis[casenum] = -1.0 / ksis[casenum];

		board.printGrid(16, 16);
		cout << endl;
		for (int i = 0; i < spcf_len; ++i)
			cout << spcf[casenum][i] << " ";
		cout << endl;
//		cout << "ksi = " << ksis[casenum] << endl;
		cout << "time elapsed: " << clock() - time_stt << "ms" << endl;
		cout << endl;
//		time_records[casenum] = clock() - time_stt;
	}

//	for (int i = 0; i < CASENUM; ++i)
//		cout << ksis[i] << " ";
//	cout << endl;

	cout << "Wrting into result_c.txt..." << endl;
	ofstream fout;
	fout.open("result_c.txt");
	for (int casenum = 0; casenum < CASENUM; ++casenum) {
		fout << betas[casenum];
		for (int i = 0; i < spcf_len; ++i)
			fout << ", " << spcf[casenum][i];
		fout << endl;
	}
	fout.close();

//	cout << "Wrting into time_records_c.txt..." << endl;
//	fout.open("time_records_c.txt");
//	for (int casenum = 0; casenum < CASENUM; ++casenum) {
//		fout << betas[casenum] << ", " << time_records[casenum] << endl;
//	}
//	fout.close();
}
// (d1) calculate alpha
void work_d1() {
	const double J = 1, Kb = 1;
	const double T_star = 2.0 * J / Kb / log(1.0 + sqrt(2.0));
	int H_CASENUM = 1;
	int CASENUM = 41;
	int N = EXP_N;
	int NN = N * N;
	const double warmup_time = 1e7;
	const double short_warmup_time = 5e6;
	const double test_time = 1e7;
	clock_t time_stt;

	// Set betas and hs for testing
	double betas[CASENUM];
	for (int i = 0; i < CASENUM; ++i)
		betas[i] = 1.0 / (T_star + 0.05 + 0.5 * (double(i) / CASENUM));

	double hs[H_CASENUM];
	hs[0] = 0.01;

	double Ms[H_CASENUM][CASENUM];
	double h, beta, M, ttp, ltp;
	for (int h_casenum = 0; h_casenum < H_CASENUM; ++h_casenum) {
		h = hs[h_casenum];
		cout << "Initializing..." << endl;
		IsingBoard board(N, betas[0], h, J, Kb, METHOD, false, true, 0);
		cout << endl << "Warming up......" << endl;
		//	Warm-up only once: make sure betas are somewhat continuous
		board.transform_for_time(warmup_time);
		for (int casenum = 0; casenum < CASENUM; ++casenum) {
			time_stt = clock();
			beta = betas[casenum];
			board.set_beta(beta);
			cout << "h: " << h << " beta: " << beta << " T: " << 1.0 / beta
					<< endl;
			// Short warm-up
			board.transform_for_time(short_warmup_time);
			double EM = 0;
			for (ttp = 0; ttp < test_time;) {
				board.transform();
				M = board.get_M();
				ttp = board.get_total_time_passed();
				ltp = board.get_last_time_passed();
				EM = EM / ttp * (ttp - ltp) + M / ttp * ltp;
			}
			Ms[h_casenum][casenum] = EM / NN;

			board.printGrid(16, 16);
			cout << endl;
			cout << "m: " << Ms[h_casenum][casenum] << endl;
			cout << "time elapsed: " << clock() - time_stt << "ms" << endl;
			cout << endl;
		}
	}

	cout << "Wrting into result_d1.txt..." << endl;
	ofstream fout;
	fout.open("result_d1.txt");
	fout << "NaN";
	for (int casenum = 0; casenum < CASENUM; ++casenum)
		fout << ", " << betas[casenum];
	fout << endl;
	for (int h_casenum = 0; h_casenum < H_CASENUM; ++h_casenum) {
		fout << hs[h_casenum];
		for (int casenum = 0; casenum < CASENUM; ++casenum)
			fout << ", " << Ms[h_casenum][casenum];
		fout << endl;
	}
	fout.close();

}
// (d2) calculate gamma
void work_d2() {
	const double h = 0, J = 1, Kb = 1;
	const double T_star = 2.0 * J / Kb / log(1.0 + sqrt(2.0));
	int CASENUM = 41;
	int N = EXP_N;
	int NN = N * N;
	const double warmup_time = 1e7;
	const double short_warmup_time = 5e6;
	const double test_time = 1e7;
	clock_t time_stt;

	// Set betas for testing
	double betas[CASENUM];
	for (int i = 0; i < CASENUM; ++i)
		betas[i] = 1.0 / (T_star + 0.05 + 0.5 * (double(i) / CASENUM));

	double Us[CASENUM];
	double Cs[CASENUM];
	cout << "Initializing..." << endl;
	IsingBoard board(N, betas[0], h, J, Kb, METHOD, true, false, 0);
	double beta, H, ttp, ltp;

	cout << endl << "Warming up......" << endl;
//	Warm-up only once: make sure betas are somewhat continuous
	board.transform_for_time(warmup_time);
	for (int casenum = 0; casenum < CASENUM; ++casenum) {
		time_stt = clock();
		beta = betas[casenum];
		board.set_beta(beta);
		cout << "beta: " << beta << " T: " << 1.0 / beta << endl;
		// Short warm-up
		board.transform_for_time(short_warmup_time);
		double EH = 0, EH2 = 0;
		for (ttp = 0; ttp < test_time;) {
			board.transform();
			H = board.get_H();
			ttp = board.get_total_time_passed();
			ltp = board.get_last_time_passed();
			EH = EH / ttp * (ttp - ltp) + H / ttp * ltp;
			EH2 = EH2 / ttp * (ttp - ltp) + H * H / ttp * ltp;
		}

		Us[casenum] = EH / NN;
		Cs[casenum] = Kb * beta * beta * (EH2 - EH * EH) / NN;

		board.printGrid(16, 16);
		cout << endl;
		cout << "U: " << Us[casenum] << endl;
		cout << "C: " << Cs[casenum] << endl;
		cout << "time elapsed: " << clock() - time_stt << "ms" << endl;
		cout << endl;
	}

	cout << "Wrting into result_d2.txt..." << endl;
	ofstream fout;
	fout.open("result_d2.txt");
	for (int casenum = 0; casenum < CASENUM; ++casenum)
		fout << betas[casenum] << ", " << Us[casenum] << ", " << Cs[casenum]
				<< endl;
	fout.close();

}
// (d3) calculate delta
void work_d3() {
	const double h = 0, J = 1, Kb = 1;
	const double T_star = 2.0 / log(1.0 + sqrt(2.0));
	const int CASENUM = 41;
	const int N = EXP_N;
	const int NN = N * N;
	const int spcf_len = N / 3 + 1;

	const double warmup_time = 1e7;
	const double short_warmup_time = 5e6;
	const double test_time = 1e7;
	clock_t time_stt;
//	int time_records[CASENUM];

	// Set betas for testing
	double betas[CASENUM];
	for (int i = 0; i < CASENUM; ++i)
		betas[i] = 1.0 / (T_star + 0.05 + 0.5 * (double(i) / CASENUM));
	cout << "Initializing..." << endl;
	IsingBoard board(N, betas[0], h, J, Kb, METHOD, false, false, spcf_len);
	double beta, ltp, ttp;

	cout << endl << "Warming up......" << endl;
//	Warm-up only once: be sure betas are somewhat continuous
	board.transform_for_time(warmup_time);
	cout << endl;

//	double ksis[CASENUM];
	double spcf[CASENUM][spcf_len];
	for (int casenum = 0; casenum < CASENUM; ++casenum) {
		time_stt = clock();
		beta = betas[casenum];
		board.set_beta(beta);
		cout << "beta: " << beta << " T: " << 1.0 / beta << endl;
		for (int j = 0; j < spcf_len; ++j)
			spcf[casenum][j] = 0;
		// Short warm-up
		board.transform_for_time(short_warmup_time);
		for (ttp = 0; ttp < test_time;) {
			board.transform();
			ttp = board.get_total_time_passed();
			ltp = board.get_last_time_passed();
			vector<double> tmp_spcf = board.get_SPCF();
			for (int j = 0; j < spcf_len; ++j)
				spcf[casenum][j] = spcf[casenum][j] / ttp * (ttp - ltp)
						+ tmp_spcf[j] / ttp * ltp;
		}
		for (int j = 0; j < spcf_len; ++j)
			spcf[casenum][j] /= NN;

		board.printGrid(16, 16);
		cout << endl;
		for (int i = 0; i < spcf_len; ++i)
			cout << spcf[casenum][i] << " ";
		cout << endl;
		cout << "time elapsed: " << clock() - time_stt << "ms" << endl;
		cout << endl;
	}

	cout << "Wrting into result_d3.txt..." << endl;
	ofstream fout;
	fout.open("result_d3.txt");
	for (int casenum = 0; casenum < CASENUM; ++casenum) {
		fout << betas[casenum];
		for (int i = 0; i < spcf_len; ++i)
			fout << ", " << spcf[casenum][i];
		fout << endl;
	}
	fout.close();
}

void Metropolis() {
	int rseed = time(0);
	cout << "random seed: " << rseed << endl;
	srand(rseed);
	cout << "method: " << METHOD_NAMES[METHOD] << endl;
	cout << "N: " << EXP_N << endl;
	cout << endl;

	clock_t time_stt = clock();
//	work_a();
//	work_b();
//	work_c();
//	work_d1();
//	work_d2();
	work_d3();
	cout << endl << "time elapsed: " << clock() - time_stt << "ms" << endl;

}

