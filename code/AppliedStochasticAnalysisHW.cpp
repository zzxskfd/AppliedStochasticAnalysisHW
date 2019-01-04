//============================================================================
// Name        : AppliedStochasticAnalysisHW.cpp
// Author      : zzx
// Version     :
// Copyright   : 
// Description :
//============================================================================

#include <iostream>
#include <fstream>
using namespace std;

extern void Metropolis();

int main() {
	// redirect cout from console to log.txt
	ofstream file;
	file.open("log.txt");
	streambuf *stream_buffer_file = file.rdbuf();
	cout.rdbuf(stream_buffer_file);
	Metropolis();
//	system("pause");
	return 0;
}
