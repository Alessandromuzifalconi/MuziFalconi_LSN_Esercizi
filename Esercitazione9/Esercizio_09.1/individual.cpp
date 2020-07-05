#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "individual.h"


using namespace std;

Individual::Individual() {
	m_N = 32;
	m_s = new int[m_N];
	for (int i = 0; i < m_N; i++) 
		m_s[i] = i;
}

Individual::Individual(int N, int *x) {
	m_N = N;
	m_s = new int[m_N];
	for (int i = 0; i < m_N; i++) 
		m_s[i] = x[i];
}

Individual :: ~Individual() {
	delete [] m_s;
}

int* Individual :: GetSequence() {
	return m_s;
}

int Individual :: GetGene(int i ) {
	return m_s[i];
}

void Individual :: SetSequence(int* x) {
	for (int i = 0;i<m_N;i++) 
		m_s[i] = x[i];
}

void Individual :: SetGene(int i, int x) {
	m_s[i] = x;
}
