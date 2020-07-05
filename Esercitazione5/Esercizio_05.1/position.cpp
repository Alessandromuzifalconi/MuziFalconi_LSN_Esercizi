#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "position.h"


using namespace std;

Position::Position() {
	m_x = 0;
	m_y = 0;
	m_z = 0;
}

Position::Position(double x,double y,double z) {
	m_x = x;
	m_y = y;
	m_z = z;
}

Position :: ~Position() {}

void Position::SetX(double x) {
	m_x = x;
	return;
}

void Position::SetY(double x) {
	m_y = x;
	return;
}	

void Position::SetZ(double x) {
	m_z = x;
	return;
}

double Position::GetX() {
	return m_x;
}

double Position::GetY() {
	return m_y;
}

double Position::GetZ() {
	return m_z;
}

double Position::GetDist() {
	return sqrt(pow(m_x, 2) + pow(m_y, 2) + pow(m_z, 2));

}
