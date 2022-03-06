

#ifndef custom_math_h
#define custom_math_h

#include <algorithm>
using std::sort;

#include <limits>
using std::numeric_limits;

#include <vector>
using std::vector;

#include <set>
using std::set;

#include <map>
using std::map;


#include <iostream>
using std::cout;
using std::endl;

#include <cmath>
#include <cstdlib>
#include <GL/glut.h>       //GLUT Library


namespace custom_math
{
	class vector_3;
	class vector_4;
	const double pi = 3.14159265358979323846;
	const double pi_half = pi/2;
	const double pi_2 = 2*pi;
	const double epsilon = 1e-6;

	double d(const double &a, const double &b);
	double d_3(const vector_3 &a, const vector_3 &b);
	double d_3_sq(const vector_3 &a, const vector_3 &b);
	double d_4(const vector_4 &a, const vector_4 &b);
};

class custom_math::vector_3
{
public:
	double x, y, z;

	vector_3(const double &src_x = 0, const double &src_y = 0, const double &src_z = 0);
	bool operator==(const vector_3 &rhs);
	bool operator!=(const vector_3 &rhs);
	void zero(void);
	void rotate_x(const double &radians);
	void rotate_y(const double &radians);
    void rotate_z(const double &radians);
    vector_3 operator+(const vector_3 &rhs);
	vector_3 operator-(const vector_3 &rhs);
	vector_3 operator*(const vector_3 &rhs);
	vector_3 operator*(const double &rhs);
	vector_3 operator/(const double &rhs);
	vector_3 &operator=(const vector_3 &rhs);
	vector_3 &operator+=(const vector_3 &rhs);
	vector_3 &operator*=(const vector_3 &rhs);
	vector_3 &operator*=(const double &rhs);
	vector_3 operator-(void);
	double length(void) const;
	vector_3 &normalize(void);
	double dot(const vector_3 &rhs) const;
	double self_dot(void) const;
	vector_3 cross(const vector_3 &rhs) const;
};


class custom_math::vector_4
{
public:
	double x, y, z, w;

	vector_4(const double &src_x = 0, const double &src_y = 0, const double &src_z = 0, const double &src_w = 0);
	void zero(void);
	vector_4 operator+(const vector_4 &rhs);
	vector_4 operator-(const vector_4 &rhs);
	vector_4 operator*(const vector_4 &rhs);
	vector_4 operator*(const double &rhs);
	vector_4 operator/(const double &rhs);
	vector_4 &operator=(const vector_4 &rhs);
	vector_4 &operator+=(const vector_4 &rhs);
	vector_4 &operator*=(const vector_4 &rhs);
	vector_4 &operator*=(const double &rhs);
	vector_4 operator-(void);
	double length(void) const;
	vector_4 &normalize(void);
	double dot(const vector_4 &rhs) const;
	double self_dot(void) const;
};



#endif

