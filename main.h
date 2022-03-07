#ifndef main_H
#define main_H

#include "uv_camera.h"
#include "custom_math.h"

#include <cstdlib>
#include <GL/glut.h>       //GLUT Library

#include <iostream>
using std::cout;
using std::endl;

#include <iomanip>
using std::setprecision;

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <sstream>
using std::ostringstream;
using std::istringstream;

#include <fstream>
using std::ofstream;
using std::ifstream;

#include <set>
using std::set;

#include <map>
using std::map;

#include <utility>
using std::pair;

#include <complex>
using namespace std;


void idle_func(void);
void init_opengl(const int &width, const int &height);
void reshape_func(int width, int height);
void display_func(void);
void keyboard_func(unsigned char key, int x, int y);
void mouse_func(int button, int state, int x, int y);
void motion_func(int x, int y);
void passive_motion_func(int x, int y);

void render_string(int x, const int y, void *font, const string &text);
void draw_objects(void);



const double speed_of_light = 1;
const double max_accel = speed_of_light * 2;
const double grav_constant = 1;
const double lj_attractive_constant = 1;
const double lj_repulsive_constant = 0.1;
const double lj_attractive_exponent = 2;
const double lj_repulsive_exponent = 4;
const size_t num_test_particles = 10000;

vector<custom_math::vector_3> julia_points;
vector<float> julia_points_mass;

vector<custom_math::vector_3> test_particle_pos;
vector<custom_math::vector_3> test_particle_vel;

vector<custom_math::vector_3> positions;

custom_math::vector_3 background_colour(1.0f, 1.0f, 1.0f);
custom_math::vector_3 control_list_colour(0.0f, 0.0f, 0.0f);

bool draw_axis = true;
bool draw_control_list = true;

uv_camera main_camera;

GLint win_id = 0;
GLint win_x = 800, win_y = 600;
float camera_w = 4;

float camera_fov = 45;
float camera_x_transform = 0;
float camera_y_transform = 0;
float u_spacer = 0.01;
float v_spacer = 0.5*u_spacer;
float w_spacer = 0.1;
float camera_near = 0.001;
float camera_far = 100;

bool lmb_down = false;
bool mmb_down = false;
bool rmb_down = false;
int mouse_x = 0;
int mouse_y = 0;






#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <ctime>
using namespace std;




float iterate_2d(vector< complex<float> >& trajectory_points,
    complex<float> Z,
    const complex<float> C,
    const short unsigned int max_iterations,
    const float threshold)
{
    trajectory_points.clear();
    trajectory_points.push_back(Z);

    for (short unsigned int i = 0; i < max_iterations; i++)
    {
        Z = Z * Z + C;
        //Z = sin(Z) + C * sin(Z);

        trajectory_points.push_back(Z);

        if (abs(Z) >= threshold)
            break;
    }

    return abs(Z);
}

int get_points(void)
{

    const unsigned short int max_iterations = 8;
    const float threshold = 4.0f;

    const float x_grid_max = 1.5;
    const float x_grid_min = -x_grid_max;
    const size_t x_res = 10;
    const complex<float> x_step_size((x_grid_max - x_grid_min) / (x_res - 1), 0);

    const float y_grid_max = 1.5;
    const float y_grid_min = -y_grid_max;
    const size_t y_res = 10;
    const complex<float> y_step_size(0, (y_grid_max - y_grid_min) / (y_res - 1));

    const complex<float> C(0.2f, 0.5f);

    complex<float> Z(x_grid_min, y_grid_min);

    vector< complex<float> > trajectory_points;

    float max_length = 0;

    for (size_t x = 0; x < x_res; x++, Z += x_step_size)
    {
        Z = complex<float>(Z.real(), y_grid_min);

        for (size_t y = 0; y < y_res; y++, Z += y_step_size)
        {
            float magnitude = iterate_2d(trajectory_points, Z, C, max_iterations, threshold);

            float length = 0.0f;

            for (size_t i = 0; i < trajectory_points.size() - 1; i++)
                length += abs(trajectory_points[i + 1] - trajectory_points[i]);

            if (length > max_length)
                    max_length = length;



            if (magnitude < threshold)
            {
                custom_math::vector_3 v;
                v.x = Z.real();
                v.y = Z.imag();
                v.z = 0;

                julia_points.push_back(v);
                julia_points_mass.push_back(length);
            }
        }
    }

    for (size_t i = 0; i < julia_points_mass.size(); i++)
        julia_points_mass[i] = 1.0f - julia_points_mass[i] / max_length;

	for (size_t i = 0; i < num_test_particles; i++)
	{
        float x = rand() / static_cast<float>(RAND_MAX);
        x *= 2;
        x -= 1;

        x *= 1.5f;

        float y = rand() / static_cast<float>(RAND_MAX);
        y *= 2;
        y -= 1;

        y *= 1.5f;

		test_particle_pos.push_back(custom_math::vector_3(x, y, 0));

        //x = rand() / static_cast<float>(RAND_MAX);
        //x *= 2;
        //x -= 1;

        //y = rand() / static_cast<float>(RAND_MAX);
        //y *= 2;
        //y -= 1;
        
        test_particle_vel.push_back(custom_math::vector_3(0, 0, 0));
	}



    return 0;
}




















#endif
