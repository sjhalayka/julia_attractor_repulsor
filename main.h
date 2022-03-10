#ifndef main_H
#define main_H

#include "uv_camera.h"
#include "custom_math.h"
#include "image.h"

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






const unsigned short int max_iterations = 8;
const float threshold = 4.0f;

const float x_grid_max = 1.5;
const float x_grid_min = -x_grid_max;
const size_t x_res = 30;
const complex<float> x_step_size((x_grid_max - x_grid_min) / (x_res - 1), 0);

const float y_grid_max = 1.5;
const float y_grid_min = -y_grid_max;
const size_t y_res = 30;
const complex<float> y_step_size(0, (y_grid_max - y_grid_min) / (y_res - 1));

const complex<float> C(0.5f, 0.5f); //  C(0.2f, 0.5f);

static const double dt = 0.01;


const double speed_of_light = 1;
const double max_accel = speed_of_light * 2;
const double grav_constant = 1;
const double lj_attractive_constant = 1;
const double lj_repulsive_constant = 0.01;
const double lj_attractive_exponent = 2;
const double lj_repulsive_exponent = 4;
const double magnetism_constant = 1;
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

bool screenshot_mode = false;

bool add_trajectory_points = true;

void idle_func(void);
void init_opengl(const int& width, const int& height);
void reshape_func(int width, int height);
void display_func(void);
void keyboard_func(unsigned char key, int x, int y);
void mouse_func(int button, int state, int x, int y);
void motion_func(int x, int y);
void passive_motion_func(int x, int y);

void render_string(int x, const int y, void* font, const string& text);
void draw_objects(void);


#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <ctime>
using namespace std;




float iterate_2d(vector< complex<float> >& trajectory_points,
    complex<float> Z,
    complex<float> C,
    const short unsigned int max_iterations,
    const float threshold)
{
	// Uncomment for Mandelbrot set
	//C = Z;
	//Z = complex<float>(0, 0);

    trajectory_points.clear();
    trajectory_points.push_back(Z);

    for (short unsigned int i = 0; i < max_iterations; i++)
    {
		Z = Z * Z + C;
        //Z = Z * Z * Z * Z + C;
        //Z = sin(Z) + C * sin(Z);

        trajectory_points.push_back(Z);

        if (abs(Z) >= threshold)
            break;
    }

    return abs(Z);
}




int get_points(void)
{
	//tga tga_texture;
	//float_grayscale luma;

	//if (false == convert_tga_to_float_grayscale("cat.tga", tga_texture, luma, true, true, true))
	//	return 0;

	//double grid_x_pos = -1.5; // Start at minimum x.
	//double grid_y_pos = 1.5;// grid_y_max; // Start at maximum y.

	//float step_size = 3.0 / static_cast<double>(luma.px);

	//// Begin march.
	//for (short unsigned int y = 0; y < luma.py; y++, grid_y_pos -= step_size, grid_x_pos = -1.5f)
	//{
	//	for (short unsigned int x = 0; x < luma.px; x++, grid_x_pos += step_size)
	//	{
	//		float pixel_datum = 0;

	//		if (luma.pixel_data[y * luma.px + x] >= 0.5)
	//			pixel_datum = 1;

	//		custom_math::vector_3 v;
	//		v.x = grid_x_pos;
	//		v.y = grid_y_pos;
	//		v.z = 0;

	//		julia_points.push_back(v);
	//		julia_points_mass.push_back(pixel_datum);
	//	}
	//}

	//for (size_t i = 0; i < num_test_particles; i++)
	//{
	//	float x = rand() / static_cast<float>(RAND_MAX);
	//	x *= 2;
	//	x -= 1;

	//	x *= 1.5f;

	//	float y = rand() / static_cast<float>(RAND_MAX);
	//	y *= 2;
	//	y -= 1;

	//	y *= 1.5f;

	//	test_particle_pos.push_back(custom_math::vector_3(x, y, 0));

	//	//x = rand() / static_cast<float>(RAND_MAX);
	//	//x *= 2;
	//	//x -= 1;

	//	//y = rand() / static_cast<float>(RAND_MAX);
	//	//y *= 2;
	//	//y -= 1;

	//	test_particle_vel.push_back(custom_math::vector_3(0, 0, 0));
	//}



	//return 0;




 

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


// TODO: fix camera bug where portrait mode crashes.
void take_screenshot(size_t num_cams_wide, const char* filename, const bool reverse_rows = false)
{
	screenshot_mode = true;

	// Set up Targa TGA image data.
	unsigned char  idlength = 0;
	unsigned char  colourmaptype = 0;
	unsigned char  datatypecode = 2;
	unsigned short int colourmaporigin = 0;
	unsigned short int colourmaplength = 0;
	unsigned char  colourmapdepth = 0;
	unsigned short int x_origin = 0;
	unsigned short int y_origin = 0;

	cout << "Image size: " << static_cast<size_t>(win_x) * num_cams_wide << "x" << static_cast<size_t>(win_y) * num_cams_wide << " pixels" << endl;

	if (static_cast<size_t>(win_x) * num_cams_wide > static_cast<unsigned short>(-1) ||
		static_cast<size_t>(win_y) * num_cams_wide > static_cast<unsigned short>(-1))
	{
		cout << "Image too large. Maximum width and height is " << static_cast<unsigned short>(-1) << endl;
		return;
	}

	unsigned short int px = win_x * static_cast<unsigned short>(num_cams_wide);
	unsigned short int py = win_y * static_cast<unsigned short>(num_cams_wide);
	unsigned char  bitsperpixel = 24;
	unsigned char  imagedescriptor = 0;
	vector<char> idstring;

	size_t num_bytes = 3 * px * py;
	vector<unsigned char> pixel_data(num_bytes);

	// Adjust some parameters for large screen format.
	bool temp_draw_control_list = draw_control_list;
	draw_control_list = false;

	vector<unsigned char> fbpixels(3 * win_x * win_y);

	const size_t total_cams = num_cams_wide * num_cams_wide;
	size_t cam_count = 0;
	// Loop through subcameras.
	for (size_t cam_num_x = 0; cam_num_x < num_cams_wide; cam_num_x++)
	{
		for (size_t cam_num_y = 0; cam_num_y < num_cams_wide; cam_num_y++)
		{
			cout << "Camera: " << cam_count + 1 << " of " << total_cams << endl;

			// Set up camera, draw, then copy the frame buffer.
			main_camera.Set_Large_Screenshot(num_cams_wide, cam_num_x, cam_num_y);
			display_func();
			glReadPixels(0, 0, win_x, win_y, GL_RGB, GL_UNSIGNED_BYTE, &fbpixels[0]);

			// Copy pixels to large image.
			for (GLint i = 0; i < win_x; i++)
			{
				for (GLint j = 0; j < win_y; j++)
				{
					size_t fb_index = 3 * (j * win_x + i);

					size_t screenshot_x = cam_num_x * win_x + i;
					size_t screenshot_y = cam_num_y * win_y + j;
					size_t screenshot_index = 3 * (screenshot_y * (win_x * num_cams_wide) + screenshot_x);

					pixel_data[screenshot_index] = fbpixels[fb_index + 2];
					pixel_data[screenshot_index + 1] = fbpixels[fb_index + 1];
					pixel_data[screenshot_index + 2] = fbpixels[fb_index];
				}
			}

			cam_count++;
		}

	}

	screenshot_mode = false;

	// Restore the parameters.
	draw_control_list = temp_draw_control_list;

	main_camera.Set();

	// Write Targa TGA file to disk.
	ofstream out(filename, ios::binary);

	if (!out.is_open())
	{
		cout << "Failed to open TGA file for writing: " << filename << endl;
		return;
	}

	out.write(reinterpret_cast<char*>(&idlength), 1);
	out.write(reinterpret_cast<char*>(&colourmaptype), 1);
	out.write(reinterpret_cast<char*>(&datatypecode), 1);
	out.write(reinterpret_cast<char*>(&colourmaporigin), 2);
	out.write(reinterpret_cast<char*>(&colourmaplength), 2);
	out.write(reinterpret_cast<char*>(&colourmapdepth), 1);
	out.write(reinterpret_cast<char*>(&x_origin), 2);
	out.write(reinterpret_cast<char*>(&y_origin), 2);
	out.write(reinterpret_cast<char*>(&px), 2);
	out.write(reinterpret_cast<char*>(&py), 2);
	out.write(reinterpret_cast<char*>(&bitsperpixel), 1);
	out.write(reinterpret_cast<char*>(&imagedescriptor), 1);

	out.write(reinterpret_cast<char*>(&pixel_data[0]), num_bytes);
}










#endif
