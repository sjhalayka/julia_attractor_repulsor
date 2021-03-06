#include "main.h"


int main(int argc, char **argv)
{
	cout << setprecision(20) << endl;

	get_points();

	glutInit(&argc, argv);
	init_opengl(win_x, win_y);
	glutReshapeFunc(reshape_func);
	glutIdleFunc(idle_func);
	glutDisplayFunc(display_func);
	glutKeyboardFunc(keyboard_func);
	glutMouseFunc(mouse_func);
	glutMotionFunc(motion_func);
	glutPassiveMotionFunc(passive_motion_func);
	//glutIgnoreKeyRepeat(1);
	glutMainLoop();
	glutDestroyWindow(win_id);

	return 0;
}

custom_math::vector_3 acceleration(const custom_math::vector_3& pos, const custom_math::vector_3& vel)
{
	custom_math::vector_3 total_accel;

	for (size_t i = 0; i < julia_points.size(); i++)
	{
		custom_math::vector_3 grav_dir = julia_points[i] - pos;

		float distance = grav_dir.length();
		grav_dir.normalize();

		if (julia_points[i].y > 0)
			grav_dir = -grav_dir;

		custom_math::vector_3 grav_accel = grav_dir * (lj_attractive_constant * julia_points_mass[i] / pow(distance, lj_attractive_exponent) - lj_repulsive_constant * julia_points_mass[i] / pow(distance, lj_repulsive_exponent));

		//custom_math::vector_3 ang_vel_dir;
		//ang_vel_dir.z = -magnetism_constant * julia_points_mass[i] / powf(distance, 2.0f);

		//custom_math::vector_3 magnus_accel = vel.cross(ang_vel_dir);

		total_accel += grav_accel;// +magnus_accel;
	}

	if (total_accel.length() > max_accel)
	{
		total_accel.normalize();
		total_accel *= max_accel;
	}

	return total_accel;
}



void proceed_Euler(custom_math::vector_3 &pos, custom_math::vector_3 &vel)
{
	custom_math::vector_3 accel = acceleration(pos, vel);

	vel += accel * dt;

	if (vel.length() > speed_of_light)
	{
		vel.normalize();
		vel *= speed_of_light;
	}

	pos += vel * dt;
}


void idle_func(void)
{
	if (add_trajectory_points)
	{
		for (size_t i = 0; i < test_particle_pos.size(); i++)
		{
			proceed_Euler(test_particle_pos[i], test_particle_vel[i]);
			positions.push_back(test_particle_pos[i]);
		}
	}

    glutPostRedisplay();
}

void init_opengl(const int &width, const int &height)
{
	win_x = width;
	win_y = height;

	if(win_x < 1)
		win_x = 1;

	if(win_y < 1)
		win_y = 1;

	glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(win_x, win_y);
	win_id = glutCreateWindow("orbit");

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glDepthMask(GL_TRUE);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	glClearColor(background_colour.x, background_colour.y, background_colour.z, 1);
	glClearDepth(1.0f);

	main_camera.Set(0, 0, camera_w, camera_fov, win_x, win_y, camera_near, camera_far);
}

void reshape_func(int width, int height)
{
	win_x = width;
	win_y = height;

	if(win_x < 1)
		win_x = 1;

	if(win_y < 1)
		win_y = 1;

	glutSetWindow(win_id);
	glutReshapeWindow(win_x, win_y);
	glViewport(0, 0, win_x, win_y);

	main_camera.Set(main_camera.u, main_camera.v, main_camera.w, main_camera.fov, win_x, win_y, camera_near, camera_far);
}

// Text drawing code originally from "GLUT Tutorial -- Bitmap Fonts and Orthogonal Projections" by A R Fernandes
void render_string(int x, const int y, void *font, const string &text)
{
	for(size_t i = 0; i < text.length(); i++)
	{
		glRasterPos2i(x, y);
		glutBitmapCharacter(font, text[i]);
		x += glutBitmapWidth(font, text[i]) + 1;
	}
}
// End text drawing code.

void draw_objects(void)
{
    glDisable(GL_LIGHTING);
    
	glPushMatrix();
  



	glEnable(GL_ALPHA);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	if (screenshot_mode)
		glPointSize(4.0);
	else
		glPointSize(1.0);
    
    glBegin(GL_POINTS);

    glColor4f(0.0, 0.0, 0.0, 1);
		
	for (size_t i = 0; i < positions.size(); i++)
	{
		if ((positions[i].x > -grid_max && positions[i].x < grid_max) &&
			(positions[i].y > -grid_max && positions[i].y < grid_max))
		{
		
		complex<float> Z = complex<float>(positions[i].x, positions[i].y);

		vector<complex<float>> trajectory_points;

		float magnitude = iterate_2d(trajectory_points, Z, C, max_iterations, threshold);

		if (1)//magnitude >= threshold)
			glVertex3f(positions[i].x, positions[i].y, positions[i].z);
		}
	}
    glEnd();
    
	glPointSize(4.0);

	glBegin(GL_POINTS);

	for (size_t i = 0; i < julia_points.size(); i++)
	{
		glColor4f(1.0f, 0.5f, 0.0f, pow(julia_points_mass[i], 10.0f));
		glVertex3f(julia_points[i].x, julia_points[i].y, julia_points[i].z);
	}

	glDisable(GL_BLEND);
	glDisable(GL_ALPHA);

	glEnd();
    
    
    
    glLineWidth(1.0f);
 //   
 //   
	//// If we do draw the axis at all, make sure not to draw its outline.
	//if(true == draw_axis)
	//{
	//	glBegin(GL_LINES);

	//	glColor3f(1, 0, 0);
	//	glVertex3f(0, 0, 0);
	//	glVertex3f(1, 0, 0);
	//	glColor3f(0, 1, 0);
	//	glVertex3f(0, 0, 0);
	//	glVertex3f(0, 1, 0);
	//	glColor3f(0, 0, 1);
	//	glVertex3f(0, 0, 0);
	//	glVertex3f(0, 0, 1);

	//	glColor3f(0.5, 0.5, 0.5);
	//	glVertex3f(0, 0, 0);
	//	glVertex3f(-1, 0, 0);
	//	glVertex3f(0, 0, 0);
	//	glVertex3f(0, -1, 0);
	//	glVertex3f(0, 0, 0);
	//	glVertex3f(0, 0, -1);

	//	glEnd();
	//}

	glPopMatrix();
}




void display_func(void)
{
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	// Draw the model's components using OpenGL/GLUT primitives.
	draw_objects();

	if(true == draw_control_list)
	{
		// Text drawing code originally from "GLUT Tutorial -- Bitmap Fonts and Orthogonal Projections" by A R Fernandes
		// http://www.lighthouse3d.com/opengl/glut/index.php?bmpfontortho
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		gluOrtho2D(0, win_x, 0, win_y);
		glScalef(1, -1, 1); // Neat. :)
		glTranslatef(0, -win_y, 0); // Neat. :)
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();

		glColor3f(control_list_colour.x, control_list_colour.y, control_list_colour.z);

		size_t break_size = 22;
		size_t start = 20;
		ostringstream oss;

		render_string(10, start, GLUT_BITMAP_HELVETICA_18, string("Mouse controls:"));
		render_string(10, start + 1*break_size, GLUT_BITMAP_HELVETICA_18, string("  LMB + drag: Rotate camera"));
		render_string(10, start + 2*break_size, GLUT_BITMAP_HELVETICA_18, string("  RMB + drag: Zoom camera"));

		render_string(10, start + 4*break_size, GLUT_BITMAP_HELVETICA_18, string("Keyboard controls:"));
        render_string(10, start + 5*break_size, GLUT_BITMAP_HELVETICA_18, string("  w: Draw axis"));
		render_string(10, start + 6*break_size, GLUT_BITMAP_HELVETICA_18, string("  e: Draw text"));
		render_string(10, start + 7*break_size, GLUT_BITMAP_HELVETICA_18, string("  u: Rotate camera +u"));
		render_string(10, start + 8*break_size, GLUT_BITMAP_HELVETICA_18, string("  i: Rotate camera -u"));
		render_string(10, start + 9*break_size, GLUT_BITMAP_HELVETICA_18, string("  o: Rotate camera +v"));
		render_string(10, start + 10*break_size, GLUT_BITMAP_HELVETICA_18, string("  p: Rotate camera -v"));


		
        custom_math::vector_3 eye = main_camera.eye;
		custom_math::vector_3 eye_norm = eye;
		eye_norm.normalize();

		oss.clear();
		oss.str("");		
		oss << "Camera position: " << eye.x << ' ' << eye.y << ' ' << eye.z;
		render_string(10, win_y - 2*break_size, GLUT_BITMAP_HELVETICA_18, oss.str());

		oss.clear();
		oss.str("");
		oss << "Camera position (normalized): " << eye_norm.x << ' ' << eye_norm.y << ' ' << eye_norm.z;
		render_string(10, win_y - break_size, GLUT_BITMAP_HELVETICA_18, oss.str());

		glPopMatrix();
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		// End text drawing code.
	}

	if (false == screenshot_mode)
	{
		glFlush();
		glutSwapBuffers();
	}
}

void keyboard_func(unsigned char key, int x, int y)
{
	switch(tolower(key))
	{
	case 'w':
		{
			draw_axis = !draw_axis;
			break;
		}
	case 'e':
		{
			draw_control_list = !draw_control_list;
			break;
		}
	case 'u':
		{
			main_camera.u -= u_spacer;
			main_camera.Set();
			break;
		}
	case 'i':
		{
			main_camera.u += u_spacer;
			main_camera.Set();
			break;
		}
	case 'o':
		{
			main_camera.v -= v_spacer;
			main_camera.Set();
			break;
		}
	case 'p':
		{
			main_camera.v += v_spacer;
			main_camera.Set();
			break;
		}
	case 'l':
	{
		add_trajectory_points = false;
		take_screenshot(8, "screenshot.tga");
		add_trajectory_points = true;
		
		break;
	}



	default:
		break;
	}
}

void mouse_func(int button, int state, int x, int y)
{
	if(GLUT_LEFT_BUTTON == button)
	{
		if(GLUT_DOWN == state)
			lmb_down = true;
		else
			lmb_down = false;
	}
	else if(GLUT_MIDDLE_BUTTON == button)
	{
		if(GLUT_DOWN == state)
			mmb_down = true;
		else
			mmb_down = false;
	}
	else if(GLUT_RIGHT_BUTTON == button)
	{
		if(GLUT_DOWN == state)
			rmb_down = true;
		else
			rmb_down = false;
	}
}

void motion_func(int x, int y)
{
	int prev_mouse_x = mouse_x;
	int prev_mouse_y = mouse_y;

	mouse_x = x;
	mouse_y = y;

	int mouse_delta_x = mouse_x - prev_mouse_x;
	int mouse_delta_y = prev_mouse_y - mouse_y;

	if(true == lmb_down && (0 != mouse_delta_x || 0 != mouse_delta_y))
	{
		main_camera.u -= static_cast<float>(mouse_delta_y)*u_spacer;
		main_camera.v += static_cast<float>(mouse_delta_x)*v_spacer;
	}
	else if(true == rmb_down && (0 != mouse_delta_y))
	{
		main_camera.w -= static_cast<float>(mouse_delta_y)*w_spacer;

		if(main_camera.w < 1.1f)
			main_camera.w = 1.1f;

	}

	main_camera.Set(); // Calculate new camera vectors.
}

void passive_motion_func(int x, int y)
{
	mouse_x = x;
	mouse_y = y;
}

