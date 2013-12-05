/* gears.c */
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <GL/glut.h>
#include "argon.hpp"

using namespace std;

#define NPARTICLES 125
#define STEP_MS 50
#define SCALE 3
static GLfloat particlePos[NPARTICLES][3];
static const GLfloat particleSize = 0.06;
static const GLfloat sphereSize = 2.3;
static const GLfloat particleColor[4] = {0.8, 0.5, 1., 0.0};
static const GLfloat outsideParticleColor[4] = {1.0, 0.0, 0.0, 0.0};
static const GLfloat sphereColor[4] = {1.0, 0.5, 0.2, 0.0};

static GLint particles[NPARTICLES];
static GLint sphere;
static GLfloat angle = 0.0;

static GLuint limit;
static GLuint count = 1;
static GLfloat pos[4] = {0.0, 0.0, 10.0, 0.0};

static GLfloat view_rotx = 90.0, view_roty = 0.0, view_rotz = 135.0;

static unsigned long int counter = 0;
static unsigned long int lastRefresh = 0;

vector<vect*> crystal;

void loadData(vector<vect*> &crystal)
{
  vect *ptr;
	char buffer[256];
	int nlCounter = 0, i,j;
	double x,y,z, trash;
	ifstream infile("avs.dat", ifstream::in);
	while(infile.good())
	{
		ptr = new vect[NPARTICLES];
		for(j = 0; j < NPARTICLES; ++j)
		{
			for(i=0; i < 256; ++i)
				buffer[i] = '\0';
			infile >> buffer;
			if(buffer[0] == '\0')
				continue;
			std::stringstream(buffer) >> ptr[j].x;
			infile >> ptr[j].y;
			infile >> ptr[j].z;
			infile >> trash >> trash >> trash; // momenta
			//cout << j << ":\t" << ptr[j].x << "\t" << ptr[j].y << "\t" << ptr[j].z << endl;
		}
		//cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		crystal.push_back(ptr);
	}
	infile.close();
}

static void draw(void)
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glPushMatrix();
  glRotatef(view_rotx, 1.0, 0.0, 0.0);
	glRotatef(view_roty, 0.0, 1.0, 0.0);
  glRotatef(view_rotz, 0.0, 0.0, 1.0);
  glCallList(sphere);
	for(int i = 0; i < NPARTICLES ; ++i)
	{
	  glPushMatrix();
    //if(i == 5)
    //  cout << i << ":::" << particlePos[i][0] << endl;
    // cout << particlePos[i][0] << "\t" << sqrt(particlePos[i][0]*particlePos[i][0] + particlePos[i][1]*particlePos[i][1] + particlePos[i][2]*particlePos[i][2]) << "\t" << sphereSize*SCALE << endl;
    // if(sqrt(particlePos[i][0]*particlePos[i][0] + particlePos[i][1]*particlePos[i][1] + particlePos[i][2]*particlePos[i][2]) > sphereSize*SCALE)
    //   glColor3f(outsideParticleColor[0], outsideParticleColor[1], outsideParticleColor[2]);
    // else
    //   glColor3f(particleColor[0], particleColor[1], particleColor[2]);
    glTranslatef(particlePos[i][0], particlePos[i][2], particlePos[i][3]);
	  glCallList(particles[i]);
	  glPopMatrix();
 	}
  glPopMatrix();
  glutSwapBuffers();
}

static void idle()
{
	unsigned long int time = glutGet(GLUT_ELAPSED_TIME);
	if(counter < crystal.size() - 4)
	{
		if(lastRefresh+STEP_MS < time)
		{
			++counter;
			lastRefresh = time;
		}
	}
	else
		counter = 0;
	for (int i=0;i<NPARTICLES;i++)
	{
	 // range -3.5...3.5
    if(counter > 0)
    {
      // particlePos[i][0] = (crystal[counter][i].x-crystal[counter-1][i].x)*SCALE;
      // particlePos[i][1] = (crystal[counter][i].y-crystal[counter-1][i].y)*SCALE;
      // particlePos[i][2] = (crystal[counter][i].z-crystal[counter-1][i].z)*SCALE;
      particlePos[i][0] = crystal[counter][i].x;
      particlePos[i][1] = crystal[counter][i].y;
      particlePos[i][2] = crystal[counter][i].z;
    }
    else
    {
      particlePos[i][0] = crystal[counter][i].x*SCALE;
      particlePos[i][1] = crystal[counter][i].y*SCALE;
      particlePos[i][2] = crystal[counter][i].z*SCALE;
    }
	}
	glutPostRedisplay();
}

static void key(unsigned char k, int x, int y)
{
  switch (k) {
  case 'z':
    view_rotz += 5.0;
    break;
  case 'Z':
    view_rotz -= 5.0;
    break;
  case 27:  /* Escape */
    exit(0);
    break;
  default:
    return;
  }
  glutPostRedisplay();
}

static void special(int k, int x, int y)
{
  switch (k) {
  case GLUT_KEY_UP:
    view_rotx += 5.0;
    break;
  case GLUT_KEY_DOWN:
    view_rotx -= 5.0;
    break;
  case GLUT_KEY_LEFT:
    view_roty += 5.0;
    break;
  case GLUT_KEY_RIGHT:
    view_roty -= 5.0;
    break;
  default:
    return;
  }
  glutPostRedisplay();
}

static void reshape(int width, int height)
{
  GLfloat h = (GLfloat) height / (GLfloat) width;

  glViewport(0, 0, (GLint) width, (GLint) height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glFrustum(-1.0, 1.0, -h, h, 5.0, 60.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(0.0, 0.0, -40.0);
}

static void init(void)
{
  glLightfv(GL_LIGHT0, GL_POSITION, pos);
  glEnable(GL_CULL_FACE);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_DEPTH_TEST);

	int t = 0;
  sphere = glGenLists(1);
  glNewList(sphere, GL_COMPILE);
  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, sphereColor);
  glTranslatef(0,0,0);
  t = 12 + (int) (sphereSize * SCALE);
  glutWireSphere(sphereSize * SCALE, t, t);
  glEndList();

  for(int i = 0; i < NPARTICLES; ++i)
  {
	  particles[i] = glGenLists(1);
	  glNewList(particles[i], GL_COMPILE);
	  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, particleColor);
    glTranslatef(particlePos[i][0], particlePos[i][1], particlePos[i][2]);
    t = 12 + (int) (particleSize * 12);
    glutSolidSphere(particleSize, t, t);
	  glEndList();
  }

  glEnable(GL_NORMALIZE);
}

void visible(int vis)
{
  if (vis == GLUT_VISIBLE)
    glutIdleFunc(idle);
  else
  {
  	cout << "NV" << endl;
    glutIdleFunc(NULL);
  }
}

main(int argc, char *argv[])
{
	loadData(crystal);

  glutInit(&argc, argv);
  glutInitWindowSize( 800, 800 );
  glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);

  glutCreateWindow("Argon");
  idle();
  init();

  glutDisplayFunc(draw);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(key);
  glutSpecialFunc(special);
  glutVisibilityFunc(visible);
  
  glutMainLoop();
  return 0;
}
