/*******************************************************************/
//OpenGL Solution Source code, V1.0 to display and manipulate a Bezier curve of arbitrary degree
//Anurag Purwar, MEC572, Spring 2015
/*******************************************************************/

#include <iostream>
#include <windows.h>
//#include <GL/glui.h>
//#include <freeglut.h>

//#include <GL/glut.h
#include <GL/freeglut.h>
#include "glui.h"
#include "point2D.h"
#include "curve.h"
#include <vector>
#include <cmath>
#include "utility.h"
    
using namespace std;

int DIV=20.0;

Curve *C;
bool Bernstein = false, Bspline = false;
Point2D point;
int p = 2;

/*Uncomment following line to get rid of console window
 and display only graphics window */
//#pragma comment( linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"" )

int width, height;

int win_id;            // each curve on each screen
int showCtrlPts=1;	   //Switch to toggle between showing ctrl pts and hiding

void mouseButtonEvents(int, int, int, int);
void mouseMotionEvents(int, int);
void mousePassiveMotion(int, int);
void keyboard(unsigned char, int, int);
void draw(int);
void display(void);
void reshape(int, int);
void shellmessages(void);
void selectMessage(int);
void init(void);



//OpenGL initialization function. DNT.
void init(void)
{
   shellmessages();
   glClearColor(.0, .0, .0, .0);
   glShadeModel(GL_FLAT);
   glEnable(GL_POINT_SMOOTH);
   glEnable(GL_LINE_SMOOTH);
   glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);	// Make round points, not square points
   glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);		// Antialias the lines

}


void shellmessages(void)
{
	cout<<"\n* Single Left Button Click  - Adds a Point\n\n"
		<<"* <Ctrl>+Single Left Button click on an existing Point - Removes that Point\n\n"
		<<"* Right Click an existing Point and move mouse - Drags that Point\n\n"
		<<"* Middle Button Click - Menu for clearing screen, etc.\n\n";
}


void display(void)
{
	
	glClear (GL_COLOR_BUFFER_BIT);
 	glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
	glColor3f (.0, .0, .0);
	glutSetWindow(win_id);
    
	draw(GL_RENDER);

	glutSwapBuffers();

}

Point2D deCasteljau(float t)
{
	int i, n, r;
	vector<Point2D> b;

	n = C->B.size()-1; //degree, one less than number of points

	for (i=0; i<C->B.size();i++)
		b.push_back(C->B[i]);

	for (r=1; r<=n; r++)
		for (i=0; i<=n-r; i++)
			{
			b[i].x = (1.0-t)*b[i].x + t*b[i+1].x;
			b[i].y = (1.0-t)*b[i].y + t*b[i+1].y;
			
			}

	return b[0];
}

//vector<double> createU(int degree, vector<Point2D> P){
//	vector<double> U;
//	int p = degree;
//	int m = P.size() + p;
//
//	for (int i = 0; i < p + 1; i++){
//		U.push_back(0.0);
//	}
//
//	double index = 1.0;
//	for (int i = p; i < m - p - 1; i++){
//		U.push_back(index / (double)(p + 1));
//		index++;
//	}
//
//	for (int i = 0; i < p + 1; i++){
//		U.push_back(1.0);
//	}
//	return U;
//}

vector<double> createU(int degree, vector<Point2D> P){
	vector<double> U;
	int p = degree;
	int m = P.size() + p;

	for (int i = 0; i <= p; i++){
		U.push_back(0.0);
	}

	double index = 1.0, tmp = 0.0, old_tmp;
	for (double j = 1; j < m - 2 * p - 1; j++){
		tmp = j / (double)(m - 2 * p - 1);
		U.push_back(tmp);

	}

	for (int i = m - p; i <= m; i++){
		U.push_back(1.0);
	}
	return U;
}

int findSpan(int degree, double u, vector<double> U){
	int m = U.size() - 1;
	int p = degree;
	int n = m - p - 1;
	if (u == U[n + 1]) return(n);
	n = U.size() - 1;
	int low = p, high = n + 1;
	int mid = (low + high) / 2;

	while (u < U[mid] || u > U[mid + 1]){
		if (u < U[mid])
			high = mid;
		else
			low = mid;

		mid = (low + high) / 2;
	}
	return (mid);
}

vector<double> BasisFunc(int span, int degree, double u, vector<double> U){
	int p = degree;
	vector<double> N(p + 1);
	vector<double> right(p + 1), left(p + 1);
	double saved, temp;
	N[0] = 1.0;

	for (int j = 1; j <= p; j++){
		left[j] = u - U[span + 1 - j];
		right[j] = U[span + j] - u;
		saved = 0.0;
		for (int r = 0; r < j; r++){
			temp = N[r] / (right[r + 1] + left[j - r]);
			N[r] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		N[j] = saved;
	}
	return N;
}

Point2D curvePoint(int degree, vector<Point2D> P, double u){
	vector<double> U, N;
	int p = degree;
	U = createU(p, P);
	int span = findSpan(p, u, U);
	N = BasisFunc(span, p, u, U);

	Point2D cPoint;
	cPoint.x = 0.0;
	cPoint.y = 0.0;
	for (int i = 0; i <= p; i++){
		if (span - p + i >= P.size()){
			break;
		}
		cPoint.x = cPoint.x + N[i] * P[span - p + i].x;
		cPoint.y = cPoint.y + N[i] * P[span - p + i].y;
	}
	return cPoint;
}

vector<Point2D> Allbernstein(Curve P, Curve Q, double u) {
	Point2D one;
	one.x = 1.0;
	one.y = 1.0;
	Point2D saved, temp;
	for (size_t i = 0; i < P.B.size(); i++)
		Q.B.push_back(one);

	for (size_t i = 1; i < P.B.size(); i++) {
		//if (i + 1 > P.B.size()) break;
		saved.x = 0.0;
		saved.y = 0.0;
		for (size_t j = 0; j < i; j++) {
			temp.x = Q.B[j].x;
			temp.y = Q.B[j].y;
			Q.B[j].x = saved.x + (1.0 - u)*temp.x;
			Q.B[j].y = saved.y + (1.0 - u)*temp.y;
			saved.x = u*temp.x;
			saved.y = u*temp.y;
		}
		Q.B[i].x = saved.x;
		Q.B[i].y = saved.y;
	}
	return Q.B;
}

Point2D PointOnBezierCurve(Curve P, Point2D t, double u) {
	Curve Q;
	Q.B = Allbernstein(P, Q, u);

	t.x = 0.0;
	t.y = 0.0;
	//C.push_back(t);


	for (size_t i = 0; i < Q.B.size(); i++) {
		//if (i + 1 > Q.B.size()|| i+1>P.B.size()) break;
		t.x = t.x + Q.B[i].x * P.B[i].x;
		t.y = t.y + Q.B[i].y * P.B[i].y;
		//C.push_back(t);
	}

	return t;
}
/* This is the function that needs to be changed for drawing anything else. */
void draw(int mode)
{
	int i;
	float t;

// The following lines clubbed together draw all the points
	glColor3f (0.0, 1.0, 1.0);
	glPointSize(5.0);
	glBegin(GL_POINTS);
	for (i=0; i<C->B.size(); i++)
	glVertex2f(C->B[i].x, C->B[i].y);
	glEnd();

// The following lines draw the boundary of the control polygon 
//	except last segment
    glLineWidth(1.5);

	if (showCtrlPts)
	{
	glColor3f (1.0, 1.0, 1.0);
	glBegin(GL_LINE_STRIP);
	for (i=0; i<C->B.size(); i++)
	glVertex2f(C->B[i].x, C->B[i].y);
	glEnd();
	}

	Point2D P; //a temporary point 

	if (C->B.size() > 0)
	{
		//DIV = (C->B.size() > 4)? C->B.size()*DIV/4:DIV;

		glColor3f (1.0, 1.0, 0.0);
		glBegin(GL_LINE_STRIP);
		for (t=0; t<1.0; t+=1.0/DIV)
			{
			P =	deCasteljau(t);
			glVertex2f(P.x, P.y);
			} //for
		glVertex2f(C->B[C->B.size()-1].x, C->B[C->B.size()-1].y);
		glEnd();
	} //if

	if (C->B.size()>0 && Bspline == true){
		glColor3f(00.294, 0.000, 0.510);
		glBegin(GL_LINE_STRIP);
		for (t = 0; t<=1.0; t += 1.0 / DIV)
		{
			P = curvePoint(p,C->B,t);
			glVertex2f(P.x, P.y);
		} //for
		
		glEnd();
	}

	
	if (Bernstein){
		glColor3f(0.0, 1.0, 0.0);

		glBegin(GL_LINE_STRIP);
		double u = 0.0;
		while (u <= 1.0) {
			 point = PointOnBezierCurve(*C, point, u);
			//for (size_t i = 0; i < C->B.size(); i++)
			glVertex2f(point.x, point.y);
			u += 1.0/DIV;
		}
		glEnd();
		glFlush();


	}

}


// When you resize/create the window, this gets called. DNT.
void reshape(int w, int h)
{
	width  = w;
	height = h;

	glViewport(0,0,width,height) ;
    glMatrixMode(GL_PROJECTION) ;
    glLoadIdentity() ;
    glOrtho(.0, w, .0, h, -1.0, 1.0) ;

}


void selectMessage(int msg)
{
  switch (msg) {
  case 1:
    if (! C->B.empty() )
		C->B.clear();
	glutPostRedisplay();
    break;
  case 2:
	glutPositionWindow(50, 50);
    glutReshapeWindow ( 600, 400 );
	break;
  case 3:
	  glutFullScreen();
	  break;
  case 4:
	  showCtrlPts = !showCtrlPts;
	  break;
  case 5:
	  Bernstein = !Bernstein;
	  break;
  case 6:
	  Bspline = !Bspline;
	  break;
  case 7:
	  p++;
	  break;
  case 8:
	  if (p>0)
	  p--;
	  break;
  case 10:
	  exit(0);
	  break;
  default:
	  break;
  }
  glutPostRedisplay();
}

/** These are the live variables passed into GLUI ***/
int   wireframe = 0;
int   segments = 8;
int   main_window;


/***************************************** myGlutIdle() ***********/

void myGlutIdle(void)
{
	/* According to the GLUT specification, the current window is
	undefined during an idle callback.  So we need to explicitly change
	it if necessary */
	if (glutGetWindow() != main_window)
		glutSetWindow(main_window);

	glutPostRedisplay();
}


void idle()
{
	glutSetWindow(main_window);
	glutPostRedisplay();
	Sleep(100);
}



int main(int argc, char** argv)
{
	Curve C1;
	C = &C1;

	

   glutInit(&argc, argv);
   glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA);
   glutInitWindowSize (600, 400);
   glutInitWindowPosition (50, 50);
   win_id = glutCreateWindow("CAGD, MEC572 Spring 2015");
   init();
   glutDisplayFunc(display);
   glutReshapeFunc(reshape);
   glutKeyboardFunc(keyboard);
   glutMouseFunc(mouseButtonEvents);
   glutMotionFunc(mouseMotionEvents);
   glutPassiveMotionFunc(mousePassiveMotion);

   /**********************************************************/
   // Glui gui //
   //main_window = glutCreateWindow("GLUI");
 /*  GLUI *glui = GLUI_Master.create_glui("GLUI");
   glui->add_checkbox("Wireframe", &wireframe);
   GLUI_Spinner *segment_spinner =
	   glui->add_spinner("Segments:", GLUI_SPINNER_INT, &segments);
   segment_spinner->set_int_limits(3, 60);

   glui->set_main_gfx_window(win_id);*/

   /* We register the idle callback with GLUI, *not* with GLUT */
   //GLUI_Master.set_glutIdleFunc(idle);

  int glut_menu = glutCreateMenu(selectMessage);
  glutAddMenuEntry("Clear Screen", 1);
  glutAddMenuEntry("Normal Window", 2);
  glutAddMenuEntry("Full Screen", 3);
  glutAddMenuEntry("Show/Hide CH", 4);
  glutAddMenuEntry("Show/Hide Bernstein", 5);
  glutAddMenuEntry("Show/Hide B-spline", 6);
  
  glutAddMenuEntry("Increase degree +1", 7);
  glutAddMenuEntry("Decrease degree -1", 8);
  glutAddMenuEntry("Exit", 10);
  glutAttachMenu(GLUT_MIDDLE_BUTTON);





   glutMainLoop();


   return 0;
}

