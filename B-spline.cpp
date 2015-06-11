// This is a template for assignment 3 of DM6122: 3D Modeling and Reconstruction.
// NTU
// August 2008
//
// Open a new project "Win32 Console Application" and add sample.c to Source Files 
#include <stdlib.h>
#include <stdio.h>
#include <GL/glut.h>
#include <GL/glu.h>
#include <GL/glui.h>
#include <math.h>
#include <malloc.h>



typedef struct {
	double x, y;      // x, y coordinates of a 2D point
} Point2d;


typedef struct {
	int      degree;   // degree of the B-spline curve
	int      cntNum;   // number of the deBoor points of the B-spline curve
	double   *knots;   // knot vector of the B-spline curve
	Point2d  *cnt;     // control points of the B-spline curve
} Bspline;


// B-spline curve
Bspline  bcr;

// The parameters are used to define a visible window in this application's World Coordinate System. 
double	winLLx = 0.0;
double  winLLy = 0.0;
double	winLen = 100.0;

// 
int    displayCP     = 1;
int    adaptivePlot  = 1;
int    samplingPnt   = 0;
int    tessNum       = 10;
double tessEps     = 2.;

//============================================================
static void Init(void)
{
	glClearColor (1.0, 1.0, 1.0, 0.0);         // set display-window color to white
    glMatrixMode (GL_PROJECTION);	           // set projection parameters
	gluOrtho2D (winLLx, winLLx+winLen, winLLy, winLLy+winLen);   // set an orthogonal projection
}

Point2d deBoor (int k, int h, int i, double t, double *knots, Point2d  *cnt){
	Point2d p;
	p.x=0; p.y=0;
	if (h==0){
		p=cnt[i];
	}
	else{
		p.x=( 1-(t-knots[i])/(knots[i+k+1-h]-knots[i]) )*deBoor(k,h-1,i-1,t,knots,cnt).x + (t-knots[i])/(knots[i+k+1-h]-knots[i])*deBoor(k,h-1,i,t,knots, cnt).x;
		p.y=( 1-(t-knots[i])/(knots[i+k+1-h]-knots[i]) )*deBoor(k,h-1,i-1,t,knots,cnt).y + (t-knots[i])/(knots[i+k+1-h]-knots[i])*deBoor(k,h-1,i,t,knots, cnt).y;
	}
	return p;
}

//============================================================
void uniformRender()
{
	// TODO:  add your own codes
	Point2d p;
	p.x =0; p.y=0;
	double t;
	glBegin(GL_LINE_STRIP);      // display the curve points
	double t0 = bcr.knots[bcr.degree]; // determine the t value of the first point
	p = deBoor(bcr.degree,bcr.degree,bcr.degree,t0,bcr.knots, bcr.cnt); // de Boor algorithm
	glVertex2f(p.x, p.y); 
	double tmax = bcr.knots[bcr.cntNum];  // determine the t value of the last point
	for (int index=1; index<tessNum-1; index++){
		t = t0 + index*(tmax-t0)/(tessNum-1); // insert points
		for(int j=bcr.degree;j<bcr.cntNum;j++){
			if (t>=bcr.knots[j] && t<=bcr.knots[j+1]){ // find the section in which t belongs to
				p = deBoor(bcr.degree,bcr.degree,j,t,bcr.knots, bcr.cnt);
				glVertex2f(p.x, p.y);
				break;
			}
		}
	}
	p = deBoor(bcr.degree,bcr.degree,bcr.cntNum-1,tmax,bcr.knots, bcr.cnt); // get the coordinate of the last point
	glVertex2f(p.x, p.y);
	glEnd();
}

//============================================================
void extractBezier (Point2d* bez, int ind)
{
	int     i, j;
	int     k;
	double  knots[50];
	Point2d cnt[30];

	k = bcr.degree;

	// copy one segment
	for (i=ind-k, j=0; i<=ind; i++) {
		cnt[j].x = bcr.cnt[i].x;
		cnt[j].y = bcr.cnt[i].y;
		j++;
	}
	for (i=ind-k, j=0; i<= ind+k+1; i++) {
		knots[j] = bcr.knots[i];
		j++;
	}

	// insert knots to make the left end be Bezier end
	while(1) {
		for (i=k-1; i>0; i--) {
			if (knots[i] < knots[k]) {
				j = i;
				break;
			}
			j = 0;
		}

		if(j==0) break;

		// update control points
		for (i=0; i<j; i++) {
			cnt[i].x = ((knots[k+1+i]-knots[k])/(knots[k+i+1]-knots[i+1]))*cnt[i].x 
					  + ((knots[k]-knots[i+1])/(knots[k+i+1]-knots[i+1]))*cnt[i+1].x;
			cnt[i].y = ((knots[k+1+i]-knots[k])/(knots[k+i+1]-knots[i+1]))*cnt[i].y 
					  + ((knots[k]-knots[i+1])/(knots[k+i+1]-knots[i+1]))*cnt[i+1].y;
		}
		// update knots
		for (i=0; i<j; i++)
			knots[i] = knots[i+1];
		knots[j] = knots[k];
	}

	// insert knots to make the right end be Bezier end
	while(1) {
		for (i=k+2; i< k+k+1; i++) {
			if (knots[i] > knots[k+1]) {
				j = i;
				break;
			}
			j = 0;
		}

		if(j==0) break;

		// update control points
		for (i=k; i>=j-k; i--) {
			cnt[i].x = ((knots[k+i]-knots[k+1])/(knots[k+i]-knots[i]))*cnt[i-1].x 
					  + ((knots[k+1]-knots[i])/(knots[k+i]-knots[i]))*cnt[i].x;
			cnt[i].y = ((knots[k+i]-knots[k+1])/(knots[k+i]-knots[i]))*cnt[i-1].y 
					  + ((knots[k+1]-knots[i])/(knots[k+i]-knots[i]))*cnt[i].y;
		}
		// update knots
		for (i=k+k+1; i>j; i--)
			knots[i] = knots[i-1];
		knots[j] = knots[k+1];
	}

	// return the Bezier control points
	for (i=0; i< bcr.cntNum; i++) {
		bez[i].x = cnt[i].x;
		bez[i].y = cnt[i].y;
	}
}

double getError(Point2d* bez, int deg){
	// bound of the flatness of the curve
	double hmax=0;
	double h=0;
	for (int i=1; i<deg; i++){
		double p0pn = sqrt( (bez[deg].x-bez[0].x)*(bez[deg].x-bez[0].x)+(bez[deg].y-bez[0].y)*(bez[deg].y-bez[0].y) );
		double crossproduct = (bez[deg].x-bez[0].x)*(bez[i].y-bez[0].y)-(bez[i].x-bez[0].x)*(bez[deg].y-bez[0].y);
		h = abs( crossproduct )/p0pn;
		if (h>hmax)
			hmax = h;
	}
	return hmax;
}

Point2d** deCasteljau (int deg, double t, Point2d* bez){
	Point2d p;
	p.x = 0; p.y =0;
	Point2d* pi = (Point2d *) malloc ((deg+1)*sizeof(Point2d));
	Point2d** pki = (Point2d **) malloc ((deg+1)*sizeof(Point2d*));
	for (int i=0; i<=deg; i++){
		pi[i] = bez[i];
		pki[0] = pi;
	}
	for (int k=1; k<=deg; k++){
		Point2d* newpi = (Point2d *) malloc ((deg+1)*sizeof(Point2d));
		for (int i=deg-k; i>=0; i--){
			p.x = (1-t)*pki[k-1][i].x + t*pki[k-1][i+1].x;
			p.y = (1-t)*pki[k-1][i].y + t*pki[k-1][i+1].y;
			newpi[i] = p;
            pki[k] = newpi;
		}
	}
	return pki;
}

void subdivision (Point2d* bez, int deg)
{
	//subdivide a bezier curve into two curves
	Point2d p;
	p.x = 0; p.y = 0;
	Point2d** pki = (Point2d **) malloc ((deg+1)*(deg+1)*sizeof(Point2d));
	Point2d* bez1 = (Point2d *) malloc ((deg+1)*sizeof(Point2d));
	Point2d* bez2 = (Point2d *) malloc ((deg+1)*sizeof(Point2d));
	double t;
	t = 0.5;
	double error = getError(bez, deg);
	// if the maximal distance is greater than a threshold, the curve is to be subdivided into two parts
	if ( error>tessEps){
		pki = deCasteljau(deg, t, bez); // de Casteljau algorithm
		for (int i=0; i<=deg; i++){
			bez1[i] = pki[i][0]; // the left sub-curve
			bez2[i] = pki[deg-i][i]; // the right sub-curve
		}
		subdivision(bez1, deg);
		subdivision(bez2, deg);
	}
	// else, draw a line from the first control point and the last one
	else{
		glVertex2f(bez[0].x, bez[0].y);
		glVertex2f(bez[deg].x, bez[deg].y);
		return;
	}
	
}

//============================================================
void plotBezier(Point2d* bez, int deg)
{
	// TODO:  add your own codes
	glBegin(GL_LINE_STRIP);
	subdivision(bez, deg);	//subdivide a bezier curve into two curves
	glEnd();
}


//============================================================
void adaptiveRender()
{
	Point2d  bez[30];  // assume the degree is not greater than 29.
	int      i;

	for (i=bcr.degree; i< bcr.cntNum; i++) {
		if (fabs(bcr.knots[i]-bcr.knots[i+1]) < 0.00001) continue;  // no segment, skip over
		extractBezier (bez, i);        // extract the i-th Bezier curve
		plotBezier(bez, bcr.degree);   // adaptively plot a Bezier curve 
	}
}



//============================================================
static void drawCurve( void )
{
	int i;
    glClear(GL_COLOR_BUFFER_BIT);	// clear display window
	glColor3f (1.0,0.0,0.0);   // set line segment color to red


// Draw the control polygon
	glColor3f(1.0, 0.0, 0.0);
	glLineWidth(3.0);
	if (displayCP != 0) {
		glBegin(GL_LINE_STRIP);      // display the control polygon
			for (i=0; i<bcr.cntNum; i++)
				glVertex2f(bcr.cnt[i].x, bcr.cnt[i].y);
		glEnd();
		
		glPointSize(6.0);            // display the control points
		glBegin(GL_POINTS);
			for (i=0; i<bcr.cntNum; i++)
				glVertex2f(bcr.cnt[i].x, bcr.cnt[i].y);
		glEnd();
	}

// Draw the curve
	glLineWidth(2.0);
	if (adaptivePlot) {  // plot adaptively
		glColor3f(0.0, 1.0, 0.0);
		adaptiveRender();
	}
	else {  // plot uniformly
		glColor3f(0.0, 0.0, 1.0);
		uniformRender();
	}


	glFlush();		    // process all openGL routines as quickly as possible	         
    glutSwapBuffers();  // swap buffers to display the current frame
}

//============================================================
static void idle( void )
{
   drawCurve();
}


//============================================================
static void hotkey(unsigned char k, int x, int y)
{
// Here we are processing keyboard events.
   switch (k) 
   {
      case 27:
		  free (bcr.cnt);
		  free (bcr.knots);
		  exit (0);
		  break;

// Toggle plotting the control polygon
	  case 'C':
	  case 'c':
		  displayCP = !displayCP;
		  break;

// Toggle sampling points
	  case 'P':
	  case 'p':
		  samplingPnt = !samplingPnt;
		  break;

// Toggle adaptive/uniform plotting
	  case 'A':
	  case 'a':
		  adaptivePlot = !adaptivePlot;
		  break;

// Increase tessellation
	  case '+':
	  case '=':
		  if (adaptivePlot) {
			  tessEps *= 0.7;
			  if (tessEps < 0.5)  tessEps = 0.01;
		  }
		  else {
			  tessNum += 1;
			  if (tessNum > 100) tessNum = 100;
		  }
		  break;

// Decrease tessellation
	  case '-':
	  case '_':
		  if (adaptivePlot) {
			  tessEps *= 1.4;
			  if (tessEps > 50)  tessEps = 50;
		  }
		  else {
			  tessNum -= 1;
			  if (tessNum < 2)
				  tessNum = 2;
		  }
		  break;
   }
}

//============================================================
void chooseWindow()
{
	int    i;
	double left, right, bottom, top;


	left = right = bcr.cnt[0].x;
	for (i=1; i< bcr.cntNum; i++) {
		if (left > bcr.cnt[i].x)  left = bcr.cnt[i].x;
		if (right < bcr.cnt[i].x) right = bcr.cnt[i].x;
	}

	bottom = top = bcr.cnt[0].y;
	for (i=1; i< bcr.cntNum; i++) {
		if (bottom > bcr.cnt[i].y)  bottom = bcr.cnt[i].y;
		if (top < bcr.cnt[i].y) top = bcr.cnt[i].y;
	}

	winLen = top-bottom;
	if (winLen < right-left) winLen = right-left;

	winLen += 100;
	winLLy = bottom - 50;
	winLLx = left - 50 ;
}



//============================================================
int readFile( char* filename )
{
	FILE *fp;
	int  i;

	if((fp = fopen(filename, "r")) == NULL) return 0;  // fail to open the file
 
	fscanf(fp, "%d%d", &(bcr.degree), &(bcr.cntNum));
	bcr.knots = (double *) malloc ((bcr.cntNum+bcr.degree+1)*sizeof(double));
	bcr.cnt = (Point2d *) malloc (bcr.cntNum*sizeof(Point2d));

	for (i=0; i<= bcr.cntNum+bcr.degree; i++)
		fscanf(fp, "%lf", &(bcr.knots[i]));

	for (i=0; i< bcr.cntNum; i++)
		fscanf(fp, "%lf%lf", &(bcr.cnt[i].x),&(bcr.cnt[i].y));
	fclose (fp);

	chooseWindow();

	return 1;
}



//============================================================
void main( int argc, char *argv[] )
{
// load the curve from a file
	char filename[20];

	printf ("\n Please enter a filename: ");
	scanf("%s",filename);
	
    if (readFile(filename)==0) 
		return;

// help information
	printf("\n\nB-spline curve plotting\n");
    printf("NTU, September 2006\n");
    printf("\n");
    printf(" ESC      - Quit program\n");
    printf("\n");
	printf(" A/a : Toggle adaptive/uniform plotting (Default adaptive)\n");
	printf(" C/c : Toggle plotting the control polygon (Default On)\n");
	printf(" P/p : Toggle sampling points (Default Off)\n");
	printf(" +   : Increase tessellation\n");
	printf(" -   : Decrease tessellation\n");
    printf("\n");

	
// set up graphics window
   glutInit(&argc, argv);                         // initialize GLUT
   glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);  // set display mode
   glutInitWindowSize (650, 650);                 // set display window width and height
   glutInitWindowPosition (100, 100);             // set top-left display window position
   glutCreateWindow ("2D B-spline curve plotting:        use +, -, c, a, p, and Esc keys.");

   Init();                        // execute initialization procedure
   glutIdleFunc(idle);            // enables us to make interaction.
   glutDisplayFunc(drawCurve);    // send graphics to display window
   glutKeyboardFunc(hotkey);

   glutMainLoop();                // display everything and wait
}
