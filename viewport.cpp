//-------------------------------------------------------------------------------
///
/// \file       viewport.cpp 
/// \author     Cem Yuksel (www.cemyuksel.com)
/// \version    2.0
/// \date       August 28, 2017
///
/// \brief Example source for CS 6620 - University of Utah.
///
//-------------------------------------------------------------------------------

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#include "scene.h"
#include "objects.h"
#include "lights.h"
#include "materials.h"
#include <stdlib.h>
#include <time.h>

#ifdef USE_GLUT
# ifdef __APPLE__
#  include <GLUT/glut.h>
# else
#  include <GL/glut.h>
# endif
#else
# include <GL/freeglut.h>
#endif

//-------------------------------------------------------------------------------

void BeginRender();	// Called to start rendering (renderer must run in a separate thread)
void StopRender();	// Called to end rendering (if it is not already finished)

extern Node rootNode;
extern Camera camera;
extern RenderImage renderImage;
extern LightList lights;

//-------------------------------------------------------------------------------

enum Mode {
	MODE_READY,			// Ready to render
	MODE_RENDERING,		// Rendering the image
	MODE_RENDER_DONE	// Rendering is finished
};

enum ViewMode
{
	VIEWMODE_OPENGL,
	VIEWMODE_IMAGE,
	VIEWMODE_Z,
};

enum MouseMode {
	MOUSEMODE_NONE,
	MOUSEMODE_DEBUG,
	MOUSEMODE_ROTATE,
};

static Mode		mode		= MODE_READY;		// Rendering mode
static ViewMode	viewMode	= VIEWMODE_OPENGL;	// Display mode
static MouseMode mouseMode	= MOUSEMODE_NONE;	// Mouse mode
static int		startTime;						// Start time of rendering
static int mouseX=0, mouseY=0;
static float viewAngle1=0, viewAngle2=0;
static GLuint viewTexture;

//-------------------------------------------------------------------------------

void GlutDisplay();
void GlutReshape(int w, int h);
void GlutIdle();
void GlutKeyboard(unsigned char key, int x, int y);
void GlutMouse(int button, int state, int x, int y);
void GlutMotion(int x, int y);

//-------------------------------------------------------------------------------

void ShowViewport()
{
	int argc = 1;
	char argstr[] = "raytrace";
	char *argv = argstr;
	glutInit(&argc,&argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );
	if (glutGet(GLUT_SCREEN_WIDTH) > 0 && glutGet(GLUT_SCREEN_HEIGHT) > 0){
		glutInitWindowPosition( (glutGet(GLUT_SCREEN_WIDTH) - camera.imgWidth)/2, (glutGet(GLUT_SCREEN_HEIGHT) - camera.imgHeight)/2 );
	}
	else glutInitWindowPosition( 50, 50 );
	glutInitWindowSize(camera.imgWidth, camera.imgHeight);

	glutCreateWindow("Ray Tracer - CS 6620");
	glutDisplayFunc(GlutDisplay);
	glutReshapeFunc(GlutReshape);
	glutIdleFunc(GlutIdle);
	glutKeyboardFunc(GlutKeyboard);
	glutMouseFunc(GlutMouse);
	glutMotionFunc(GlutMotion);

	glClearColor(0,0,0,0);

	glPointSize(3.0);
	glEnable( GL_CULL_FACE );

	float zero[] = {0,0,0,0};
	glLightModelfv( GL_LIGHT_MODEL_AMBIENT, zero );

	glEnable(GL_NORMALIZE);

	glLineWidth(2);

	glGenTextures(1,&viewTexture);
	glBindTexture(GL_TEXTURE_2D, viewTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	glutMainLoop();
}

//-------------------------------------------------------------------------------

void GlutReshape(int w, int h)
{
	if( w != camera.imgWidth || h != camera.imgHeight ) {
		glutReshapeWindow( camera.imgWidth, camera.imgHeight);
	} else {
		glViewport( 0, 0, w, h );

		glMatrixMode( GL_PROJECTION );
		glLoadIdentity();
		float r = (float) w / float (h);
		gluPerspective( camera.fov, r, 0.02, 1000.0);

		glMatrixMode( GL_MODELVIEW );
		glLoadIdentity();
	}
}

//-------------------------------------------------------------------------------

void DrawNode( Node *node )
{
	glPushMatrix();

	const Material *mtl = node->GetMaterial();
	if ( mtl ) mtl->SetViewportMaterial();

	Matrix3 tm = node->GetTransform();
	Point3 p = node->GetPosition();
	float m[16] = { tm[0],tm[1],tm[2],0, tm[3],tm[4],tm[5],0, tm[6],tm[7],tm[8],0, p.x,p.y,p.z,1 };
	glMultMatrixf( m );

	Object *obj = node->GetNodeObj();
	if ( obj ) obj->ViewportDisplay(mtl);

	for ( int i=0; i<node->GetNumChild(); i++ ) {
		DrawNode( node->GetChild(i) );
	}

	glPopMatrix();
}

//-------------------------------------------------------------------------------

void DrawScene()
{
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );

	glEnable( GL_LIGHTING );
	glEnable( GL_DEPTH_TEST );

	glPushMatrix();
	Point3 p = camera.pos;
	Point3 t = camera.pos + camera.dir;
	Point3 u = camera.up;
	gluLookAt( p.x, p.y, p.z,  t.x, t.y, t.z,  u.x, u.y, u.z );

	glRotatef( viewAngle1, 1, 0, 0 );
	glRotatef( viewAngle2, 0, 0, 1 );

	if ( lights.size() > 0 ) {
		for ( unsigned int i=0; i<lights.size(); i++ ) {
			lights[i]->SetViewportLight(i);
		}
	} else {
		float white[] = {1,1,1,1};
		float black[] = {0,0,0,0};
		Point4 p(camera.pos, 1);
		glEnable ( GL_LIGHT0 );
		glLightfv( GL_LIGHT0, GL_AMBIENT,  black );
		glLightfv( GL_LIGHT0, GL_DIFFUSE,  white );
		glLightfv( GL_LIGHT0, GL_SPECULAR, white );
		glLightfv( GL_LIGHT0, GL_POSITION, &p.x );
	}

	DrawNode(&rootNode);

	glPopMatrix();

	glDisable( GL_DEPTH_TEST );
	glDisable( GL_LIGHTING );
}

//-------------------------------------------------------------------------------

void DrawImage( void *data, GLenum type, GLenum format )
{
	glBindTexture(GL_TEXTURE_2D, viewTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, renderImage.GetWidth(), renderImage.GetHeight(), 0, format, type, data); 

	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );

	glEnable(GL_TEXTURE_2D);

	glMatrixMode( GL_PROJECTION );
	glPushMatrix();
	glLoadIdentity();

	glColor3f(1,1,1);
	glBegin(GL_QUADS);
	glTexCoord2f(0,1);
	glVertex2f(-1,-1);
	glTexCoord2f(1,1);
	glVertex2f(1,-1);
	glTexCoord2f(1,0);
	glVertex2f(1,1);
	glTexCoord2f(0,0);
	glVertex2f(-1,1);
	glEnd();

	glPopMatrix();
	glMatrixMode( GL_MODELVIEW );

	glDisable(GL_TEXTURE_2D);
}

//-------------------------------------------------------------------------------

void DrawProgressBar(float done)
{
	glMatrixMode( GL_PROJECTION );
	glPushMatrix();
	glLoadIdentity();

	glBegin(GL_LINES);
	glColor3f(1,1,1);
	glVertex2f(-1,-1);
	glVertex2f(done*2-1,-1);
	glColor3f(0,0,0);
	glVertex2f(done*2-1,-1);
	glVertex2f(1,-1);
	glEnd();

	glPopMatrix();
	glMatrixMode( GL_MODELVIEW );
}

//-------------------------------------------------------------------------------

void DrawRenderProgressBar()
{
	int rp = renderImage.GetNumRenderedPixels();
	int np = renderImage.GetWidth() * renderImage.GetHeight();
	if ( rp >= np ) return;
	float done = (float) rp / (float) np;
	DrawProgressBar(done);
}

//-------------------------------------------------------------------------------

void GlutDisplay()
{
	switch ( viewMode ) {
	case VIEWMODE_OPENGL:
		DrawScene();
		break;
	case VIEWMODE_IMAGE:
		//DrawImage( renderImage.GetPixels(), GL_UNSIGNED_BYTE, GL_RGB );
		glDrawPixels( renderImage.GetWidth(), renderImage.GetHeight(), GL_RGB, GL_UNSIGNED_BYTE, renderImage.GetPixels() );
		DrawRenderProgressBar();
		break;
	case VIEWMODE_Z:
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );
		if ( ! renderImage.GetZBufferImage() ) renderImage.ComputeZBufferImage();
		//DrawImage( renderImage.GetZBufferImage(), GL_UNSIGNED_BYTE, GL_LUMINANCE );
		glDrawPixels( renderImage.GetWidth(), renderImage.GetHeight(), GL_LUMINANCE, GL_UNSIGNED_BYTE, renderImage.GetZBufferImage() );
		break;
	}

	glutSwapBuffers();
}

//-------------------------------------------------------------------------------

void GlutIdle()
{
	static int lastRenderedPixels = 0;
	if ( mode == MODE_RENDERING ) {
		int nrp = renderImage.GetNumRenderedPixels();
		if ( lastRenderedPixels != nrp ) {
			lastRenderedPixels = nrp;
			if ( renderImage.IsRenderDone() ) {
				mode = MODE_RENDER_DONE;
				int endTime = (int) time(NULL);
				int t = endTime - startTime;
				int h = t / 3600;
				int m = (t % 3600) / 60;
				int s = t % 60;
				printf("\nRender time is %d:%02d:%02d.\n",h,m,s);
			}
			glutPostRedisplay();
		}
	}
}

//-------------------------------------------------------------------------------

void GlutKeyboard(unsigned char key, int x, int y)
{
	switch ( key ) {
	case 27:	// ESC
		exit(0);
		break;
	case ' ':
		switch ( mode ) {
		case MODE_READY: 
			mode = MODE_RENDERING;
			viewMode = VIEWMODE_IMAGE;
			DrawScene();
			glReadPixels( 0, 0, renderImage.GetWidth(), renderImage.GetHeight(), GL_RGB, GL_UNSIGNED_BYTE, renderImage.GetPixels() );
			{
				Color24 *c = renderImage.GetPixels();
				for ( int y0=0, y1=renderImage.GetHeight()-1; y0<y1; y0++, y1-- ) {
					int i0 = y0 * renderImage.GetWidth();
					int i1 = y1 * renderImage.GetWidth();
					for ( int x=0; x<renderImage.GetWidth(); x++, i0++, i1++ ) {
						Color24 t=c[i0]; c[i0]=c[i1]; c[i1]=t;
					}
				}
			}
			startTime = (int) time(NULL);
			BeginRender();
			break;
		case MODE_RENDERING:
			mode = MODE_READY;
			StopRender();
			glutPostRedisplay();
			break;
		case MODE_RENDER_DONE: 
			mode = MODE_READY;
			viewMode = VIEWMODE_OPENGL;
			glutPostRedisplay();
			break;
		}
		break;
	case '1':
		viewAngle1 = viewAngle2 = 0;
		viewMode = VIEWMODE_OPENGL;
		glutPostRedisplay();
		break;
	case '2':
		viewMode = VIEWMODE_IMAGE;
		glutPostRedisplay();
		break;
	case '3':
		viewMode = VIEWMODE_Z;
		glutPostRedisplay();
		break;
	}
}

//-------------------------------------------------------------------------------

void PrintPixelData(int x, int y)
{
	if ( x < renderImage.GetWidth() && y < renderImage.GetHeight() ) {
		Color24 *colors = renderImage.GetPixels();
		float *zbuffer = renderImage.GetZBuffer();
		int i = (renderImage.GetHeight() - y - 1 ) *renderImage.GetWidth() + x;
		printf("Pixel [ %d, %d ] Color24: %d, %d, %d   Z: %f\n", x, y, colors[i].r, colors[i].g, colors[i].b, zbuffer[i] );
	} else {
		printf("-- Invalid pixel (%d,%d) --\n",x,y);
	}
}

//-------------------------------------------------------------------------------

void GlutMouse(int button, int state, int x, int y)
{
	if ( state == GLUT_UP ) {
		mouseMode = MOUSEMODE_NONE;
	} else {
		switch ( button ) {
			case GLUT_LEFT_BUTTON:
				mouseMode = MOUSEMODE_DEBUG;
				PrintPixelData(x,y);
				break;
			case GLUT_RIGHT_BUTTON:
				mouseMode = MOUSEMODE_ROTATE;
				mouseX = x;
				mouseY = y;
				break;
		}
	}
}

//-------------------------------------------------------------------------------

void GlutMotion(int x, int y)
{
	switch ( mouseMode ) {
		case MOUSEMODE_NONE:
				// Do nothing.
				break;
		case MOUSEMODE_DEBUG:
			PrintPixelData(x,y);
			break;
		case GLUT_RIGHT_BUTTON:
			viewAngle1 -= 0.2f * ( mouseY - y );
			viewAngle2 -= 0.2f * ( mouseX - x );
			mouseX = x;
			mouseY = y;
			glutPostRedisplay();
			break;
	}
}

//-------------------------------------------------------------------------------
// Viewport Methods for various classes
//-------------------------------------------------------------------------------
bool Sphere::IntersectRay(const Ray &ray, HitInfo &hInfo, int hitSide) const
{
	Point3 o = Point3(0, 0, 0);
	float r = 1;
	float b = ray.dir % ray.p;
	float a = ray.dir % ray.dir;
	float c = ray.p % ray.p - r * r;
	float d = b * b - a * c;
	if (d < 0) {
		return false;
	}

	d = sqrtf(d);
	float t0 = (-ray.dir % ray.p - d) / a;
	float t1 = (-ray.dir % ray.p + d) / a;
	if (t0 < 0 && t1 < 0) {
		return false;
	}
	
	if (t0 == t1 || (t0 > 0 && t1 > 0)) {
		hInfo.z = min(t0, t1);	
	} else {
		hInfo.z = max(t0, t1);
	}
	hInfo.p = ray.p + ray.dir * hInfo.z;
	hInfo.N = (hInfo.p - o).GetNormalized();
	return true;
}

void Sphere::ViewportDisplay(const Material *mtl) const
{
	static GLUquadric *q = NULL;
	if ( q == NULL ) {
		q = gluNewQuadric();
	}
	gluSphere(q,1,50,50);
}

HitInfo cast(Ray ro, Node *n = &rootNode)
{
	assert(!isnan(ro.dir.x));
  // Compute the ray in parent space.
  Ray r = Ray(ro);
  r = n->ToNodeCoords(r);
  // Renormalize since scaling may denormalize the direction.
  r.Normalize();

  HitInfo h = HitInfo();
  Object *o = n->GetNodeObj();
  if (o) {
    if (o->IntersectRay(r, h)) {
      // h.z is set in IntersectRay.
      h.node = n;
    }
  }

  for (int i = 0; i < n->GetNumChild(); i++) {
    HitInfo hc = cast(r, n->GetChild(i));
    if (hc.z < h.z) {
      h = hc;
    }
  }

  if (h.node) {
    // Compute the intersection and normal in parent space.
    n->FromNodeCoords(h);
    // Compute the distance in parent space.
    h.z = (h.p - ro.p).Length();
  }
  return h;
}

Point3 reflect(Point3 d, Point3 n)
{
	// Shirley 10.6
	// Construction of Shirley 9.7
	// Assume normalized vectors l and n.
	// vproj(l, n) = sproj(l, n)n = (|l|cos(t))n = cos(t)n
	return d - 2 * n * (d % n);
}

Point3 refract(Point3 d, Point3 n, float e)
{
	// Shirley 10.7
	float k = 1 - e * e * (1 - (d % n) * (d % n));
	if (k < 0) {  
		return Point3(0, 0, 0);
	}
	return e * (d - n * (d % n)) - n * sqrt(k);
}

float GenLight::Shadow(Ray ray, float t_max)
{
	// TODO: t_max for point lights was modified to be |L| instead of 1.
	// Figure out what exactly 1 corresponds to.
	float e = 0.01;
	Ray r = Ray(ray.p, ray.dir.GetNormalized());
	r.p = r.p + e * r.dir;
	HitInfo h = cast(r);
	if (h.node && h.z <= t_max) {
		return 0;
	}
	return 1;
}

Color MtlBlinn::Shade(const Ray &ray, const HitInfo &hInfo, const LightList &lights, int bounceCount) const
{
	int bounceMax = 3;
	Color color = Color::Black();
	if (!hInfo.node || bounceCount > bounceMax) {
    return Color::Black();
  }

	float e = 0.001;

	for (int i = 0; i < lights.size(); i++) {
    Light *light = lights[i];
		Color li = light->Illuminate(hInfo.p, hInfo.N);
    if (!light->IsAmbient()) {
      float s = 0;
      Point3 l = light->Direction(hInfo.p).GetNormalized() * -1;
			float lambertian = max(l % hInfo.N, 0);
      if (lambertian > 0) {
        // Blinn-Phong
        // https://en.wikipedia.org/wiki/Blinn–Phong_shading_model
        Point3 v = (ray.p - hInfo.p).GetNormalized();
        Point3 h = (l + v).GetNormalized();
        float sa = max(h % hInfo.N, 0);
        s = powf(sa, glossiness);
				color += specular * lambertian * s * li;
      }
			color += diffuse * lambertian * li;
    } else {
			color +=  diffuse * li;
		}
	}

	if (reflection != Color::Black()) {
		// Shirley 10.6
		Point3 r = reflect(ray.dir, hInfo.N);
		Ray rn = Ray(hInfo.p + e * r, r);
		HitInfo hn = cast(rn);
		if (hn.node) {
			// c *= (Color::White() - reflection);?
			color += reflection * specular * hn.node->GetMaterial()->Shade(rn, hn, lights, bounceCount + 1);
		}
	}

	if (refraction != Color::Black()) {
		// Shirley 10.7
		// Split the recursion path for dielectrics.
		Point3 t;
		Color k;
		float c;
		float n = 1 / ior;
		float R = 0;
		if (ray.dir % hInfo.N < 0) {
			// "On the outisde looking in".
			t = refract(ray.dir, hInfo.N, n);
			k = Color::White();
			c = -ray.dir % hInfo.N;
		} else {
			t = refract(ray.dir, -hInfo.N, 1 / n);
			k.r = expf(-absorption.r * e);
			k.g = expf(-absorption.g * e);
			k.b = expf(-absorption.b * e);
			if (!t.IsZero()) {
				c = t % hInfo.N;
			} else {
				// TODO: Consider returning early after refactoring cast -> check node -> shade.
				R = 1;
			}
		}

		Color cRefract = Color::Black();
		if (!t.IsZero()) {
			Ray r = Ray(hInfo.p + e * t, t);
			HitInfo h = cast(r);
			if (h.node) {
				cRefract = h.node->GetMaterial()->Shade(r, h, lights, bounceCount + 1);
			}
		}

		Color cReflect = Color::Black();
		Point3 rd = reflect(ray.dir, hInfo.N);
		Ray r = Ray(hInfo.p + e * rd, rd);
		HitInfo h = cast(r);
		if (h.node) {
			cReflect = h.node->GetMaterial()->Shade(r, h, lights, bounceCount + 1);
		}
		
		if (R == 0) {
			float Ro = ((n - 1) * (n - 1)) / ((n + 1) * (n + 1));
			R = Ro + (1 - Ro) * powf(1 - c, 5);
		}
		color += k * (R * cReflect + (1 - R) * cRefract * refraction);
	}

	return color;
}

void MtlBlinn::SetViewportMaterial(int subMtlID) const
{
	ColorA c;
	c = ColorA(diffuse);
	glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, &c.r );
	c = ColorA(specular);
	glMaterialfv( GL_FRONT, GL_SPECULAR, &c.r );
	glMaterialf( GL_FRONT, GL_SHININESS, glossiness*1.5f );
}

void GenLight::SetViewportParam(int lightID, ColorA ambient, ColorA intensity, Point4 pos ) const
{
	glEnable ( GL_LIGHT0 + lightID );
	glLightfv( GL_LIGHT0 + lightID, GL_AMBIENT,  &ambient.r );
	glLightfv( GL_LIGHT0 + lightID, GL_DIFFUSE,  &intensity.r );
	glLightfv( GL_LIGHT0 + lightID, GL_SPECULAR, &intensity.r );
	glLightfv( GL_LIGHT0 + lightID, GL_POSITION, &pos.x );
}
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
