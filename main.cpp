#define USE_GLUT

#include <iostream>
#include "viewport.cpp"
#include "xmlload.cpp"

Camera camera;
Node rootNode;
RenderImage renderImage;

void print(Point3 p, const char *name) {
  std::cout << name << ": " << p.x << ", " << p.y << ", " << p.z << std::endl;
}

float dtor(float d) {
  float pi = 3.141592653589793238462643383279502884197169;
  return d * (pi / 180);
}

HitInfo cast(Ray r, Node *n)
{
  Point3 po = Point3(r.p);
  HitInfo h = HitInfo();
  Object *o = n->GetNodeObj();
  if (o) {
    r = n->ToNodeCoords(r);
    // Renormalize since scaling may denormalize the direction.
    r.Normalize();
    if (o->IntersectRay(r, h)) {
      // h.z is set in IntersectRay.
      h.node = n;
    }
  }

  for (int i = 0; i < n->GetNumChild(); i++) {
    HitInfo hc = cast(r, n->GetChild(i));
    if (hc.node && hc.z < h.z) {
        h = hc;
    }
	}

  if (o && h.node) {
    // The intersection in object space.
    Point3 i = r.p + r.dir * h.z;
    // The intersection in parent space.
    i = h.node->TransformFrom(i);
    // Compute z in parent space.
    h.z = (i - po).Length();
  }
  return h;
}

void render(int i, int j, Ray r)
{
  HitInfo h = cast(r, &rootNode);
  
  renderImage.setZBufferPixel(i, j, h.z);

  Color24 c;
  int v = h.node ? 255 : 0;
  c.r = v;
  c.g = v;
  c.b = v;
  renderImage.setRenderedPixel(i, j, c);
}

void BeginRender()
{
  std::cout << "begin render" << std::endl;
  
  for (int i = 0; i < camera.imgWidth; i++) {
    for (int j = 0; j < camera.imgHeight; j++) {
      // Useful reference.
      // https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-generating-camera-rays/generating-camera-rays

      // Construct the orthographic view volume.
      // Shirley 7.2
      float b = -1, f = -1, l = -1;
      float t =  1, n =  1, r =  1;

      // Construct a coordinate system (orthonormal frame).
      // Shirley 7.2.1
      Point3 w = Point3(camera.dir).GetNormalized() * -1;
      Point3 u = camera.up.Cross(w).GetNormalized();
      Point3 v = w.Cross(u);

      // Compute for aspect ratio.
      // Shirley 7.5
      float ar = camera.imgWidth > camera.imgHeight
        ? float(camera.imgWidth) / camera.imgHeight
        : float(camera.imgHeight) / camera.imgWidth;
      
      // Compute field of view.
      // Shirley 7.5
      float tf = tan(dtor(camera.fov / 2)) * abs(n);

      // Compute screen coordinates in camera space.
      // Shirley 10.2
      Point3 sc = Point3(
        (l + (r - l) * ((i + 0.5) / camera.imgWidth)) * ar * tf,
        // Corrected to invert y from [1, -1] to [-1, 1].
        (-b - (t - b) * ((j + 0.5) / camera.imgHeight)) * tf,
        -1
      );

      // Compute screen position in world space.
      // Shirley 10.2
      Point3 sw = camera.pos + sc.x * u + sc.y * v + sc.z * w;

      Point3 rd = (sw - camera.pos).GetNormalized();
      render(i, j, Ray(camera.pos, rd));
    }
  }
}

void StopRender()
{
  // TODO
  std::cout << "stop render" << std::endl;
}

int main()
{
  std::cout << "tracer" << std::endl;
  LoadScene("scene.xml");
  ShowViewport();
  return 0;
}