/**
 * TODO
 * - Use valgrind to inspect for memory leaks.
 */

#define USE_GLUT

#include <iostream>
#include "viewport.cpp"
#include "xmlload.cpp"

Camera camera;
LightList lights;
MaterialList materials;
ObjectList objects;
Node rootNode;
RenderImage renderImage;

void print(Point3 p, const char *name) {
  std::cout << name << ": " << p.x << ", " << p.y << ", " << p.z << std::endl;
}

float dtor(float d) {
  float pi = 3.141592653589793238462643383279502884197169;
  return d * (pi / 180);
}

HitInfo cast(Ray ro, Node *n)
{
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

float map(float input, float inputStart, float inputEnd, float outputStart, float outputEnd)
{
  float inputRange = inputEnd - inputStart;
  float outputRange = outputEnd - outputStart;

  return (input - inputStart) * outputRange / inputRange + outputStart;
}

Color24 normalMap(Point3 n)
{
  // https://en.wikipedia.org/wiki/Normal_mapping
  return Color24(
    map(n.x, -1, 1, 0, 255),
    map(n.y, -1, 1, 0, 255),
    map(n.z, 0, 1, 128, 255)
  );
}

Color24 shade(Ray &r, HitInfo &h)
{
  if (!h.node) {
    return Color24::Black();
  }
  // return normalMap(h.N);
  return Color24(h.node->GetMaterial()->Shade(r, h, lights, 1));
}

void render(int i, int j, Ray r)
{
  HitInfo h = cast(r, &rootNode);
  
  renderImage.setZBufferPixel(i, j, h.z);
  
  Color24 c = shade(r, h);

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
        (b + (t - b) * ((j + 0.5) / camera.imgHeight)) * tf,
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