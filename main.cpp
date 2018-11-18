#define USE_GLUT

#include <iostream>
#include "viewport.cpp"
#include "xmlload.cpp"

Camera camera;
Node rootNode;
RenderImage renderImage;

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

Color24 render(Ray r)
{
  HitInfo h = cast(r, &rootNode);
  Color24 c;
  int v = h.node ? 255 : 0;
  c.r = v;
  c.g = v;
  c.b = v;
  return c;
}

void BeginRender()
{
  std::cout << "begin render" << std::endl;
  
  for (int x = 0; x < camera.imgWidth; x++) {
    for (int y = 0; y < camera.imgHeight; y++) {
      // TODO: Account for FOV and aspect ratio.
      // https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-generating-camera-rays/generating-camera-rays
      float xw = 2 * float(x) / camera.imgWidth - 1;
      float yw = 2 * float(y) / camera.imgHeight - 1;
      Point3 pd = Point3(xw, yw, camera.dir.z);
      pd.Normalize();
      Ray r = Ray(camera.pos, pd);
      Color24 c = render(r);
      renderImage.setRenderedPixel(x, y, c);
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