/**
 * TODO
 * - Use valgrind to inspect for memory leaks.
 */

#define USE_GLUT

#include <iostream>
#include <thread>
#include "viewport.cpp"
#include "xmlload.cpp"

Camera camera;
LightList lights;
MaterialList materials;
ObjectList objects;
ObjFileList objList;
Node rootNode;
RenderImage renderImage;
TexturedColor background;
TexturedColor environment;
TextureList textureList;

float dtor(float d) {
  float pi = 3.141592653589793238462643383279502884197169;
  return d * (pi / 180);
}

float map(float input, float inputStart, float inputEnd, float outputStart, float outputEnd)
{
  float inputRange = inputEnd - inputStart;
  float outputRange = outputEnd - outputStart;

  return (input - inputStart) * outputRange / inputRange + outputStart;
}

Color24 normalMap(HitInfo &h)
{
  if (!h.node) {
    return Color24::Black();
  }
  // https://en.wikipedia.org/wiki/Normal_mapping
  return Color24(
    map(h.N.x, -1, 1, 0, 255),
    map(h.N.y, -1, 1, 0, 255),
    map(h.N.z, 0, 1, 128, 255)
  );
}

Color shade(Ray &r, HitInfo &h)
{
  Color c = Color::Black();
  if (!h.node) {
    return c;
  }
  return h.node->GetMaterial()->Shade(r, h, lights, 0);
}

Color color(int i, int j, Ray r)
{
  HitInfo h = cast(r);
  renderImage.setZBufferPixel(i, j, h.z);
  Color c = shade(r, h);
  // Color c = normalMap(h);
  return c;
}

void render(int w1, int w2, int h1, int h2)
{
  for (int i = w1; i < w2; i++) {
    for (int j = h1; j < h2; j++) {
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
      float arw = camera.imgWidth > camera.imgHeight
        ? float(camera.imgWidth) / camera.imgHeight
        : 1;
      float arh = camera.imgHeight > camera.imgWidth
        ? float(camera.imgHeight) / camera.imgWidth
        : 1;
      
      // Compute field of view.
      // Shirley 7.5
      float tf = tan(dtor(camera.fov / 2)) * abs(n);

      // Supersample anti-aliasing
      Color c = Color::Black();
      int AA = 1;
      for (int a1 = 0; a1 < AA; a1++) {
        for (int a2 = 0; a2 < AA; a2++) {
          Point2 o = Point2(a1, a2) / float(AA) - 0.5;

          // Compute screen coordinates in camera space.
          // Shirley 10.2
          Point3 sc = Point3(
            (l + (r - l) * ((i + 0.5 + o.x) / camera.imgWidth)) * arw * tf,
            (b + (t - b) * ((j + 0.5 + o.y) / camera.imgHeight)) * arh * tf,
            -1
          );

          // Compute screen position in world space.
          // Shirley 10.2
          Point3 sw = camera.pos + sc.x * u + sc.y * v + sc.z * w;
          Point3 rd = (sw - camera.pos).GetNormalized();

          c += color(i, j, Ray(camera.pos, rd));
        }
      }
      renderImage.setRenderedPixel(i, j, Color24(c / float(AA * AA)));
    }
  }
}

void BeginRender()
{
  std::cout << "begin render" << std::endl;
  
  int n = 2;
  std::thread threads[n * n];
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      threads[i * n + j] = std::thread(
        render,
        (camera.imgWidth / n) * i,
        (camera.imgWidth / n) * (i + 1),
        (camera.imgHeight / n) * j,
        (camera.imgHeight / n) * (j + 1)
      );
    }
  }
  for (int i = 0; i < n * n; i++) {
    threads[i].join();
  }
  renderImage.SaveImage("image.png");
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