#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>

using namespace std;
using glm::ivec2;
using glm::mat3;
using glm::mat4;
using glm::vec2;
using glm::vec3;
using glm::vec4;

SDL_Event event;

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
// #define FOCAL_LENGTH (SCREEN_HEIGHT)
#define FULLSCREEN_MODE false

/* ----------------------------------------------------------------------------*/
/* VARIABLES                                                                   */
vec3 black(0, 0, 0);
vec3 white(1, 1, 1);
vec4 cameraPos(0, 0, -3.001, 1);
vector<Triangle> triangles;
mat4 R;
float theta_x = 0, theta_y = 0, theta_z = 0;
float depthBuffer[SCREEN_HEIGHT * SCREEN_WIDTH];
int FOCAL_LENGTH = SCREEN_HEIGHT;

/* ----------------------------------------------------------------------------*/
/* STRUCTS                                                                     */
struct Pixel
{
  int x;
  int y;
  float zinv;
};

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void LoadRotationMatrix();

void Draw(screen *screen);

void VertexShader(const vec4 &v, glm::ivec2 &p);
void VertexShader(const vec4 &v, Pixel &p);

void Interpolate(ivec2 a, ivec2 b, vector<ivec2> &result);
void Interpolate(Pixel a, Pixel b, vector<Pixel> &result);

void DrawPolygonEdges(screen *screen, const vector<vec4> &vertices);

void DrawLineSDL(screen *screen, ivec2 a, ivec2 b, vec3 color);
void DrawLineSDL(screen *screen, Pixel a, Pixel b, vec3 color);

void ComputePolygonRows(const vector<ivec2> &vertexPixels, vector<ivec2> &leftPixels, vector<ivec2> &rightPixels);
void ComputePolygonRows(const vector<Pixel> &vertexPixels, vector<Pixel> &leftPixels, vector<Pixel> &rightPixels);

void DrawPolygonRows(screen *screen, const vector<ivec2> &leftPixels, const vector<ivec2> &rightPixels, vec3 color);
void DrawPolygonRows(screen *screen, const vector<Pixel> &leftPixels, const vector<Pixel> &rightPixels, vec3 color);

void DrawPolygon(screen *screen, const vector<vec4> &vertices, vec3 color);
void DrawPolygonPixel(screen *screen, const vector<vec4> &vertices, vec3 color);

void TestPolyRows();

int main(int argc, char *argv[])
{

  screen *screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE);
  LoadTestModel(triangles);
  LoadRotationMatrix();

  while (Update())
  {
    Draw(screen);
    SDL_Renderframe(screen);
  }

  SDL_SaveImage(screen, "screenshot.bmp");

  KillSDL(screen);
  return 0;
}

/*Place your drawing here*/
void Draw(screen *screen)
{
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height * screen->width * sizeof(uint32_t));
  memset(depthBuffer, 0, screen->height * screen->width * sizeof(float));
  for (uint32_t i = 0; i < triangles.size(); ++i)
  {
    vector<vec4> vertices(3);
    vertices[0] = triangles[i].v0;
    vertices[1] = triangles[i].v1;
    vertices[2] = triangles[i].v2;
    DrawPolygonPixel(screen, vertices, triangles[i].color);
  }
}
/*Place updates of parameters here*/
bool Update()
{
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2 - t);
  t = t2;

  SDL_Event e;
  while (SDL_PollEvent(&e))
  {
    if (e.type == SDL_QUIT)
    {
      return false;
    }
    else if (e.type == SDL_KEYDOWN)
    {
      int key_code = e.key.keysym.sym;
      switch (key_code)
      {
      case SDLK_UP:
        FOCAL_LENGTH += 5;
        break;
      case SDLK_DOWN:
        FOCAL_LENGTH -= 5;
        break;
      case SDLK_LEFT:
        /* Move camera left */
        break;
      case SDLK_RIGHT:
        /* Move camera right */
        break;
      case SDLK_ESCAPE:
        /* Move camera quit */
        return false;
      }
    }
  }
  return true;
}

void DrawPolygon(screen *screen, const vector<vec4> &vertices, vec3 color)
{
  int V = vertices.size();
  vector<ivec2> vertexPixels(V);
  for (int i = 0; i < V; ++i)
    VertexShader(vertices[i], vertexPixels[i]);
  vector<ivec2> leftPixels;
  vector<ivec2> rightPixels;
  ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
  DrawPolygonRows(screen, leftPixels, rightPixels, color);
}

void DrawPolygonPixel(screen *screen, const vector<vec4> &vertices, vec3 color)
{
  int V = vertices.size();
  vector<Pixel> vertexPixels(V);
  for (int i = 0; i < V; ++i)
    VertexShader(vertices[i], vertexPixels[i]);
  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;
  ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
  DrawPolygonRows(screen, leftPixels, rightPixels, color);
}

void DrawPolygonRows(screen *screen, const vector<ivec2> &leftPixels,
                     const vector<ivec2> &rightPixels, vec3 color)
{
  for (int i = 0; i < leftPixels.size(); ++i)
    DrawLineSDL(screen, leftPixels[i], rightPixels[i], color);
}
void DrawPolygonRows(screen *screen, const vector<Pixel> &leftPixels, const vector<Pixel> &rightPixels, vec3 color)
{
  for (int i = 0; i < leftPixels.size(); ++i)
    DrawLineSDL(screen, leftPixels[i], rightPixels[i], color);
}

void ComputePolygonRows(const vector<ivec2> &vertexPixels, vector<ivec2> &leftPixels, vector<ivec2> &rightPixels)
{
  int vp_size = vertexPixels.size();
  //STEP 1: MAX AND MIN Y VALUE
  int max_y = INT_MIN, min_y = INT_MAX;
  for (int i = 0; i < vp_size; ++i)
  {
    if (max_y < vertexPixels[i].y)
      max_y = vertexPixels[i].y;
    if (min_y > vertexPixels[i].y)
      min_y = vertexPixels[i].y;
  }
  //STEP 2: RESIZE LEFT-RIGHT PIXELS
  int num_rows = max_y - min_y + 1;
  leftPixels.resize(num_rows);
  rightPixels.resize(num_rows);
  //STEP 3: x in left to really large, right to really small
  for (int i = 0; i < num_rows; ++i)
  {
    leftPixels[i].x = INT_MAX;
    rightPixels[i].x = INT_MIN;
  }
  //STEP 4: ITERATE AND UPDATE VALUES
  for (int i = 0; i < vp_size; ++i)
  {
    int j = (i + 1) % vp_size;
    ivec2 a = vertexPixels[i];
    ivec2 b = vertexPixels[j];
    ivec2 delta = glm::abs(a - b);
    int pixels = glm::max(delta.x, delta.y) + 1;
    vector<ivec2> line(pixels);
    Interpolate(a, b, line);
    for (int k = 0; k < pixels; ++k)
    {
      int c = line[k].y - min_y;
      if (leftPixels[c].x > line[k].x)
        leftPixels[c] = line[k];
      if (rightPixels[c].x < line[k].x)
        rightPixels[c] = line[k];
    }
  }
}

void ComputePolygonRows(const vector<Pixel> &vertexPixels, vector<Pixel> &leftPixels, vector<Pixel> &rightPixels)
{
  int vp_size = vertexPixels.size();
  //STEP 1: MAX AND MIN Y VALUE
  int max_y = INT_MIN, min_y = INT_MAX;
  for (int i = 0; i < vp_size; ++i)
  {
    if (max_y < vertexPixels[i].y)
      max_y = vertexPixels[i].y;
    if (min_y > vertexPixels[i].y)
      min_y = vertexPixels[i].y;
  }
  //STEP 2: RESIZE LEFT-RIGHT PIXELS
  int num_rows = max_y - min_y + 1;
  leftPixels.resize(num_rows);
  rightPixels.resize(num_rows);
  //STEP 3: x in left to really large, right to really small
  for (int i = 0; i < num_rows; ++i)
  {
    leftPixels[i].x = INT_MAX;
    rightPixels[i].x = INT_MIN;
  }
  //STEP 4: ITERATE AND UPDATE VALUES
  for (int i = 0; i < vp_size; ++i)
  {
    int j = (i + 1) % vp_size;
    Pixel a = vertexPixels[i];
    Pixel b = vertexPixels[j];
    ivec2 delta = glm::abs(ivec2(a.x, a.y) - ivec2(b.x, b.y));
    int pixels = glm::max(delta.x, delta.y) + 1;
    vector<Pixel> line(pixels);
    Interpolate(a, b, line);
    for (int k = 0; k < pixels; ++k)
    {
      int c = glm::abs(line[k].y - min_y);
      if (leftPixels[c].x > line[k].x)
        leftPixels[c] = line[k];
      if (rightPixels[c].x < line[k].x)
        rightPixels[c] = line[k];
    }
  }
}

void DrawPolygonEdges(screen *screen, const vector<vec4> &vertices)
{
  int V = vertices.size();
  vector<ivec2> projectedVertices(V);
  for (int i = 0; i < V; ++i)
    VertexShader(vertices[i], projectedVertices[i]);
  for (int i = 0; i < V; ++i)
  {
    int j = (i + 1) % V;
    DrawLineSDL(screen, projectedVertices[i], projectedVertices[j], white);
  }
}

void DrawLineSDL(screen *screen, ivec2 a, ivec2 b, vec3 color)
{
  ivec2 delta = glm::abs(a - b);
  int pixels = glm::max(delta.x, delta.y) + 1;
  vector<ivec2> line(pixels);
  Interpolate(a, b, line);
  for (int i = 0; i < pixels; ++i)
    PutPixelSDL(screen, line[i].x, line[i].y, color);
}

void DrawLineSDL(screen *screen, Pixel a, Pixel b, vec3 color)
{
  ivec2 delta = glm::abs(ivec2(a.x, a.y) - ivec2(b.x, b.y));
  int pixels = glm::max(delta.x, delta.y) + 1;
  vector<Pixel> line(pixels);
  Interpolate(a, b, line);
  for (int i = 0; i < pixels; ++i)
  {
    std::cout << "START" << std::endl;
    int index = line[i].y * SCREEN_WIDTH + line[i].x;
    std::cout << index << std::endl;
    if (index > SCREEN_HEIGHT * SCREEN_WIDTH - 1 || index < 0) continue;
    if (depthBuffer[index] < line[i].zinv)
    {
      depthBuffer[index] = line[i].zinv;
      std::cout << "BETWEEN" << std::endl;
      PutPixelSDL(screen, line[i].x, line[i].y, color);
    }
    std::cout << "END" << std::endl;
  }
}

void VertexShader(const vec4 &v, ivec2 &p)
{
  vec4 _p = R * (v - cameraPos);
  p.x = (FOCAL_LENGTH * (_p.x / _p.z)) + (SCREEN_WIDTH / 2);
  p.y = (FOCAL_LENGTH * (_p.y / _p.z)) + (SCREEN_HEIGHT / 2);
}

void VertexShader(const vec4 &v, Pixel &p)
{
  vec4 _p = R * (v - cameraPos);
  p.zinv = 1 / _p.z; 
  p.x = (FOCAL_LENGTH * (_p.x / _p.z)) + (SCREEN_WIDTH / 2);
  p.y = (FOCAL_LENGTH * (_p.y / _p.z)) + (SCREEN_HEIGHT / 2);
}

void Interpolate(ivec2 a, ivec2 b, vector<ivec2> &result)
{
  int N = result.size();
  vec2 step = vec2(b - a) / float(max(N - 1, 1));
  vec2 current(a);
  for (unsigned int i = 0; i < N; ++i)
  {
    result[i] = round(current);
    current += step;
  }
}

void Interpolate(Pixel a, Pixel b, vector<Pixel> &result)
{
  int N = result.size();
  float div = float(max(N - 1, 1));
  float step_x = (b.x - a.x) / div;
  float step_y = (b.y - a.y) / div;
  float step_z = (b.zinv - a.zinv) / div;
  vec3 step(step_x, step_y, step_z);
  vec3 current(float(a.x), float(a.y), float(a.zinv));
  for (unsigned int i = 0; i < N; ++i)
  {
    result[i].x = round(current.x);
    result[i].y = round(current.y);
    result[i].zinv = current.z;
    current += step;
  }
}

void LoadRotationMatrix()
{
  R[0][0] = cos(theta_y) * cos(theta_z);
  R[0][1] = -1 * cos(theta_x) * sin(theta_z) + (sin(theta_x) * sin(theta_y) * cos(theta_z));
  R[0][2] = sin(theta_x) * sin(theta_z) + cos(theta_x) * sin(theta_y) * cos(theta_z);

  R[1][0] = cos(theta_y) * sin(theta_z);
  R[1][1] = cos(theta_x) * cos(theta_z) + sin(theta_x) * sin(theta_y) * sin(theta_z);
  R[1][2] = -1 * sin(theta_x) * cos(theta_z) + cos(theta_x) * sin(theta_y) * sin(theta_z);

  R[2][0] = -1 * sin(theta_y);
  R[2][1] = sin(theta_x) * cos(theta_y);
  R[2][2] = cos(theta_x) * cos(theta_y);

  R[3][0] = 0;
  R[3][1] = 0;
  R[3][2] = 0;

  R[0][3] = 1;
  R[1][3] = 1;
  R[2][3] = 1;
  R[3][3] = 1;
}

void TestPolyRows()
{
  vector<ivec2> vertexPixels(3);
  vertexPixels[0] = ivec2(10, 5);
  vertexPixels[1] = ivec2(5, 10);
  vertexPixels[2] = ivec2(15, 15);
  vector<ivec2> leftPixels;
  vector<ivec2> rightPixels;
  ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
  for (int row = 0; row < leftPixels.size(); ++row)
  {
    cout << "Start: ("
         << leftPixels[row].x << ","
         << leftPixels[row].y << "). "
         << "End: ("
         << rightPixels[row].x << ","
         << rightPixels[row].y << "). " << endl;
  }
}