#include <iostream>
#include <chrono>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <tbb/tbb.h>
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

#define SCREEN_WIDTH 1024
#define SCREEN_HEIGHT 1024

/* ----------------------------------------------------------------------------*/
/* VARIABLES                                                                   */
float theta_x = 0, theta_y = 0, theta_z = 0;
float depthBuffer[SCREEN_HEIGHT * SCREEN_WIDTH];
int FOCAL_LENGTH = SCREEN_HEIGHT;
mat4 R;
vector<Triangle> triangles;
vec3 black(0, 0, 0);
vec3 white(1, 1, 1);
vec4 cameraPos(0, 0, -3.001, 1);
vec3 lightPos(0, -0.5, -0.7);
vec3 lightPower = 14.1f * white;
vec3 indirect_light = 0.5f * white;
bool change = true;
/* ----------------------------------------------------------------------------*/
/* STRUCTS and DEFINITIONS                                                                    */

#define FULLSCREEN_MODE false

typedef std::chrono::high_resolution_clock Clock;

struct Pixel
{
  int x;
  int y;
  float zinv;
  vec4 pos3d;
};

struct Vertex
{
  vec4 position;
};

//1024 X 1024 OPTIMISATION
//INITIAL ~430 ms
//STAGE 1 O3 ~53 ms
//STAGE 2 PARALLELISE COMPUTE ROWS - NOT SUCCESSFUL ~53 ms
//STAGE 3 PARALLELISE DRAW LOSSLY SUCCESSFUL ~34 ms

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void LoadRotationMatrix();

void Draw(screen *screen);

void VertexShader(const Vertex &v, Pixel &p);

void PixelShader(screen *screen, const Pixel &p, vec3 currentReflectance, vec4 currentNormal);

void Interpolate(Pixel a, Pixel b, vector<Pixel> &result);

void DrawLineSDL(screen *screen, Pixel a, Pixel b, vec3 color, vec4 normal);

void ComputePolygonRows(const vector<Pixel> &vertexPixels, vector<Pixel> &leftPixels, vector<Pixel> &rightPixels);

void DrawPolygonRows(screen *screen, const vector<Pixel> &leftPixels, const vector<Pixel> &rightPixels, vec3 color, vec4 normal);

void DrawPolygonPixel(screen *screen, const vector<Vertex> &vertices, vec3 color, vec4 normal);

int main(int argc, char *argv[])
{

  screen *screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE);
  LoadTestModel(triangles);
  LoadRotationMatrix();

  while (Update())
  {
    if (!change)
      continue;
    auto t1 = Clock::now();
    Draw(screen);
    auto t2 = Clock::now();
    std::cout << "Delta t2-t1: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
              << " milliseconds" << std::endl;
    SDL_Renderframe(screen);
    change = false;
  }

  SDL_SaveImage(screen, "screenshot.bmp");

  KillSDL(screen);
  return 0;
}

/*Place your drawing here*/
void Draw(screen *screen)
{
  /* Clear buffers */
  memset(screen->buffer, 0, screen->height * screen->width * sizeof(uint32_t));
  memset(depthBuffer, 0, screen->height * screen->width * sizeof(float));
  tbb::parallel_for(tbb::blocked_range<int>(0, triangles.size()), [&](tbb::blocked_range<int> r) {
    for (uint32_t i = r.begin(); i < r.end(); ++i)
    {
      vector<Vertex> vertices(3);
      vertices[0].position = triangles[i].v0;
      vertices[1].position = triangles[i].v1;
      vertices[2].position = triangles[i].v2;
      DrawPolygonPixel(screen, vertices, triangles[i].color, triangles[i].normal);
    }
  });
}

void DrawPolygonPixel(screen *screen, const vector<Vertex> &vertices, vec3 color, vec4 normal)
{
  int V = vertices.size();
  vector<Pixel> vertexPixels(V);
  for (int i = 0; i < V; ++i)
    VertexShader(vertices[i], vertexPixels[i]);
  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;
  ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
  DrawPolygonRows(screen, leftPixels, rightPixels, color, normal);
}

void DrawPolygonRows(screen *screen, const vector<Pixel> &leftPixels, const vector<Pixel> &rightPixels, vec3 color, vec4 normal)
{
  for (int i = 0; i < leftPixels.size(); ++i)
    DrawLineSDL(screen, leftPixels[i], rightPixels[i], color, normal);
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
  //

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

void DrawLineSDL(screen *screen, Pixel a, Pixel b, vec3 color, vec4 normal)
{
  ivec2 delta = glm::abs(ivec2(a.x, a.y) - ivec2(b.x, b.y));
  int pixels = glm::max(delta.x, delta.y) + 1;
  vector<Pixel> line(pixels);
  Interpolate(a, b, line);
  for (int i = 0; i < pixels; ++i)
  {
    PixelShader(screen, line[i], color, normal);
  }
}

void VertexShader(const Vertex &v, Pixel &p)
{
  vec4 _p = R * (v.position - cameraPos);
  p.zinv = 1 / _p.z;
  p.x = int(FOCAL_LENGTH * _p.x * p.zinv) + (SCREEN_WIDTH / 2);
  p.y = int(FOCAL_LENGTH * _p.y * p.zinv) + (SCREEN_HEIGHT / 2);
  p.pos3d = v.position;
}

void PixelShader(screen *screen, const Pixel &p, vec3 currentReflectance, vec4 currentNormal)
{
  int index = p.y * SCREEN_WIDTH + p.x;
  if (index > SCREEN_HEIGHT * SCREEN_WIDTH - 1 || index < 0)
    return;
  if (p.zinv > depthBuffer[index])
  {
    depthBuffer[index] = p.zinv;
    vec3 color = black;
    float r = glm::distance(lightPos, vec3(p.pos3d));
    vec3 direction = vec3(lightPos - vec3(p.pos3d));
    float rn = glm::dot(glm::normalize(direction), glm::normalize(vec3(currentNormal)));
    vec3 direct_light = (lightPower * glm::max(rn, (float)0.0)) / (4 * M_PI * r * r);
    color = currentReflectance * (direct_light + indirect_light);
    PutPixelSDL(screen, p.x, p.y, color);
  }
}

void Interpolate(Pixel a, Pixel b, vector<Pixel> &result)
{
  int N = result.size();
  float div = float(max(N - 1, 1));
  float step_x = (b.x - a.x) / div;
  float step_y = (b.y - a.y) / div;
  float step_z = (b.zinv - a.zinv) / div;
  vec4 pos_current = a.pos3d * a.zinv;
  vec4 pos_step = ((b.pos3d * b.zinv) - (a.pos3d * a.zinv)) / div;
  vec3 step(step_x, step_y, step_z);
  vec3 current(float(a.x), float(a.y), float(a.zinv));
  for (unsigned int i = 0; i < N; ++i)
  {
    result[i].x = round(current.x);
    result[i].y = round(current.y);
    result[i].zinv = current.z;
    result[i].pos3d = pos_current / result[i].zinv;
    pos_current += pos_step;
    current += step;
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
        cameraPos.z += 0.1f;
        change = true;
        break;
      case SDLK_DOWN:
        cameraPos.z -= 0.1f;
        change = true;
        break;
      case SDLK_w:
        lightPos.z += 0.1f;
        change = true;
        break;
      case SDLK_s:
        lightPos.z -= 0.1f;
        change = true;
        break;
      case SDLK_a:
        lightPos.x -= 0.1f;
        change = true;
        break;
      case SDLK_d:
        lightPos.x += 0.1f;
        change = true;
        break;
      case SDLK_ESCAPE:
        /* Move camera quit */
        return false;
      }
    }
  }
  return true;
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

// void TestPolyRows()
// {
//   vector<ivec2> vertexPixels(3);
//   vertexPixels[0] = ivec2(10, 5);
//   vertexPixels[1] = ivec2(5, 10);
//   vertexPixels[2] = ivec2(15, 15);
//   vector<ivec2> leftPixels;
//   vector<ivec2> rightPixels;
//   ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
//   for (int row = 0; row < leftPixels.size(); ++row)
//   {
//     cout << "Start: ("
//          << leftPixels[row].x << ","
//          << leftPixels[row].y << "). "
//          << "End: ("
//          << rightPixels[row].x << ","
//          << rightPixels[row].y << "). " << endl;
//   }
// }