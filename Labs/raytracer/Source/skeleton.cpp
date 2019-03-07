#include <stdint.h>
#include <iostream>
#include <tbb/tbb.h>
#include <chrono>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"

using namespace std;
using glm::mat3;
using glm::mat4;
using glm::vec3;
using glm::vec4;

#define SCREEN_WIDTH 1024
#define SCREEN_HEIGHT 1024
#define FULLSCREEN_MODE true
#define STEP (4)
#define LIGHT_COUNT 4

//02 SPEED BEFORE EXTENSION
//95.87% in Closest
//03.25% in DirectLight
//01.27% in Draw

//03 SPEED BEFORE EXTENSION
//95.09% in Closest
//03.60% in DirectLight
//01.70% in Draw

//03 SPEED STAGE 1
//94.34% in Closest
//03.41% in DirectLight
//02.27% in Draw

//CHANGES
//STAGE 1 AUTO for triangles, shortened code, ternary at the end ~ 805 milliseconds
//STAGE 2 parallel_for
//STAGE 2 TRIAL 1: PARALLELIZED Y ~ 210 milliseconds --SUCCESSFUL
//STAGE 2 TRIAL 2: PARALLELIZED X ~ 220 milliseconds --UNSUCCESSFUL
//STAGE 3 SHADOW RETURN TRUE IF HIT ~ 203 milliseconds
//STAGE 4 REPLACED UNNECESSARY COMPUTATION WITH AUTO PREPARED ~190 milliseconds
//STAGE 5 HORIZONTAL INTERPOLATION ~ 80 milliseconds
//STAGE 5 TRIAL 2 VERTICAL DOESNT STACK WITH PARALLELISM
//STAGE 6 HORIZONTAL INTERPOLATION WITH Parallelism ~ 57 milliseconds

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */
bool changed = true;
float focal_length = SCREEN_HEIGHT;
float theta_x = 0.0, theta_y = 0.0, theta_z = 0.0;
vec3 black(0, 0, 0);
vec3 white(1, 1, 1);
vec4 cameraPos(0, 0, -3, 1.0);

vec4 lightPos(0, -0.5, -0.7, 1.0);
vector<vec4> lights;
vec3 lightColor = 14.f * white;
vec3 indirectLight = 0.5f * white;

vec4 zetward(0, 0, 0.1, 1);
vec4 xerward(0.1, 0, 0, 1);

vector<Triangle> triangles;
mat4 R;
/* ----------------------------------------------------------------------------*/
/* DEFINITIONS                                                                 */

typedef std::chrono::high_resolution_clock Clock;

class Intersection
{
public:
  vec4 position;
  float distance;
  int index;

  Intersection()
  {
    position = vec4(0, 0, 0, 1);
    distance = std::numeric_limits<float>::max();
    index = -1;
  }

  Intersection(vec3 p, float d, int i)
  {
    position = vec4(p, 1);
    distance = d;
    index = i;
  }
};

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen *screen);
bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle> &triangles, Intersection &closest);
bool ShadowIntersection(vec4 start, vec4 dir, const vector<Triangle> &triangles, Intersection &closest, vec4 im, vec4 light_pos);
void LoadRotationMatrix();
vec3 DirectLight(const Intersection &i);
float RightMost(const vector<Triangle> &triangles);
float LeftMost(const vector<Triangle> &triangles);
vector<vec3> InterpolateLine(uint32_t a, uint32_t b);
void Interpolate(vec3 a, vec3 b, vector<vec3> &result);
vec3 GetPixel(uint32_t pixel);
void init_cameras();

std::ostream &operator<<(std::ostream &os, vec4 const &v)
{
  return os << "(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w
            << ")";
}

std::ostream &operator<<(std::ostream &os, vec3 const &v)
{
  return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}

std::ostream &operator<<(std::ostream &os, mat4 const &m)
{
  glm::mat4 mt = transpose(m);
  return os << mt[0] << endl
            << mt[1] << endl
            << mt[2] << endl
            << mt[3];
}

int main(int argc, char *argv[])
{

  screen *screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE);

  LoadTestModel(triangles);
  LoadRotationMatrix();

  while (Update())
  {
    if (!changed)
      continue;
    auto t1 = Clock::now();
    Draw(screen);
    auto t2 = Clock::now();
    std::cout << "Delta t2-t1: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
              << " milliseconds" << std::endl;
    SDL_Renderframe(screen);
    changed = false;
  }

  SDL_SaveImage(screen, "screenshot.bmp");

  KillSDL(screen);
  return 0;
}

/*Place your drawing here*/
void Draw(screen *screen)
{
  init_cameras();
  memset(screen->buffer, 0, screen->height * screen->width * sizeof(uint32_t));

  tbb::parallel_for(tbb::blocked_range<int>(0, SCREEN_HEIGHT), [&](tbb::blocked_range<int> r) {
    for (int y = r.begin(); y < r.end(); y++)
    {
      for (int x = 0; x < SCREEN_WIDTH; x += STEP)
      {
        vec4 dir = vec4(x - (SCREEN_WIDTH / 2), y - (SCREEN_HEIGHT / 2), focal_length, 1);
        dir = R * dir;
        Intersection closest;
        vec3 colour(0.0, 0.0, 0.0);
        if (ClosestIntersection(cameraPos, dir, triangles, closest))
        {
          colour = (DirectLight(closest) + indirectLight) * triangles[closest.index].color;
        }
        PutPixelSDL(screen, x, y, colour);
      }
    }
    for (int y = r.begin(); y < r.end(); y++)
    {
      for (int x = 0; x < SCREEN_WIDTH; x += STEP)
      {
        uint32_t pixel_left = screen->buffer[y * SCREEN_WIDTH + x];
        uint32_t pixel_right = pixel_left;
        if (x + STEP < SCREEN_WIDTH)
          pixel_right = screen->buffer[y * SCREEN_WIDTH + x + STEP];

        vector<vec3> interpolated_pixels = InterpolateLine(pixel_left, pixel_right);
        for (int z = 1; z < STEP; ++z)
          PutPixelSDL(screen, x + z, y, interpolated_pixels[z]);
      }
    }
  });

  // for (int i = 0; i < 25; ++i)
  // Stencilize(screen);
}

vector<vec3> InterpolateLine(uint32_t a, uint32_t b)
{
  vector<vec3> interpolated_pixels(STEP);
  vec3 color_a = GetPixel(a);
  vec3 color_b = GetPixel(b);
  Interpolate(color_a, color_b, interpolated_pixels);
  return interpolated_pixels;
}

vec3 GetPixel(uint32_t pixel)
{
  float r = ((float)((pixel >> 16) & 0xFF)) / 255.0f;
  float g = ((float)((pixel >> 8) & 0xFF)) / 255.0f;
  float b = ((float)((pixel)&0xFF)) / 255.0f;
  return vec3(r, g, b);
}

void Interpolate(vec3 a, vec3 b, vector<vec3> &result)
{
  uint32_t len = result.size();

  if (len == 1)
  {
    result[0] = a;
    return;
  }

  float step_x = (b.x - a.x) / (len - 1);
  float step_y = (b.y - a.y) / (len - 1);
  float step_z = (b.z - a.z) / (len - 1);

  for (int i = 0; i < len; ++i)
  {
    result[i].x = a.x + (step_x * i);
    result[i].y = a.y + (step_y * i);
    result[i].z = a.z + (step_z * i);
  }
}

bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle> &triangles, Intersection &closest)
{
  auto lv3d = length(vec3(dir));
  auto nv3d = -vec3(dir);
  for (int i = 0; i < triangles.size(); ++i)
  {
    auto triangle = triangles[i];
    vec3 b = vec3(start - triangle.v0);
    mat3 A(nv3d, triangle.e1, triangle.e2);
    vec3 x = glm::inverse(A) * b;
    float dist = x.x / lv3d;
    if (x.x >= 0.0f &&
        x.y >= 0.0f &&
        x.z >= 0.0f &&
        x.y + x.z <= 1 &&
        closest.distance > dist)
      closest = Intersection(vec3(start + x.x * dir), dist, i);
  }

  return (closest.index == -1) ? false : true;
}

bool ShadowIntersection(vec4 start, vec4 dir, const vector<Triangle> &triangles, Intersection &closest, vec4 im, vec4 light_pos)
{

  auto lsim = length(light_pos - im);
  auto lv3d = length(vec3(dir));
  auto nv3d = -vec3(dir);
  for (int i = 0; i < triangles.size(); ++i)
  {
    auto triangle = triangles[i];
    vec3 b = vec3(start - triangle.v0);
    mat3 A(nv3d, triangle.e1, triangle.e2);
    vec3 x = glm::inverse(A) * b;
    float dist = x.x / lv3d;
    if (x.x >= 0.0f &&
        x.y >= 0.0f &&
        x.z >= 0.0f &&
        x.y + x.z <= 1 &&
        closest.distance > dist)
    {
      closest.distance = dist;
      if ((length((start + x.x * dir) - im) < lsim))
        return true;
    }
  }

  return false;
}

// bool TriangleIntersection(vec4 start, vec4 dir, const vector<Triangle> &triangles, Intersection &closest, vec4 im, int t_index, )
// {
//   auto lsim = length(light_pos- im);
//   auto lv3d = length(vec3(dir));
//   auto nv3d = -vec3(dir);

//   auto triangle = triangles[t_index];
//   vec3 b = vec3(start - triangle.v0);
//   mat3 A(nv3d, triangle.e1, triangle.e2);
//   vec3 x = glm::inverse(A) * b;
//   float dist = x.x / lv3d;
//   if (x.x >= 0.0f &&
//       x.y >= 0.0f &&
//       x.z >= 0.0f &&
//       x.y + x.z <= 1 )
//   {
//     return true;
//   }
//   return false;
// }

vec3 DirectLight(const Intersection &i)
{
  Intersection closest;
  vec3 r_v = vec3(lightPos - i.position);

  float count = lights.size();
  if (ShadowIntersection(i.position + 0.001f * triangles[i.index].normal, vec4(r_v, 1), triangles, closest, i.position, lightPos))
  {
    count = 0;
    for (int it = 0; it < lights.size(); it++)
    {
      Intersection c = closest;
      if (!ShadowIntersection(i.position + 0.001f * triangles[i.index].normal, vec4(vec3(lights[it] - i.position), 1), triangles, c, i.position, lights[it]))
      {
        count++;
      }
    }
  }
  float r = sqrt(pow(r_v.x, 2.0) + pow(r_v.y, 2.0) + pow(r_v.z, 2.0));
  r_v = glm::normalize(r_v);
  vec3 n_u = glm::normalize(vec3(triangles[i.index].normal));
  return ((lightColor / (12.56f * r * r)) * sqrt(glm::dot(r_v, n_u) * glm::dot(r_v, n_u))) * (count / lights.size());
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

/*Place updates of parameters here*/
bool Update()
{
  // static int t = SDL_GetTicks();
  /* Compute frame time */
  // int t2 = SDL_GetTicks();
  // float dt = float(t2 - t);
  // t = t2;

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
        cameraPos += R * zetward;
        changed = true;
        break;
      case SDLK_DOWN:
        cameraPos -= R * zetward;
        changed = true;
        break;
      case SDLK_LEFT:
        cameraPos -= R * xerward;
        changed = true;
        break;
      case SDLK_RIGHT:
        cameraPos += R * xerward;
        changed = true;
        break;
      case SDLK_w:
        theta_x -= 0.1f;
        LoadRotationMatrix();
        changed = true;
        break;
      case SDLK_s:
        theta_x += 0.1f;
        LoadRotationMatrix();
        changed = true;
        break;
      case SDLK_a:
        theta_y += 0.1f;
        LoadRotationMatrix();
        changed = true;
        break;
      case SDLK_d:
        theta_y -= 0.1f;
        LoadRotationMatrix();
        changed = true;
        break;
      case SDLK_KP_4:
        lightPos.x -= 0.1f;
        changed = true;
        break;
      case SDLK_KP_6:
        lightPos.x += 0.1f;
        changed = true;
        break;
      case SDLK_KP_8:
        lightPos.y -= 0.1f;
        changed = true;
        break;
      case SDLK_KP_2:
        lightPos.y += 0.1f;
        changed = true;
        break;
      case SDLK_KP_5:
        lightPos.z -= 0.1f;
        changed = true;
        break;
      case SDLK_ESCAPE:
        /* Move camera quit */
        return false;
      }
    }
  }
  return true;
}

void init_cameras()
{
  lights.clear();
  for (int i = -LIGHT_COUNT / 2; i < LIGHT_COUNT / 2; i++)
  {
    for (int j = -LIGHT_COUNT / 2; j < LIGHT_COUNT / 2; j++)
    {
      lights.push_back(vec4(lightPos.x + 0.1f * i, lightPos.y, lightPos.z + 0.1f * j, 1));
    }
  }
}

// void Stencilize(screen *screen)
// {
//   uint32_t * new_buffer = (uint32_t *) malloc(sizeof(uint32_t) * SCREEN_WIDTH * SCREEN_HEIGHT);

//   for(int x = 0; x < SCREEN_WIDTH; x++) {
//     new_buffer[x] = screen->buffer[x];
//     new_buffer[(SCREEN_HEIGHT - 1) * SCREEN_WIDTH + x] = screen->buffer[(SCREEN_HEIGHT - 1) * SCREEN_WIDTH + x];
//   }

//   for (int y = 0; y < SCREEN_HEIGHT; y++) {
//     new_buffer[SCREEN_WIDTH * y] = screen->buffer[SCREEN_WIDTH * y];
//     new_buffer[SCREEN_WIDTH * (y + 1) - 1] = screen->buffer[SCREEN_WIDTH * (y + 1) - 1];
//   }

//   for (int y = 1; y < SCREEN_HEIGHT - 1; y++)
//   {
//     for (int x = 1; x < SCREEN_WIDTH - 1; x++)
//     {
//       new_buffer[y * SCREEN_WIDTH + x] = WeightColor(screen->buffer[y * SCREEN_WIDTH + x], 0.2);
//       new_buffer[y * SCREEN_WIDTH + x] += WeightColor(screen->buffer[y * SCREEN_WIDTH + x - 1], 0.1);
//       new_buffer[y * SCREEN_WIDTH + x] += WeightColor(screen->buffer[y * SCREEN_WIDTH + x + 1], 0.1);
//       new_buffer[y * SCREEN_WIDTH + x] += WeightColor(screen->buffer[(y - 1) * SCREEN_WIDTH + x], 0.1);
//       new_buffer[y * SCREEN_WIDTH + x] += WeightColor(screen->buffer[(y - 1) * SCREEN_WIDTH + x - 1], 0.1);
//       new_buffer[y * SCREEN_WIDTH + x] += WeightColor(screen->buffer[(y - 1) * SCREEN_WIDTH + x + 1], 0.1);
//       new_buffer[y * SCREEN_WIDTH + x] += WeightColor(screen->buffer[(y + 1) * SCREEN_WIDTH + x], 0.1);
//       new_buffer[y * SCREEN_WIDTH + x] += WeightColor(screen->buffer[(y + 1) * SCREEN_WIDTH + x - 1], 0.1);
//       new_buffer[y * SCREEN_WIDTH + x] += WeightColor(screen->buffer[(y + 1) * SCREEN_WIDTH + x + 1], 0.1);
//     }
//   }
//   memcpy(screen->buffer, new_buffer, sizeof(uint32_t) * SCREEN_WIDTH * SCREEN_HEIGHT);
// }

// uint32_t WeightColor(uint32_t color, float weight) {
//   float red = weight * ((float)((color >> 16) & 0xFF) / 255.0f);
//   float green = weight * ((float)((color >> 8) & 0xFF) / 255.0f);
//   float blue = weight * ((float)((color) & 0xFF) / 255.0f);
//   uint32_t r = uint32_t( glm::clamp( 255*red, 0.f, 255.f ) );
//   uint32_t g = uint32_t( glm::clamp( 255*green, 0.f, 255.f ) );
//   uint32_t b = uint32_t( glm::clamp( 255*blue, 0.f, 255.f ) );
//   return (128<<24) + (r<<16) + (g<<8) + b;
// }