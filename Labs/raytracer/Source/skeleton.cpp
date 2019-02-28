#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>

using namespace std;
using glm::mat3;
using glm::mat4;
using glm::vec3;
using glm::vec4;

#define SCREEN_WIDTH 500
#define SCREEN_HEIGHT 500
#define FULLSCREEN_MODE true

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */
float focal_length = SCREEN_HEIGHT;
float theta_x = 0.0, theta_y = 0.0, theta_z = 0.0;
vec3 black(0, 0, 0);
vec3 white(1, 1, 1);
vec4 cameraPos(0, 0, -3, 1.0);

vec4 lightPos(0, -0.5, -0.7, 1.0);
vec3 lightColor = 14.f * white;
vec3 indirectLight = 0.5f * white;

vec4 zetward(0, 0, 0.1, 1);
vec4 xerward(0.1, 0, 0, 1);

vector<Triangle> triangles;
mat4 R;
/* ----------------------------------------------------------------------------*/
/* DEFINITIONS                                                                 */
struct Intersection
{
  vec4 position;
  float distance;
  int triangleIndex;
};

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen *screen);
bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle> &triangles, Intersection &closest);
void LoadRotationMatrix();
vec3 DirectLight(const Intersection &i);

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

  for (int y = 0; y < SCREEN_HEIGHT; ++y)
  {
    for (int x = 0; x < SCREEN_WIDTH; ++x)
    {
      vec4 dir = vec4(x - (SCREEN_WIDTH / 2), y - (SCREEN_HEIGHT / 2), focal_length, 1);

      dir = R * dir;

      Intersection closest;
      vec3 colour(0.0, 0.0, 0.0);
      if (ClosestIntersection(cameraPos, dir, triangles, closest))
      {
        colour = (DirectLight(closest) + indirectLight) * triangles[closest.triangleIndex].color;
      }
      PutPixelSDL(screen, x, y, colour);
    }
  }
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
        break;
      case SDLK_DOWN:
        cameraPos -= R * zetward;
        break;
      case SDLK_LEFT:
        cameraPos -= R * xerward;
        break;
      case SDLK_RIGHT:
        cameraPos += R * xerward;
        break;
      case SDLK_w:
        theta_x -= 0.1f;
        LoadRotationMatrix();
        break;
      case SDLK_s:
        theta_x += 0.1f;
        LoadRotationMatrix();
        break;
      case SDLK_a:
        theta_y += 0.1f;
        LoadRotationMatrix();
        break;
      case SDLK_d:
        theta_y -= 0.1f;
        LoadRotationMatrix();
        break;
      case SDLK_KP_4:
        lightPos.x -= 0.1f;
        break;
      case SDLK_KP_6:
        lightPos.x += 0.1f;
        break;
      case SDLK_KP_8:
        lightPos.y -= 0.1f;
        break;
      case SDLK_KP_2:
        lightPos.y += 0.1f;
        break;
      case SDLK_KP_5:
        lightPos.z -= 0.1f;
        break;
      case SDLK_ESCAPE:
        /* Move camera quit */
        return false;
      }
    }
  }
  return true;
}

bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle> &triangles, Intersection &closest)
{
  closest.position = vec4(0, 0, 0, 1);
  closest.distance = std::numeric_limits<float>::max();
  closest.triangleIndex = -1;

  //dir = vec4(glm::normalize(vec3(dir)), 1.0f);

  for (unsigned int i = 0; i < triangles.size(); ++i)
  {
    Triangle t = triangles[i];
    vec4 v0 = t.v0;

    vec3 b = vec3(start - v0);
    mat3 A(-vec3(dir), t.e1, t.e2);

    vec3 x = glm::inverse(A) * b;

    float t_dist = x.x;
    float u = x.y;
    float v = x.z;

    if (t_dist >= 0.0 && closest.distance > (t_dist / length(vec3(dir))) && u >= 0.0 && v >= 0.0 && u + v <= 1)
    {
      closest.distance = t_dist / length(vec3(dir));
      closest.position = start + t_dist * dir;
      closest.position.w = 1;
      closest.triangleIndex = i;
    }
  }

  if (closest.triangleIndex == -1)
    return false;

  return true;
}

vec3 DirectLight(const Intersection &i)
{
  Intersection closest;
  vec3 r_v = vec3(lightPos - i.position);
  if (ClosestIntersection(i.position + 0.001f * triangles[i.triangleIndex].normal, vec4(r_v, 1), triangles, closest) && (length(closest.position - i.position) < length(lightPos - i.position)))
  {
    return vec3(0, 0, 0);
  }
  float r = sqrt(pow(r_v.x, 2.0) + pow(r_v.y, 2.0) + pow(r_v.z, 2.0));
  r_v = glm::normalize(r_v);
  vec3 n_u = glm::normalize(vec3(triangles[i.triangleIndex].normal));
  return (lightColor / (4.0f * 3.14f * r * r)) * sqrt(glm::dot(r_v, n_u) * glm::dot(r_v, n_u)); //max(glm::dot(r_v, n_u), 0.0f);
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