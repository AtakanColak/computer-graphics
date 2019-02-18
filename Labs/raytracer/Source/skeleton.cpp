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

#define SCREEN_WIDTH 100
#define SCREEN_HEIGHT 100
#define FULLSCREEN_MODE true

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */
float focal_length = SCREEN_HEIGHT;
float yaw = 0;
float xaw;
mat4 R;
vec3 white(0, 0, 0);
vec4 cameraPos(0, 0, -3, 1.0);
vector<Triangle> triangles;

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





int main(int argc, char *argv[])
{

  screen *screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE);


  LoadTestModel(triangles);

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
      Intersection closest;
      vec3 colour(0.0, 0.0, 0.0);
      if (ClosestIntersection(cameraPos, dir, triangles, closest))
      {
        colour = triangles[closest.triangleIndex].color;
      }
      PutPixelSDL(screen, x, y, colour);
    }
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
        cameraPos.z += 0.1;
        break;
      case SDLK_DOWN:
        cameraPos.z -= 0.1;
        break;
      case SDLK_LEFT:
        cameraPos.x -= 0.1;
        break;
      case SDLK_RIGHT:
        cameraPos.x += 0.1;
        break;
      case SDLK_w:
        //Rotate(0, 1);
        break;
      case SDLK_s:
        //Rotate(0, -1);
        break;
      case SDLK_a:
        //Rotate(-1, 0);
        break;
      case SDLK_d:
        //Rotate(1, 0);
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

  for (int i = 0; i < triangles.size(); ++i)
  {
    Triangle t = triangles[i];
    vec4 v0 = t.v0;
    vec4 v1 = t.v1;
    vec4 v2 = t.v2;

    vec3 e1 = vec3(v1.x - v0.x, v1.y - v0.y, v1.z - v0.z);
    vec3 e2 = vec3(v2.x - v0.x, v2.y - v0.y, v2.z - v0.z);

    vec3 b = vec3(start.x - v0.x, start.y - v0.y, start.z - v0.z);

    vec3 d(dir.x, dir.y, dir.z);

    mat3 A(-d, e1, e2);

    vec3 x = glm::inverse(A) * b;

    float t_dist = x.x;
    float u = x.y;
    float v = x.z;

    if (closest.distance > t_dist && u >= 0 && v >= 0 && u + v <= 1)
    {
      closest.position = vec4(x.x, x.y, x.z, 1);
      closest.distance = t_dist;
      closest.triangleIndex = i;
    }
  }

  if (closest.triangleIndex == -1)
    return false;

  return true;
}