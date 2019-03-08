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

SDL_Event event;

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FOCAL_LENGTH (SCREEN_HEIGHT)
#define FULLSCREEN_MODE false

/* ----------------------------------------------------------------------------*/
/* VARIABLES                                                                   */
vec3 black(0, 0, 0);
vec3 white(1, 1, 1);
vec4 cameraPos(0,0,-3.001,1);
vector<Triangle> triangles;
mat4 R;
float theta_x = 0.1, theta_y = 0.1, theta_z = 0.1;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen *screen);
void VertexShader(const vec4 &v, glm::ivec2 &p);
void LoadRotationMatrix();

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

  for (uint32_t i = 0; i < triangles.size(); ++i)
  {
    vector<vec4> vertices(3);
    vertices[0] = triangles[i].v0;
    vertices[1] = triangles[i].v1;
    vertices[2] = triangles[i].v2;
    for (int v = 0; v < 3; ++v)
    {
      glm::ivec2 projPos;
      VertexShader(vertices[v], projPos);
      // if(projPos.x < 0 || projPos.x >= SCREEN_WIDTH) continue;
      // if(projPos.y < 0 || projPos.y >= SCREEN_HEIGHT) continue;
      PutPixelSDL(screen, projPos.x, projPos.y, white);
      std::cout << "(" << projPos.x << "," << projPos.y << ")" << std::endl;
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
        /* Move camera forward */
        break;
      case SDLK_DOWN:
        /* Move camera backwards */
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
void VertexShader(const vec4 &v, glm::ivec2 &p)
{
  vec4 _p = R * (v - cameraPos);
  p.x = (FOCAL_LENGTH * (_p.x / _p.z)) + (SCREEN_WIDTH / 2);
  p.y = (FOCAL_LENGTH * (_p.y / _p.z)) + (SCREEN_HEIGHT / 2);
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