#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <stdint.h>

using namespace std;
using glm::vec3;
using glm::mat3;

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FOCAL_LENGTH (SCREEN_HEIGHT / 2)
#define FULLSCREEN_MODE false


/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */
int t;
vector<vec3> stars(1000);


/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw(screen* screen);
void Interpolate(float a, float b, vector<float>& result);
void Interpolate(vec3 a, vec3 b, vector<vec3>& result);
void Test_Interpolate();

int main( int argc, char* argv[] ) {

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  t = SDL_GetTicks();	/*Set start value for timer.*/

  //Test_Interpolate();



  for(int i = 0; i < 1000; ++i) {
    float val_x = float(rand()) / float(RAND_MAX);
    float val_y = float(rand()) / float(RAND_MAX);
    float val_z = float(rand()) / float(RAND_MAX);
    bool sign_x = (float(rand()) / float(RAND_MAX)) >= 0.5;
    bool sign_y = (float(rand()) / float(RAND_MAX)) >= 0.5;
    val_x *= (sign_x * -2) + 1;
    val_y *= (sign_y * -2) + 1;
    stars[i].x = val_x;
    stars[i].y = val_y;
    stars[i].z = val_z;
  }

  while( NoQuitMessageSDL() )
    {
      Draw(screen);
      Update();
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}

/*Place your drawing here*/
void Draw(screen* screen){
  /* Clear buffer */
  //memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  for (int i = 0; i < SCREEN_HEIGHT; ++i) {
    for (int j = 0; j < SCREEN_WIDTH; ++j) {
      uint32_t rgb = screen->buffer[i * SCREEN_WIDTH + j];
      uint32_t alpha = rgb & 0xFF000000;
      uint32_t red = (rgb >> 16) & 0xFF;
      uint32_t green = (rgb >> 8) & 0xFF;
      uint32_t blue = rgb & 0xFF;
      red *= 0.8;
      green *= 0.8;
      blue *= 0.8;
      screen->buffer[i * SCREEN_WIDTH + j] = alpha + (red << 16) + (green << 8) + blue;
    }
  }

  vec3 white(1.0, 1.0, 1.0);

  for (size_t s = 0; s < stars.size(); ++s) {
    float u = FOCAL_LENGTH * (stars[s].x / stars[s].z) + (SCREEN_WIDTH / 2);
    float v = FOCAL_LENGTH * (stars[s].y / stars[s].z) + (SCREEN_HEIGHT / 2);
    vec3 color = (0.2f * white) / (stars[s].z * stars[s].z);
    PutPixelSDL(screen, u, v, color);
  }
}

/*Place updates of parameters here*/
void Update(){
  static int t = SDL_GetTicks();
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;

  static float v = 0.001;

  /* Update variables*/
  for (size_t s = 0; s < stars.size(); ++s) {
    stars[s].z -= v * dt;
    if (stars[s].z <= 0)
      stars[s].z += 1;
    if (stars[s].z > 1)
      stars[s].z -= 1;
  }
}

void Interpolate(float a, float b, vector<float>& result) {
  uint32_t len = result.size();

  if (len == 1) {
    result[0] = a;
    return;
  }

  float step = (b - a) / (len - 1);
  for (int i = 0; i < len; ++i)
    result[i] = a + (step * i);
}

void Interpolate(vec3 a, vec3 b, vector<vec3>& result) {
  uint32_t len = result.size();

  if (len == 1) {
    result[0] = a;
    return;
  }

  float step_x = (b.x - a.x) / (len - 1);
  float step_y = (b.y - a.y) / (len - 1);
  float step_z = (b.z - a.z) / (len - 1);

  for (int i = 0; i < len; ++i) {
    result[i].x = a.x + (step_x * i);
    result[i].y = a.y + (step_y * i);
    result[i].z = a.z + (step_z * i);
  }
}

void Test_Interpolate() {
  vector<float> result(10);
  Interpolate(5,14,result);
  for (int i = 0; i < result.size(); ++i)
    std::cout << result[i] << " ";
  std::cout << std::endl;

  vector<vec3> result2(4);
  vec3 a(1,4,9.2);
  vec3 b(4,1,9.8);
  Interpolate(a,b, result2);
  for(int i = 0; i < result2.size(); ++i) {
    std::cout << "( " << result2[i].x << ", " << result2[i].y << ", " << result2[i].z << " ) ";
  }
}

//Drawing the light spectrum.

// void Draw(screen* screen){
//   /* Clear buffer */
//   memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
//
//   vec3 topLeft(1,0,0);    // red
//   vec3 topRight(0,0,1);   // blue
//   vec3 bottomRight(0,1,0);// green
//   vec3 bottomLeft(1,1,0); // yellow
//
//   vector<vec3> leftSide( SCREEN_HEIGHT );
//   vector<vec3> rightSide( SCREEN_HEIGHT );
//   Interpolate(topLeft, bottomLeft, leftSide);
//   Interpolate(topRight, bottomRight, rightSide);
//
//   for(int i=0; i < SCREEN_HEIGHT; i++)
//     {
//       vector<vec3> line ( SCREEN_WIDTH );
//       Interpolate(leftSide[i], rightSide[i], line);
//       for (int j = 0; j < SCREEN_WIDTH; ++j)
//         PutPixelSDL(screen, j, i, line[j]);
//     }
// }
