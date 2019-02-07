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
#define FULLSCREEN_MODE false


/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */
int t;

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
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  vec3 topLeft(1,0,0);    // red
  vec3 topRight(0,0,1);   // blue
  vec3 bottomRight(0,1,0);// green
  vec3 bottomLeft(1,1,0); // yellow

  vector<vec3> leftSide( SCREEN_HEIGHT );
  vector<vec3> rightSide( SCREEN_HEIGHT );
  Interpolate(topLeft, bottomLeft, leftSide);
  Interpolate(topRight, bottomRight, rightSide);

  for(int i=0; i < SCREEN_HEIGHT; i++)
    {
      vector<vec3> line ( SCREEN_WIDTH );
      Interpolate(leftSide[i], rightSide[i], line);
      for (int j = 0; j < SCREEN_WIDTH; ++j)
        PutPixelSDL(screen, j, i, line[j]);
    }
}

/*Place updates of parameters here*/
void Update(){
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;
  /*Good idea to remove this*/
  //std::cout << "Render time: " << dt << " ms." << std::endl;
  /* Update variables*/
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
