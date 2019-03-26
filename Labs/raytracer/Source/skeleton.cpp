#include <stdint.h>
#include <iostream>
#include <tbb/tbb.h>
#include <chrono>
#include <random>
#include <glm/glm.hpp>
#include <glm/gtx/vector_angle.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"

using namespace std;
using glm::mat3;
using glm::mat4;
using glm::vec3;
using glm::vec4;

#define SCREEN_WIDTH 512
#define SCREEN_HEIGHT 512
#define FULLSCREEN_MODE false
#define STEP 1
#define LIGHT_COUNT 1
#define AA 1
#define BOUNCE 0
// #define BOUNCE_POWER 2
#define MTF90THETA (false)
#define BRDF (false)
#define BRDF_RANGE 5.01f
#define DINIT 0.9999f
#define DEATH 0.0f

//CHANGES
//TIMES FOR 1024 X 1024
//STAGE 1 AUTO for triangles, shortened code, ternary at the end                ~  800 milliseconds
//STAGE 2 parallel_for
//STAGE 2 TRIAL 1: PARALLELIZED Y                                               ~  210 milliseconds
//STAGE 2 TRIAL 2: PARALLELIZED X                                               ~  220 milliseconds
//STAGE 3 SHADOW RETURN TRUE IF HIT                                             ~  200 milliseconds
//STAGE 4 REPLACED UNNECESSARY COMPUTATION WITH AUTO PREPARED                   ~  190 milliseconds
//STAGE 5 HORIZONTAL INTERPOLATION                                              ~   80 milliseconds
//STAGE 5 TRIAL 2 VERTICAL DOESNT STACK WITH PARALLELISM
//STAGE 6 HORIZONTAL INTERPOLATION WITH Parallelism                             ~   60 milliseconds
//STAGE 7 SOFT SHADOWS LIGHTCOUNT = 3 (27 light sources)                        ~  730 milliseconds
//STAGE 8 ANTI ALIZING at AA 3                                                  ~ 6400 milliseconds
//STAGE 9 LIGHT BOUNCES W MTFTHETA90 PROTOCOL INTRODUCED                        ~69500 milliseconds (bounce light count aa is 3, step is 1)
//STAGE 10 BOUNCE REWORKED INTO PROBABILISTIC
//SAMPLING ON A SPHERE
//DIFFUSE
//PROBABILITY BOUNCE

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */
bool changed = true;
float focal_length = SCREEN_HEIGHT;
float theta_x = 0.0, theta_y = 0.0, theta_z = 0.0;
mat4 R;
vec4 zetward(0, 0, 0.1, 1);
vec4 xerward(0.1, 0, 0, 1);
vec4 cameraPos(0, 0, -3, 1.0);
vec4 lightPos1(0, -0.5, -0.7, 1.0);
vector<vec4> light_1;
vector<Triangle> triangles;
vec3 black(0, 0, 0), white(1, 1, 1);
vec3 lightColor = 20.f * white;
vec3 indirectLight = 0.5f * white;

/* ----------------------------------------------------------------------------*/
/* DEFINITIONS                                                                 */

typedef std::chrono::high_resolution_clock Clock;

class Intersection
{
  public:
	vec4 position;
	float distance;
	int index;
	bool sphere;

	Intersection()
	{
		position = vec4(0, 0, 0, 1);
		distance = std::numeric_limits<float>::max();
		index = -1;
		sphere = false;
	}

	Intersection(vec3 p, float d, int i)
	{
		position = vec4(p, 1);
		distance = d;
		index = i;
		sphere = false;
	}
};

class Sphere
{
  public:
	vec4 center;
	vec4 color;
	float radius;
	float radius2;
	Sphere()
	{
		center = vec4(0, -0.3, 0, 1.0);
		radius = 0.3f;
		radius2 = radius * radius;
		color = vec4(1, 0, 0, 0);
	}
};

vector<Sphere> spheres;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen *screen);

void LoadRotationMatrix();
// vec3 DirectLight(const Intersection &i, int bounce, vec4 I);
vector<vec3> InterpolateLine(uint32_t a, uint32_t b);
void Interpolate(vec3 a, vec3 b, vector<vec3> &result);
vec3 GetPixel(uint32_t pixel);
void init_light_1();
vec4 IntersectionNormal(const Intersection &i);
vec4 IntersectionColor(const Intersection &i);
bool TriangleIntersection(vec3 start, vec3 dir, Intersection &in);
bool SphereIntersection(vec3 start, vec3 dir, Intersection &in);
bool GetIntersection(vec3 start, vec3 dir, Intersection &in);
vec3 LightRay(Intersection &i, vec4 I, int b, float d);

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
	init_light_1();
	Sphere s = Sphere();
	spheres.push_back(s);
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

	memset(screen->buffer, 0, screen->height * screen->width * sizeof(uint32_t));

	tbb::parallel_for(tbb::blocked_range<int>(0, SCREEN_HEIGHT), [&](tbb::blocked_range<int> r) {
		for (int y = r.begin(); y < r.end(); y++)
		{
			for (int x = 0; x < SCREEN_WIDTH; x += STEP)
			{
				vec3 colour(0.0, 0.0, 0.0);
				float count = 0;
				for (float i = 0; i < AA; i += 1.0f)
				{
					for (float j = 0; j < AA; j += 1.0f)
					{
						vec4 dir = vec4((x - (SCREEN_WIDTH / 2)) + (i / AA), (y - (SCREEN_HEIGHT / 2)) + (j / AA), focal_length, 1);
						dir = R * dir;
						Intersection closest;
						if (GetIntersection(vec3(cameraPos), vec3(dir), closest))
						{
							// std::cout << "Sphere hit" << std::endl;
							colour += LightRay(closest, dir, BOUNCE, DINIT) + indirectLight; //(DirectLight(closest, BOUNCE, dir) + indirectLight);
							// if (!MTF90THETA)
							colour *= vec3(IntersectionColor(closest));
							count += 1.0f;
						}
						// if (TriangleIntersection(vec3(cameraPos), vec3(dir), closest))
						// {
						//   colour += triangles[closest.index].color;//LightRay(closest, dir, BOUNCE, DINIT) + indirectLight; //(DirectLight(closest, BOUNCE, dir) + indirectLight);
						//   // if (!MTF90THETA)
						//   //   colour *= vec3(IntersectionColor(closest));
						//   count += 1.0f;
						// }
					}
				}
				if (colour != vec3(0, 0, 0))
					colour /= count;
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

bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1)
{
	float discr = b * b - 4 * a * c;
	if (discr < 0)
		return false;
	else if (discr == 0)
		x0 = x1 = -0.5 * b / a;
	else
	{
		float q = (b > 0) ? -0.5 * (b + sqrt(discr)) : -0.5 * (b - sqrt(discr));
		x0 = q / a;
		x1 = c / q;
	}
	if (x0 > x1)
		std::swap(x0, x1);

	return true;
}

bool TriangleIntersection(vec3 start, vec3 dir, Intersection &in)
{
	auto lv3d = length(dir);
	auto nv3d = -dir;
	for (int i = 0; i < triangles.size(); ++i)
	{
		auto triangle = triangles[i];
		vec3 b = start - vec3(triangle.v0);
		mat3 A(nv3d, triangle.e1, triangle.e2);
		vec3 x = glm::inverse(A) * b;
		float dist = x.x; // / lv3d;

		if (x.x >= 0.0f &&
			x.y >= 0.0f &&
			x.z >= 0.0f &&
			x.y + x.z <= 1 &&
			in.distance > dist)
		{
			in = Intersection(vec3(start + x.x * dir), dist, i);
			in.sphere = false;
			// std::cout << "triangle t" << x.x << std::endl;
		}
	}

	return (in.index == -1) ? false : true;
}

bool SphereIntersection(vec3 start, vec3 dir, Intersection &in)
{
	for (int i = 0; i < spheres.size(); ++i)
	{
		auto sphere = spheres[i];
		float t0, t1;
		// vec3 ndir = glm::normalize(dir);

#if 0
		// geometric solution
		vec3 L = vec3(sphere.center) - start;
		float tca = glm::dot(L,dir);
		// if (tca < 0) return false;
		float d2 = glm::dot(L,L) - tca * tca;
		if (d2 > sphere.radius2) return false;
		float thc = sqrt(sphere.radius2 - d2);
		t0 = tca - thc;
		t1 = tca + thc;
#else
		// analytic solution
		vec3 L = start - vec3(sphere.center);
		float a = 1;
		float b = 2 * glm::dot(dir, L);
		float c = glm::dot(L, L) - sphere.radius2;
		if (!solveQuadratic(a, b, c, t0, t1))
			continue;
#endif
		if (t0 > t1)
			std::swap(t0, t1);
		if (t0 < 0)
			t0 = t1;

		// if (ndir.x * ndir.x + ndir.y * ndir.y >= 0.1f) continue;

		float dist = t0; // / length(dir);
		if (t0 >= 0 && in.distance > dist)
		{
			// std::cout << "sphere t" << t0 << std::endl;
			in.distance = dist;
			in.position = vec4(start + t0 * dir, 1);
			in.index = i;
			in.sphere = true;
		}
		// vec3 nhit = glm::normalize(phit - vec3(sphere.center));
	}

	return (in.index == -1) ? false : true;
}

bool GetIntersection(vec3 start, vec3 dir, Intersection &in)
{
	bool t = TriangleIntersection(start, glm::normalize(dir), in);
	bool s = SphereIntersection(start, glm::normalize(dir), in);
	return t || s;
}

vec3 point_light(vec3 color, vec3 r_v, vec3 n_u, float r)
{
	return ((color / (12.56f * r * r)) * sqrt(glm::dot(r_v, n_u) * glm::dot(r_v, n_u)));
}

vec4 IntersectionColor(const Intersection &i)
{

	if (i.sphere)
	{
		vec4 N = IntersectionNormal(i);
		float weight = 0.5f / light_1.size();
		vec4 color_to_return = 0.5f * weight * spheres[i.index].color;
		// vec4 lightdir = vec4(glm::normalize(vec3(light_1[0] - i.position)), 1);
		// float x = glm::angle(N, lightdir);
		for (int l = 0; l < light_1.size(); ++l)
		{
			vec4 lightdir = vec4(glm::normalize(vec3(light_1[l] - i.position)), 1);
			float x = glm::angle(N, lightdir);
			color_to_return += weight * (1.57f - x) * spheres[i.index].color;
		}
		// if (x > 0) return spheres[i.index].color * 0.8f;
		return color_to_return; //0.5f * (1.0f - x) * spheres[i.index].color + 0.5f *  spheres[i.index].color;
	}
	else
		return vec4(triangles[i.index].color, 0);
}

vec4 IntersectionNormal(const Intersection &i)
{
	if (i.sphere)
	{

		return vec4(glm::normalize(vec3(i.position) - vec3(spheres[i.index].center)), 1);
	}
	else
		return vec4(glm::normalize(vec3(triangles[i.index].normal)), 1);
}

bool TriangleShadowIntersection(vec4 start, vec4 dir, vec4 im, vec4 light_pos)
{
	Intersection closest;
	auto lsim = length(light_pos - im);
	auto lv3d = length(vec3(dir));
	auto nv3d = -vec3(dir);
	closest.distance = std::numeric_limits<float>::max();
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

			closest.sphere = false;
			if ((length((start + x.x * dir) - im) < lsim))
				return true;
		}
	}
	return false;
}

bool SphereShadowIntersection(vec4 start, vec4 dir, vec4 im, vec4 light_pos, bool isSphere)
{
	Intersection closest;
	auto lsim = length(light_pos - im);
	// auto lv3d = length(vec3(dir));
	// auto nv3d = -vec3(dir);
	closest.distance = std::numeric_limits<float>::max();

	for (int i = 0; i < spheres.size(); ++i)
	{
		auto sphere = spheres[i];
		float t0, t1;
		vec3 L = vec3(start) - vec3(sphere.center);
		float a = 1;
		float b = 2 * glm::dot(vec3(dir), L);
		float c = glm::dot(L, L) - sphere.radius2;
		if (!solveQuadratic(a, b, c, t0, t1))
			continue;

		if (t0 > t1)
			std::swap(t0, t1);
		if (t0 < 0)
			t0 = t1;

		// if (ndir.x * ndir.x + ndir.y * ndir.y >= 0.1f) continue;

		float dist = t0; // / length(dir);
		if (t0 >= 0 && closest.distance > dist)
		{
			// std::cout << "sphere t" << t0 << std::endl;
			closest.distance = dist;
			closest.sphere = true;
			isSphere = true;
			// std::cout << length((start + t0 * dir) - im) << std::endl;
			// std::cout << lsim << std::endl;
			if (length((start + t0 * dir) - im) < lsim) {
				// std::cout << "TRUE WORKS" << std::endl;
				return true;
				}
		}
		// vec3 nhit = glm::normalize(phit - vec3(sphere.center));
	}

	return false;
}

// ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
// ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
// ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
// THE RECURSIVE LIGHT RAY FUNCTION
// ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
// ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
// ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//i is the intersection, I is ray direction, b is bounce count, d is death chance
vec3 LightRay(Intersection &i, vec4 I, int b, float d)
{
	//1. DECIDE IF THIS RAY DIED
	std::random_device r;
	std::default_random_engine e(r());
	std::uniform_real_distribution<float> uniform_dist(0, 1);
	vec3 color = black;
	bool died = d < uniform_dist(e);
	if (died)
	{
		// std::cout << "DIED" << std::endl;
		return color;
	}
	//2. IF NOT DEAD, CALCULATE HOW MUCH LIGHT THE INTERSECTION GETS
	vec3 N = vec3(IntersectionNormal(i));
	vec4 start = i.position + 0.001f * vec4(N, 0);
	float weight = 1.0f / (light_1.size());

	for (int l = 0; l < light_1.size(); ++l)
	{
		vec3 rHat = vec3(light_1[l] - start);
		// if (GetIntersection(vec3(start), rHat, buf))
//!i.sphere && TriangleShadowIntersection(start, vec4(rHat, 1), i.position, light_1[l]) ||
bool sphere = false;
		if ( SphereShadowIntersection(start, vec4(rHat, 1), i.position, light_1[l], sphere))
		{
			// if (length(buf.position - i.position) < lsim)
			if(sphere)
			continue;
			// std::cout << "i dist is " << i.distance << std::endl;
			// std::cout << "cam dist is " << length( - start) << std::endl;
			// if (i.distance < (length(vec3(light_1[l]) - vec3(start))) / length(rHat))
		}
		float r = sqrt(pow(rHat.x, 2.0f) + pow(rHat.y, 2.0f) + pow(rHat.z, 2.0f));
		rHat = glm::normalize(rHat);
		color += weight * point_light(lightColor, rHat, N, r);
	}

	return color;
	if (b == 0)
		return (color + color);
	//3. MAIN BOUNCE
	float new_death = d - DEATH;
	vec4 reflecting = vec4(vec3(I) - 2 * (glm::dot(N, vec3(I))) * N, 1);
	Intersection reflection;
	if (GetIntersection(vec3(start), vec3(reflecting), reflection))
	{
		vec3 reflectionColor = vec3(IntersectionColor(reflection));
		vec3 reflected = weight * reflectionColor * LightRay(reflection, reflecting, b - 1, new_death);
		if (!BRDF)
			reflected *= (light_1.size() - 1) * weight;
		color += reflected;
		if (MTF90THETA)
			color *= reflectionColor;
	}
	//4. BRDF
	if (BRDF)
	{
		for (int l = 0; l < LIGHT_COUNT; ++l)
		{
			vec4 deviation(BRDF_RANGE * uniform_dist(e), BRDF_RANGE * uniform_dist(e), BRDF_RANGE * uniform_dist(e), 0);
			deviation += reflecting;
			if (GetIntersection(vec3(start), vec3(deviation), reflection))
				color += LIGHT_COUNT * weight * vec3(IntersectionColor(reflection)) * LightRay(reflection, deviation, b - 1, new_death);
		}
	}

	return color;
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
				for (unsigned int i = 0; i < light_1.size(); i++)
					light_1[i].x -= 0.1f;
				changed = true;
				break;
			case SDLK_KP_6:
				for (unsigned int i = 0; i < light_1.size(); i++)
					light_1[i].x += 0.1f;
				changed = true;
				break;
			case SDLK_KP_8:
				for (unsigned int i = 0; i < light_1.size(); i++)
					light_1[i].y -= 0.1f;
				changed = true;
				break;
			case SDLK_KP_2:
				for (unsigned int i = 0; i < light_1.size(); i++)
					light_1[i].y += 0.1f;
				changed = true;
				break;
			case SDLK_KP_5:
				for (unsigned int i = 0; i < light_1.size(); i++)
					light_1[i].z -= 0.1f;
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

void init_light_1()
{
	//vec4 lightPos1(0, -0.5, -0.7, 1.0);
	light_1.clear();
	// light_1.push_back(vec4(0, -0.5, -0.7, 1.0));
	// light_1.push_back(vec4(0.01, -0.5, -0.7, 1.0));
	//-0.5 + 0.5
	// float w = 0.2f;
	// float it = w / LIGHT_COUNT;
	for (int i = 0; i < LIGHT_COUNT; i++)
		for (int j = 0; j < LIGHT_COUNT; j++)
			for (int z = 0; z < LIGHT_COUNT; z++)
				light_1.push_back(vec4(0.5 + i * 0.05, -0.5 + z * 0.05, -0.70 + j * 0.05, 1.0));
}

// bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle> &triangles, Intersection &closest);
// bool ShadowIntersection(vec4 start, vec4 dir, const vector<Triangle> &triangles, vec4 im, vec4 light_pos);
// bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle> &triangles, Intersection &closest)
// {
// 	auto lv3d = length(vec3(dir));
// 	auto nv3d = -vec3(dir);
// 	for (int i = 0; i < triangles.size(); ++i)
// 	{
// 		auto triangle = triangles[i];
// 		vec3 b = vec3(start - triangle.v0);
// 		mat3 A(nv3d, triangle.e1, triangle.e2);
// 		vec3 x = glm::inverse(A) * b;
// 		float dist = x.x / lv3d;
// 		if (x.x >= 0.0f &&
// 			x.y >= 0.0f &&
// 			x.z >= 0.0f &&
// 			x.y + x.z <= 1 &&
// 			closest.distance > dist)
// 			closest = Intersection(vec3(start + x.x * dir), dist, i);
// 		closest.sphere = false;
// 	}

// 	for (int i = 0; i < spheres.size(); ++i)
// 	{
// 		auto sphere = spheres[i];
// 		float t0, t1;
// 		vec3 L = vec3(start - sphere.center);
// 		float a = glm::dot(vec3(dir), vec3(dir));
// 		float b = glm::dot(vec3(dir), L);
// 		float c = glm::dot(L, L) - sphere.radius2;
// 		if (!solveQuadratic(a, b, c, t0, t1))
// 			continue;
// 		vec3 phit = vec3(start) + glm::normalize(vec3(dir)) * t0;
// 		float dist = length(phit - vec3(start));
// 		if (closest.distance > dist)
// 		{
// 			std::cout << "Sphere hit" << std::endl;
// 			closest.position = vec4(phit, 1);
// 			closest.distance = dist;
// 			closest.index = i;
// 			closest.sphere = true;
// 		}
// 		// vec3 nhit = glm::normalize(phit - vec3(sphere.center));
// 	}

// 	return (closest.index == -1) ? false : true;
// }

// bool ShadowIntersection(vec4 start, vec4 dir, const vector<Triangle> &triangles, vec4 im, vec4 light_pos)
// {
// 	Intersection closest;
// 	auto lsim = length(light_pos - im);
// 	auto lv3d = length(vec3(dir));
// 	auto nv3d = -vec3(dir);
// 	closest.distance = 100.f;
// 	for (int i = 0; i < triangles.size(); ++i)
// 	{
// 		auto triangle = triangles[i];
// 		vec3 b = vec3(start - triangle.v0);
// 		mat3 A(nv3d, triangle.e1, triangle.e2);
// 		vec3 x = glm::inverse(A) * b;
// 		float dist = x.x / lv3d;
// 		if (x.x >= 0.0f &&
// 			x.y >= 0.0f &&
// 			x.z >= 0.0f &&
// 			x.y + x.z <= 1 &&
// 			closest.distance > dist)
// 		{
// 			closest.distance = dist;
// 			closest.sphere = false;
// 			if ((length((start + x.x * dir) - im) < lsim))
// 				return true;
// 		}
// 	}
// 	for (int i = 0; i < spheres.size(); ++i)
// 	{
// 		auto sphere = spheres[i];
// 		float t0, t1;
// 		vec3 L = vec3(start - sphere.center);
// 		float tca = glm::dot(L, vec3(dir));
// 		float d2 = glm::dot(L, L) - tca * tca;
// 		if (d2 > sphere.radius2)
// 			continue;
// 		float thc = sqrt(sphere.radius2 - d2);
// 		t0 = tca - thc;
// 		t1 = tca + thc;
// 		if (t0 > t1)
// 			std::swap(t0, t1);
// 		if (t0 < 0)
// 		{
// 			t0 = t1;
// 			if (t0 < 0)
// 				continue;
// 		}
// 		vec3 phit = vec3(start) + vec3(dir) * t0;
// 		float dist = length(phit - vec3(start));
// 		if (closest.distance > dist)
// 		{
// 			closest.distance = dist;
// 			closest.sphere = true;
// 			if (length(phit - vec3(im)) < lsim)
// 				return true;
// 		}
// 	}
// 	return false;
// }

// vec3 DirectLight(const Intersection &i, int bounce, vec4 I)
// {
//   vec3 direct_light = black;
//   vec3 n_u = glm::normalize(vec3(triangles[i.index].normal));
//   for (int l = 0; l < light_1.size(); ++l)
//   {
//     vec3 r_v = vec3(light_1[l] - i.position);
//     float count = 1.0f / light_1.size();
//     if (ShadowIntersection(i.position + 0.001f * triangles[i.index].normal, vec4(r_v, 1), triangles, i.position, light_1[l]))
//       count = 0.0f;
//     float r = sqrt(pow(r_v.x, 2.0f) + pow(r_v.y, 2.0f) + pow(r_v.z, 2.0f));
//     r_v = glm::normalize(r_v);
//     direct_light += point_light(lightColor, r_v, n_u, r, count);
//   }
//   //BOUNCE
//   if (bounce > 0)
//   {
//     float Bp = pow(float(BOUNCE), float(BOUNCE_POWER));
//     float bp = pow(bounce, float(BOUNCE_POWER));

//     float weight = (bp / Bp);
//     if (BRDF)
//       weight /= (BOUNCE * BOUNCE);

//     vec4 ref_dir = vec4(vec3(I) - 2 * (glm::dot(n_u, vec3(I))) * n_u, 1);
//     Intersection closest;
//     if (ClosestIntersection(i.position + 0.001f * triangles[i.index].normal, ref_dir, triangles, closest))
//     {
//       vec3 reflected_light = weight * (DirectLight(closest, bounce - 1, ref_dir)) * triangles[closest.index].color;
//       direct_light += reflected_light;
//       if (MTF90THETA)
//         direct_light *= triangles[closest.index].color;
//     }

//     if (BRDF)
//     {
//       for (int j = 0; j < BOUNCE * BOUNCE; ++j)
//       {
//         std::random_device r;
//         std::default_random_engine e(r());
//         std::uniform_real_distribution<float> uniform_dist(-0.1f, 0.1f);
//         float rand_x = uniform_dist(e);
//         float rand_y = uniform_dist(e);
//         float rand_z = uniform_dist(e);
//         vec4 rand_dev(rand_x, rand_y, rand_z, 0);
//         rand_dev += ref_dir;
//         Intersection closest_dev;
//         if (ClosestIntersection(i.position + 0.001f * triangles[i.index].normal, rand_dev, triangles, closest_dev))
//         {
//           vec3 reflected_light = weight * (DirectLight(closest_dev, bounce - 1, rand_dev)) * triangles[closest_dev.index].color;
//           direct_light += reflected_light;
//         }
//       }
//     }
//   }

//   return direct_light;
// }