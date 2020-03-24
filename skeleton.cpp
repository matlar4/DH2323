// Introduction lab that covers:
// * C++
// * SDL
// * 2D graphics
// * Plotting pixels
// * Video memory
// * Color representation
// * Linear interpolation
// * glm::vec3 and std::vector

#include "SDL.h"
#include <iostream>
#include <glm/glm.hpp>
#include <vector>
#include "SDLauxiliary.h"

using namespace std;
using glm::vec3;

// --------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 640;
const int SCREEN_HEIGHT = 480;
SDL_Surface* screen;

// --------------------------------------------------------
// FUNCTION DECLARATIONS

void Draw();
void interpolate(vec3 a, vec3 b, vector<vec3>& result);
void Update();

// --------------------------------------------------------
// FUNCTION DEFINITIONS


vec3 topLeft(1, 0, 0);
vec3 topRight(0, 0, 1);
vec3 bottomRight(0, 1, 0); //example code was wrong here
vec3 bottomLeft(1, 1, 0);
vector<vec3> leftSide(SCREEN_HEIGHT);
vector<vec3> rightSide(SCREEN_HEIGHT);

vector<vec3> stars(1000);
float v = 0.0005;


int t;


int main( int argc, char* argv[] )
{	

	interpolate(topLeft, bottomLeft, leftSide);
	interpolate(topRight, bottomRight, rightSide);

	
	for (int i = 0; i < stars.size(); i++) {
		float rZ = float(rand()) / float(RAND_MAX);
		float rX = -1 + float(rand()) / float(RAND_MAX / 2);
		float rY = -1 + float(rand()) / float(RAND_MAX / 2);
		stars[i] = vec3(rX, rY, rZ);
	}
	

	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();
	while( NoQuitMessageSDL() )
	{
		Update();
		Draw();
	}
	SDL_SaveBMP( screen, "screenshot.bmp" );
	return 0;
}

void Draw()
{	
	float f = SCREEN_HEIGHT / 2;
	float w2 = SCREEN_WIDTH / 2;

	//uncomment for interpolation
	/*for( int y=0; y<SCREEN_HEIGHT; ++y )
	{
		vector<vec3> row(SCREEN_WIDTH);
		interpolate(leftSide[y], rightSide[y], row);
		for( int x=0; x<SCREEN_WIDTH; ++x )
		{
			PutPixelSDL( screen, x, y, row[x] );
		}
	}*/

	SDL_FillRect(screen, 0, 0);
	if (SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);
	for (size_t s = 0; s < stars.size(); ++s) {
		vec3 star = stars[s];
		//project screen
		float u = (f * (star.x / star.z)) + (SCREEN_WIDTH / 2);
		float v = (f * (star.y / star.z)) + (SCREEN_HEIGHT / 2);
		
		vec3 color = 0.2f * vec3(1, 1, 1) / (star.z * star.z);

		PutPixelSDL(screen, u, v, color);
	}

	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

void interpolate(vec3 a, vec3 b, vector<vec3>& result) {
	if (result.size() == 1) {
		result[0] = a;
		return;
	}
	float deltaX = (b.x - a.x) / (result.size() -1);
	float deltaY = (b.y - a.y) / (result.size() -1);
	float deltaZ = (b.z - a.z ) / (result.size() -1);
	for (float i = 0; i < result.size(); i++) {
		float x = a.x + i * deltaX;
		float y = a.y + i * deltaY;
		float z = a.z + i * deltaZ;
		result[i] = vec3(x,y,z);
	}
}

void Update() {
	int t2 = SDL_GetTicks();
	float dt = float(t2 - t);
	t = t2;

	for (size_t s = 0; s < stars.size(); ++s) {
		//update stars
		stars[s].z = stars[s].z - v * dt;

		if (stars[s].z < 0) {
			stars[s].z += 1;
		}
		if (stars[s].z > 1) {
			stars[s].z -= 1;
		}
	}
}