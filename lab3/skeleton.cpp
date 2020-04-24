#include <iostream>
#include <glm/glm.hpp>
#include "SDL.h"
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::ivec2;
using glm::mat3;
using glm::vec2;

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;
vector<Triangle> triangles;
float f =SCREEN_HEIGHT; //assuming 90 degree viewing angle
vec3 cameraPos(0, 0, -3.001);
mat3 R = mat3(vec3(1, 0, 0), vec3(0, 1, 0), vec3(0, 0, 1));
float yaw = 0.0f;

// ----------------------------------------------------------------------------
// FUNCTIONS

void Update();
void Draw();
void VertexShader(const vec3& v, ivec2& p);
void DrawPolygonEdges(const vector<vec3>& vertices);
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result);
void DrawLineSDL(SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color);
void ComputePolygonRows(const vector<ivec2>& vertexpixels, vector<ivec2>& leftPixels, vector<ivec2>& rightpixels);

int main( int argc, char* argv[] )
{
	LoadTestModel( triangles );
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.
	vector<ivec2> vertexPixels(3);
	vertexPixels[0] = ivec2(10, 5);
	vertexPixels[1] = ivec2(5, 10); 
	vertexPixels[2] = ivec2(15, 15);
	vector<ivec2> leftPixels;
	vector<ivec2> rightPixels; 
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	for (int row = 0; row < leftPixels.size(); ++row){
		cout << "Start: ("
			<< leftPixels[row].x << ","
			<< leftPixels[row].y << ")."
			<< "End: ("
			<< rightPixels[row].x << ","
			<< rightPixels[row].y << ")." << endl;
		}

	while( NoQuitMessageSDL() )
	{
		Update();
		Draw();
	}

	SDL_SaveBMP( screen, "screenshot.bmp" );
	return 0;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	//cout << "Render time: " << dt << " ms." << endl;

	vec3 forward(R[2][0], R[2][1], R[2][2]);
	vec3 right(R[0][0], R[0][1], R[0][2]);

	Uint8* keystate = SDL_GetKeyState(0);

	if( keystate[SDLK_UP] )
		cameraPos.z += 0.01f ;

	if( keystate[SDLK_DOWN] )
		cameraPos.z -= 0.01f ;

	if( keystate[SDLK_RIGHT] )
		yaw = yaw - 0.01;
		R = mat3(vec3(glm::cos(yaw), float(0), -glm::sin(yaw)), vec3(0, 1, 0), vec3(glm::sin(yaw), float(0), glm::cos(yaw)));
		cameraPos = cameraPos + 0.01f * right;
		

	if( keystate[SDLK_LEFT] )
		yaw = yaw + 0.01;
		R = mat3(vec3(glm::cos(yaw), float(0), -glm::sin(yaw)), vec3(0, 1, 0), vec3(glm::sin(yaw), float(0), glm::cos(yaw)));
		cameraPos = cameraPos - 0.01f * right;

	if( keystate[SDLK_RSHIFT] )
		;

	if( keystate[SDLK_RCTRL] )
		;

	if( keystate[SDLK_w] )
		;

	if( keystate[SDLK_s] )
		;

	if( keystate[SDLK_d] )
		;

	if( keystate[SDLK_a] )
		;

	if( keystate[SDLK_e] )
		;

	if( keystate[SDLK_q] )
		;
}

void Draw()
{
	SDL_FillRect( screen, 0, 0 );

	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);
	
	for( int i=0; i<triangles.size(); ++i )
	{
		vector<vec3> vertices(3);

		vertices[0] = triangles[i].v0;
		vertices[1] = triangles[i].v1;
		vertices[2] = triangles[i].v2;

		// Add drawing
		DrawPolygonEdges(vertices);
	}
	
	if ( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

void VertexShader(const vec3& v, ivec2& p) {
	vec3 P = (v - cameraPos) * R;
	p.x = (f * P.x/P.z) + (SCREEN_WIDTH / 2);
	p.y = (f * P.y / P.z) + (SCREEN_HEIGHT /2);

}

void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result) {
	int N = result.size();
	vec2 step = vec2(b - a) / float(glm::max(N - 1, 1));
	vec2 current(a); 
	for (int i = 0; i < N; ++i) {
		result[i] = current;
		current += step;
	}

}

void DrawPolygonEdges(const vector<vec3>& vertices) {
	int delta = vertices.size();
	vector<ivec2> projectedVertices(delta);
	for (int i = 0; i < delta; ++i) {
		VertexShader(vertices[i], projectedVertices[i]);
	}
	//draw edges between vetrices 
	for (int i = 0; i < delta; ++i) {
		int j = (i + 1) % delta; // The next vertex
		vec3 color(1, 1, 1);
		DrawLineSDL(screen, projectedVertices[i], projectedVertices[j], color);
	}
}

void DrawLineSDL(SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color) {
	ivec2 delta = glm::abs(a - b);
	int pixels = glm::max(delta.x, delta.y) + 1;

	vector<ivec2> line(pixels);
	Interpolate(a, b, line);

	for (int i = 0; i < pixels; ++i) {
		ivec2 point = line[i];
		PutPixelSDL(surface, point.x, point.y, color);
	}
}

void ComputePolygonRows(const vector<ivec2>& vertexpixels, vector<ivec2>& leftPixels, vector<ivec2>& rightpixels) {
	// find max and min y-value of polygon
	//compute number of rows
	int min = +numeric_limits<int>::max();
	int max = -numeric_limits<int>::max();
	for (int i = 0; i < vertexpixels.size(); ++i) {
		if (vertexpixels[i].y > max){
			max = vertexpixels[i].y;
		}
		if (vertexpixels[i].y < min) {
			min = vertexpixels[i].y;
		}
	}
	int rows = (max - min) + 1;
	cout << rows << endl;
	//resize leftPixels, rightPixels
	leftPixels.resize(rows);
	rightpixels.resize(rows);

	//initialize x-coordinates in both arrays
	for (int i = 0; i < rows; ++i) {
		leftPixels[i].x = +numeric_limits<int>::max();
		rightpixels[i].x = -numeric_limits<int>::max();
	}

	//loop through, interpolate, update values
	for (int i = 0; i < vertexpixels.size(); ++i) {
		int j = (i + 1) % vertexpixels.size();
		ivec2 delta = glm::abs(vertexpixels[i] - vertexpixels[j]);
		int pixels = delta.x + 1;
		vector<ivec2> result(pixels);
		Interpolate(vertexpixels[i], vertexpixels[j], result);
		for (int row = 0; row < pixels; ++row) {
			cout << result[row].x << "  " << result[row].y << endl;
			if (rightpixels[i].x < result[row].x) {
				rightpixels[i] = result[row];
			}
			if (leftPixels[i].x > result[row].x) {
				leftPixels[i] = result[row];
			}
		}
	}
}