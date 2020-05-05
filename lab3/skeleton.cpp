#define _USE_MATH_DEFINES
#include <iostream>
#include <glm/glm.hpp>
#include "SDL.h"
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <math.h>

using namespace std;
using glm::vec3;
using glm::ivec2;
using glm::mat3;
using glm::vec2;

struct Pixel{
	int x;
	int y;
	float zinv;
	vec3 illumination;
};

struct Vertex{
	vec3 position;
	vec3 normal;
	vec3 reflectance;
};

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
vec3 currentcolor;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
vec3 lightPos(0, -0.5, -0.7); 
vec3 lightPower = 1.1f*vec3(1, 1, 1); 
vec3 indirectLightPowerPerArea = 0.5f*vec3(1, 1, 1);
vec3 reflectance = vec3(0.1, 0.1, 0.1);
// ----------------------------------------------------------------------------
// FUNCTIONS

void Update();
void Draw();
void VertexShader(const Vertex& v, Pixel& p);
void DrawPolygonEdges(const vector<vec3>& vertices);
//void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result);
void DrawLineSDL(SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color);
void ComputePolygonRows(const vector<Pixel>& vertexpixels, vector<Pixel>& leftPixels, vector<Pixel>& rightpixels);
void DrawRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels);
void DrawPolygon(const vector<Vertex>& vertices);
void Interpolate(Pixel a, Pixel b, vector<Pixel>& result);
void PixelShader(const Pixel& p);

int main( int argc, char* argv[] )
{
	LoadTestModel( triangles );
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.
	

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
	cout << "Render time: " << dt << " ms." << endl;

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

	for(int y=0;y<SCREEN_HEIGHT;++y){
		for(int x=0; x<SCREEN_WIDTH;++x){
			depthBuffer[y][x] = 0;
		}
	}
	
	for( int i=0; i<triangles.size(); ++i )
	{
		currentcolor = triangles[i].color;
		vector<Vertex> vertices(3);

		vertices[0].position = triangles[i].v0;
		vertices[1].position = triangles[i].v1;
		vertices[2].position = triangles[i].v2;
		vertices[0].normal = triangles[i].normal;
		vertices[1].normal = triangles[i].normal;
		vertices[2].normal = triangles[i].normal;

		// Add drawing
		DrawPolygon(vertices);
	}
	
	if ( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

void VertexShader(const Vertex& v, Pixel& p) {
	vec3 P = (v.position - cameraPos) * R;
	p.x = (f * P.x/P.z) + (SCREEN_WIDTH / 2);
	p.y = (f * P.y / P.z) + (SCREEN_HEIGHT /2);
	p.zinv = 1.0f / P.z;
	float r = glm::length(lightPos - v.position);
	vec3 r_hat = glm::normalize(lightPos - v.position);
	vec3 n_hat = glm::normalize(v.normal);
	vec3 D =lightPower * (glm::max(glm::dot(r_hat, n_hat), float(0))) / (4.0f * float(M_PI) * r * r);
	p.illumination = reflectance * (D + indirectLightPowerPerArea);

}

/*void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result) {
	int N = result.size();
	vec2 step = vec2(b - a) / float(glm::max(N - 1, 1));
	vec2 current(a); 
	for (int i = 0; i < N; ++i) {
		result[i] = current;
		current += step;
	}

}*/

void Interpolate(Pixel a, Pixel b, vector<Pixel>& result){
	int N = result.size();
	vec3 step;
	step.x = float(b.x - a.x) / float(glm::max(N - 1, 1));
	step.y = float(b.y - a.y) / float(glm::max(N - 1, 1));
	step.z = float(b.zinv - a.zinv) / float(glm::max(N - 1, 1));
	vec3 il = b.illumination - a.illumination / float(glm::max(N - 1, 1));
	Pixel current(a);
	for (int i = 0; i < N; ++i) {
		result[i] = current;
		current.x = a.x + (i+1)*step.x;
		current.y = a.y + (i+1)*step.y;
		current.zinv = a.zinv + (i+1)* step.z;

		current.illumination += il;
		
	}
}

void PixelShader(const Pixel& p) {
	int x = p.x;
	int y = p.y;
	if (p.zinv > depthBuffer[y][x]) {
		depthBuffer[y][x] = p.zinv;
		PutPixelSDL(screen, x, y, currentcolor * p.illumination);
	}
}

/*void DrawPolygonEdges(const vector<vec3>& vertices) {
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
}*/

/*void DrawLineSDL(SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color) {
	ivec2 delta = glm::abs(a - b);
	int pixels = glm::max(delta.x, delta.y) + 1;

	vector<ivec2> line(pixels);
	Interpolate(a, b, line);

	for (int i = 0; i < pixels; ++i) {
		ivec2 point = line[i];
		PutPixelSDL(surface, point.x, point.y, color);
	}
}*/

void ComputePolygonRows(const vector<Pixel>& vertexpixels, vector<Pixel>& leftPixels, vector<Pixel>& rightpixels) {
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
	//resize leftPixels, rightPixels
	leftPixels.resize(rows);
	rightpixels.resize(rows);

	//initialize x-coordinates in both arrays
	for (int i = 0; i < rows; ++i) {
		leftPixels[i].x = +numeric_limits<int>::max();
		leftPixels[i].y = min + i;
		rightpixels[i].x = -numeric_limits<int>::max();
		rightpixels[i].y = min + i;
	}

	//loop through, interpolate, update values
	for (int i = 0; i < vertexpixels.size(); ++i) {
		int j = (i + 1) % vertexpixels.size();
		ivec2 delta;
		delta.x = glm::abs(vertexpixels[i].x - vertexpixels[j].x);
		delta.y = glm::abs(vertexpixels[i].y - vertexpixels[j].y);
		int pixels = glm::max(delta.x,delta.y) + 1;
		vector<Pixel> line(pixels);
		Interpolate(vertexpixels[i], vertexpixels[j], line);
		for (int l = 0; l < line.size(); ++l) {
			for (int c = 0; c < rows; ++c) {
				if (leftPixels[c].y == line[l].y) {
					if (rightpixels[c].x <= line[l].x) {
						rightpixels[c].x = line[l].x;
						rightpixels[c].zinv = line[l].zinv;
						rightpixels[c].illumination = line[l].illumination;
					}
					if (leftPixels[c].x >= line[l].x) {
						leftPixels[c].x = line[l].x;
						leftPixels[c].zinv = line[l].zinv;
						leftPixels[c].illumination = line[l].illumination;
					}
				}
			}
		}
	}
}

void DrawRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels){
	int pixels = leftPixels.size();
	for (int i = 0; i < leftPixels.size(); ++i) {
		vector<Pixel> row(pixels);
		Interpolate(leftPixels[i], rightPixels[i], row);
		for (int j = 0; j < row.size(); ++j) {
			PixelShader(row[j]);
		}
	}
}

void DrawPolygon(const vector<Vertex>& vertices){
	int V = vertices.size();
	vector<Pixel> vertexPixels(V);
	for(int i=0;i<V;++i){
		VertexShader(vertices[i], vertexPixels[i]);
	}
	vector<Pixel> leftPixels;
	vector<Pixel> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	DrawRows(leftPixels, rightPixels);
}
