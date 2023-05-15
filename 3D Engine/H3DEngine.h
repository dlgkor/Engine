#pragma once
#include<vector>
#include<math.h>
#include<algorithm>
#include<fstream>
#include<strstream>

struct vec3d {
	float x, y, z;
};

struct triangle {
	vec3d p[3];
	float shadow;
};

struct mesh {
	std::vector<triangle> tris;

	bool loadFromObjectFil(std::string sFilename) {
		std::ifstream f(sFilename);
		if (!f.is_open())
			return false;

		std::vector<vec3d> verts;

		while (!f.eof()) {
			char line[128];
			f.getline(line, 128);

			std::strstream s;
			s << line;

			char junk;

			if (line[0] == 'v') {
				vec3d v;
				s >> junk >> v.x >> v.y >> v.z;
				verts.push_back(v);
			}

			if (line[0] == 'f') {
				int f[3];
				s >> junk >> f[0] >> f[1] >> f[2];
				tris.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });

			}
			
		}

		return true;
	}
};

struct mat4x4 {
	float m[4][4] = { 0 };
	//[열][행] = [x][y] = [세로묶음][가로묶음]
};

struct camera {
	vec3d position;
	float pitch;
	float yaw;
};

float dotProduct(vec3d& i, vec3d& j) {
	return i.x * j.x + i.y * j.y + i.z * j.z;
}

vec3d crossProduct(vec3d& i, vec3d& j) {
	vec3d tmp = { 
		i.y * j.z - i.z * j.y, 
		i.z * j.x - i.x * j.z, 
		i.x * j.y - i.y * j.x 
	};
	return tmp;
}

void normalize(vec3d& v) {
	float length = sqrtf((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
	v.x /= length; v.y /= length; v.z /= length;
	return;
}

vec3d operator*(const mat4x4& m, const vec3d& v)
{
	vec3d nVector;
	nVector.x = m.m[0][0] * v.x + m.m[1][0] * v.y + m.m[2][0] * v.z + m.m[3][0] * 1;
	nVector.y = m.m[0][1] * v.x + m.m[1][1] * v.y + m.m[2][1] * v.z + m.m[3][1] * 1;
	nVector.z = m.m[0][2] * v.x + m.m[1][2] * v.y + m.m[2][2] * v.z + m.m[3][2] * 1;

	return nVector;
}

vec3d operator-(vec3d& i, vec3d& j) {
	vec3d tmp = { i.x - j.x, i.y - j.y, i.z - j.z };
	return tmp;
}


mat4x4 FPSViewRH(vec3d eye, float pitch, float yaw)
{
	// I assume the values are already converted to radians.
	float cosPitch = cos(pitch);
	float sinPitch = sin(pitch);
	float cosYaw = cos(yaw);
	float sinYaw = sin(yaw);

	vec3d xaxis = { cosYaw, 0, -sinYaw };
	vec3d yaxis = { sinYaw * sinPitch, cosPitch, cosYaw * sinPitch };
	vec3d zaxis = { sinYaw * cosPitch, -sinPitch, cosPitch * cosYaw };

	// Create a 4x4 view matrix from the right, up, forward and eye position vectors
	mat4x4 viewMatrix;
	viewMatrix.m[0][0] = xaxis.x; viewMatrix.m[0][1] = yaxis.x; viewMatrix.m[0][2] = zaxis.x;
	viewMatrix.m[1][0] = xaxis.y; viewMatrix.m[1][1] = yaxis.y; viewMatrix.m[1][2] = zaxis.y;
	viewMatrix.m[2][0] = xaxis.z; viewMatrix.m[2][1] = yaxis.z; viewMatrix.m[2][2] = zaxis.z;
	viewMatrix.m[3][0] = -dotProduct(xaxis, eye); viewMatrix.m[3][1] = -dotProduct(yaxis, eye); viewMatrix.m[3][2] = -dotProduct(zaxis, eye);
	viewMatrix.m[3][3] = 1;

	return viewMatrix;
}

class Engine3D {
private:
	mesh meshCube;
	mat4x4 matProj;
	
	float ScreenHeight;
	float ScreenWidth;

	float fNear;
	float fFar;
	float fFov;
	float fAspectRatio;
	float fFovRad;


	float fTheta;

	HWND hWnd;
	camera Camera;

	void MultiplyMatrixVector(vec3d& i, vec3d& o, mat4x4& m) {
		//지금 vector에서 W의 값을 1이라 가정하고 계산중이다
		o.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + m.m[3][0];
		o.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + m.m[3][1];
		o.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + m.m[3][2];
		float w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + m.m[3][3];
		if (w != 0.0f) {
			o.x /= w; o.y /= w; o.z /= w;
		}
	}
public:
	Engine3D(HWND h):hWnd(h) {	
		Camera.position = { 0.0f, 0.0f, -9.5f};
		Camera.pitch = 0.0f;
		Camera.yaw = 0.0f;
	}
	bool OnUserCreate() {
		RECT rt;
		GetWindowRect(hWnd, &rt);
		ScreenHeight = abs(rt.top - rt.bottom);
		ScreenWidth = abs(rt.left - rt.right);

		/*
		meshCube.tris = {
		{ 0.0f, 0.0f, 0.0f,    0.0f, 1.0f, 0.0f,    1.0f, 1.0f, 0.0f },
		{ 0.0f, 0.0f, 0.0f,    1.0f, 1.0f, 0.0f,    1.0f, 0.0f, 0.0f },

		// EAST                                                      
		{ 1.0f, 0.0f, 0.0f,    1.0f, 1.0f, 0.0f,    1.0f, 1.0f, 1.0f },
		{ 1.0f, 0.0f, 0.0f,    1.0f, 1.0f, 1.0f,    1.0f, 0.0f, 1.0f },

		// NORTH                                                     
		{ 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f },
		{ 1.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 1.0f },

		// WEST                                                      
		{ 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 0.0f },
		{ 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 0.0f,    0.0f, 0.0f, 0.0f },

		// TOP                                                       
		{ 0.0f, 1.0f, 0.0f,    0.0f, 1.0f, 1.0f,    1.0f, 1.0f, 1.0f },
		{ 0.0f, 1.0f, 0.0f,    1.0f, 1.0f, 1.0f,    1.0f, 1.0f, 0.0f },

		// BOTTOM                                                    
		{ 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f },
		{ 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f,    1.0f, 0.0f, 0.0f },

		};
		*/
		meshCube.loadFromObjectFil("tongs.obj");

		//행렬연산
		fNear = 0.1f;
		fFar = 1000.0f;
		fFov = 90.0f;
		fAspectRatio = ScreenHeight / ScreenWidth; //height / width

		fFovRad = 1.0f / tanf(fFov * 0.5f / 180.0f * 3.14159f);

		//matProj.m[0][0] = fAspectRatio * fFovRad;
		//matProj.m[1][1] = fFovRad;
		matProj.m[0][0] = fAspectRatio * fNear; 
		matProj.m[1][1] = fNear;
		//확대 할꺼면 같은 값을 곱해줘야 함.
		//ex Width를 각각 곱하기
		//비율은 어차피 유지된다.
		//javidx9에 나온 행렬에서 시야각과 절두체 윗면까지의 거리 사이의 관계가 명시되어 있지 않으므로 잘못된 공식이다.
		//둘 사이에 연관성을 만들거나, 절두체 윗면까지의 거리로 통일하자. 나는 후술한 해결책을 사용했다.

		matProj.m[2][2] = fFar / (fFar - fNear);
		matProj.m[3][2] = (-fFar * fNear) / (fFar - fNear);
		matProj.m[2][3] = 1.0f;
		matProj.m[3][3] = 0.0f;

		return true;
	}
	void Render(HDC hdc, float fElapsedTime) {
		mat4x4 matRotZ, matRotX;
		fTheta += 1.0f * fElapsedTime;

		// Rotation Z
		matRotZ.m[0][0] = cosf(fTheta);
		matRotZ.m[0][1] = sinf(fTheta);
		matRotZ.m[1][0] = -sinf(fTheta);
		matRotZ.m[1][1] = cosf(fTheta);
		matRotZ.m[2][2] = 1;
		matRotZ.m[3][3] = 1;

		// Rotation X
		matRotX.m[0][0] = 1;
		matRotX.m[1][1] = cosf(fTheta * 0.5f);
		matRotX.m[1][2] = sinf(fTheta * 0.5f);
		matRotX.m[2][1] = -sinf(fTheta * 0.5f);
		matRotX.m[2][2] = cosf(fTheta * 0.5f);
		matRotX.m[3][3] = 1;

		
		//Camera.position.x += fElapsedTime * 1;
		//Camera.position.y += fElapsedTime * 1;
		//Camera.yaw -= fElapsedTime * 0.1f;

		//We need Clipping
		

		mat4x4 viewmatrix = FPSViewRH(Camera.position, Camera.pitch, Camera.yaw);
		std::vector<triangle> tristoRoaster;

		//Draw triangle
		for (auto tri : meshCube.tris) {
			triangle triProjected, triTranslated, triRotatedZ, triRotatedZX;

			// Rotate in Z-Axis
			MultiplyMatrixVector(tri.p[0], triRotatedZ.p[0], matRotZ);
			MultiplyMatrixVector(tri.p[1], triRotatedZ.p[1], matRotZ);
			MultiplyMatrixVector(tri.p[2], triRotatedZ.p[2], matRotZ);

			// Rotate in X-Axis
			MultiplyMatrixVector(triRotatedZ.p[0], triRotatedZX.p[0], matRotX);
			MultiplyMatrixVector(triRotatedZ.p[1], triRotatedZX.p[1], matRotX);
			MultiplyMatrixVector(triRotatedZ.p[2], triRotatedZX.p[2], matRotX);
			
			for (int i = 0; i < 3; i++) {
				triTranslated.p[i] = viewmatrix * triRotatedZX.p[i];
				//triTranslated.p[i] = viewmatrix * tri.p[i];
			}

			int maxWidth, maxHeight;
			bool isOnBoundary = true;

			for (int i = 0; i < 3; i++) {
				maxWidth = (triTranslated.p[i].z / fNear) * ScreenWidth / 2.0f;
				maxHeight = maxWidth * fAspectRatio / 2.0f;

				if(triTranslated.p[i].z < fFar && triTranslated.p[i].z > fNear){}
				else {
					isOnBoundary = false;
				}
				
				if (abs(triTranslated.p[i].x) < maxWidth && abs(triTranslated.p[i].y) < maxHeight){}
				else {
					isOnBoundary = false;
				}
			}

			if (isOnBoundary == false) {
				continue;
			}

			vec3d vecA = triTranslated.p[1] - triTranslated.p[0];
			vec3d vecB = triTranslated.p[2] - triTranslated.p[0];
			vec3d crossVec3d = crossProduct(vecA, vecB);
			normalize(crossVec3d);

			vec3d cameraLine = triTranslated.p[0] - Camera.position;

			if (dotProduct(cameraLine, crossVec3d) > 0.0f) {
				continue;
			}

			// Project triangles from 3D --> 2D
			MultiplyMatrixVector(triTranslated.p[0], triProjected.p[0], matProj);
			MultiplyMatrixVector(triTranslated.p[1], triProjected.p[1], matProj);
			MultiplyMatrixVector(triTranslated.p[2], triProjected.p[2], matProj);

			vec3d light_direction = { 0.0f, 0.0f, -1.0f }; //빛이 튕겨져 나올 때의 벡터를 정의한다.
			normalize(light_direction);
			float dp = max(dotProduct(crossVec3d, light_direction), 0.0f);

			triProjected.shadow = dp;

			tristoRoaster.push_back(triProjected);
		}

		std::sort(tristoRoaster.begin(), tristoRoaster.end(), [](triangle& t1, triangle& t2) {
			float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
			float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
			return z1 > z2;
		});

		for (auto tri : tristoRoaster) {
			HBRUSH oldBrush, DLGBrush;
			HPEN oldPen, DLGPen;

			int rN = 255*tri.shadow;
			int gN = 153*tri.shadow;
			int bN = 51*tri.shadow;
			/*
			74 168 216
			128 0 128
			*/

			DLGBrush = CreateSolidBrush(RGB(rN, gN, bN));
			oldBrush = (HBRUSH)SelectObject(hdc, DLGBrush);

			DLGPen = CreatePen(PS_SOLID, 1, RGB(rN, gN, bN));
			oldPen = (HPEN)SelectObject(hdc, DLGPen);

			POINT pList[3];
			float multiply;
			multiply = ScreenWidth;
			for (int i = 0; i < 3; i++) {
				pList[i].x = tri.p[i].x * multiply;
				pList[i].y = tri.p[i].y * multiply;
			}
			Polygon(hdc, pList, 3);

			DeleteObject(SelectObject(hdc, oldBrush));
			DeleteObject(SelectObject(hdc, oldPen));
		}
	}
	~Engine3D() {}
};