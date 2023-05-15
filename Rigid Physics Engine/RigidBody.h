#pragma once
#include<math.h>
#include<iostream>
#include<SDL2/SDL.h>

class matrix2x2{
public:
    float m[2][2];
};

class vector2d {
public:
	float x, y;
public:
	vector2d() :x(0.0f), y(0.0f) {}
	vector2d(float _x, float _y) :x(_x), y(_y) {}

	float Size() {
		float size = sqrt((x * x) + (y * y));
		return size;
	}

	vector2d UnitVector() {
		return vector2d(x / Size(), y / Size());
	}

	vector2d operator+(const vector2d& a) const{
		vector2d vec(x + a.x, y + a.y);
		return vec;
	}
	vector2d operator-(const vector2d& a) const{
		vector2d vec(x - a.x, y - a.y);
		return vec;
	}
	vector2d operator*(const float& a) const{
		vector2d vec(x * a, y * a);
		return vec;
	}

	vector2d operator*(const matrix2x2& matrix) const{
		vector2d vec(x * matrix.m[0][0] + y * matrix.m[1][0], x * matrix.m[0][1] + y * matrix.m[1][1]);
		//matrix m[0][] i hat
		//matrix m[1][] j hat
		return vec;
	}

    float dot(const vector2d& a){
        return x*a.x + y*a.y;
    }
	float cross(const vector2d& a){
        return x*a.y - y*a.x;
    }        
};

class Vertex{
public:
	vector2d point;
	Vertex* nextVertex;
};

class RigidBody{
public:
    Vertex* vertex; //allcate and deallocate
    int vertexCount;
    vector2d* sidenormal;
    //relative vecter that sets the original point by center of mass
    //rotation applied

    vector2d centerOfMass; 
    //translation applied
    //related to collision

    vector2d velocity;
	double mass;
	double inverseMass;

    double angular_velocity; 
	double inertia_tensor;
	double inverse_inertensor;

	float maxLength; //used to set circle boundary to check collision arbitarily

	bool pin;


	vector2d accumulatedImpulse;
	vector2d accumulatedRotateImpulse;

public:
	RigidBody(vector2d* _vertex, int _vertexCount, 
	vector2d _centerOfMass, double _mass, double _inertia_tensor) : 
	vertexCount(_vertexCount), centerOfMass(_centerOfMass), mass(_mass), inertia_tensor(_inertia_tensor) {
		vertex = new Vertex;
		vertex->point = _vertex[0];
		Vertex* CurrentVertex = vertex;
		maxLength = CurrentVertex->point.Size();
		for(int i=1;i<vertexCount;i++){
			CurrentVertex->nextVertex = new Vertex;
			CurrentVertex = CurrentVertex->nextVertex;
			CurrentVertex->point = _vertex[i];
			
			maxLength = std::max(maxLength,CurrentVertex->point.Size());
			
		}
		CurrentVertex->nextVertex = vertex;

		
		inertia_tensor = 0.0;
		CurrentVertex = vertex;
		for(int i=0;i<vertexCount;i++){
			inertia_tensor += pow(CurrentVertex->point.Size(),2) * mass / (vertexCount);
			CurrentVertex = CurrentVertex->nextVertex;
		}
		


		sidenormal = new vector2d[vertexCount];
		UpdateNormalVector();

		velocity = {0,0};
		angular_velocity = 0.0;

		inverseMass = 1.0/mass;
		inverse_inertensor = 1.0/inertia_tensor;

		pin = false;

	}
	void UpdateNormalVector(){
		//Update normal unit vector of each line
		matrix2x2 rotation;
		rotation.m[0][0] = 0;
		rotation.m[0][1] = 1;
		rotation.m[1][0] = -1;
		rotation.m[1][1] = 0;

		Vertex* CurrentVertex = vertex;
		for(int i=0;i<vertexCount;i++){
			Vertex* NextVertex = CurrentVertex->nextVertex;
			sidenormal[i] = (NextVertex->point-CurrentVertex->point) * rotation;
			sidenormal[i] = sidenormal[i].UnitVector();
			CurrentVertex = NextVertex;
		}
		
	}
	void Render(SDL_Renderer* render){
		Vertex* CurrentVertex = vertex;
		for(int i=0;i<vertexCount;i++){
			Vertex* NextVertex = CurrentVertex->nextVertex;
			SDL_RenderDrawLine(render, CurrentVertex->point.x + centerOfMass.x, CurrentVertex->point.y + centerOfMass.y, NextVertex->point.x + centerOfMass.x, NextVertex->point.y + centerOfMass.y);
			CurrentVertex = NextVertex;
		}
		return;
	}

	float CalculateRotateImpulse(vector2d actionPoint, vector2d Impulse){
		vector2d r = actionPoint - centerOfMass;
		return r.cross(Impulse);
	}

	void AddImpulse(vector2d F, float Time){
		if(pin)
			return;

		velocity = velocity + F * Time * inverseMass;
	}
	void AddRotationImpulse(double impulseSize, float Time){
		if(pin)
			return;
		angular_velocity += impulseSize * Time * inverse_inertensor;
	}
	void Rotate(float deltaTime){
		if(pin)
			return;

		matrix2x2 rotation;
		rotation.m[0][0] = cos(angular_velocity * deltaTime);
		rotation.m[0][1] = sin(angular_velocity * deltaTime);
		rotation.m[1][0] = -sin(angular_velocity * deltaTime);
		rotation.m[1][1] = cos(angular_velocity * deltaTime);

		Vertex* CurrentVertex = vertex;
		for(int i=0;i<vertexCount;i++){
			CurrentVertex->point = CurrentVertex->point * rotation;
			CurrentVertex = CurrentVertex->nextVertex;
		}

		UpdateNormalVector();
	}

	void MoveCenter(vector2d displacement){
		if(pin)
			return;
		
		centerOfMass = centerOfMass + displacement;
	}
	void Update(float deltaTime){
		MoveCenter(velocity * deltaTime);
		Rotate(deltaTime);
	}
	void isFixed(bool _pin){
		pin = _pin;
	}
	~RigidBody(){
		Vertex* CurrentVertex = vertex;
		for(int i=0;i<vertexCount;i++){
			Vertex* NextVertex = CurrentVertex->nextVertex;
			delete CurrentVertex;
			CurrentVertex = NextVertex;
		}
		delete[] sidenormal;
	}
};

class CollisionData{
public:
	int ContactPointCount;
	float penetration_depth;
	vector2d contact_normal;
	vector2d contact_point;	
};

class CollisionDataList{
public:
	CollisionData data[8];
	int dataCount;
public:
	CollisionDataList() : dataCount(0){}
};

class RigidBodyControl{
private:
	RigidBody* RigidBodyList[100];
	int RigidBodyCount;
public:
	RigidBodyControl():RigidBodyCount(0){
	}
	int AddRigidBody(vector2d* _vertex, int _vertexCount, 
	vector2d _centerOfMass, double _mass, double _inertia_tensor){
		RigidBodyList[RigidBodyCount++] = new RigidBody(_vertex, _vertexCount, _centerOfMass,  _mass, _inertia_tensor);
		return RigidBodyCount-1;
		//Add RigidBody and return key value
	}
	RigidBody* GetRigidBody(int key){
		if(key >= RigidBodyCount){
			return nullptr;
		}
		return RigidBodyList[key];
	}
	bool SimpleCollisionCheck(const RigidBody& rigid1, const RigidBody& rigid2){
		//Check By Circle
		float centerDistance = (rigid1.centerOfMass - rigid2.centerOfMass).Size();
		if(centerDistance > std::max(rigid1.maxLength,rigid2.maxLength)){
			if(centerDistance <= rigid1.maxLength + rigid2.maxLength)
				return true;
			else
				return false;
		}
		else{
			return true;
		}

		return false;
	}
	CollisionDataList GetCollisionData(const RigidBody& rigid1, const RigidBody& rigid2){
		//find inner vertex
		//return 1~4 collisiondata
		CollisionDataList collisiondata;

		Vertex* Rigid1CurrentVertex = rigid1.vertex;
		for(int i=0;i<rigid1.vertexCount;i++){
			CollisionData tmpcollisiondata;
			tmpcollisiondata = GetPointCollisionData(Rigid1CurrentVertex, rigid1.centerOfMass, rigid2);
			if(tmpcollisiondata.ContactPointCount == 0){
				Rigid1CurrentVertex = Rigid1CurrentVertex->nextVertex;
				continue;
			}

			collisiondata.data[collisiondata.dataCount++] = tmpcollisiondata;
			Rigid1CurrentVertex = Rigid1CurrentVertex->nextVertex;
		}	

		return collisiondata;
	}

	CollisionData GetPointCollisionData(const Vertex* vertex, const vector2d com, const RigidBody& rigid){
		CollisionData collisiondata;
		collisiondata.ContactPointCount = 1;
		collisiondata.penetration_depth = 100000.0;

		Vertex* RigidCurrentVertex = rigid.vertex;
		for(int j=0;j<rigid.vertexCount;j++){
			float dotproduct = rigid.sidenormal[j].dot(vertex->point + com - (RigidCurrentVertex->point + rigid.centerOfMass));
			if(dotproduct > 0){
				collisiondata.ContactPointCount = 0;
				//point is outside the polygon
				break;
			}
			if(-dotproduct < collisiondata.penetration_depth){
				//renovate penetration depth and contact normal of the point
				collisiondata.penetration_depth = -dotproduct;
				collisiondata.contact_normal = rigid.sidenormal[j];
				//find contact vertex by penetration depth
				collisiondata.contact_point = vertex->point + collisiondata.contact_normal * collisiondata.penetration_depth + com;
			}

			RigidCurrentVertex = RigidCurrentVertex->nextVertex;
		}

		return collisiondata;
	}

	void Update(float deltaTime){
		//Update and Correction * 15

		for(int i=0;i<RigidBodyCount;i++){

			RigidBodyList[i]->AddImpulse(vector2d(0.0, RigidBodyList[i]->mass * 10.0), deltaTime);
			//Gravity

			RigidBodyList[i]->Update(deltaTime);
			//Update position

			//Contact point = penetration depth max point
		}

		for(int i=0;i<15;i++){
			//Correction 15 times
			for(int j=0;j<RigidBodyCount;j++){
				PenetrationRevise(j,deltaTime);
			}
		}
	}
	void PenetrationRevise(int i, float deltaTime){
		for(int j=0;j<RigidBodyCount;j++){
			if(j==i)
				continue;
			if(!SimpleCollisionCheck(*RigidBodyList[i],*RigidBodyList[j]))
				continue;	
			
			CollisionDataList collisiondata[2];
			collisiondata[0] = GetCollisionData(*RigidBodyList[i],*RigidBodyList[j]);
			collisiondata[1] = GetCollisionData(*RigidBodyList[j],*RigidBodyList[i]);
			for(int k=0;k<collisiondata[1].dataCount;k++){
				collisiondata[1].data[k].contact_normal = collisiondata[1].data[k].contact_normal * (-1); //normalize
			}
			
			float max_penetration_depth = 0.0f;
			vector2d contact_normal;
			for(int k=0;k<2;k++){
				for(int l=0;l<collisiondata[k].dataCount;l++){
					//max penetration depth means contact point
					max_penetration_depth = std::max(max_penetration_depth,collisiondata[k].data[l].penetration_depth);
					contact_normal = collisiondata[k].data[l].contact_normal;
				}
			}

			if(RigidBodyList[i]->pin){
				RigidBodyList[j]->MoveCenter(contact_normal * max_penetration_depth * (-1));
			}
			else if(RigidBodyList[j]->pin){
				RigidBodyList[i]->MoveCenter(contact_normal * max_penetration_depth);
			}
			if(!RigidBodyList[i]->pin && !RigidBodyList[j]->pin){
				RigidBodyList[i]->MoveCenter(contact_normal * max_penetration_depth * 
				(RigidBodyList[j]->mass / (RigidBodyList[i]->mass + RigidBodyList[j]->mass)));
				RigidBodyList[j]->MoveCenter(contact_normal * max_penetration_depth * 
				(RigidBodyList[i]->mass / (RigidBodyList[i]->mass + RigidBodyList[j]->mass)) * (-1));
			}
			
			for(int k=0;k<2;k++){
				for(int l=0;l<collisiondata[k].dataCount;l++){
					float dImpulse = 0.0f;
					float dRImpulsei = 0.0f;
					float dRImpulsej = 0.0f;

					float mue = 0.0f;
					float friction_Impulse = 0.0f;

					matrix2x2 rotation;
					rotation.m[0][0] = 0;
					rotation.m[0][1] = 1;
					rotation.m[1][0] = -1;
					rotation.m[1][1] = 0;

					if(collisiondata[k].data[l].penetration_depth < max_penetration_depth - 0.000001)
						continue;
					
					float invMass = 1.0 / (RigidBodyList[i]->inverseMass + RigidBodyList[j]->inverseMass + 
						pow(collisiondata[k].data[l].contact_normal.dot((collisiondata[k].data[l].contact_point - RigidBodyList[i]->centerOfMass) * rotation),2) * RigidBodyList[i]->inverse_inertensor
						+ pow(collisiondata[k].data[l].contact_normal.dot((collisiondata[k].data[l].contact_point - RigidBodyList[j]->centerOfMass) * rotation),2) * RigidBodyList[j]->inverse_inertensor);
					
					float invRMass = 1.0 / (RigidBodyList[i]->inverse_inertensor + RigidBodyList[j]->inverse_inertensor);

					vector2d iRvelocity = (collisiondata[k].data[l].contact_point - RigidBodyList[i]->centerOfMass) * rotation * RigidBodyList[i]->angular_velocity;
					vector2d jRvelocity = (collisiondata[k].data[l].contact_point - RigidBodyList[j]->centerOfMass) * rotation * RigidBodyList[j]->angular_velocity;
					
					
					dImpulse = (-1) * collisiondata[k].data[l].contact_normal.dot(RigidBodyList[i]->velocity + iRvelocity - RigidBodyList[j]->velocity - jRvelocity) * invMass * 1.5; 
					
					if(dImpulse < 0)
						dImpulse = 0.0;
					
					dRImpulsei = RigidBodyList[i]->CalculateRotateImpulse(collisiondata[k].data[l].contact_point, collisiondata[k].data[l].contact_normal * dImpulse);
					dRImpulsej = RigidBodyList[j]->CalculateRotateImpulse(collisiondata[k].data[l].contact_point, collisiondata[k].data[l].contact_normal * -dImpulse);

					mue = (collisiondata[k].data[l].contact_normal * rotation).dot((RigidBodyList[i]->velocity + iRvelocity - RigidBodyList[j]->velocity - jRvelocity));

					if(abs(mue) > 0.8)
						mue = (mue / abs(mue)) * 0.8f;

					friction_Impulse = mue * abs(dImpulse);

					RigidBodyList[i]->AddImpulse(collisiondata[k].data[l].contact_normal * (dImpulse), deltaTime);
					RigidBodyList[i]->AddRotationImpulse(dRImpulsei, deltaTime);

					RigidBodyList[i]->AddImpulse(collisiondata[k].data[l].contact_normal * rotation * (-friction_Impulse), deltaTime);
					dRImpulsei = RigidBodyList[i]->CalculateRotateImpulse(collisiondata[k].data[l].contact_point, collisiondata[k].data[l].contact_normal * rotation * (-friction_Impulse));
					RigidBodyList[i]->AddRotationImpulse(dRImpulsei, deltaTime);
					


					RigidBodyList[j]->AddImpulse(collisiondata[k].data[l].contact_normal * (-dImpulse), deltaTime);
					RigidBodyList[j]->AddRotationImpulse(dRImpulsej, deltaTime);
					
					RigidBodyList[j]->AddImpulse(collisiondata[k].data[l].contact_normal * rotation * (friction_Impulse), deltaTime);
					dRImpulsej = RigidBodyList[j]->CalculateRotateImpulse(collisiondata[k].data[l].contact_point, collisiondata[k].data[l].contact_normal * rotation * (friction_Impulse));
					RigidBodyList[j]->AddRotationImpulse(dRImpulsej, deltaTime);
				}


			}		
		
		}
	}

	void Render(SDL_Renderer* render){
		for(int i=0;i<RigidBodyCount;i++){
			RigidBodyList[i]->Render(render);
		}
		return;
	}

	~RigidBodyControl(){
		for(int i=0;i<RigidBodyCount;i++){
			delete RigidBodyList[i];
		}
	}
};