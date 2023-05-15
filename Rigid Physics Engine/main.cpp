#include<iostream>
#include<string>
#include<SDL2/SDL.h>
#include<SDL2/SDL_image.h>
#include"RigidBody.h"

namespace DLG{
    class SDL2App{
    private:
        int WIDTH, HEIGHT;
        SDL_Window* window;
    public:
        SDL_Renderer* renderer;
    public:
        SDL2App(int, int);
        int InitApp();
        void PrepareScreen();
        void PrepareScreen(int r, int g, int b);
        void UpdateScreen();
        SDL_Texture* GetTexture(std::string, std::string);
        void DestroyApp();
    };
}

DLG::SDL2App::SDL2App(int width, int height){
    WIDTH = width;
    HEIGHT = height;
    window = nullptr;
    renderer = nullptr;
}

int DLG::SDL2App::InitApp(){
    SDL_Init(SDL_INIT_VIDEO);

    window = SDL_CreateWindow("Hello SDL world", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WIDTH, HEIGHT, SDL_WINDOW_ALLOW_HIGHDPI);
    if(window==NULL){
        std::cout<<"Could not create Window: " << SDL_GetError() <<std::endl;
        return -1;
    }

    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED|SDL_TEXTUREACCESS_TARGET);
    if(renderer==NULL){
        std::cout<<"failed to create renderer"<<std::endl;
        return -1;
    }

    return 0;
}

void DLG::SDL2App::PrepareScreen(){
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, SDL_ALPHA_OPAQUE);
    SDL_RenderClear(renderer);
}

void DLG::SDL2App::PrepareScreen(int r, int g, int b){
    SDL_SetRenderDrawColor(renderer, r, g, b, SDL_ALPHA_OPAQUE);
    SDL_RenderClear(renderer);
}

void DLG::SDL2App::UpdateScreen(){
    SDL_RenderPresent(renderer);
}

SDL_Texture* DLG::SDL2App::GetTexture(std::string filename, std::string suffix){
    SDL_Surface *imgSurface = IMG_Load((filename+suffix).c_str());
    if(imgSurface==NULL){
        std::cout<<"failed to load image"<<std::endl;
        return NULL;
    }
    SDL_Texture *imgTexture = SDL_CreateTextureFromSurface(renderer, imgSurface);
	SDL_FreeSurface(imgSurface);
    return imgTexture;
}
/*
suffix ex: .jpg, .png, .bmp
*/

void DLG::SDL2App::DestroyApp(){
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    return;
}


class HRectangle{
public:
    vector2d point[4];
};

class HTriangle{
public:
    vector2d point[3];
};

HRectangle GetRectangle(float width, float height){
    HRectangle rect;
    rect.point[0] = {-width/2,height/2};
    rect.point[1] = {width/2,height/2};
    rect.point[2] = {width/2,-height/2};
    rect.point[3] = {-width/2,-height/2};
    return rect;
}

float GetRectInertiaTensor(float width, float height, float mass){
    return mass * (width*width + height*height) / 12;
}

HTriangle GetTriangle(vector2d p1, vector2d p2, vector2d p3){
    HTriangle tri;
    tri.point[0] = p1;
    tri.point[1] = p2;
    tri.point[2] = p3;
    return tri;
}



int main(int argc, char* argv[]){
    DLG::SDL2App app(640, 720);
    app.InitApp();

    SDL_Event windowEvent;
    bool quitapp = false;

    RigidBodyControl rbc;

    /*
    rbc.AddRigidBody(GetRectangle(70,50).point,4,{450,340},40.0,
    GetRectInertiaTensor(70,50,40));
    rbc.GetRigidBody(0)->AddImpulse({-1000,0},1.0);

    rbc.AddRigidBody(GetRectangle(161,50).point,4,{360,305},60.0,
    GetRectInertiaTensor(160,50,60));
    rbc.AddRigidBody(GetRectangle(160,100).point,4,{280,105},200.0,
    GetRectInertiaTensor(160,100,200));
    rbc.AddRigidBody(GetRectangle(20,70).point,4,{280,500},60.0,
    GetRectInertiaTensor(20,70,60));
    */

    rbc.AddRigidBody(GetRectangle(800,70).point,4,{400,700},1000.0,
    GetRectInertiaTensor(600,70,1000));
    rbc.GetRigidBody(0)->isFixed(true);
    
    int width;
    int height;

    while(!quitapp){
        if(SDL_PollEvent(&windowEvent)){
            switch (windowEvent.type)
            {
            case SDL_QUIT:
                quitapp = true;
                /* code */
                break;
            case SDL_MOUSEBUTTONDOWN:
                width = rand() % 100 + 20;
                height = rand() % 100 + 20;
                rbc.AddRigidBody(GetRectangle(width,height).point,4,{(float)windowEvent.button.x,(float)windowEvent.button.y},width*height,
                GetRectInertiaTensor(width,height,width*height));
                break;
            default:
                break;
            }
        }
        app.PrepareScreen(150,150,150);
        SDL_SetRenderDrawColor(app.renderer, 242, 242, 242, 255);

        rbc.Update(0.01);
        rbc.Render(app.renderer);
        app.UpdateScreen();

        SDL_Delay(10);
    }

    app.DestroyApp();
    SDL_Quit();
}