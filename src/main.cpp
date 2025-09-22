#include <iostream>
#include <iomanip>

#include <cmath>

#include <array>
#include <string>

#include<chrono>
#include<thread>

#include <random>

#include <Particle.h>
#include <Collider.h>
#include <FluidSimulation.h>
#include<Block.h>

#include <raylib.h>


const double PARTICLE_RADIUS = 1.;

const int SCREEN_HEIGHT = 700;
const int SCREEN_WIDTH = 900;

double getRoundedRandom(double min, double max, int decimals)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dist(min,max);

    double scale = std::pow(10.0,decimals);
    return std::round(dist(gen) * scale)/scale;
}
double roundDouble(double num,uint8_t decimals)
{
    double scale = std::pow(10.0,decimals);
    return std::round(num*scale) /scale;
}

std::vector<Particle> populateParticles(int numberOfParticles, std::vector<Particle> particles){
    double baseSpacing = PARTICLE_RADIUS * 2.;
    int side = std::cbrt(numberOfParticles);
    int count = 0;
    for(int x = 0 ; x < side; x++){
        for(int y = 0; y < side ; y++){
            for(int z = 0; z < side; z++){
                float3 pos = {x*baseSpacing,y*baseSpacing,z*baseSpacing};
                float3 vel = {getRoundedRandom(-1.,1.,1),getRoundedRandom(-1.,1.,1),getRoundedRandom(-1.,1.,1)};
                float3 acc = {0.,0.,0.};

                particles.emplace_back(pos,vel,acc,count);
                count++;                
            }
        }
    }
    return particles;
}


int main (int argc, char** argv)
{
    std::vector<Particle> particles;
    particles.reserve(500);
    populateParticles(500,particles);

    float3 minBounds = {0,0,0};
    float3 maxBounds = {50,50,50};

    Collider boxCollider(minBounds,maxBounds);

    FluidSimulation simulation(particles,boxCollider);
    simulation.setTargetPhysicsRate(60.0f);
    
    InitWindow(SCREEN_WIDTH,SCREEN_HEIGHT,"simulation");
    SetTargetFPS(60);

    Camera3D camera;
    camera.position = {30.0f, 30.0f, 30.0f};
    camera.target = {0.0f, 0.0f, 0.0f};
    camera.up = {0.0f, 1.0f, 0.0f};
    camera.fovy = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    while(!WindowShouldClose()){
        simulation.update();
        BeginDrawing();
        ClearBackground(RAYWHITE);

        BeginMode3D(camera);

        // Draw bounds
        DrawCubeWires({25, 25, 25}, 50, 50, 50, GRAY);

        // Draw particles
        for (const auto& p : particles) {
            auto pos = p.getPosition();
            DrawSphere({(float)pos.x, (float)pos.y, (float)pos.z}, 0.2f, BLUE);
        }

        EndMode3D();

        DrawFPS(10, 10);
        EndDrawing();
    }
    CloseWindow();

    return 0;
}