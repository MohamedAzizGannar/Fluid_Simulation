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

void populateParticles(int numberOfParticles, std::vector<Particle>& particles){
    double baseSpacing = 2.;
    int side = std::cbrt(numberOfParticles);
    int count = 0;
    float maxX = 0;
    float maxY = 0;

    for(int x = 0 ; x < side; x++){
        for(int y = 0; y < side ; y++){
            for(int z = 0; z < side; z++){
                float3 pos = {x*baseSpacing,y*baseSpacing,z*baseSpacing};
                float3 vel = {getRoundedRandom(-1.,1.,1),getRoundedRandom(-1.,1.,1),getRoundedRandom(-1.,1.,1)};
                float3 acc = {0.,0.,0.};
                if(x*baseSpacing > maxX) maxX = x*baseSpacing;
                if(y*baseSpacing > maxY) maxY = y*baseSpacing;
                
                particles.emplace_back(pos,vel,acc,count);
                count++;                
            }
        }
    }
    std::cout<<maxX<<"\n "<<maxY<<"\n";
}


int main (int argc, char** argv)
{
    std::vector<Particle> particles;
    particles.reserve(1500);
    populateParticles(1500,particles);

    float3 minBounds = {0,0,0};
    float3 maxBounds = {20,20,20};

    Collider boxCollider(minBounds,maxBounds);

    FluidSimulation simulation(particles,boxCollider);
    simulation.setTargetPhysicsRate(120.0f);
    
    InitWindow(SCREEN_WIDTH,SCREEN_HEIGHT,"simulation");
    SetTargetFPS(120);

    double multiplicationFactorX = 900. / 20.;
    double multiplicationFactorY = 700. / 20.;


    while(!WindowShouldClose()){
        simulation.update();
        BeginDrawing();
        ClearBackground(RAYWHITE);


        // Draw bounds
        DrawRectangleLines(0,0,900,700,RED);

        // Draw particles
        for(int i = 0; i < simulation.getParticles().size();i ++){
    
            auto pos = simulation.getParticles()[i].getPosition();
            double x = pos.x * multiplicationFactorX;
            double y = pos.y * multiplicationFactorY;

            auto vel = simulation.getParticles()[i].getVelocity();

            double velx = vel.x ;
            double vely = vel.y ;
            DrawCircle(x,y,5.f,BLUE);




        }
        DrawFPS(10, 10);
        EndDrawing();
    }
    CloseWindow();

    return 0;
}