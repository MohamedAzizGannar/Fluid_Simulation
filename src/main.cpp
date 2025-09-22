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
const double PARTICLE_RADIUS = 1.;
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

std::vector<Particle> populateParticles(int rows, int cols , int depth){
    std::vector<Particle> particles;
    double baseSpacing = PARTICLE_RADIUS * 2.;
    int count = 0;
    for(int i = 0 ; i < rows; i++){
        for(int j = 0; j < cols ; j++){
            for(int k = 0; k < depth; k++){
                float3 pos = {(double)i + baseSpacing,(double)j + baseSpacing,(double)k + baseSpacing};
                float3 vel = {getRoundedRandom(-1.,1.,1),getRoundedRandom(-1.,1.,1),getRoundedRandom(-1.,1.,1)};
                float3 acc = {0.,0.,0.};

                Particle newParticle(pos,vel,acc,count);
                count++;

                particles.push_back(newParticle);
                
            }
        }
    }
    return particles;
}


int main (int argc, char** argv)
{
    std::vector<Particle> particles = populateParticles(10,10,10);

    float3 minBounds = {-1,-1,-1};
    float3 maxBounds = {12,12,12};

    Collider boxCollider(minBounds,maxBounds);

    FluidSimulation simulation(particles,boxCollider);
    simulation.setTargetPhysicsRate(60.0f);

    int frames = 2;
    int currentFrame = 0;
    auto startTime = std::chrono::high_resolution_clock::now();
    while (currentFrame < frames)
    {
        simulation.update();
        std::cout<<"Particle 0 : \n"<<
                    "Position : ("<<simulation.getParticles()[0].getPosition().x<<","<<simulation.getParticles()[0].getPosition().y<<","<<simulation.getParticles()[0].getPosition().z<<")\n"<<
                    "Velocity : ("<<simulation.getParticles()[0].getVelocity().x<<","<<simulation.getParticles()[0].getVelocity().y<<","<<simulation.getParticles()[0].getVelocity().z<<")\n"<<
                    "Density : "<<simulation.getParticles()[0].getDensity()<<"\n"<<
                    "Pressure : "<<simulation.getParticles()[0].getDensity()<<"\n\n";

        currentFrame ++;
    }
    auto endTime = std::chrono::high_resolution_clock::now();
    auto simulationTime = std::chrono::duration<double>(endTime-startTime).count();
    auto frameTime = simulationTime / frames;

    std::cout<<"Total Simulation Time for "<<frames<<"Iterations : "<<simulationTime<<"\n";
    std::cout<<"Average Frame Time for "<<frames<<"Iterations : "<<frameTime<<"\n";
    std::cout<<"Average FPS "<<frames<<"Iterations : "<<1/frameTime<<"\n";


    return 0;
}