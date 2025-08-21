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
bool running = true;
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

std::vector<Particle> initialiseParticlesArray(int rows,int cols,int depth)
{
    int counter = 0;
    std::vector<Particle> particles;
    for(int i = 0; i< rows;i++){
        for(int j = 0; j < cols; j++){
            for(int k = 0; k < depth; k++){
                const float3 pos = {static_cast<double>(i)*2. + getRoundedRandom(-0.1, 0.1, 2),
                             static_cast<double>(j)*2. + getRoundedRandom(-0.1, 0.1, 2),
                             static_cast<double>(k)*2. + getRoundedRandom(-0.1, 0.1, 2)};
                const float3 cst = {0,0,0};
                const float3 vel = {static_cast<double>(i) + getRoundedRandom(-0.1, 0.1, 2),
                             static_cast<double>(j) + getRoundedRandom(-0.1, 0.1, 2),
                             static_cast<double>(k) + getRoundedRandom(-0.1, 0.1, 2)};

                Particle newParticle = Particle(pos,vel,cst,counter);
                counter++;
                particles.push_back(newParticle);
            }
        }
    }
    return particles;
}

float3 minBounds = {-0.5, -0.5, -0.5};
float3 maxBounds = {50.5, 50.5, 50.5};

Collider boxCollider = Collider(minBounds,maxBounds);
std::vector<Particle> particles = initialiseParticlesArray(20,20,20);

void runSimulation(double durationSeconds){
    FluidSimulation sim(particles,boxCollider);
    auto testStart = std::chrono::high_resolution_clock::now();
    auto lastPrint = testStart;
    while(true){
        auto now = std::chrono::high_resolution_clock::now();

        auto elapsed = std::chrono::duration<double>(now - testStart).count();

        if (elapsed >= durationSeconds) break;

        sim.update(8);

         auto printElapsed = std::chrono::duration<double>(now- lastPrint).count();
        if(printElapsed >1.0){
            std::cout<<1000/printElapsed<<std::endl;
            lastPrint = now;
        }

    }

}

int main (int argc, char** argv)
{
    unsigned int frames = 60;
    unsigned int repetitions = 10;
    double targetFPS = 60.;
    double targetFrameTime = 1/targetFPS;


    runSimulation(5);

    return 0;
}