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
void simulationLoop(std::vector<Particle>& particles, Collider& worldBounds,unsigned int counter)
{
    FluidSimulation simulation(particles,worldBounds);
    auto lastTime = std::chrono::high_resolution_clock::now();
    unsigned int frameCounter = 0;
    while (frameCounter < counter ){
        int index = 1;

        
        auto currentTime = std::chrono::high_resolution_clock::now();
        double deltaTime = std::chrono::duration<double>(currentTime - lastTime).count();
        lastTime = currentTime;

        simulation.update(deltaTime);

        std::cout<<"COUNTER"<<frameCounter<<std::endl;
        frameCounter ++;
        std::this_thread::sleep_for(std::chrono::milliseconds(32));

    }
}



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
                std::array<double, 3> pos = {static_cast<double>(i)*2. + getRoundedRandom(-0.1, 0.1, 2),
                             static_cast<double>(j)*2. + getRoundedRandom(-0.1, 0.1, 2),
                             static_cast<double>(k)*2. + getRoundedRandom(-0.1, 0.1, 2)};
                std::array <double, 3> cst = {0,0,0};
                std::array <double, 3> vel = {static_cast<double>(i) + getRoundedRandom(-0.1, 0.1, 2),
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

int main (int argc, char** argv)
{
    std::vector<Particle> particles = initialiseParticlesArray(7,7,7);
    std::array<double,3> minBounds = {-0.5, -0.5, -0.5};
    std::array<double,3> maxBounds = {1000.5, 1000.5, 1000.5};

    Collider boxCollider = Collider(minBounds,maxBounds);
    unsigned int iterations = 100;
    auto beginTime = std::chrono::high_resolution_clock::now();
    simulationLoop(particles,boxCollider,iterations);
    auto endTime = std::chrono::high_resolution_clock::now();
    auto runTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - beginTime);

    std::cout<<runTime.count()<<"ms RUNTIME\n";

 

  
    return 0;
}