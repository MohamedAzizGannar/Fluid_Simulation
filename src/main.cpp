#include <iostream>
#include <iomanip>

#include <cmath>

#include <array>
#include <string>

#include<chrono>
#include<thread>

#include <random>
#include <Particle.h>



const double CORE_RADIUS = 2.0;

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
    std::vector<Particle> particles;
    for(int i = 0; i< rows;i++){
        for(int j = 0; j < cols; j++){
            for(int k = 0; k < depth; k++){
                std::array<double, 3> pos = {static_cast<double>(i),static_cast<double>(j),static_cast<double>(k)};
                std::array <double, 3> cst = {0,0,0};
                Particle newParticle = Particle(pos,cst,cst);
                particles.push_back(newParticle);
            }
        }
    }
    return particles;
}

int main (int argc, char** argv)
{
    std::cout<<"Running\n";
    std::vector<Particle> particles = initialiseParticlesArray(5,2,3);
    
   
    const auto interval_ms = std::chrono::milliseconds(1800);

    auto nextFrame = std::chrono::steady_clock::now();

    int counter = 1;

    while(true){
        int index = 1;
        std::cout<<std::setw(10)<<"Counter : "<<counter<<std::endl;
        for(Particle& particle : particles)
        {
            
            particle.updateDensity(particles,CORE_RADIUS);
            particle.updatePressure(particles,CORE_RADIUS);

            auto [xPressure, yPressure, zPressure] = particle.calculatePressureForce(particles,CORE_RADIUS);
            auto [xViscosity, yViscosity, zViscosity] = particle.calculateViscosityForce(particles,CORE_RADIUS);
            auto [x, y, z] = particle.getPosition();


            std::string valuesData = "Particle " +std::to_string(index)+ " Position : "+ "(" + std::to_string(x) + ", " + std::to_string(y) + ", " +std::to_string(z) + ")" + " || Pressure : " + std::to_string(particle.getPressure()) + " || Density : " + std::to_string(particle.getDensity()) + "\n";
            std::cout<<valuesData;

            std::string pressureForceData = "PressureForce : (" + std::to_string(xPressure) + ", " + std::to_string(yPressure) + ", " +std::to_string(zPressure) + ")\n";
            std::string viscosityForceData = "ViscosityForce : (" + std::to_string(xViscosity) + ", " + std::to_string(yViscosity) + ", " +std::to_string(zViscosity) + ")\n";
            std::cout<<pressureForceData<<viscosityForceData;
            std::cout<<std::setw(50)<<std::setfill('-')<<std::endl;
            index++;
        }

        std::cout<<std::endl;

        counter ++;
        

        if(counter>1)
        {
            break;
        }

        nextFrame += interval_ms;
        std::this_thread::sleep_until(nextFrame);  
    }
    return 0;
}