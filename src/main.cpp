#include <iostream>
#include <iomanip>

#include <cmath>

#include <array>
#include <string>

#include<chrono>
#include<thread>

#include <random>
#include <Particle.h>

const int rows = 2;
const int cols = 5;

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

std::vector<Particle> initialiseParticlesArray(int arraySize)
{
    std::vector<Particle> particles(arraySize);
    for(int i = 0; i< rows;i++){
        for(int j = 0; j < cols; j++){
            std::array<double, 3> pos = {static_cast<double>(i),static_cast<double>(j),0.};
            std::array <double, 3> cst = {0,0,0};
            Particle newParticle = Particle(pos,cst,cst);
            particles.push_back(newParticle);
        }
    }
    return particles;
}

int main (int argc, char** argv)
{
    std::cout<<"Running\n";
    std::vector<Particle> particles = initialiseParticlesArray(10);
    
   
    const auto interval_ms = std::chrono::milliseconds(900);

    auto nextFrame = std::chrono::steady_clock::now();

    int counter = 1;

    while(true){

        std::cout<<std::setw(5)<<"Counter: "<<counter<<std::endl;
        for(int index = 0; index < particles.size(); index ++){
            auto [x,y,z] = particles[index].calculatePressureForce(particles,CORE_RADIUS);
            auto pressure = particles[index].calculatePressure(particles,CORE_RADIUS);
            std::string output = "Pressure :" + std::to_string(roundDouble(pressure,2)) + "Pressure Force Vector: (" + std::to_string(roundDouble(x,2)) + ", " + std::to_string(roundDouble(y,2)) + ", "+ std::to_string(roundDouble(z,2))+")\n";
            std::cout<<std::setw(20)<<output;
            
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