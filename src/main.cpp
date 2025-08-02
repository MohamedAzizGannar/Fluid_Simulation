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

std::vector<Particle> initialiseParticlesArray(int rows,int cols)
{
    std::vector<Particle> particles;
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
    std::vector<Particle> particles = initialiseParticlesArray(10,10);
    
   
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

            std::string data = "Particle " +std::to_string(index)+ " : Pressure : " + std::to_string(particle.getPressure()) + "|| Density : " + std::to_string(particle.getDensity()) + "\n";
            std::cout<<data;
            index++;
        }

        std::cout<<std::endl;

        counter ++;
        

        if(counter>10)
        {
            break;
        }

        nextFrame += interval_ms;
        std::this_thread::sleep_until(nextFrame);  
    }
    return 0;
}