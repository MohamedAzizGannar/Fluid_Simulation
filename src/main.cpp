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
                std::array <double, 3> vel = {1,0.5,2};

                Particle newParticle = Particle(pos,vel,cst);
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
    std::array<double,3> minBounds = {-0.5, -0.5, -0.5};
    std::array<double,3> maxBounds = {10.5, 10.5, 10.5};

    Collider boxCollider = Collider(minBounds,maxBounds);
   
    const auto interval_ms = std::chrono::milliseconds(450);

    auto nextFrame = std::chrono::steady_clock::now();

    int counter = 1;

    while(true){
        int index = 1;
        std::cout<<std::setw(10)<<"Counter : "<<counter<<std::endl;
        for(Particle& particle : particles)
        {
            std::cout<<std::boolalpha;
            std::cout<<boxCollider.resolveSphereAABBCollision(particle)<<std::endl;

            particle.updatePosition();

            auto [x, y, z] = particle.getPosition();
            auto [vx, vy, vz] = particle.getVelocity();

            std::string valuesData = "Particle " +std::to_string(index)+ " Position : "+ "(" + std::to_string(x) + ", " + std::to_string(y) + ", " +std::to_string(z) + ")"
            +" Velocity : "+ "(" + std::to_string(vx) + ", " + std::to_string(vy) + ", " +std::to_string(vz) + ")";
            
            std::cout<<valuesData<<std::endl;
            index++;
        }

        std::cout<<std::endl;

        counter ++;
        

        if(counter>50)
        {
            break;
        }

        nextFrame += interval_ms;
        std::this_thread::sleep_until(nextFrame);  
    }
    return 0;
}