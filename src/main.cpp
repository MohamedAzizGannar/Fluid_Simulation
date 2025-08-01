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

double getRoundedRandom(double min, double max, int decimals)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dist(min,max);

    double scale = std::pow(10.0,decimals);
    return std::round(dist(gen) * scale)/scale;
}

std::vector<Particle> initialiseParticlesArray(int arraySize)
{
    std::vector<Particle> particles(arraySize);
    for(int i = 0; i< particles.size();i++){
        //Initialising Positions in a grid pattern
        double xPosition = (i) % rows;

        double yPosition = std::floor((i)/rows);


        double newPosition[3] = {xPosition,yPosition,0.};

        //Adding a random velocity to each particle
        

        double xVelocity = getRoundedRandom(-3.0,3.0,1);
        double yVelocity = getRoundedRandom(-3.0,3.0,1);
        double zVelocity = getRoundedRandom(-3.0,3.0,1);


        double newVelocity[3] = {xVelocity,yVelocity,zVelocity};

        particles[i].setPosition(newPosition);
        particles[i].setVelocity(newVelocity);
    }
    return particles;


}

int main (int argc, char** argv){
    
    std::vector<Particle> particles = initialiseParticlesArray(10);

    
    
    
     

    int i = 0;
    while( i < particles.size()){
        for(int row = 0; row < rows ; row++){
            while(particles[i].getPosition()[1] == row){
                std::cout<<std::setw(2)<<"("<<particles[i].getPosition()[0]<<", "<<particles[i].getPosition()[1]<<", "<<particles[i].getPosition()[2]<<")";
                ++i;
            }
            std::cout<<std::endl;
        }
    }

    
   
    const auto interval_ms = std::chrono::milliseconds(900);

    auto nextFrame = std::chrono::steady_clock::now();

    int counter = 1;

    while(true){

        std::cout<<std::setw(5)<<"Counter: "<<counter<<std::endl;
        for(int i = 0; i < particles.size(); i++){
            particles[i].updatePosition();
            auto [x,y,z] = particles[i].getPosition();
            std::cout<<std::setw(5)<<"("<<x<<","<<y<<", "<<z<<")";
            if((i+1)%cols == 0)
            {
                std::cout<<std::endl;
            }
        }
        std::cout<<std::endl;

        counter ++;
        

        if(counter>6)
        {
            break;
        }

        nextFrame += interval_ms;
        std::this_thread::sleep_until(nextFrame);  
    }

    return 0;
}