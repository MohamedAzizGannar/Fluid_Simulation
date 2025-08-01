#include <iostream>
#include <iomanip>

#include <cmath>

#include <array>
#include <string>

#include<chrono>
#include<thread>

#include <random>
#include <Particle.h>
double getRoundedRandom(double min, double max, int decimals)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dist(min,max);

    double scale = std::pow(10.0,decimals);
    return std::round(dist(gen) * scale)/scale;
}

int main (int argc, char** argv){
    
    std::array<Particle,10> particles;

    int rows = 2;
    int cols = 5;
    for(int i = 0; i< particles.size();i++){
        //Initialising Positions in a grid pattern
        double xPosition = (i) % 10;

        double yPosition = std::floor((i)/10);

        double newPosition[2] = {xPosition,yPosition};

        //Adding a random velocity to each particle
        

        double xVelocity = getRoundedRandom(-3.0,3.0,1);
        double yVelocity = getRoundedRandom(-3.0,3.0,1);

        double newVelocity[2] = {xVelocity,yVelocity};

        particles[i].setPosition(newPosition);
        particles[i].setVelocity(newVelocity);
    }

    
     

    int i = 0;
    while( i < particles.size()){
        for(int row = 0; row < rows ; row++){
            while(particles[i].getPosition()[1] == row){
                std::cout<<std::setw(2)<<"("<<particles[i].getPosition()[0]<<","<<particles[i].getPosition()[1]<<")";
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
            auto [x,y] = particles[i].getPosition();
            std::cout<<std::setw(5)<<"("<<x<<","<<y<<")";
            if((i+1)%5 == 0)
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