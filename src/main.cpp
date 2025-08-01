#include <iostream>
#include <cmath>
#include <array>
#include <string>
#include <iomanip>

#include <Particle.h>


int main (int argc, char** argv){
    
    std::array<Particle,50> particles;

    int rows = 5;
    int cols = 10;
    for(int i = 0; i< particles.size();i++){
        double xPosition = (i) % 10;

        double yPosition = std::floor((i)/10);
        std::cout<<"y = "<<yPosition<<"\n";

        double newPosition[2] = {xPosition,yPosition};

        particles[i].setPosition(newPosition);
    }


    for(auto particle : particles)
    {
        std::array<double,2> pos = particle.getPosition();

        std::cout<<pos[0]<<","<<pos[1]<<"\n";
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

    std::cout<<"AAAAA"<<std::endl;

    return 0;
}