#include <Particle.h>
#include <array>

Particle::Particle(std::array<double,2> initialPosition,std::array<double,2> initialVelocity,std::array<double,2> intialAcceleration,const double initialRadius){
    for(int i =0; i < 2; i++){
        position[i] = initialPosition[i];
        velocity[i] = initialVelocity[i];
        acceleration[i] = intialAcceleration[i];
    }
    radius = initialRadius;
}
//Default Constructor
Particle::Particle()
{
    for(int i =0; i < 2; i++){
        position[i] = 0;
        velocity[i] = 0;
        acceleration[i] = 0;
    }
    radius = 1.0;
}

void Particle::setPosition(double newPosition[2]){
    position[0] = newPosition[0];
    position[1] = newPosition[1];
}
void Particle::setVelocity(double newVelocity[2]){
    velocity[0] = newVelocity[0];
    velocity[1] = newVelocity[1];
}
void Particle::setAcceleration(double newAcceleration[2]){
    acceleration[0] = newAcceleration[0];
    acceleration[1] = newAcceleration[1];
}


std::array<double,2> Particle::getPosition(){
    return position;
}

std::array<double,2> Particle::getAcceleration(){
    return acceleration;
}
std::array<double,2> Particle::getVelocity(){
    return velocity;
}



void Particle::updatePosition()
{
    position[0] += velocity[0];
    position[1] += velocity[1];
}

void Particle::updateVelocity()
{
    velocity[0] += acceleration[0];
    velocity[1] += acceleration[1];
}