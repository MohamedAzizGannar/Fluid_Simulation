#include <Particle.h>
#include <array>

Particle::Particle(std::array<double,3> initialPosition,std::array<double,3> initialVelocity,std::array<double,3> intialAcceleration){
    for(int i =0; i < 3; i++){
        position[i] = initialPosition[i];
        velocity[i] = initialVelocity[i];
        acceleration[i] = intialAcceleration[i];
    }
}
//Default Constructor
Particle::Particle()
{
    for(int i =0; i < 3; i++){
        position[i] = 0;
        velocity[i] = 0;
        acceleration[i] = 0;
    }
}

void Particle::setPosition(double newPosition[3]){
    position[0] = newPosition[0];
    position[1] = newPosition[1];
    position[2] = newPosition[2];

}
void Particle::setVelocity(double newVelocity[3]){
    velocity[0] = newVelocity[0];
    velocity[1] = newVelocity[1];
    velocity[2] = newVelocity[2];

}
void Particle::setAcceleration(double newAcceleration[3]){
    acceleration[0] = newAcceleration[0];
    acceleration[1] = newAcceleration[1];
    acceleration[2] = newAcceleration[2];

}


std::array<double,3> Particle::getPosition(){
    return position;
}

std::array<double,3> Particle::getAcceleration(){
    return acceleration;
}
std::array<double,3> Particle::getVelocity(){
    return velocity;
}



void Particle::updatePosition()
{
    position[0] += velocity[0];
    position[1] += velocity[1];
    position[2] += velocity[2];

}

void Particle::updateVelocity()
{
    velocity[0] += acceleration[0];
    velocity[1] += acceleration[1];
    velocity[2] += acceleration[2];

    
}