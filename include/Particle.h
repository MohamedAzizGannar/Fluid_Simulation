#ifndef PARTICLE_H
#define PARTICLE_H

#include <array>
class Particle{

    private:
    std::array<double,2> position;
    std::array<double,2> velocity;
    std::array<double,2> acceleration;


    double radius;

    public:
    //Constructor
    Particle(std::array<double,2> initialPosition,std::array<double,2> initialVelocity,std::array<double,2> intialAcceleration,const double radius);

    //Default Constructor (Array Initialization)
    Particle();

    //Getters
    std::array<double,2> getPosition();
    std::array<double,2> getVelocity();
    std::array<double,2> getAcceleration();
    double getRadius();

    //Setters
    void setPosition(double newPosition[2]);
    void setVelocity(double newVelocity[2]);
    void setAcceleration(double newAcceleration[2]);

    //Update Position 
    void updatePosition();

    //Update Velocity

    void updateVelocity();

    

};
#endif