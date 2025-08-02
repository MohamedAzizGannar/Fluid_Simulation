#ifndef PARTICLE_H
#define PARTICLE_H

#include <array>
#include <vector>
class Particle{

    private:
    std::array<double,3> position;
    std::array<double,3> velocity;
    std::array<double,3> acceleration;


    const double radius = 1.0;

    const double mass = 1.0;
    double pressure;

    public:
    //Constructor
    Particle(std::array<double,3> initialPosition,std::array<double,3> initialVelocity,std::array<double,3> intialAcceleration);

    //Default Constructor (Array Initialization)
    Particle();

    //Getters
    std::array<double,3> getPosition();
    std::array<double,3> getAcceleration();
    std::array<double,3> getVelocity();
    double getMass();
    double getRadius();

    //Setters
    void setPosition(double newPosition[3]);
    void setVelocity(double newVelocity[3]);
    void setAcceleration(double newAcceleration[3]);

    //Update Position 
    void updatePosition();

    //Update Velocity
    void updateVelocity();

    //SPH Functions

    double calculateDensity(std::vector<Particle> particles, Particle consideredParticle);

    double calculatePressure(std::vector<Particle> particles, Particle consideredParticle);

    std::array<double,3> calculatePressureForce(std::vector<Particle> particles, double coreRadius,Particle consideredParticle);


    

};
#endif