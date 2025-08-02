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

    const double mass = 2600.0;

    double density;
    double pressure;
    double viscosity;

    public:
    //Constructor
    Particle(std::array<double,3> initialPosition,std::array<double,3> initialVelocity,std::array<double,3> intialAcceleration);

    //Default Constructor (Array Initialization)
    Particle();

    //Getters
    std::array<double,3> getPosition() const;
    std::array<double,3> getAcceleration()const;
    std::array<double,3> getVelocity()const;
    double getMass()const;
    double getRadius()const;
    double getPressure()const;
    double getViscosity()const;
    double getDensity()const;


    //Setters
    void setPosition(double newPosition[3]);
    void setVelocity(double newVelocity[3]);
    void setAcceleration(double newAcceleration[3]);

    //Update Position 
    void updatePosition();

    //Update Velocity
    void updateVelocity();

    //SPH Functions

    double calculateDensity(const std::vector<Particle>& particles, double coreRadius);

    double calculatePressure(const std::vector<Particle>& particles, double coreRadius);

    std::array<double,3> calculatePressureForce(const std::vector<Particle>& particles, double coreRadius);

    std::array<double,3> calculateViscosityForce(std::vector<Particle>& particles, double coreRadius);
};
#endif