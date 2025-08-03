#ifndef FLUIDSIMULATION_H
#define FLUIDSIMULATION_H

#include <iostream>
#include <iomanip>

#include <cmath>

#include <array>
#include <string>

#include<chrono>
#include<thread>

#include <random>
#include <vector>

#include <Particle.h>
#include <Collider.h>
class FluidSimulation{
    private:
    std::vector<Particle> particles;
    Collider worldBounds;


    const double fixedTimeStep = 1.0 / 120.0;
    double maxAccumulator = 0.25;
    double accumulator = 0.0;

    const double DAMPING = 0.99;

    public:
    FluidSimulation(const std::vector<Particle>& particles,const Collider& worldBounds);
    void update(double deltaTime);
    inline bool isValidVector(const std::array<double, 3>& vec) const;

    std::vector<Particle> getParticles(){return particles;}
    private:
    void updateDensityAndPressure();

    void applyForces();
    void integrateVel(double dt);
    void predictPositions(double dt);
    void resolveCollisions(double dt);
    void updatePositions();
    void physicsOperations(double dt);
    void applyDamping();


 
};

#endif