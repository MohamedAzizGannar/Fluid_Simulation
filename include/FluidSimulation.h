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
#include <Grid.h>   

class FluidSimulation{
    private:
    std::vector<Particle> particles;
    Collider worldBounds;
    Grid grid;

    float targetPhysicsRate;
    float fixedTimeStep;
    float timeAccumulator; 
    int maxPhysicsStepsPerFrame;
    std::chrono::high_resolution_clock::time_point lastFrameTime;


    const double DAMPING = 0.99;

    public:
    FluidSimulation(const std::vector<Particle>& particles,const Collider& worldBounds);
    inline bool isValidVector(const float3& vec) const;
    std::vector<Particle> getParticles(){return particles;}
    void setTargetPhysicsRate(float hz);
    void update();

    Grid getGrid();
    private:

    void updateDensityAndPressure();
    void applyForces();
    void integrateVel(double dt);
    void predictPositions(double dt);
    void resolveCollisions(double dt);
    void updatePositions();
    void applyDamping();

    void updatePhysics(float dt);


 
};

#endif