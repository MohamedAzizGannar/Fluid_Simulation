#include <iostream>
#include <iomanip>

#include <cmath>

#include <array>
#include <string>

#include<chrono>
#include<thread>

#include <random>

#include <Particle.h>
#include <Collider.h>
#include <FluidSimulation.h>
#include <Grid.h>


FluidSimulation::FluidSimulation(const std::vector<Particle>& particles,const Collider& worldBounds): particles(particles),worldBounds(worldBounds){
    grid = Grid(particles[0].getRadius());
}
Grid FluidSimulation::getGrid(){return grid;}
void FluidSimulation::update(double deltaTime)
{
    deltaTime = std::min(deltaTime,maxAccumulator);
    accumulator += deltaTime;
    while(accumulator >= fixedTimeStep){
        physicsOperations(fixedTimeStep);
        accumulator -= fixedTimeStep;
    }
} 
void FluidSimulation::updateDensityAndPressure(){
    for(auto& particle : particles){
        std::vector<int> neighborIDs = getGrid().getNeighbors(particle.getPosition());
        std::vector<Particle> effectiveParticles = {};
        for(int neighborID : neighborIDs){effectiveParticles.push_back(particles[neighborID]);}
        particle.updateDensity(effectiveParticles);
    }
    
    for(auto& particle : particles){

        std::vector<int> neighborIDs = getGrid().getNeighbors(particle.getPosition());
        std::vector<Particle> effectiveParticles = {};
        for(int neighborID : neighborIDs){effectiveParticles.push_back(particles[neighborID]);}
        particle.updatePressure(effectiveParticles);
    }
}
void FluidSimulation::applyForces(){
    updateDensityAndPressure();
    for(auto& particle:particles){

        std::vector<int> neighborIDs = getGrid().getNeighbors(particle.getPosition());
        std::vector<Particle> effectiveParticles = {};
        for(int neighborID : neighborIDs){effectiveParticles.push_back(particles[neighborID]);}
        particle.applyForcesOptimised(effectiveParticles);
    }
}

void FluidSimulation::integrateVel(double dt){
    for (auto& particle : particles){
        auto velocity = particle.getVelocity();
        const auto& acceleration = particle.getAcceleration();
        
        velocity += acceleration * dt;
        
        particle.setVelocity(velocity);
    }
}
inline bool FluidSimulation::isValidVector(const float3& vec) const {
    return std::isfinite(vec.x) && std::isfinite(vec.y) && std::isfinite(vec.z);
}
void FluidSimulation::predictPositions(double dt){
    for(auto& particle : particles){
        const auto& position = particle.getPosition();
        const auto& velocity = particle.getVelocity();

        if ( !isValidVector(position) || !isValidVector(velocity)){
            particle.setPosition({0.0, 0.0, 0.0});
            particle.setVelocity({0.0, 0.0, 0.0});
            continue;
        }

        

        float3 predictedPosition = position + velocity * dt;
       
        particle.setPredictedPosition(predictedPosition);
    }
}
void FluidSimulation::resolveCollisions(double dt){
    const uint8_t collisionPasses = 2;

    for(int pass = 0;pass < collisionPasses; pass++){
        for(auto& particle : particles){
            worldBounds.resolveSphereAABBCollision(particle);
            
        }
    }
}
void FluidSimulation::updatePositions(){
    for(auto& particle : particles){
        particle.setPosition(particle.getPredictedPosition());
    }
}
void FluidSimulation::applyDamping(){
    for (auto& particle : particles) {
        auto velocity = particle.getVelocity();
        velocity *= DAMPING;
        particle.setVelocity(velocity);
    }
}
void FluidSimulation::physicsOperations(double dt){

    applyForces();
    integrateVel(dt);
    predictPositions(dt);
    resolveCollisions(dt);
    updatePositions();
}

