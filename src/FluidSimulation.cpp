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

inline bool FluidSimulation::isValidVector(const float3& vec) const {
    return std::isfinite(vec.x) && std::isfinite(vec.y) && std::isfinite(vec.z);
}
FluidSimulation::FluidSimulation(const std::vector<Particle>& particles,const Collider& worldBounds): 
                    particles(particles),
                    worldBounds(worldBounds),
                    targetPhysicsRate(60.0f), 
                    fixedTimeStep(1.0f / 60.0f),
                    timeAccumulator(0.0f),
                    maxPhysicsStepsPerFrame(1){
    grid = Grid(2.0);
}
Grid FluidSimulation::getGrid(){return grid;}
void FluidSimulation::setTargetPhysicsRate(float hz){
    targetPhysicsRate = hz;
    fixedTimeStep = 1/ targetPhysicsRate;
}
void FluidSimulation::update(){
    auto currentTime = std::chrono::high_resolution_clock::now();
    float deltaTime = std::chrono::duration<float>(currentTime - lastFrameTime).count();
    lastFrameTime = currentTime;
    deltaTime = std::min(deltaTime,0.25f);
    timeAccumulator += deltaTime;
    int physicsSteps = 0;
    while (timeAccumulator >= fixedTimeStep && physicsSteps < maxPhysicsStepsPerFrame) {
        updatePhysics(fixedTimeStep);
        timeAccumulator -= fixedTimeStep;
        physicsSteps++;
    }
    
}
void FluidSimulation::updateDensityAndPressure(){
    grid.rebuild(particles);
    std::vector<int> neighborIDs;

    for(auto& particle : particles){
        grid.getNeighbors(particle.getPosition(),neighborIDs);
        std::vector<Particle> effectiveParticles = {};
        for(int neighborID : neighborIDs){effectiveParticles.push_back(particles[neighborID]);}
        particle.updateDensity(effectiveParticles);
        particle.updatePressure(effectiveParticles);
    }
}
void FluidSimulation::applyForces(){
    updateDensityAndPressure();
    std::vector<int> reusableNeighborsId;
    reusableNeighborsId.reserve(50);
    for(auto& particle:particles){
        std::vector<Particle> effectiveParticles = {};
        for(int neighborID : reusableNeighborsId){effectiveParticles.push_back(particles[neighborID]);}
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

void FluidSimulation::updatePhysics(float dt){
    applyForces();
    integrateVel(dt);
    predictPositions(dt);
    resolveCollisions(dt);
    updatePositions();

}