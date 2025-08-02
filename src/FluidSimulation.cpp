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


FluidSimulation::FluidSimulation(const std::vector<Particle>& particles,const Collider& worldBounds): particles(particles),worldBounds(worldBounds){}
void FluidSimulation::update(double deltaTime)
{
    deltaTime = std::min(deltaTime,maxAccumulator);
    accumulator += deltaTime;
    while(accumulator >= fixedTimeStep){
        physicsOperations(fixedTimeStep);
        accumulator -= fixedTimeStep;
    }

    for(auto& particle:particles){
        auto [x,y,z] = particle.getPosition();
        std::cout<<"Position : ("<<x<<","<<y<<","<<z<<")"<<std::endl;
    }
    std::cout<<std::endl;
} 
void FluidSimulation::updateDensityAndPressure(){
    for(auto& particle : particles){
        particle.updateDensity(particles);
    }
    
    for(auto& particle : particles){
        particle.updatePressure(particles);
    }
}
void FluidSimulation::applyForces(){
    updateDensityAndPressure();



    for(auto& particle:particles){
        particle.applyForces(particles);
    }
}

void FluidSimulation::integrateVel(double dt){
    for (auto& particle : particles){
        auto velocity = particle.getVelocity();
        auto acceleration = particle.getAcceleration();
        
        for(int i = 0; i < 3; i ++){
            velocity[i] += acceleration[i] * dt;
        }
        particle.setVelocity(velocity);
    }
}
void FluidSimulation::predictPositions(double dt){
    for(auto& particle : particles){
        auto position = particle.getPosition();
        auto velocity = particle.getVelocity();

      

        for(int i = 0; i < 3; i++) {
            if(!std::isfinite(position[i]) || !std::isfinite(velocity[i])) {
                
                return; 
            }
        }

        auto prevPosition = particle.getPosition();
        std::array<double,3> predictedPosition = {
            position[0] + velocity[0] * dt,
            position[1] + velocity[1]* dt,
            position[2] + velocity[2]* dt
        };
       
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
        velocity[0] *= DAMPING;
        velocity[1] *= DAMPING;
        velocity[2] *= DAMPING;
        particle.setVelocity(velocity);
    }
}
void FluidSimulation::physicsOperations(double dt){

    applyForces();
    integrateVel(dt);
    predictPositions(dt);
    resolveCollisions(dt);
    updatePositions();
    applyDamping();
}

