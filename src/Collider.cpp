#include <Collider.h>
#include <iostream>
#include <array>
#include <cmath>

double dotProductofVect(std::array<double,3> V1, std::array<double,3> V2){
    return V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2] ;
}
std::array<double,3> substractArrays(std::array<double,3> v1, std::array<double,3> v2){
    return {v1[0]-v2[0],v1[1]-v2[1],v1[2]-v2[2]};
}
std::array<double,3> addArrays(std::array<double,3> v1, std::array<double,3> v2){
    return {v1[0]+v2[0],v1[1]+v2[1],v1[2]+v2[2]};
}
std::array<double,3> scaleArray(std::array<double,3> v1,double a){
    return {v1[0]*a,v1[1]*a,v1[2]*a};
}


Collider::Collider(const std::array<double,3>& minBounds,const std::array<double,3>& maxBounds):maxBounds(maxBounds),minBounds(minBounds){}

void Collider::resolveSphereAABBCollision(Particle& particle)const {

    auto predictedPosition = particle.getPredictedPosition();
    auto velocity = particle.getVelocity();
    const double radius = particle.getRadius();

    bool collisionOccured = false;

    for(int i = 0; i < 3; i ++){
        double minEdge = minBounds[i] + radius;
        double maxEdge = maxBounds[i] - radius;
        if(predictedPosition[i] < minEdge){
            predictedPosition[i] = minEdge;
            if(velocity[i] < 0)
            {
                velocity[i] = - velocity[i] * RESTITUTION;
                velocity[i] *= (1.0  - FRICTION);
            }
            collisionOccured = true;

        }
        else if (predictedPosition[i] > maxEdge){
            predictedPosition[i] = maxEdge;
            if(velocity[i] > 0){
                velocity[i] = -velocity[i] * RESTITUTION;
                velocity[i] *= (1.0 - FRICTION);
            }
            collisionOccured = true;
        }
    }
    if(collisionOccured){
        particle.setPredictedPosition(predictedPosition);
        particle.setVelocity(velocity);
    }
}

void Collider::resolveSphereCollision(Particle& p1, Particle& p2)const {
    const double radiusSum = p1.getRadius() + p2.getRadius();
    const auto pos1 = p1.getPosition();
    const auto pos2 = p2.getPosition();
    
    // Calculate squared distance
    double distSq = 0;
    std::array<double,3> normal{0,0,0};
    for(int i=0; i<3; i++) {
        normal[i] = pos1[i] - pos2[i];
        distSq += normal[i]*normal[i];
    }
    
    // Check collision
    if(distSq >= radiusSum*radiusSum) return;
    
    // Normalize collision normal
    const double dist = sqrt(distSq);
    for(int i=0; i<3; i++) normal[i] /= dist;

    // Calculate relative velocity
    const auto vel1 = p1.getVelocity();
    const auto vel2 = p2.getVelocity();
    double velAlongNormal = 0;
    for(int i=0; i<3; i++) 
        velAlongNormal += (vel1[i]-vel2[i])*normal[i];

    // Only resolve if moving toward each other
    if(velAlongNormal > 0) return;

    // Calculate impulse 
    const double restitution = 0.9; // Elasticity coefficient
    double impulse = -(1 + restitution) * velAlongNormal / 2.0;

    // Apply impulse
    std::array<double,3> newVel1 = vel1;
    std::array<double,3> newVel2 = vel2;
    for(int i=0; i<3; i++) {
        newVel1[i] += impulse * normal[i];
        newVel2[i] -= impulse * normal[i];
    }

    // Position correction to prevent sticking
    const double overlap = (radiusSum - dist) * 0.5;
    std::array<double,3> correctedPos1 = pos1;
    std::array<double,3> correctedPos2 = pos2;
    for(int i=0; i<3; i++) {
        correctedPos1[i] += overlap * normal[i];
        correctedPos2[i] -= overlap * normal[i];
    }

    // Commit changes
    p1.setVelocity(newVel1);
    p2.setVelocity(newVel2);
    p1.setPosition(correctedPos1);
    p2.setPosition(correctedPos2);
}