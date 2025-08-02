#include <Collider.h>
#include <iostream>
#include <array>


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

void Collider::resolveSphereCollision(Particle& p1, Particle& p2)const{
    auto pos1 = p1.getPredictedPosition();
    auto pos2 = p2.getPredictedPosition();

    auto vel1 = p1.getVelocity();
    auto vel2 = p2.getVelocity();


    double radius1 = p1.getRadius();
    double radius2 = p2.getRadius();


    std::array<double,3> vect = {0.0,0.0,0.0};
    double distanceSquared = 0;
    for(int i = 0; i < 3 ; i ++){
        vect[i] = pos1[i] - pos2[i];
        distanceSquared += vect[i] * vect[i];
    }
    double distance = std::sqrt(distanceSquared);

    if(distance < radius1 + radius2){

    }

}