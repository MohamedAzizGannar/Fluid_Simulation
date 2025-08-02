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
        std::array<double,3> newV1 = {0,0,0};
        std::array<double,3> newV2 = {0,0,0};

        std::array<double,3> differenceOfVect = substractArrays(vel2,vel1);
        double dotProduct = dotProductofVect(differenceOfVect,vect);
        double coef = - dotProduct/distanceSquared;

        newV1 = addArrays(vel1,scaleArray(vect,coef));
        newV2 = addArrays(vel2,scaleArray(vect,coef));

        p1.setVelocity(newV1);
        p2.setVelocity(newV2);

    }

}