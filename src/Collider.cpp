#include <Collider.h>
#include <iostream>
#include <array>


Collider::Collider(const std::array<double,3>& minBounds,const std::array<double,3>& maxBounds):maxBounds(maxBounds),minBounds(minBounds){}

bool Collider::resolveSphereAABBCollision(Particle& particle)const {

    auto position = particle.getPosition();
    auto velocity = particle.getVelocity();
    const double radius = particle.getRadius();

    for(int i = 0; i < 3; i ++){
        double minEdge = minBounds[i] + radius;
        double maxEdge = maxBounds[i] - radius;

        if(position[i] < minEdge){
            position[i] = minEdge;
            if(velocity[i] < 0)
            {
                velocity[i] = - velocity[i] * RESTITUTION;
                velocity[i] *= (1.0  - FRICTION);
            }
            particle.setPosition(position);
            particle.setVelocity(velocity);
            return true;
            
        }
        else if (position[i] > maxEdge){
            position[i] = maxEdge;
            if(velocity[i] > 0){
                velocity[i] = -velocity[i] * RESTITUTION;
                velocity[i] *= (1.0 - FRICTION);
            }
            particle.setPosition(position);
            particle.setVelocity(velocity);
            return true;
        }
    }
    return false;
}
