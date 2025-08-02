#ifndef COLLIDER_H
#define COLLIDER_H
#include <array>
#include <Particle.h>

class Collider{
    private:

    const std::array<double,3> minBounds;
    const std::array<double,3> maxBounds;

    const double FRICTION = 0.1;
    const double RESTITUTION = 1.0;

    public:

    Collider(const std::array<double,3>& minBounds,const std::array<double,3>& maxBounds);


    void resolveSphereAABBCollision(Particle& particle)const;

};

#endif