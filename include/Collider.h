#ifndef COLLIDER_H
#define COLLIDER_H
#include <array>
#include <Particle.h>

class Collider{
    private:

    const std::array<double,3> minBounds;
    const std::array<double,3> maxBounds;

    const double FRICTION = 0.2;
    const double RESTITUTION = 0.8;

    public:

    Collider(const std::array<double,3>& minBounds,const std::array<double,3>& maxBounds);


    void resolveSphereAABBCollision(Particle& particle)const;
    void resolveSphereCollision(Particle& p1, Particle& p2)const;

};

#endif