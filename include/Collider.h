#ifndef COLLIDER_H
#define COLLIDER_H
#include <array>
#include <Particle.h>

class Collider{
    private:

    const float3 minBounds;
    const float3 maxBounds;

    const double FRICTION = 0.1;
    const double RESTITUTION = 0.9;

    public:

    Collider(const float3& minBounds,const float3& maxBounds);


    void resolveSphereAABBCollision(Particle& particle)const;
    void resolveSphereCollision(Particle& p1, Particle& p2)const;

};

#endif