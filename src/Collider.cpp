
#include <Collider.h>
#include <iostream>
#include <array>
#include <cmath>



Collider::Collider(const float3& minBounds,const float3& maxBounds):maxBounds(maxBounds),minBounds(minBounds){}

void Collider::resolveSphereAABBCollision(Particle& particle)const {

    auto predictedPosition = particle.getPredictedPosition();
    auto velocity = particle.getVelocity();
    const double radius = particle.getRadius();

    const double dampingFactor = RESTITUTION * (1.0 - FRICTION);
    
    bool collisionOccurred = false;

    const float3 minEdges = minBounds+ radius;
    
    const float3 maxEdges = maxBounds - radius;

    for (int i = 0; i < 3; ++i) {
        if (predictedPosition[i] < minEdges[i]) {
            predictedPosition[i] = minEdges[i];
            if (velocity[i] < 0) {
                velocity[i] = -velocity[i] * dampingFactor;
            }
            collisionOccurred = true;
        } else if (predictedPosition[i] > maxEdges[i]) {
            predictedPosition[i] = maxEdges[i];
            if (velocity[i] > 0) {
                velocity[i] = -velocity[i] * dampingFactor;
            }
            collisionOccurred = true;
        }
    }

    if (collisionOccurred) {
        particle.setPredictedPosition(predictedPosition);
        particle.setVelocity(velocity);
    }
}

void Collider::resolveSphereCollision(Particle& p1, Particle& p2)const {
    const double radiusSum = p1.getRadius() + p2.getRadius();
    const auto& pos1 = p1.getPosition();
    const auto& pos2 = p2.getPosition();
    
    const float3 collision = pos1 - pos2;
    const double distSq = collision.lengthSQR();
    const double radiusSumSq = radiusSum * radiusSum;
    
    if (distSq >= radiusSumSq || distSq < 1e-12) return;
    
    const double dist = std::sqrt(distSq);
    const double invDist = 1.0 / dist;
    const float3 normal = collision*invDist;

    const auto& vel1 = p1.getVelocity();
    const auto& vel2 = p2.getVelocity();
    const float3 relativeVelocity = vel1 - vel2;
    const double velAlongNormal = relativeVelocity.dot(normal);

    if (velAlongNormal > 0) return;

    const double impulseMagnitude = -(1.0 + RESTITUTION) * velAlongNormal * 0.5;
    const float3 impulse = normal * impulseMagnitude;

    const float3 newVel1 = vel1 + impulse;
    const float3 newVel2 = vel2 - impulse;

    const double overlap = radiusSum - dist;
    const double correctionMagnitude = overlap * 0.5;
    const float3 correction = normal * correctionMagnitude;
    
    const float3 correctedPos1 = pos1 + correction;
    const float3 correctedPos2 = pos2 - correction;

    p1.setVelocity(newVel1);
    p2.setVelocity(newVel2);
    p1.setPosition(correctedPos1);
    p2.setPosition(correctedPos2);
}