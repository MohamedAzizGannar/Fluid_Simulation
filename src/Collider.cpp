#include <Collider.h>
#include <iostream>
#include <array>
#include <cmath>

inline double dotProduct(const std::array<double,3>& v1, const std::array<double,3>& v2) {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

inline std::array<double,3> subtract(const std::array<double,3>& v1, const std::array<double,3>& v2) {
    return {v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]};
}

inline std::array<double,3> add(const std::array<double,3>& v1, const std::array<double,3>& v2) {
    return {v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2]};
}

inline std::array<double,3> scale(const std::array<double,3>& v, double scalar) {
    return {v[0]*scalar, v[1]*scalar, v[2]*scalar};
}

inline double lengthSquared(const std::array<double,3>& v) {
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

inline double length(const std::array<double,3>& v) {
    return std::sqrt(lengthSquared(v));
}


Collider::Collider(const std::array<double,3>& minBounds,const std::array<double,3>& maxBounds):maxBounds(maxBounds),minBounds(minBounds){}

void Collider::resolveSphereAABBCollision(Particle& particle)const {

    auto predictedPosition = particle.getPredictedPosition();
    auto velocity = particle.getVelocity();
    const double radius = particle.getRadius();

    const double dampingFactor = RESTITUTION * (1.0 - FRICTION);
    
    bool collisionOccurred = false;

    const std::array<double,3> minEdges = {
        minBounds[0] + radius,
        minBounds[1] + radius,
        minBounds[2] + radius
    };
    
    const std::array<double,3> maxEdges = {
        maxBounds[0] - radius,
        maxBounds[1] - radius,
        maxBounds[2] - radius
    };

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
    
    const std::array<double,3> collision = subtract(pos1, pos2);
    const double distSq = lengthSquared(collision);
    const double radiusSumSq = radiusSum * radiusSum;
    
    if (distSq >= radiusSumSq || distSq < 1e-12) return;
    
    const double dist = std::sqrt(distSq);
    const double invDist = 1.0 / dist;
    const std::array<double,3> normal = scale(collision, invDist);

    const auto& vel1 = p1.getVelocity();
    const auto& vel2 = p2.getVelocity();
    const std::array<double,3> relativeVelocity = subtract(vel1, vel2);
    const double velAlongNormal = dotProduct(relativeVelocity, normal);

    if (velAlongNormal > 0) return;

    constexpr double restitution = 0.9;
    const double impulseMagnitude = -(1.0 + restitution) * velAlongNormal * 0.5;
    const std::array<double,3> impulse = scale(normal, impulseMagnitude);

    const std::array<double,3> newVel1 = add(vel1, impulse);
    const std::array<double,3> newVel2 = subtract(vel2, impulse);

    const double overlap = radiusSum - dist;
    const double correctionMagnitude = overlap * 0.5;
    const std::array<double,3> correction = scale(normal, correctionMagnitude);
    
    const std::array<double,3> correctedPos1 = add(pos1, correction);
    const std::array<double,3> correctedPos2 = subtract(pos2, correction);

    p1.setVelocity(newVel1);
    p2.setVelocity(newVel2);
    p1.setPosition(correctedPos1);
    p2.setPosition(correctedPos2);
}