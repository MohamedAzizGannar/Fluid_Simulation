#include <Particle.h>
#include <iostream>
#include <array>
#include <vector>
#include <cmath>

const double TENSION_COEFFICIENT = 7.28;
const double GAS_CONSTANT = 300.0;
const double REST_DENSITY = 1000.0;
const double GRAVITATIONAL_CONSTANT = 6.674e-11;
const double CORE_RADIUS = 1.0;
const double CORE_RADIUS9 = std::pow(CORE_RADIUS, 9);
const double CORE_RADIUS6 = std::pow(CORE_RADIUS,6);
const double CORE_RADIUS2 = std::pow(CORE_RADIUS, 2);

const double MU = 0.00089;


const double WPOLY6_COEFF = 315.0 / (64.0 * M_PI * CORE_RADIUS9);
const double GRAD_WPOLY6_COEFF = -945.0 / (32.0 * M_PI * CORE_RADIUS9);
const double GRAD_WSPIKY_COEFF = -45.0 / (M_PI * CORE_RADIUS6);
const double LAPLACIAN_VISC_COEFF = 45.0 / (M_PI * CORE_RADIUS6);


double Wpoly6(double distanceSq){
    if( distanceSq > CORE_RADIUS)return 0.0;
    const double coefficient = 315.0/(64.0*M_PI*CORE_RADIUS9);
    const double diff = CORE_RADIUS2 - distanceSq;
    return coefficient * diff * diff * diff;
}
float3 gradientWpoly6(const float3& vect){
    const double distanceSqrd = vect.lengthSQR();
    if(distanceSqrd > CORE_RADIUS2) return {0.0,0.0,0.0};
    const double secondMember = CORE_RADIUS2 -distanceSqrd;
    return vect * GRAD_WPOLY6_COEFF*secondMember*secondMember;
}
double magnitudeGradientWpoly6(const float3& vect){
    const double distanceSqrd = vect.lengthSQR();
    const double distance = vect.length();

    if(distanceSqrd > CORE_RADIUS2)return 0.0;
    const double coefficient = - 945.0 / (32.0 * M_PI * CORE_RADIUS9);
    const double secondMember = std::pow(CORE_RADIUS2 - distanceSqrd , 2);
    return coefficient * secondMember * distance;
}
double  laplacianWpoly6(const float3& vect){
    const double distanceSqrd = vect.lengthSQR();
    if(distanceSqrd > CORE_RADIUS2) return 0.0;
    const double secondMember = (CORE_RADIUS2 - distanceSqrd);
    const double thirdMember = (3*distanceSqrd - CORE_RADIUS2) ;
    return GRAD_WPOLY6_COEFF * secondMember * thirdMember;
}
float3 gradientWspiky( const float3& vect){
    const double distanceSqrd = vect.lengthSQR();
    const double distance = vect.length();

    if(distanceSqrd > CORE_RADIUS2 || distance < 1e-12){return {0.0,0.0,0.0};}    
    const double secondMember = CORE_RADIUS - distance;
    const double coefficient = GRAD_WSPIKY_COEFF * 3 * secondMember * secondMember / distance;
    const double factor = coefficient * secondMember;
    return vect*factor;
}
double laplacianViscosityKernel(const double& distance){
    
    
    if(distance >= CORE_RADIUS2 || distance < 0) return 0.0;
    return LAPLACIAN_VISC_COEFF * (CORE_RADIUS - distance);
}

Particle::Particle(float3 initialPosition,float3 initialVelocity,float3 intialAcceleration, int id):id(id){

    position = float3(initialPosition.x,initialPosition.y,initialPosition.z);
    velocity = float3(initialVelocity.x,initialVelocity.y,initialVelocity.z);
    acceleration = float3(intialAcceleration.x,intialAcceleration.y,intialAcceleration.z);

}

void Particle::setPosition(float3  newPosition){
    position = newPosition;
}
void Particle::setPredictedPosition(float3  newPredictedPosition){
    predictedPosition = newPredictedPosition;
}
void Particle::setVelocity(float3  newVelocity){
    velocity = newVelocity;
}
void Particle::setAcceleration(float3  newAcceleration){
    acceleration = newAcceleration;
}
void Particle::setDensity(double newDensity){density = newDensity;}
void Particle::setPressure(double newPressure){pressure = newPressure;}
void Particle::setId(int i){id = i;}



float3 Particle::getPosition()const{
    return position;
}
float3 Particle::getPredictedPosition()const{
    return predictedPosition;
}

float3 Particle::getAcceleration()const{
    return acceleration;
}
float3 Particle::getVelocity()const{
    return velocity;
}
double Particle::getMass()const{
    return mass;
}
double Particle::getRadius()const{
    return radius;
}
double Particle::getPressure()const{return pressure;}
double Particle::getDensity()const{return density;}
int Particle::getId()const{return id;}



void Particle::updatePosition(){
    position += velocity;
}

void Particle::updateVelocity(){
    velocity += acceleration;
}
void Particle::updateDensity(const std::vector<Particle>& particles){
    setDensity(calculateDensity(particles));
}

void Particle::updatePressure(const std::vector<Particle>& particles){
    setPressure(calculatePressure(particles));
}

//SPH FUNCTIONS
double Particle::calculateDensity(const std::vector<Particle>& particles)
{
    double density = 0.0;
    const auto& myPos = getPosition();

    for( const Particle& neighbor : particles){
        if (neighbor.getId( ) == this->getId()) continue;
        const auto& neighborPos = neighbor.getPosition();
        float3 distanceVector= myPos - neighborPos;
        double distanceSqrd = distanceVector.lengthSQR();
        if(distanceSqrd>CORE_RADIUS)continue;
        density+= neighbor.getMass() * Wpoly6(distanceSqrd);
    }
    return std::max(density,1e-12);
}
double Particle::calculatePressure(const std::vector<Particle>& particles){

    double pressure = GAS_CONSTANT * (calculateDensity(particles) - REST_DENSITY);
    return std::max(0.0,pressure);

}
float3 Particle::calculateGravity(){
    return {0.0,-9.81,0.0};
}


double calculateCurvature(double laplacian,const float3& gradient){
    const double magnitudeSq = gradient.lengthSQR();
    if(magnitudeSq < 1e-12)return 0.0;
    const double K = - laplacian / sqrt(magnitudeSq);
    return K;
}


ColorFieldProperties Particle::calculateColorFieldProperties(const std::vector<Particle>& particles) {
    ColorFieldProperties props = {0.0, {0.0, 0.0, 0.0}, 0.0};
    const auto& myPos = getPosition();
    
    for (const auto& neighbor : particles) {
        if (neighbor.getId() == id) continue;
        
        const auto& neighborPos = neighbor.getPosition();
        const float3 vect = myPos - neighborPos;
        
        const double distanceSq = vect.lengthSQR();
        if (distanceSq >= CORE_RADIUS2) continue;
        
        const double densityJ = neighbor.getDensity();
        if (densityJ < 1e-6) continue;
        
        const double coefficient = neighbor.getMass() / densityJ;
        
        const double poly6Value = Wpoly6(distanceSq);
        const float3 gradientValue = gradientWpoly6(vect);
        const double laplacianValue = laplacianWpoly6(vect);
        
        props.field += poly6Value * coefficient;
        props.gradient += gradientValue*coefficient;
        props.laplacian += laplacianValue * coefficient;
    }
    
    return props;
}
void Particle::applyForcesOptimised(const std::vector<Particle>& particles){
    float3 pressureForce = float3();
    float3 viscosityForce = float3();
    ColorFieldProperties colorProps = {0.0,float3(),0.0};

    const auto& myPos = getPosition();
    const auto& myVel = getVelocity();
    const double myPressure = getPressure();
    const double myDensity = getDensity();

    if(myDensity<1e-6){
        setAcceleration({0.0,0.0,0.0});
        return;
    }

    for(const auto& neighbor : particles){
        if(neighbor.getId() == getId())continue;

        const auto& neighborPos = neighbor.getPosition();
        const float3 vect = myPos - neighborPos;
        const double distanceSqrd = vect.lengthSQR();

        if(distanceSqrd >= CORE_RADIUS2) continue;

        const double distance = sqrt(distanceSqrd);
        const double neighborPressure = neighbor.getPressure();
        const double neighborDensity = neighbor.getDensity();
        const double neighborMass = neighbor.getMass();

        if(density < 1e-16)continue;

        const double coefficient = neighborMass / neighborDensity;

        //PRESSURE FORCE CALCULATION
        if(distance < 0.0001){
            const auto gradW = gradientWspiky(vect);
            const double avgPressure = (myPressure + neighborPressure ) * 0.5;
            const double pressureCoefficient = -avgPressure * coefficient;
            pressureForce += gradW*avgPressure;
        }

        //VISCOSITY FORCE CALCULATION
        const auto& neighborVel = neighbor.getVelocity();

        const float3 velocityDifference = myVel - neighborVel;

        const double laplacienViscosity = laplacianViscosityKernel(distance);
        const double viscosityCoefficient = MU * neighborMass * laplacienViscosity / neighborDensity;

        viscosityForce = viscosityForce + velocityDifference*viscosityCoefficient;


        const double poly6Value = Wpoly6(distanceSqrd);
        const auto gradientValue = gradientWpoly6(vect);
        const double laplacianValue = laplacianWpoly6(vect);

        colorProps.field += poly6Value * coefficient;
        colorProps.gradient = colorProps.gradient + gradientValue*coefficient;
        colorProps.laplacian += laplacianValue * coefficient;

        const double curvature = calculateCurvature(colorProps.laplacian, colorProps.gradient);
        const double tensionCoefficient = -TENSION_COEFFICIENT * curvature;
        const float3 surfaceTensionForce = colorProps.gradient*tensionCoefficient;

        const auto gravity = calculateGravity();


        const float3 totalForce = pressureForce + viscosityForce + surfaceTensionForce + gravity;

        const double invDensity = 1.0 / myDensity;
        setAcceleration(totalForce*invDensity);
    }

}

