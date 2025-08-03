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


double findVector3LengthSquared(const std::array<double,3>& vect){
    auto [x,y,z] = vect;
    return x*x + y*y + z*z;
}
double findVector3Length(const std::array<double,3>& vect){
    auto [x,y,z] = vect;
    return sqrt(x*x + y*y + z*z);
}
std::array<double,3> normalizeVector3(const std::array<double,3>& vect){
    auto [x,y,z] = vect;
    double lenSqr = findVector3LengthSquared(vect);
    if (lenSqr == 0 )return{0,0,0};
    const double invLen = 1.0 / std::sqrt(lenSqr);
    return {
        invLen * vect[0],
        invLen * vect[1],
        invLen * vect[2],
    };
}
double Wpoly6(double distanceSq){
    if( distanceSq > CORE_RADIUS)return 0.0;
    const double coefficient = 315.0/(64.0*M_PI*CORE_RADIUS9);
    const double diff = CORE_RADIUS2 - distanceSq;
    return coefficient * diff * diff * diff;
}
std::array<double,3> gradientWpoly6(const std::array<double,3>& vect){
    const double distanceSqrd = findVector3LengthSquared(vect);
    if(distanceSqrd > CORE_RADIUS2) return {0.0,0.0,0.0};
    const double secondMember = CORE_RADIUS2 -distanceSqrd;
    return {
        GRAD_WPOLY6_COEFF*secondMember*secondMember*vect[0],
        GRAD_WPOLY6_COEFF*secondMember*secondMember*vect[1],
        GRAD_WPOLY6_COEFF*secondMember*secondMember*vect[2]
    };
}
double magnitudeGradientWpoly6(const std::array<double,3>& vect){
    const double distanceSqrd = findVector3LengthSquared(vect);
    const double distance = findVector3Length(vect);

    if(distanceSqrd > CORE_RADIUS2)return 0.0;
    const double coefficient = - 945.0 / (32.0 * M_PI * CORE_RADIUS9);
    const double secondMember = std::pow(CORE_RADIUS2 - distanceSqrd , 2);
    return coefficient * secondMember * distance;
}
double  laplacianWpoly6(const std::array<double,3>& vect){
    const double distanceSqrd = findVector3LengthSquared(vect);
    if(distanceSqrd > CORE_RADIUS2) return 0.0;
    const double secondMember = (CORE_RADIUS2 - distanceSqrd);
    const double thirdMember = (3*distanceSqrd - CORE_RADIUS2) ;
    return GRAD_WPOLY6_COEFF * secondMember * thirdMember;
}
std::array<double,3> gradientWspiky( const std::array<double,3>& vect){
    const double distanceSqrd = findVector3LengthSquared(vect);
    const double distance = findVector3Length(vect);

    if(distanceSqrd > CORE_RADIUS2 || distance < 1e-12){return {0.0,0.0,0.0};}    
    const double secondMember = CORE_RADIUS - distance;
    const double coefficient = GRAD_WSPIKY_COEFF * 3 * secondMember * secondMember / distance;
    
    return {
        coefficient * secondMember * vect[0],
        coefficient * secondMember * vect[1],
        coefficient * secondMember * vect[2]
    };
}
double laplacianViscosityKernel(const double& distance){
    
    
    if(distance >= CORE_RADIUS2 || distance < 0) return 0.0;
    return LAPLACIAN_VISC_COEFF * (CORE_RADIUS - distance);
}

Particle::Particle(std::array<double,3> initialPosition,std::array<double,3> initialVelocity,std::array<double,3> intialAcceleration, int id):id(id){
    for(int i =0; i < 3; i++){
        position[i] = initialPosition[i];
        velocity[i] = initialVelocity[i];
        acceleration[i] = intialAcceleration[i];
    }
}

void Particle::setPosition(std::array<double,3>  newPosition){
    position[0] = newPosition[0];
    position[1] = newPosition[1];
    position[2] = newPosition[2];
}
void Particle::setPredictedPosition(std::array<double,3>  newPredictedPosition){
    predictedPosition[0] = newPredictedPosition[0];
    predictedPosition[1] = newPredictedPosition[1];
    predictedPosition[2] = newPredictedPosition[2];
}
void Particle::setVelocity(std::array<double,3>  newVelocity){
    velocity[0] = newVelocity[0];
    velocity[1] = newVelocity[1];
    velocity[2] = newVelocity[2];
}
void Particle::setAcceleration(std::array<double,3>  newAcceleration){
    acceleration[0] = newAcceleration[0];
    acceleration[1] = newAcceleration[1];
    acceleration[2] = newAcceleration[2];
}
void Particle::setDensity(double newDensity){density = newDensity;}
void Particle::setPressure(double newPressure){pressure = newPressure;}
void Particle::setId(int i){id = i;}



std::array<double,3> Particle::getPosition()const{
    return position;
}
std::array<double,3> Particle::getPredictedPosition()const{
    return predictedPosition;
}

std::array<double,3> Particle::getAcceleration()const{
    return acceleration;
}
std::array<double,3> Particle::getVelocity()const{
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
    position[0] += velocity[0];
    position[1] += velocity[1];
    position[2] += velocity[2];
}

void Particle::updateVelocity(){
    velocity[0] += acceleration[0];
    velocity[1] += acceleration[1];
    velocity[2] += acceleration[2];
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
    const auto& pi = getPosition();

    for( const Particle& neighbor : particles){
        if (neighbor.getId( ) == this->getId()) continue;
        const auto& pj = neighbor.getPosition();
        std::array<double,3> distanceVector={
            pi[0] - pj[0],
            pi[1] - pj[1],
            pi[2] - pj[2]
        };
        
        double distanceSqrd = findVector3LengthSquared(distanceVector);
        if(distanceSqrd>CORE_RADIUS)continue;
        density+= neighbor.getMass() * Wpoly6(distanceSqrd);
    }
    return std::max(density,1e-12);
}
double Particle::calculatePressure(const std::vector<Particle>& particles){

    double pressure = GAS_CONSTANT * (calculateDensity(particles) - REST_DENSITY);
    return std::max(0.0,pressure);

}
std::array<double,3> Particle::calculateGravity(){
    return {0.0,-9.81,0.0};
}


double calculateCurvature(double laplacian,const std::array<double,3>& gradient){
    const double magnitudeSq = findVector3LengthSquared(gradient);
    if(magnitudeSq < 1e-12)return 0.0;
    const double K = - laplacian / sqrt(magnitudeSq);
    return K;
}


ColorFieldProperties Particle::calculateColorFieldProperties(const std::vector<Particle>& particles) {
    ColorFieldProperties props = {0.0, {0.0, 0.0, 0.0}, 0.0};
    const auto& pos = getPosition();
    
    for (const auto& pj : particles) {
        if (pj.getId() == id) continue;
        
        const auto& posj = pj.getPosition();
        const std::array<double,3> vect = {
            pos[0] - posj[0], 
            pos[1] - posj[1], 
            pos[2] - posj[2]
        };
        
        const double distanceSq = findVector3LengthSquared(vect);
        if (distanceSq >= CORE_RADIUS2) continue;
        
        const double densityJ = pj.getDensity();
        if (densityJ < 1e-6) continue;
        
        const double coefficient = pj.getMass() / densityJ;
        
        const double poly6Value = Wpoly6(distanceSq);
        const auto gradientValue = gradientWpoly6(vect);
        const double laplacianValue = laplacianWpoly6(vect);
        
        props.field += poly6Value * coefficient;
        props.gradient[0] += gradientValue[0] * coefficient;
        props.gradient[1] += gradientValue[1] * coefficient;
        props.gradient[2] += gradientValue[2] * coefficient;
        props.laplacian += laplacianValue * coefficient;
    }
    
    return props;
}
void Particle::applyForcesOptimised(const std::vector<Particle>& particles){
    std::array<double,3> pressureForce = {0.0,0.0,0.0};
    std::array<double,3> viscosityForce = {0.0,0.0,0.0};
    ColorFieldProperties colorProps = {0.0,{0.0,0.0,0.0},0.0};

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
        const std::array<double,3> vect = {
            myPos[0] - neighborPos[0],
            myPos[1] - neighborPos[1],
            myPos[2] - neighborPos[2]
        };
        const double distanceSqrd = findVector3LengthSquared(vect);

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

            pressureForce[0] += gradW[0] * pressureCoefficient;
            pressureForce[1] += gradW[1] * pressureCoefficient;
            pressureForce[2] += gradW[2] * pressureCoefficient;
        }

        //VISCOSITY FORCE CALCULATION
        const auto& neighborVel = neighbor.getVelocity();

        const std::array<double,3> velocityDifference = {
            myVel[0] - neighborVel[0],
            myVel[1] - neighborVel[1],
            myVel[2] - neighborVel[2]
        };

        const double laplacienViscosity = laplacianViscosityKernel(distance);
        const double viscosityCoefficient = MU * neighborMass * laplacienViscosity / neighborDensity;

        viscosityForce[0] += viscosityCoefficient * velocityDifference[0];
        viscosityForce[1] += viscosityCoefficient * velocityDifference[1];
        viscosityForce[2] += viscosityCoefficient * velocityDifference[2];


        const double poly6Value = Wpoly6(distanceSqrd);
        const auto gradientValue = gradientWpoly6(vect);
        const double laplacianValue = laplacianWpoly6(vect);

        colorProps.field += poly6Value * coefficient;
        colorProps.gradient[0] += gradientValue[0] * coefficient;
        colorProps.gradient[1] += gradientValue[1] * coefficient;
        colorProps.gradient[2] += gradientValue[2] * coefficient;
        colorProps.laplacian += laplacianValue * coefficient;

        const double curvature = calculateCurvature(colorProps.laplacian, colorProps.gradient);
        const double tensionCoefficient = -TENSION_COEFFICIENT * curvature;
        const std::array<double,3> surfaceTensionForce = {
            tensionCoefficient * colorProps.gradient[0],
            tensionCoefficient * colorProps.gradient[1],
            tensionCoefficient * colorProps.gradient[2]
        };

        const auto gravity = calculateGravity();


        const std::array<double,3> totalForce = {
            viscosityForce[0] + pressureForce[0] + surfaceTensionForce[0],
            viscosityForce[1] + pressureForce[1] + surfaceTensionForce[1]+gravity[1],
            viscosityForce[2] + pressureForce[2] + surfaceTensionForce[2],
        };

        const double invDensity = 1.0 / myDensity;
        setAcceleration({
            totalForce[0] * invDensity,
            totalForce[1] * invDensity,
            totalForce[2] * invDensity
        });
    }

}

