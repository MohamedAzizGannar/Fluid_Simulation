#include <Particle.h>
#include <iostream>
#include <array>
#include <vector>
#include <cmath>


double findVector3Length(std::array<double,3> vect){
    auto [x,y,z] = vect;
    return std::sqrt(x*x + y*y + z*z);
}
std::array<double,3> normalizeVector3(std::array<double,3> vect){
    auto [x,y,z] = vect;
    double len = findVector3Length(vect);
    if (len == 0 )return{0,0,0};
    return {x/len,y/len,z/len};
}
double Wpoly6(double distance, double coreRadius)
{
    if(distance < 0.0 || distance > coreRadius)return 0.0;
    double k = 315.0/(64.0*M_PI*std::pow(coreRadius,9));
    double secondMember = std::pow((coreRadius*coreRadius - distance*distance),3);
    return k * secondMember;
}
double Wspike(double distance, double coreRadius)
{
    if(distance<0 || distance > coreRadius){return 0.0;}
    const double coefficient = -45.0/(M_PI * std::pow(coreRadius,6));
    return coefficient * std::pow(coreRadius-distance,2);
}
double laplacianViscosityKernel(double distance,double coreRadius)
{
    if(distance < 0 || distance > coreRadius) return 0.0;
    const double coefficient = 45.0/(M_PI * std::pow(coreRadius,6));
    return coefficient * (coreRadius - distance);
}

Particle::Particle(std::array<double,3> initialPosition,std::array<double,3> initialVelocity,std::array<double,3> intialAcceleration){
    for(int i =0; i < 3; i++){
        position[i] = initialPosition[i];
        velocity[i] = initialVelocity[i];
        acceleration[i] = intialAcceleration[i];
    }
}
//Default Constructor
Particle::Particle()
{
    for(int i =0; i < 3; i++){
        position[i] = 0;
        velocity[i] = 0;
        acceleration[i] = 0;
    }
}

void Particle::setPosition(double newPosition[3]){
    position[0] = newPosition[0];
    position[1] = newPosition[1];
    position[2] = newPosition[2];

}
void Particle::setVelocity(double newVelocity[3]){
    velocity[0] = newVelocity[0];
    velocity[1] = newVelocity[1];
    velocity[2] = newVelocity[2];

}
void Particle::setAcceleration(double newAcceleration[3]){
    acceleration[0] = newAcceleration[0];
    acceleration[1] = newAcceleration[1];
    acceleration[2] = newAcceleration[2];

}


std::array<double,3> Particle::getPosition()const{
    return position;
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
double Particle::getViscosity()const{return viscosity;}
double Particle::getDensity()const{return density;}


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

//SPH FUNCTIONS
double Particle::calculateDensity(const std::vector<Particle>& particles, double coreRadius)
{
    double density = 0.0;
    for( const Particle& neighbor : particles){
        std::array<double,3> distanceVector;
        auto pi = this->getPosition();
        auto pj = neighbor.getPosition();
        for(int i =0 ; i < 3; i++)
        {
            distanceVector[i] = pi[i] - pj[i];
        }
        double distance = findVector3Length(distanceVector);
        if(distance>coreRadius)continue;
        density+= neighbor.getMass() * Wpoly6(distance,coreRadius);
    }
    this->density = density;
    return density < 1e-12 ? 1e-12:density;
}
double Particle::calculatePressure(const std::vector<Particle>& particles,double coreRadius){

    double GAS_CONSTANT = 1000.0;
    double REST_DENSITY = 1000.0;

    double pressure = GAS_CONSTANT * (this->calculateDensity(particles,coreRadius) - REST_DENSITY);
    pressure = std::max(0.0,pressure);
    this -> pressure = pressure;
    return pressure;
}
std::array<double,3> Particle::calculatePressureForce(const std::vector<Particle>& particles, double coreRadius){

    std::array<double,3> pressureForce = {0.0,0.0,0.0};
    double pressureI = this->getPressure();


    for(const Particle& neighbor: particles){
        if(&neighbor == this) continue;
        std::array<double,3> particlePosition = neighbor.getPosition();

        std::array<double,3> direction = {position[0]-particlePosition[0],position[1]-particlePosition[1],position[2]-particlePosition[2]};

        double distance = findVector3Length(direction);
        std::cout<<"Distance = "<< distance<<std::endl;
        if(distance<= 0.0001 || distance>=coreRadius) continue;

        std::array<double,3> normalizedDirection = normalizeVector3(direction);

        double particleMass = neighbor.getMass();

        std::array<double,3> gradW = {0.0,0.0,0.0};
        for(int i = 0; i < gradW.size(); i++)
        {
            gradW[i] = normalizedDirection[i] * Wspike(distance, coreRadius);
            std::cout<<"GradW at index"<<i<<" : "<<gradW[i]<<"\n";
        }
        double pressureJ = neighbor.getPressure();
        double densityJ = neighbor.getDensity();
        if(densityJ < 1e-12) continue;
        std::cout<<"Densityj = "<< densityJ<<std::endl;

        double avgPressure = (pressureI + pressureJ) / 2.0;
        std::cout<<"PressureI = "<< pressureI<<std::endl;
        std::cout<<"pressureJ = "<< pressureJ<<std::endl;

        double coefficient = -particleMass * avgPressure / densityJ;
        std::cout<<"Coefficient = "<< coefficient<<std::endl;

        for(int i = 0; i < pressureForce.size(); i++)
        {
            pressureForce[i] += gradW[i] * coefficient;
            if(std::isnan(pressureForce[i]) || std::isinf(pressureForce[i])) pressureForce[i] = 0.0;
            
            
        }
        
    }
    return pressureForce;
}


std::array<double,3> Particle::calculateViscosityForce(std::vector<Particle>& particles, double coreRadius){

    std::array<double,3> viscosityForce = {0.0,0.0,0.0};
    const double mu = 0.00089;

    for(int i =0 ; i < particles.size(); i++){

        std::array<double,3> particlePosition = particles[i].getPosition();

        std::array<double,3> direction = {position[0]-particlePosition[0],position[1]-particlePosition[1],position[2]-particlePosition[2]};
        double distance = findVector3Length(direction);
        std::array<double,3> normalizedDirection = normalizeVector3(direction);

        double particleMass = particles[i].getMass();


        if(distance < coreRadius && distance > 0){
            std::array<double,3> velocityDifference= {0,0,0};
            for(int j = 0; j < velocityDifference.size();j++){
                velocityDifference[j] = particles[i].getVelocity()[j] - this->getVelocity()[j];
            }
            double laplacianW = Wspike(distance,coreRadius);
            for(int j = 0; j <viscosityForce.size(); j++){
                viscosityForce[j] += mu*particleMass*velocityDifference[j]*laplacianW;
            }
        }
    }
    return viscosityForce;
}