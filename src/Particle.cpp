#include <Particle.h>
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
    double k = 315/64*M_PI*std::pow(coreRadius,9);
    double secondMember = 0<= distance && distance <= coreRadius ? std::pow((coreRadius*coreRadius - distance*distance),3):0;
    return k * secondMember;
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


std::array<double,3> Particle::getPosition(){
    return position;
}

std::array<double,3> Particle::getAcceleration(){
    return acceleration;
}
std::array<double,3> Particle::getVelocity(){
    return velocity;
}
double Particle::getMass(){
    return mass;
}
double Particle::getRadius(){
    return radius;
}


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
double Particle::calculateDensity(std::vector<Particle> particles,Particle consideredParticle)
{
    double density = 0;
    for(int i = 0 ; i < particles.size(); i++){
        std::array<double,3> particlePosition = particles[i].getPosition();
        double distSquared = std::pow(particlePosition[0]-consideredParticle.getPosition()[0],2) + std::pow(particlePosition[1]-consideredParticle.getPosition()[1],2) + std::pow(particlePosition[2]-consideredParticle.getPosition()[3],2);
        double dist = std::sqrt(distSquared);
        double particleMass = particles[i].getMass();
        /*pressure += particleMass*densitySmoothingKernel(dist,h)
        With: 
            -h = core radius
            -densitySmoothingKernel(r-rj,h) the smoothing function associated to density giving us the density
            We will be using Wpoly6 for density
        */
    }
    return density;
}
double Particle::calculatePressure(std::vector<Particle> particles, Particle consideredParticle)
{
    double GAS_CONSTANT = 9e4 - 1e7;
    double REST_DENSITY = 1000;
    double pressure = GAS_CONSTANT * (calculateDensity(particles,consideredParticle)-REST_DENSITY);
    return pressure;
}
std::array<double,3> Particle::calculatePressureForce(std::vector<Particle> particles, double coreRadius,Particle consideredParticle)
{
    std::array<double,3> pressure = {0,0,0};
    for(int i = 0 ; i < particles.size(); i++){
        std::array<double,3> particlePosition = particles[i].getPosition();

        std::array<double,3> direction = {position[0]-particlePosition[0],position[1]-particlePosition[1],position[2]-particlePosition[2]};
        double distance = findVector3Length(direction);
        std::array<double,3> normalizedDirection = normalizeVector3(direction);

        double particleMass = particles[i].getMass();

        std::array<double,3> gradW = {0,0,0};
        for(int i = 0; i < gradW.size(); i++)
        {
            gradW[i] = normalizedDirection[i] * Wpoly6(distance, coreRadius);
        }
        double coefficient = (-particleMass) * (calculatePressure(particles,consideredParticle) + calculatePressure(particles,particles[i]))/2 * calculateDensity(particles,particles[i]);
        for(int i = 0; i < pressure.size(); i++)
        {
            pressure[i] += gradW[i] * coefficient;
        }
        /* 
        With:
            -pressureSmoothingKernel being the smoothing function used for pressure: in this case Wspiky
            -h =  core radius
            -gradW = the gradient of the smoothing kernel
        This will be used to calculate the pressure Force Vector with this formula

        Fpressure = - Sum(massj*(Pressurei + Pressurej)/(2*ρj)*gradW(ri-rj,h))
        With:
            -Pressurei = pressure at the considered point
            -Pressurej = pressure at the inflicting point
            -Massj = mass of inflicting particle
            -ρj = density at inflicting point
            -ri = position of the considered particle
            -rj = position of the inflicting particle
            -h = core radius
        */
    }
    return pressure;
}
