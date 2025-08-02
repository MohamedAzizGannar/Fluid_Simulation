#include <Particle.h>
#include <iostream>
#include <array>
#include <vector>
#include <cmath>

const double TENSION_COEFFICIENT = 7.28;
const double GAS_CONSTANT = 300.0;
const double REST_DENSITY = 1000.0;
const double GRAVITATIONAL_CONSTANT = 6.674e-11;
const double CORE_RADIUS = 1.5;
const double CORE_RADIUS9 = std::pow(CORE_RADIUS, 9);
const double MU = 0.00089;


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
double Wpoly6(double distance){
    if(distance < 0.0 || distance > CORE_RADIUS)return 0.0;
    const double coefficient = 315.0/(64.0*M_PI*std::pow(CORE_RADIUS,9));
    const double secondMember = std::pow((CORE_RADIUS*CORE_RADIUS - distance*distance),3);
    return coefficient * secondMember;
}
std::array<double,3> gradientWpoly6(const std::array<double,3> vect){
    const double distance = findVector3Length(vect);
    if(distance< 0.0 || distance > CORE_RADIUS) return {0.0,0.0,0.0};
    const double coefficient = - 945.0 / (32.0 * M_PI * std::pow(CORE_RADIUS,9));
    const double secondMember = std::pow(CORE_RADIUS*CORE_RADIUS - distance * distance , 2);
    return {
        coefficient*secondMember*vect[0]/distance,
        coefficient*secondMember*vect[1]/distance,
        coefficient*secondMember*vect[2]/distance
    };
}
double magnitudeGradientWpoly6(const std::array<double,3>& vect){
    const double distance = findVector3Length(vect);
    if(distance < 0 || distance > CORE_RADIUS)return 0.0;
    const double coefficient = - 945.0 / (32.0 * M_PI * std::pow(CORE_RADIUS,9));
    const double secondMember = std::pow(CORE_RADIUS*CORE_RADIUS - distance * distance , 2);
    return coefficient * secondMember * distance;
}
double  laplacianWpoly6(const std::array<double,3>& vect){
    const double distance = findVector3Length(vect);
    if(distance< 0.0 || distance > CORE_RADIUS) return 0.0;
    const double coefficient = - 945.0 / (32.0 * M_PI * std::pow(CORE_RADIUS,9));
    const double secondMember = (CORE_RADIUS*CORE_RADIUS - distance * distance);
    const double thirdMember = (3*distance*distance - CORE_RADIUS * CORE_RADIUS) ;
    return coefficient * secondMember * thirdMember;
}
std::array<double,3> gradientWspiky( const std::array<double,3>& vect){
    const double distance = findVector3Length(vect);
    if(distance<0 || distance > CORE_RADIUS){return {0.0,0.0,0.0};}
    const double coefficient = -45.0/(M_PI * std::pow(CORE_RADIUS,6));
    const double secondMember =std::pow(CORE_RADIUS-distance,2); 
    return {
        coefficient * secondMember * vect[0]/distance,
        coefficient * secondMember * vect[1]/distance,
        coefficient * secondMember * vect[2]/distance
    };
}
double laplacianViscosityKernel(const std::array<double,3>& vect){
    const double distance = findVector3Length(vect);
    if(distance < 0 || distance > CORE_RADIUS) return 0.0;
    const double coefficient = 45.0/(M_PI * std::pow(CORE_RADIUS,6));
    return coefficient * (distance - CORE_RADIUS);
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
    double newDensity = this->calculateDensity(particles);
    this->setDensity(newDensity);
}

void Particle::updatePressure(const std::vector<Particle>& particles){
    double newPressure = this->calculatePressure(particles);
    this->setPressure(newPressure);
}

//SPH FUNCTIONS
double Particle::calculateDensity(const std::vector<Particle>& particles)
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
        if(distance>CORE_RADIUS)continue;
        density+= neighbor.getMass() * Wpoly6(distance);
    }
    return density < 1e-12 ? 1e-12:density;
}
double Particle::calculatePressure(const std::vector<Particle>& particles){

 

    double pressure = GAS_CONSTANT * (this->calculateDensity(particles) - REST_DENSITY);
    pressure = std::max(0.0,pressure);
    return pressure;
}
std::array<double,3> Particle::calculatePressureForce(const std::vector<Particle>& particles){

    std::array<double,3> pressureForce = {0.0,0.0,0.0};
    double pressureI = this->getPressure();


    for(const Particle& neighbor: particles){
        if(&neighbor == this) continue;
        std::array<double,3> particlePosition = neighbor.getPosition();

        std::array<double,3> direction = {position[0]-particlePosition[0],position[1]-particlePosition[1],position[2]-particlePosition[2]};

        double distance = findVector3Length(direction);
        if(distance<= 0.0001 || distance>=CORE_RADIUS) continue;

        std::array<double,3> normalizedDirection = normalizeVector3(direction);

        double particleMass = neighbor.getMass();

        std::array<double,3> gradW = {0.0,0.0,0.0};
        std::array<double,3> gradientWSpikyVal = gradientWspiky(normalizedDirection);
        for(int i = 0; i < gradW.size(); i++)
        {
            gradW[i] = gradientWSpikyVal[i];
        }
        double pressureJ = neighbor.getPressure();
        double densityJ = neighbor.getDensity();
        if(densityJ < 1e-12) continue;

        double avgPressure = (pressureI + pressureJ) / 2.0;

        double coefficient = -particleMass * avgPressure / densityJ;

        for(int i = 0; i < pressureForce.size(); i++)
        {
            pressureForce[i] += gradW[i] * coefficient;
            if(std::isnan(pressureForce[i]) || std::isinf(pressureForce[i])) pressureForce[i] = 0.0;
        }
    }
    return pressureForce;
}


std::array<double,3> Particle::calculateViscosityForce(const std::vector<Particle>& particles){

    std::array<double,3> viscosityForce = {0.0,0.0,0.0};

    for(const Particle& neighbor : particles){

        std::array<double,3> particlePosition = neighbor.getPosition();

        std::array<double,3> direction = {position[0]-particlePosition[0],position[1]-particlePosition[1],position[2]-particlePosition[2]};
        double distance = findVector3Length(direction);
        std::array<double,3> normalizedDirection = normalizeVector3(direction);

        double particleMass = neighbor.getMass();


        if(distance < CORE_RADIUS && distance > 0){
            std::array<double,3> velocityDifference= {0,0,0};
            for(int j = 0; j < velocityDifference.size();j++){
                velocityDifference[j] = neighbor.getVelocity()[j] - this->getVelocity()[j];
            }
            double densityJ = neighbor.getDensity();
            double laplacianW = laplacianViscosityKernel(direction);
            for(int j = 0; j <viscosityForce.size(); j++){
                viscosityForce[j] += MU*particleMass*velocityDifference[j]/densityJ*laplacianW;
            }
        }
    }
    return viscosityForce;
}
std::array<double,3> Particle::calculateGravitationalPull(const std::vector<Particle>& particles){
    std::array<double,3> gravityForce = {0.0,0.0,0.0};
    for(const auto& neighbor : particles)
    {
        if (&neighbor == this) continue;
        const std::array<double,3> neighborPosition = neighbor.getPosition();
        const double neighborMass = neighbor.getMass();
        std::array<double,3> vect;
        for(int i = 0; i < 3 ; i ++)
        {
            vect[i] = neighborPosition[i] - this->getPosition()[i];
        }
        const double distance = findVector3Length(vect);
        if (distance < 1e-6) continue;
        std::array<double,3> unitDirection;
        for (int i = 0; i < 3; i++) {
            unitDirection[i] = vect[i] / distance;
        }
        const double gravitationalPullCoefficient = GRAVITATIONAL_CONSTANT * neighborMass * this->getMass() / (distance * distance);
        for(int i = 0; i < 3 ; i ++)
        {
            gravityForce[i] += gravitationalPullCoefficient * unitDirection[i];
        }
    }
    return gravityForce;
}
std::array<double,3> Particle::calculateGravity(){
    return {0.0,-9.81,0.0};
}

double Particle::calculateSmoothedColorField( const std::vector<Particle>& particles){
    double sum = 0.0;
    const std::array<double,3> pos = this ->getPosition();
    for(const auto& pj : particles){

        const std::array<double,3> posj = pj.getPosition();
        const std::array<double,3> vect = {pos[0]-posj[0],pos[1]-posj[1],pos[2]-posj[2]};
        const double distance = findVector3Length(vect);
        const double massJ = pj.getMass();
        const double densityJ = pj.getDensity();
        const double poly6Value = Wpoly6(distance);
        if(densityJ < 1e-6)continue;
        const double coefficient = massJ / densityJ;
        sum += poly6Value * coefficient;
    }
    return sum;
}

std::array<double,3> Particle::calculateColorFieldGradient( const std::vector<Particle>& particles){
    std::array<double,3> sum = {0.0,0.0,0.0};
    const std::array<double,3> pos = this ->getPosition();
    for(const auto& pj : particles){

        const std::array<double,3> posj = pj.getPosition();
        const std::array<double,3> vect = {pos[0]-posj[0],pos[1]-posj[1],pos[2]-posj[2]};
        const double massJ = pj.getMass();
        const double densityJ = pj.getDensity();
        const std::array<double,3> Wpoly6Gradient = gradientWpoly6(vect);
        if(densityJ < 1e-6)continue;
        const double coefficient = massJ/densityJ;
        for(int i = 0; i < sum.size(); i++){
            sum[i] +=   Wpoly6Gradient[i] * coefficient;
        }
    }
    return sum;
}
double Particle::calculateColorFieldLaplacian(  const std::vector<Particle>& particles){
    double sum = 0.0;
    const std::array<double,3> pos = this ->getPosition();
    for(const auto& pj : particles){
        const std::array<double,3> posj = pj.getPosition();
        const std::array<double,3> vect = {pos[0]-posj[0],pos[1]-posj[1],pos[2]-posj[2]};
        const double massJ = pj.getMass();
        const double densityJ = pj.getDensity();
        const double Wpoly6Laplacian = laplacianWpoly6(vect);
        if(densityJ < 1e-6)continue;
        const double coefficient = massJ/densityJ; 
        sum += massJ/densityJ * Wpoly6Laplacian;
    }
    return sum;
}
double calculateCurvature(double laplacian,const std::array<double,3>& gradient){
    const double magnitude = findVector3Length(gradient);
    if(magnitude < 1e-6)return 0.0;
    const double K = - laplacian / magnitude;
    return K;
}

std::array<double,3> Particle::calculateSurfaceTensionForce(  const std::vector<Particle>& particles){
    const std::array<double,3> pos = this ->getPosition();
    std::array<double,3> surfaceTensionForce = {0.0,0.0,0.0};


    const std::array<double,3> gradient = calculateColorFieldGradient(particles);
    const double laplacian = calculateColorFieldLaplacian(particles);
    const double curvature = calculateCurvature(laplacian,gradient);
    for(int i = 0 ; i < surfaceTensionForce.size(); i++){
        surfaceTensionForce[i] += -TENSION_COEFFICIENT * curvature * gradient[i];
    }
    
    
    return surfaceTensionForce;
}

void Particle::applyForces(const std::vector<Particle>& particles){
    std::array<double,3> totalForce = {0.0,0.0,0.0};

    const std::array<double,3> pressureForce = calculatePressureForce( particles);
    const std::array<double,3> viscosityForce = calculateViscosityForce( particles);
    const std::array<double,3> gravitationalForce = calculateGravitationalPull( particles);
    const std::array<double,3> gravity = calculateGravity();

    
    for(int i = 0; i < 3; i ++){
        totalForce[i] = pressureForce[i] + viscosityForce[i] + gravitationalForce[i] + gravity[i];
    }
    std::array<double,3> newAcceleration = {0.0, 0.0, 0.0}; 
    const double density = this->getDensity();
    if (density > 1e-6){
        for(int i = 0 ; i < 3; i++){
        newAcceleration[i] = totalForce[i] / density;
        }
    }
    this->setAcceleration(newAcceleration);
}