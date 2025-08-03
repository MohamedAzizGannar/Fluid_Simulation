#ifndef PARTICLE_H
#define PARTICLE_H

#include <array>
#include <vector>
struct ColorFieldProperties {
    double field;
    std::array<double,3> gradient;
    double laplacian;
};
class Particle{

    private:
    std::array<double,3> position;
    std::array<double,3> predictedPosition;

    std::array<double,3> velocity;
    std::array<double,3> acceleration;


    const double radius = 0.5;

    const double mass = 1.0;

    unsigned int id;

    double density;
    double pressure;

    public:
    //Constructor
    Particle(std::array<double,3> initialPosition,std::array<double,3> initialVelocity,std::array<double,3> intialAcceleration, int id);

    //Default Constructor (Array Initialization)
    Particle();

    ColorFieldProperties calculateColorFieldProperties(const std::vector<Particle>& particles);

    //Getters
    std::array<double,3> getPosition() const;
    std::array<double,3> getPredictedPosition() const;

    std::array<double,3> getAcceleration()const;
    std::array<double,3> getVelocity()const;
    double getMass()const;
    double getRadius()const;
    double getPressure()const;
    double getDensity()const;
    int getId()const;


    //Setters
    void setPosition(std::array<double,3>  newPosition);
    void setPredictedPosition(std::array<double,3>  newPredictedPosition);
    void setVelocity(std::array<double,3> newVelocity);
    void setAcceleration(std::array<double,3> newAcceleration);
    void setDensity(double newDensity);
    void setPressure(double newPressure);
    void setId(int i);


    //Update Position 
    void updatePosition();

    //Update Velocity
    void updateVelocity();
    void updatePressure(const std::vector<Particle>& particles);
    void updateDensity(const std::vector<Particle>& particles);

    //SPH Functions

    double calculateDensity(const std::vector<Particle>& particles);

    double calculatePressure(const std::vector<Particle>& particles);

    std::array<double,3> calculatePressureForce(const std::vector<Particle>& particles);

    std::array<double,3> calculateViscosityForce(const std::vector<Particle>& particles);
    
    std::array<double,3> calculateGravitationalPull(const std::vector<Particle>& particles);
    std::array<double,3> calculateGravity();

    double calculateSmoothedColorField( const std::vector<Particle>& particles);
    std::array<double,3> calculateColorFieldGradient( const std::vector<Particle>& particles);
    double calculateColorFieldLaplacian( const std::vector<Particle>& particles);
    std::array<double,3> calculateSurfaceTensionForce(  const std::vector<Particle>& particles);

    void applyForces(const std::vector<Particle>& particles);
};
#endif