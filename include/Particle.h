#ifndef PARTICLE_H
#define PARTICLE_H

#include <array>
#include <vector>

struct float3{
    double x;
    double y;
    double z;

    float3 (double _x, double _y, double _z): x(_x),y(_y),z(_z){}
    float3():x(0),y(0),z(0){}
    double& operator[](int i) {
        return (&x)[i];
    }
    
    const double& operator[](int i) const {
        return (&x)[i];
    }

    float3 operator+(const float3& other)const{
        return(float3(x+other.x,y+other.y,z+other.z));
    }
    float3 operator+(const float& other)const{
        return(float3(x+other,y+other,z+other));
    }
    float3 operator-(const float& other)const{
        return(float3(x-other,y-other,z-other));
    }
    float3 operator-(const float3& other)const{
        return(float3(x-other.x,y-other.y,z-other.z));
    }
    
    float3 operator*(float scalar) const {
        return float3(x * scalar, y * scalar, z * scalar);
    }
    
    float3 operator/(float scalar) const {
        return float3(x / scalar, y / scalar, z / scalar);
    }
    float3& operator+=(const float3& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }
    
    float3& operator-=(const float3& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }
    
    float3& operator*=(float factor) {
        x *= factor;
        y *= factor;
        z *= factor;
        return *this;
    }
    
    float3& operator/=(float factor) {
        x /= factor;
        y /= factor;
        z /= factor;
        return *this;
    }

    float dot(const float3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }
    
    
    
    double lengthSQR()const{
        return (x*x + y*y + z*z);
    }
    double length()const{
        return sqrt(x*x + y*y + z*z);
    }
    float3 normalize() const{
        double invLen = 1/length();
        
        return invLen>0? *this*invLen :float3();
    }

};
struct ColorFieldProperties {
    double field;
    float3 gradient;
    double laplacian;
};
class Particle{

    private:
    float3 position;
    float3 predictedPosition;

    float3 velocity;
    float3 acceleration;


    const double radius = 0.5;

    const double mass = 1.0;

    unsigned int id;

    double density;
    double pressure;

    public:
    //Constructor
    Particle(float3 initialPosition,float3 initialVelocity,float3 intialAcceleration, int id);

    //Default Constructor (Array Initialization)
    Particle();

    ColorFieldProperties calculateColorFieldProperties(const std::vector<Particle>& particles);

    //Getters
    float3 getPosition() const;
    float3 getPredictedPosition() const;

    float3 getAcceleration()const;
    float3 getVelocity()const;
    double getMass()const;
    double getRadius()const;
    double getPressure()const;
    double getDensity()const;
    int getId()const;


    //Setters
    void setPosition(float3  newPosition);
    void setPredictedPosition(float3  newPredictedPosition);
    void setVelocity(float3 newVelocity);
    void setAcceleration(float3 newAcceleration);
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

    float3 calculatePressureForce(const std::vector<Particle>& particles);

    float3 calculateViscosityForce(const std::vector<Particle>& particles);
    
    float3 calculateGravitationalPull(const std::vector<Particle>& particles);
    float3 calculateGravity();

    double calculateSmoothedColorField( const std::vector<Particle>& particles);
    float3 calculateColorFieldGradient( const std::vector<Particle>& particles);
    double calculateColorFieldLaplacian( const std::vector<Particle>& particles);
    float3 calculateSurfaceTensionForce(  const std::vector<Particle>& particles);

    void applyForces(const std::vector<Particle>& particles);
    void applyForcesOptimised(const std::vector<Particle>& particles);

};
#endif