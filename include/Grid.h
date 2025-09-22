#ifndef GRID_H
#define GRID_H

#include <array>
#include <vector>

#include <unordered_map>
#include <Particle.h>


class Grid{
    private:
    double cellSize;
    std::unordered_map<int,std::vector<int>> cellHashMap; // int = Hashed Value of cell Coordinates || std::vector<int> IDs of the contained Particles

    public:
    Grid();
    Grid( double cellSize);
    void setCellSize(double newCellSize);
    double getCellSize();


    void clear();
    void insertParticle(int particleID, const float3& particlePosition);
    std::array<int,3> getCellCoordinates( const float3& particlePosition);
    void getNeighbors( const float3& particlePosition, std::vector<int>& outNeighbors);
    void rebuild(const std::vector<Particle>& particles);
    void printGridStats() const ;
    private:
    uint64_t hashCellCoordinates(const std::array<int,3>& cellCoordinates);


};

#endif