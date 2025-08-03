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
    Grid(double cellSize);
    void setCellSize(double newCellSize);
    double getCellSize();


    void clear();
    void insertParticle(int particleID, const std::array<double,3>& particlePosition);
    std::array<int,3> getCellCoordinates( const std::array<double,3>& particlePosition);
    std::vector<int> getNeighbors( const std::array<double,3>& particlePosition);
    void rebuild(const std::vector<Particle>& particles);

    private:
    int hashCellCoordinates(const std::array<int,3>& cellCoordinates);


};

#endif