#include <Grid.h>

#include <array>
#include<iostream>
#include <vector>

#include <unordered_map>
#include <Particle.h>

Grid::Grid(double cellSize):cellSize(cellSize){
}
Grid::Grid(){
    cellSize = 1.0;
}

void Grid::setCellSize(double newCellSize){
    cellSize = newCellSize;
}
double Grid::getCellSize(){return cellSize;}

void Grid::insertParticle(int particleID,const float3& particlePosition){

    auto [cellX, cellY, cellZ] = getCellCoordinates(particlePosition);
    
    uint64_t cellId = hashCellCoordinates({cellX, cellY, cellZ});

    auto it = cellHashMap.find(cellId);
    if(it != cellHashMap.end()){
        it->second.push_back(particleID);
    }
    else{
        std::vector<int> newCell;
        newCell.reserve(8); 
        newCell.push_back(particleID);
        cellHashMap[cellId] = std::move(newCell);
    }

}

uint64_t Grid::hashCellCoordinates(const std::array<int,3>& cellCoordinates){
    const uint64_t p1 = 0x9e3779b97f4a7c15ULL;
    const uint64_t p2 = 0x3243f6a8885a308dULL;
    const uint64_t p3 = 0x243f6a8885a308d3ULL;

    auto [x, y, z] = cellCoordinates;
    
    uint64_t hash = 0;
    hash ^= static_cast<uint64_t>(x) * p1;
    hash ^= static_cast<uint64_t>(y) * p2;
    hash ^= static_cast<uint64_t>(z) * p3;
    
    hash ^= hash >> 32;
    hash *= 0x9e3779b97f4a7c15ULL;
    hash ^= hash >> 32;
    
    return hash;
}

std::array<int,3> Grid::getCellCoordinates( const float3& particlePosition){
    auto [x,y,z] = particlePosition;
    int cellX = std::floor(x/cellSize);
    int cellY = std::floor(y/cellSize);
    int cellZ = std::floor(z/cellSize);
    
    return {cellX,cellY,cellZ};
}

void Grid::getNeighbors(const float3& particlePosition, std::vector<int>& outNeighbors){
    outNeighbors.clear();
    outNeighbors.reserve(64);
    auto [cellX,cellY,cellZ] = getCellCoordinates(particlePosition);
    for(int dx = -1; dx <= 1; dx++){
        for(int dy = -1; dy <= 1; dy++){
            for(int dz = -1; dz <= 1; dz++){
                int currCellHash = hashCellCoordinates({
                    cellX + dx,
                    cellY + dy,
                    cellZ + dz
                });
                auto it = cellHashMap.find(currCellHash);
                if( it != cellHashMap.end()){
                    const auto& cellParticles = it ->second;
                    outNeighbors.insert(outNeighbors.end(),cellParticles.begin(),cellParticles.end());
                }
            }
        }
    }
}

void Grid::rebuild(const std::vector<Particle>& particles){
    clear();
    size_t estimatedCells = particles.size() / 4; // Rough estimate
    cellHashMap.reserve(estimatedCells);
    for(const auto& particle : particles){
        auto position = particle.getPosition();
        insertParticle(particle.getId(),position);
    }
}

void Grid::clear(){
    cellHashMap.clear();
}

void Grid::printGridStats() const {
    if(cellHashMap.empty()){
        std::cout << "Grid is empty\n";
        return;
    }
    
    size_t totalParticles = 0;
    size_t maxParticlesPerCell = 0;
    size_t minParticlesPerCell = SIZE_MAX;
    
    for(const auto& [cellId, particles] : cellHashMap){
        totalParticles += particles.size();
        maxParticlesPerCell = std::max(maxParticlesPerCell, particles.size());
        minParticlesPerCell = std::min(minParticlesPerCell, particles.size());
    }
    
    double avgParticlesPerCell = static_cast<double>(totalParticles) / cellHashMap.size();
    
    std::cout << "Grid Statistics:\n";
    std::cout << "  Total cells: " << cellHashMap.size() << "\n";
    std::cout << "  Total particles: " << totalParticles << "\n";
    std::cout << "  Avg particles per cell: " << avgParticlesPerCell << "\n";
    std::cout << "  Min particles per cell: " << minParticlesPerCell << "\n";
    std::cout << "  Max particles per cell: " << maxParticlesPerCell << "\n";
    std::cout << "  Cell size: " << cellSize << "\n";
}