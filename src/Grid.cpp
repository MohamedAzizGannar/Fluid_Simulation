#include <Grid.h>

#include <array>
#include <vector>

#include <unordered_map>
#include <Particle.h>

Grid::Grid(int cellSize):cellSize(cellSize){
}
Grid::Grid(){
    cellSize = 1.0;
}

void Grid::setCellSize(double newCellSize){
    cellSize = newCellSize;
}
double Grid::getCellSize(){return cellSize;}

void Grid::insertParticle(int particleID,const float3& particlePosition){

    auto [cellX,cellY,cellZ] = getCellCoordinates(particlePosition);
    
    int cellId = hashCellCoordinates({cellX,cellY,cellZ});

    auto it = cellHashMap.find(cellId);
    if(it != cellHashMap.end()){
        it->second.push_back(particleID);
    }
    else{
        cellHashMap[cellId] = std::vector<int>{particleID};
    }
}

int Grid::hashCellCoordinates(const std::array<int,3>& cellCoordinates){
    const int p1 = 73856093;
    const int p2 = 19349663;
    const int p3 = 83492791;

    auto [x,y,z] = cellCoordinates;

    return (x*p1 + y*p2 + z*p3);
}

std::array<int,3> Grid::getCellCoordinates( const float3& particlePosition){
    auto [x,y,z] = particlePosition;
    int cellX = std::floor(x/cellSize);
    int cellY = std::floor(y/cellSize);
    int cellZ = std::floor(z/cellSize);
    
    return {cellX,cellY,cellZ};
}

std::vector<int> Grid::getNeighbors(const float3& particlePosition){
    auto [cellX,cellY,cellZ] = getCellCoordinates(particlePosition);
    std::vector<int> neighborParticles = {};
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
                    neighborParticles.insert(neighborParticles.end(),cellParticles.end(),cellParticles.begin());
                }
            }
        }
    }
    return neighborParticles;
}

void Grid::rebuild(const std::vector<Particle>& particles){
    for(const auto& particle : particles){
        clear();
        auto position = particle.getPosition();
        insertParticle(particle.getId(),position);
    }
}

void Grid::clear(){
    cellHashMap.clear();
}