#ifndef LEGALIZER_H
#define LEGALIZER_H

#include <chrono>
#include <limits>
#include <utility>
#include <string>
#include <vector>

struct Instance {
    std::string name;
    long long x = 0, y = 0, w = 0, h = 0;
    std::string orient = "R0";
    std::string type;
};

struct Interval {
    int l = 0;
    int r = 0; // exclusive
};

struct Row {
    std::vector<Interval> freeSites;
};

struct Candidate {
    int row = -1;
    int site = -1;
    double cost = std::numeric_limits<double>::infinity();
    double displacement = std::numeric_limits<double>::infinity();
};

class Legalizer {
public:
    bool readInput(const std::string &path);
    bool legalize(double alpha, double threshold);
    bool writeOutput(const std::string &path) const;

private:
    long long dbuPerMicron = 0;
    long long dieLLX = 0, dieLLY = 0, dieURX = 0, dieURY = 0;
    long long siteW = 0, siteH = 0;
    int numSites = 0;
    int numRows = 0;
    std::vector<Instance> cells;
    std::vector<Instance> obstacles;
    std::vector<Row> rows;
    std::vector<long long> placementX;
    std::vector<long long> placementY;
    std::vector<int> placementRow;
    std::vector<int> placementSite;
    std::vector<int> cellWidthSites;
    std::vector<int> cellHeightRows;
    std::vector<int> targetRows;
    std::vector<int> targetSites;
    std::vector<int> siteOwner;
    std::vector<long long> bestPlacementX;
    std::vector<long long> bestPlacementY;
    std::vector<int> bestPlacementRow;
    std::vector<int> bestPlacementSite;
    std::vector<int> bestSiteOwner;

    double densityThreshold = 1.0; // Normalized ratio in [0, 1].
    double bestQuality = std::numeric_limits<double>::infinity();
    long long gridSize = 0;
    long long gridLLX = 0, gridLLY = 0, gridURX = 0, gridURY = 0;
    int gridCols = 0;
    int gridRows = 0;
    std::vector<double> gridAreaUsed;
    std::vector<double> bestGridAreaUsed;
    std::vector<unsigned char> gridExcluded;
    bool runtimeLimitHit = false;
    std::chrono::steady_clock::time_point runtimeDeadline;

    int targetRow(const Instance &cell) const;
    int targetSite(const Instance &cell, int widthSites) const;
    void buildRows();
    void buildSiteOwnerGrid();
    void buildDensityGrid();
    bool snapInitialLegalize(double alpha);
    bool abacusInitialLegalize(double alpha);
    std::vector<std::pair<int, int>> packAbacusSubrow(const std::vector<int> &cellIdxs,
                                                      const Interval &interval) const;
    double packedSubrowCost(const std::vector<std::pair<int, int>> &packed,
                            int row, double alpha) const;
    std::vector<int> candidateSitesInRow(int row, int target, int widthSites) const;
    bool canPlaceMultiRow(int row, int site, int widthSites, int heightRows) const;
    Candidate findBestCandidate(int cellIdx) const;
    bool ownerCanPlace(int cellIdx, int row, int site) const;
    void setOwner(int cellIdx, int row, int site, int value);
    void removePlacedCell(int cellIdx);
    void restorePlacedCell(int cellIdx, int row, int site);
    double displacementCost(int cellIdx, long long x, long long y) const;
    double placementCost(int cellIdx, long long x, long long y, double alpha) const;
    double averageDisplacement() const;
    double densityOverflowRatio() const;
    double qualityScore(double alpha) const;
    bool runtimeLimitReached();
    void saveBestPlacement(double quality);
    bool updateBestPlacement(double alpha);
    void restoreBestPlacement();
    bool tryMoveCell(int cellIdx, int rowRadius, int siteRadius, double alpha);
    void detailPlace(double alpha);
    void refineCells(double alpha);
    void refinePairSwaps(double alpha);
    void refineSmallReinsert(double alpha);
    void refineDisplacementSwaps(double alpha);
    void refineRowSliding(double alpha);
    void refineBoundaryTouches(double alpha);
    double estimateDensityPenalty(long long x, long long y, long long w, long long h) const;
    void occupy(int row, int site, int widthSites, int heightRows);
    void addCellToDensity(long long x, long long y, long long w, long long h);
    void addCellToDensityDelta(long long x, long long y, long long w, long long h, double sign);
};

#endif
