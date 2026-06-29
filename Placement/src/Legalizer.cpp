#include "Legalizer.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

using namespace std;

namespace {

constexpr int kRuntimeLimitMinutes = 25;
constexpr size_t kTetrisOnlyCellThreshold = 1000000;
constexpr double kDisplacementNormFactor = 18.2;

long long ceilDiv(long long a, long long b) {
    if (b < 0) {
        a = -a;
        b = -b;
    }
    if (a >= 0) return (a + b - 1) / b;
    return -((-a) / b);
}

long long floorDiv(long long a, long long b) {
    if (b < 0) {
        a = -a;
        b = -b;
    }
    if (a >= 0) return a / b;
    return -((-a + b - 1) / b);
}

int gridLineLowerBound(long long value, long long low, long long high, long long step, int bins) {
    if (value <= low) return 0;
    if (value >= high) return bins;
    return static_cast<int>(min(static_cast<long long>(bins), ceilDiv(value - low, step)));
}

int gridLineUpperBound(long long value, long long low, long long high, long long step, int bins) {
    if (value < low) return 0;
    if (value >= high) return bins + 1;
    return static_cast<int>(min(static_cast<long long>(bins + 1), floorDiv(value - low, step) + 1));
}

bool rectOverlap(long long ax1, long long ay1, long long ax2, long long ay2,
                 long long bx1, long long by1, long long bx2, long long by2) {
    return ax1 < bx2 && bx1 < ax2 && ay1 < by2 && by1 < ay2;
}

vector<string> splitTokens(const string &line) {
    stringstream ss(line);
    vector<string> tokens;
    string token;
    while (ss >> token) tokens.push_back(token);
    return tokens;
}

bool isCommentStart(const string &token) {
    return token == "//" || (!token.empty() && token[0] == '#');
}

int columnIndex(const unordered_map<string, int> &columns, const string &name) {
    auto it = columns.find(name);
    return it == columns.end() ? -1 : it->second;
}

string normalizedOrient(const string &orient) {
    return (orient == "R0" || orient == "MX") ? orient : "R0";
}

bool parseIntegerToken(const string &token, long long &value) {
    try {
        size_t parsed = 0;
        value = stoll(token, &parsed);
        return parsed == token.size();
    } catch (const exception &) {
        return false;
    }
}

void printStage(const string &stage) {
    cout << "[Stage] " << stage << '\n';
}

struct AbacusCluster {
    int x = 0;
    int width = 0;
    int weight = 0;
    long long q = 0;
    vector<int> cells;
};

struct AbacusSubrowState {
    int row = -1;
    Interval interval;
    int usedSites = 0;
    vector<int> cells;
    double cost = 0.0;
};

struct TetrisRowState {
    size_t interval = 0;
    int cursor = 0;
};

} // namespace

bool Legalizer::readInput(const string &path) {
    ifstream fin(path);
    if (!fin) {
        cerr << "Cannot open input file: " << path << '\n';
        return false;
    }

    dbuPerMicron = 0;
    dieLLX = dieLLY = dieURX = dieURY = 0;
    siteW = siteH = 0;
    cells.clear();
    obstacles.clear();

    string line;
    vector<string> header;
    unordered_map<string, int> columns;
    int lineNo = 0;
    while (getline(fin, line)) {
        ++lineNo;
        vector<string> parts = splitTokens(line);
        if (parts.empty() || isCommentStart(parts[0])) continue;
        const string &key = parts[0];
        if (key == "DBU_Per_Micron") {
            if (parts.size() < 2 || !parseIntegerToken(parts[1], dbuPerMicron)) {
                cerr << "Invalid DBU_Per_Micron on line " << lineNo << ".\n";
                return false;
            }
        } else if (key == "DieArea_LL") {
            if (parts.size() < 3 || !parseIntegerToken(parts[1], dieLLX) ||
                !parseIntegerToken(parts[2], dieLLY)) {
                cerr << "Invalid DieArea_LL on line " << lineNo << ".\n";
                return false;
            }
        } else if (key == "DieArea_UR") {
            if (parts.size() < 3 || !parseIntegerToken(parts[1], dieURX) ||
                !parseIntegerToken(parts[2], dieURY)) {
                cerr << "Invalid DieArea_UR on line " << lineNo << ".\n";
                return false;
            }
        } else if (key == "Site_Width") {
            if (parts.size() < 2 || !parseIntegerToken(parts[1], siteW)) {
                cerr << "Invalid Site_Width on line " << lineNo << ".\n";
                return false;
            }
        } else if (key == "Site_Height") {
            if (parts.size() < 2 || !parseIntegerToken(parts[1], siteH)) {
                cerr << "Invalid Site_Height on line " << lineNo << ".\n";
                return false;
            }
        } else if (key == "Name") {
            header = parts;
            for (int i = 0; i < static_cast<int>(header.size()); ++i) {
                columns[header[i]] = i;
            }
            for (const string &required : {"Name", "LLX", "LLY", "Width", "Height", "Type"}) {
                if (columnIndex(columns, required) < 0) {
                    cerr << "Missing input column: " << required << '\n';
                    return false;
                }
            }
            break;
        }
    }

    if (header.empty()) {
        cerr << "Missing input table header.\n";
        return false;
    }

    while (getline(fin, line)) {
        ++lineNo;
        vector<string> parts = splitTokens(line);
        if (parts.empty() || isCommentStart(parts[0])) continue;
        if (parts.size() < 6) {
            cerr << "Invalid instance row on line " << lineNo << ".\n";
            return false;
        }

        Instance inst;
        auto requiredToken = [&](const string &column, string &value) -> bool {
            int idx = columnIndex(columns, column);
            if (idx < 0 || idx >= static_cast<int>(parts.size())) return false;
            value = parts[idx];
            return true;
        };
        auto requiredInteger = [&](const string &column, long long &value) -> bool {
            string token;
            return requiredToken(column, token) && parseIntegerToken(token, value);
        };

        if (!requiredToken("Name", inst.name) ||
            !requiredInteger("LLX", inst.x) ||
            !requiredInteger("LLY", inst.y) ||
            !requiredInteger("Width", inst.w) ||
            !requiredInteger("Height", inst.h)) {
            cerr << "Invalid instance geometry on line " << lineNo << ".\n";
            return false;
        }

        int orientCol = columnIndex(columns, "Orient");
        int typeCol = columnIndex(columns, "Type");
        if (orientCol >= 0 && typeCol > orientCol &&
            parts.size() == header.size() - 1) {
            inst.orient = "R0";
            inst.type = parts[typeCol - 1];
        } else if (parts.size() == header.size()) {
            if (orientCol >= 0) inst.orient = parts[orientCol];
            inst.type = parts[typeCol];
        } else if (parts.size() == header.size() + 1 && orientCol < 0) {
            inst.orient = parts[typeCol];
            inst.type = parts[typeCol + 1];
        } else {
            inst.orient = parts.size() >= 7 ? parts[5] : "R0";
            inst.type = parts.back();
        }
        inst.orient = normalizedOrient(inst.orient);

        if (inst.type.empty()) {
            cerr << "Missing instance type on line " << lineNo << ".\n";
            return false;
        }
        if (inst.type == "CELL") {
            cells.push_back(inst);
        } else {
            obstacles.push_back(inst);
        }
    }

    if (dbuPerMicron <= 0 || siteW <= 0 || siteH <= 0 || dieURX <= dieLLX || dieURY <= dieLLY) {
        cerr << "Invalid input geometry.\n";
        return false;
    }
    return true;
}

bool Legalizer::legalize(double alpha, double threshold) {
    // Start with fast Tetris placement, then keep the best density-aware refinement.
    printStage("Initializing legalization");
    runtimeLimitHit = false;
    runtimeDeadline = chrono::steady_clock::now() + chrono::minutes(kRuntimeLimitMinutes);
    if (!isfinite(threshold) || threshold < 0.0 || threshold > 100.0) {
        cerr << "Threshold must be in [0, 100].\n";
        return false;
    }
    densityThreshold = threshold / 100.0;

    numSites = static_cast<int>((dieURX - dieLLX) / siteW);
    numRows = static_cast<int>((dieURY - dieLLY) / siteH);
    if (numSites <= 0 || numRows <= 0) {
        cerr << "No legal site rows in die.\n";
        return false;
    }

    printStage("Building rows and density grid");
    buildRows();
    buildSiteOwnerGrid();
    buildDensityGrid();

    placementX.assign(cells.size(), 0);
    placementY.assign(cells.size(), 0);
    placementRow.assign(cells.size(), -1);
    placementSite.assign(cells.size(), -1);
    cellWidthSites.assign(cells.size(), 0);
    cellHeightRows.assign(cells.size(), 0);

    bool haveInitial = false;
    bool tetrisOnly = cells.size() > kTetrisOnlyCellThreshold;
    printStage("Running Tetris initial legalization");
    if (snapInitialLegalize(alpha)) {
        saveBestPlacement(qualityScore(alpha));
        haveInitial = true;
    }

    if (tetrisOnly) {
        if (!haveInitial) return false;
        printStage("Skipping Abacus initial legalization for large Tetris-only placement");
        restoreBestPlacement();
        return true;
    }

    printStage("Preparing Abacus initial legalization");
    buildRows();
    buildSiteOwnerGrid();
    gridAreaUsed.assign(static_cast<size_t>(gridCols) * gridRows, 0.0);
    placementX.assign(cells.size(), 0);
    placementY.assign(cells.size(), 0);
    placementRow.assign(cells.size(), -1);
    placementSite.assign(cells.size(), -1);

    printStage("Running Abacus initial legalization");
    if (abacusInitialLegalize(alpha)) {
        if (haveInitial) {
            updateBestPlacement(alpha);
        } else {
            saveBestPlacement(qualityScore(alpha));
            haveInitial = true;
        }
    }
    if (!haveInitial) return false;

    restoreBestPlacement();

    printStage("Running detail placement");
    detailPlace(alpha);
    printStage("Refining adjacent row cells");
    refineCells(alpha);
    printStage("Refining pair swaps");
    refinePairSwaps(alpha);
    printStage("Refining small-case reinsertion");
    refineSmallReinsert(alpha);
    printStage("Refining displacement swaps");
    refineDisplacementSwaps(alpha);
    printStage("Refining row sliding");
    refineRowSliding(alpha);
    for (int pass = 0; pass < 6; ++pass) {
        printStage("Refining density boundary touches pass " + to_string(pass + 1));
        refineBoundaryTouches(alpha);
    }

    printStage("Restoring best placement");
    updateBestPlacement(alpha);
    restoreBestPlacement();

    if (runtimeLimitHit) {
        cerr << "Runtime limit reached; outputting best placement found.\n";
    }

    return true;
}

bool Legalizer::writeOutput(const string &path) const {
    ofstream fout(path);
    if (!fout) {
        cerr << "Cannot open output file: " << path << '\n';
        return false;
    }
    fout << fixed << setprecision(6);
    for (size_t i = 0; i < cells.size(); ++i) {
        long long originX = placementX[i];
        long long originY = placementY[i];
        if (cells[i].orient == "MX") {
            originY += cells[i].h;
        }
        fout << "place_cell -inst_name " << cells[i].name
             << " -orient " << cells[i].orient
             << " -origin {" << static_cast<double>(originX) / dbuPerMicron
             << ' ' << static_cast<double>(originY) / dbuPerMicron << "}\n";
    }
    return true;
}

int Legalizer::targetRow(const Instance &cell) const {
    long long rel = cell.y - dieLLY;
    long long rounded = floorDiv(rel + siteH / 2, siteH);
    int heightRows = max(1, static_cast<int>(ceilDiv(cell.h, siteH)));
    return static_cast<int>(max(0LL, min(static_cast<long long>(numRows - heightRows), rounded)));
}

int Legalizer::targetSite(const Instance &cell, int widthSites) const {
    long long rel = cell.x - dieLLX;
    long long rounded = floorDiv(rel + siteW / 2, siteW);
    return static_cast<int>(max(0LL, min(static_cast<long long>(numSites - widthSites), rounded)));
}

void Legalizer::buildRows() {
    rows.assign(numRows, Row());
    vector<vector<pair<int, int>>> blocked(numRows);

    for (const Instance &obs : obstacles) {
        long long ox1 = max(obs.x, dieLLX);
        long long oy1 = max(obs.y, dieLLY);
        long long ox2 = min(obs.x + obs.w, dieURX);
        long long oy2 = min(obs.y + obs.h, dieURY);
        if (ox1 >= ox2 || oy1 >= oy2) continue;

        int rowL = static_cast<int>(max(0LL, floorDiv(oy1 - dieLLY, siteH)));
        int rowR = static_cast<int>(min(static_cast<long long>(numRows), ceilDiv(oy2 - dieLLY, siteH)));
        int siteL = static_cast<int>(max(0LL, floorDiv(ox1 - dieLLX, siteW)));
        int siteR = static_cast<int>(min(static_cast<long long>(numSites), ceilDiv(ox2 - dieLLX, siteW)));
        if (siteL >= siteR) continue;
        for (int r = rowL; r < rowR; ++r) {
            long long rowY1 = dieLLY + static_cast<long long>(r) * siteH;
            long long rowY2 = rowY1 + siteH;
            if (rectOverlap(dieLLX + static_cast<long long>(siteL) * siteW, rowY1,
                            dieLLX + static_cast<long long>(siteR) * siteW, rowY2,
                            ox1, oy1, ox2, oy2)) {
                blocked[r].push_back({siteL, siteR});
            }
        }
    }

    for (int r = 0; r < numRows; ++r) {
        auto &b = blocked[r];
        sort(b.begin(), b.end());
        int cur = 0;
        for (auto [l, rr] : b) {
            l = max(l, 0);
            rr = min(rr, numSites);
            if (l > cur) rows[r].freeSites.push_back({cur, l});
            cur = max(cur, rr);
        }
        if (cur < numSites) rows[r].freeSites.push_back({cur, numSites});
    }
}

void Legalizer::buildSiteOwnerGrid() {
    siteOwner.assign(static_cast<size_t>(numRows) * numSites, -2);
    for (int r = 0; r < numRows; ++r) {
        for (const auto &seg : rows[r].freeSites) {
            for (int s = seg.l; s < seg.r; ++s) {
                siteOwner[static_cast<size_t>(r) * numSites + s] = -1;
            }
        }
    }
}

void Legalizer::buildDensityGrid() {
    gridSize = max(1LL, 10 * dbuPerMicron);
    gridLLX = dieLLX;
    gridLLY = dieLLY;
    gridURX = dieURX;
    gridURY = dieURY;

    long long dieW = dieURX - dieLLX;
    long long dieH = dieURY - dieLLY;
    constexpr double kEdgeBlockageSpanRatio = 0.999;
    for (const Instance &obs : obstacles) {
        if (obs.type != "BLOCKAGE") continue;
        long long x1 = max(obs.x, dieLLX);
        long long y1 = max(obs.y, dieLLY);
        long long x2 = min(obs.x + obs.w, dieURX);
        long long y2 = min(obs.y + obs.h, dieURY);
        if (x1 >= x2 || y1 >= y2) continue;

        bool spansDieWidth = static_cast<double>(x2 - x1) / dieW >= kEdgeBlockageSpanRatio;
        bool spansDieHeight = static_cast<double>(y2 - y1) / dieH >= kEdgeBlockageSpanRatio;
        if (y1 <= dieLLY && spansDieWidth) gridLLY = max(gridLLY, y2);
        if (y2 >= dieURY && spansDieWidth) gridURY = min(gridURY, y1);
        if (x1 <= dieLLX && spansDieHeight) gridLLX = max(gridLLX, x2);
        if (x2 >= dieURX && spansDieHeight) gridURX = min(gridURX, x1);
    }
    if (gridLLX >= gridURX || gridLLY >= gridURY) {
        gridLLX = dieLLX;
        gridLLY = dieLLY;
        gridURX = dieURX;
        gridURY = dieURY;
    }

    gridCols = static_cast<int>(ceilDiv(gridURX - gridLLX, gridSize));
    gridRows = static_cast<int>(ceilDiv(gridURY - gridLLY, gridSize));
    gridAreaUsed.assign(static_cast<size_t>(gridCols) * gridRows, 0.0);
    gridExcluded.assign(static_cast<size_t>(gridCols) * gridRows, 0);
}

bool Legalizer::snapInitialLegalize(double alpha) {
    (void)alpha;

    vector<int> order(cells.size());
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        const Instance &cell = cells[i];
        int widthSites = max(1, static_cast<int>(ceilDiv(cell.w, siteW)));
        int heightRows = max(1, static_cast<int>(ceilDiv(cell.h, siteH)));
        cellWidthSites[i] = widthSites;
        cellHeightRows[i] = heightRows;
        if (widthSites > numSites || heightRows > numRows) {
            cerr << "Cell is larger than the legal area: " << cell.name << '\n';
            return false;
        }
        if (heightRows != 1) return false;
        order[i] = i;
    }

    targetRows.resize(cells.size());
    targetSites.resize(cells.size());
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        targetRows[i] = targetRow(cells[i]);
        targetSites[i] = targetSite(cells[i], cellWidthSites[i]);
    }

    stable_sort(order.begin(), order.end(), [&](int a, int b) {
        if (targetSites[a] != targetSites[b]) return targetSites[a] < targetSites[b];
        if (targetRows[a] != targetRows[b]) return targetRows[a] < targetRows[b];
        if (cells[a].x != cells[b].x) return cells[a].x < cells[b].x;
        return a < b;
    });

    auto resetPlacement = [&]() {
        buildSiteOwnerGrid();
        gridAreaUsed.assign(static_cast<size_t>(gridCols) * gridRows, 0.0);
        placementX.assign(cells.size(), 0);
        placementY.assign(cells.size(), 0);
        placementRow.assign(cells.size(), -1);
        placementSite.assign(cells.size(), -1);
    };

    auto initRowStates = [&]() {
        vector<TetrisRowState> states(numRows);
        for (int r = 0; r < numRows; ++r) {
            states[r].interval = 0;
            states[r].cursor = rows[r].freeSites.empty() ? numSites : rows[r].freeSites.front().l;
        }
        return states;
    };

    auto maxRemainingRun = [&](const TetrisRowState &state, int row) {
        int best = 0;
        int cursor = state.cursor;
        for (size_t i = state.interval; i < rows[row].freeSites.size(); ++i) {
            const Interval &seg = rows[row].freeSites[i];
            int start = max(cursor, seg.l);
            if (start < seg.r) best = max(best, seg.r - start);
            cursor = 0;
        }
        return best;
    };

    auto firstTetrisFit = [&](const TetrisRowState &state, int row, int target, int widthSites,
                              bool preserveTargetGaps) {
        int bestSite = -1;
        int bestDistance = numeric_limits<int>::max();
        int cursor = state.cursor;
        for (size_t i = state.interval; i < rows[row].freeSites.size(); ++i) {
            const Interval &seg = rows[row].freeSites[i];
            int start = max(cursor, seg.l);
            if (start + widthSites <= seg.r) {
                int site = start;
                if (preserveTargetGaps) {
                    int latest = seg.r - widthSites;
                    site = max(start, min(target, latest));
                }
                int distance = abs(site - target);
                if (!preserveTargetGaps) return site;
                if (distance < bestDistance ||
                    (distance == bestDistance && (bestSite < 0 || site < bestSite))) {
                    bestDistance = distance;
                    bestSite = site;
                }
            }
            cursor = 0;
        }
        return bestSite;
    };

    auto commitTetrisFit = [&](TetrisRowState &state, int row, int site, int widthSites) {
        while (state.interval < rows[row].freeSites.size()) {
            const Interval &seg = rows[row].freeSites[state.interval];
            if (site >= seg.l && site + widthSites <= seg.r) {
                state.cursor = site + widthSites;
                while (state.interval < rows[row].freeSites.size() &&
                       state.cursor >= rows[row].freeSites[state.interval].r) {
                    ++state.interval;
                    if (state.interval < rows[row].freeSites.size()) {
                        state.cursor = max(state.cursor, rows[row].freeSites[state.interval].l);
                    }
                }
                return true;
            }
            ++state.interval;
            if (state.interval < rows[row].freeSites.size()) {
                state.cursor = rows[row].freeSites[state.interval].l;
            }
        }
        return false;
    };

    auto runTetris = [&](bool preserveTargetGaps) {
        resetPlacement();
        vector<TetrisRowState> rowStates = initRowStates();
        vector<int> rowMaxRun(numRows, 0);
        set<int> activeRows;
        for (int r = 0; r < numRows; ++r) {
            rowMaxRun[r] = maxRemainingRun(rowStates[r], r);
            if (rowMaxRun[r] > 0) activeRows.insert(r);
        }

        for (int idx : order) {
            int targetR = targetRows[idx];
            int targetS = targetSites[idx];
            int widthSites = cellWidthSites[idx];
            int bestRow = -1;
            int bestSite = -1;
            int bestRowDistance = numeric_limits<int>::max();
            int bestSiteDistance = numeric_limits<int>::max();
            double bestDisp = numeric_limits<double>::infinity();

            auto upper = activeRows.lower_bound(targetR);
            auto lower = upper;
            bool hasLower = false;
            if (lower != activeRows.begin()) {
                --lower;
                hasLower = true;
            }
            bool hasUpper = upper != activeRows.end();

            while (hasLower || hasUpper) {
                bool takeUpper = false;
                if (!hasLower) {
                    takeUpper = true;
                } else if (!hasUpper) {
                    takeUpper = false;
                } else {
                    takeUpper = abs(*upper - targetR) <= abs(*lower - targetR);
                }

                int row = takeUpper ? *upper : *lower;
                int rowDistance = abs(row - targetR);
                double minRowDisp = static_cast<double>(rowDistance) * siteH / dbuPerMicron;
                if (bestRow >= 0 && minRowDisp > bestDisp + 1e-9) break;

                if (takeUpper) {
                    ++upper;
                    hasUpper = upper != activeRows.end();
                } else if (lower == activeRows.begin()) {
                    hasLower = false;
                } else {
                    --lower;
                }

                if (rowMaxRun[row] < widthSites) continue;
                int site = firstTetrisFit(rowStates[row], row, targetS, widthSites, preserveTargetGaps);
                if (site < 0) continue;

                long long x = dieLLX + static_cast<long long>(site) * siteW;
                long long y = dieLLY + static_cast<long long>(row) * siteH;
                double disp = displacementCost(idx, x, y);
                int siteDistance = abs(site - targetS);
                if (disp < bestDisp - 1e-9 ||
                    (fabs(disp - bestDisp) < 1e-9 && rowDistance < bestRowDistance) ||
                    (fabs(disp - bestDisp) < 1e-9 && rowDistance == bestRowDistance &&
                     siteDistance < bestSiteDistance) ||
                    (fabs(disp - bestDisp) < 1e-9 && rowDistance == bestRowDistance &&
                     siteDistance == bestSiteDistance && row < bestRow)) {
                    bestDisp = disp;
                    bestRow = row;
                    bestSite = site;
                    bestRowDistance = rowDistance;
                    bestSiteDistance = siteDistance;
                }
            }

            if (bestRow < 0) return false;
            if (!commitTetrisFit(rowStates[bestRow], bestRow, bestSite, widthSites)) return false;
            restorePlacedCell(idx, bestRow, bestSite);

            rowMaxRun[bestRow] = maxRemainingRun(rowStates[bestRow], bestRow);
            if (rowMaxRun[bestRow] <= 0) activeRows.erase(bestRow);
        }
        return true;
    };

    if (runTetris(true)) return true;
    return runTetris(false);

}

bool Legalizer::abacusInitialLegalize(double alpha) {
    vector<int> order(cells.size());
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        const Instance &cell = cells[i];
        int widthSites = max(1, static_cast<int>(ceilDiv(cell.w, siteW)));
        int heightRows = max(1, static_cast<int>(ceilDiv(cell.h, siteH)));
        cellWidthSites[i] = widthSites;
        cellHeightRows[i] = heightRows;
        if (widthSites > numSites || heightRows > numRows) {
            cerr << "Cell is larger than the legal area: " << cell.name << '\n';
            return false;
        }
        order[i] = i;
    }

    targetRows.resize(cells.size());
    targetSites.resize(cells.size());
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        targetRows[i] = targetRow(cells[i]);
        targetSites[i] = targetSite(cells[i], cellWidthSites[i]);
    }

    stable_sort(order.begin(), order.end(), [&](int a, int b) {
        int ra = targetRows[a];
        int rb = targetRows[b];
        if (ra != rb) return ra < rb;
        if (cells[a].x != cells[b].x) return cells[a].x < cells[b].x;
        return cells[a].name < cells[b].name;
    });

    vector<int> singleRowCells;
    singleRowCells.reserve(cells.size());
    for (int idx : order) {
        if (cellHeightRows[idx] == 1) {
            singleRowCells.push_back(idx);
            continue;
        }

        Candidate cand = findBestCandidate(idx);
        if (cand.row < 0) {
            cerr << "Cannot legalize cell: " << cells[idx].name << '\n';
            return false;
        }

        occupy(cand.row, cand.site, cellWidthSites[idx], cellHeightRows[idx]);
        long long x = dieLLX + static_cast<long long>(cand.site) * siteW;
        long long y = dieLLY + static_cast<long long>(cand.row) * siteH;
        placementX[idx] = x;
        placementY[idx] = y;
        placementRow[idx] = cand.row;
        placementSite[idx] = cand.site;
        setOwner(idx, cand.row, cand.site, idx);
        addCellToDensity(x, y, cells[idx].w, cells[idx].h);
    }

    vector<AbacusSubrowState> subrows;
    vector<vector<int>> rowSubrows(numRows);
    for (int r = 0; r < numRows; ++r) {
        for (const Interval &seg : rows[r].freeSites) {
            if (seg.r <= seg.l) continue;
            int id = static_cast<int>(subrows.size());
            AbacusSubrowState state;
            state.row = r;
            state.interval = seg;
            subrows.push_back(std::move(state));
            rowSubrows[r].push_back(id);
        }
    }

    auto rebuildSingleRowByBacktracking = [&]() -> bool {
        constexpr int kBacktrackCellLimit = 128;
        constexpr size_t kBacktrackNodeLimit = 2000000;
        if (singleRowCells.empty()) return true;
        if (singleRowCells.size() > kBacktrackCellLimit || subrows.empty()) return false;

        vector<int> remaining(subrows.size(), 0);
        vector<vector<int>> assigned(subrows.size());
        for (int i = 0; i < static_cast<int>(subrows.size()); ++i) {
            AbacusSubrowState &state = subrows[i];
            state.cells.clear();
            state.usedSites = 0;
            state.cost = 0.0;
            remaining[i] = state.interval.r - state.interval.l;
        }

        vector<vector<int>> candidates(cells.size());
        for (int idx : singleRowCells) {
            int widthSites = cellWidthSites[idx];
            for (int subrowId = 0; subrowId < static_cast<int>(subrows.size()); ++subrowId) {
                if (remaining[subrowId] >= widthSites) candidates[idx].push_back(subrowId);
            }
            if (candidates[idx].empty()) return false;

            stable_sort(candidates[idx].begin(), candidates[idx].end(), [&](int a, int b) {
                const AbacusSubrowState &sa = subrows[a];
                const AbacusSubrowState &sb = subrows[b];
                int rowDistA = abs(sa.row - targetRows[idx]);
                int rowDistB = abs(sb.row - targetRows[idx]);
                if (rowDistA != rowDistB) return rowDistA < rowDistB;

                int siteA = max(sa.interval.l, min(targetSites[idx], sa.interval.r - widthSites));
                int siteB = max(sb.interval.l, min(targetSites[idx], sb.interval.r - widthSites));
                long long xa = dieLLX + static_cast<long long>(siteA) * siteW;
                long long ya = dieLLY + static_cast<long long>(sa.row) * siteH;
                long long xb = dieLLX + static_cast<long long>(siteB) * siteW;
                long long yb = dieLLY + static_cast<long long>(sb.row) * siteH;
                double dispA = displacementCost(idx, xa, ya);
                double dispB = displacementCost(idx, xb, yb);
                if (fabs(dispA - dispB) > 1e-9) return dispA < dispB;

                int slackA = remaining[a] - widthSites;
                int slackB = remaining[b] - widthSites;
                if (slackA != slackB) return slackA < slackB;
                return a < b;
            });
        }

        vector<int> order = singleRowCells;
        stable_sort(order.begin(), order.end(), [&](int a, int b) {
            if (cellWidthSites[a] != cellWidthSites[b]) return cellWidthSites[a] > cellWidthSites[b];
            if (candidates[a].size() != candidates[b].size()) return candidates[a].size() < candidates[b].size();
            int rowDistA = abs(subrows[candidates[a].front()].row - targetRows[a]);
            int rowDistB = abs(subrows[candidates[b].front()].row - targetRows[b]);
            if (rowDistA != rowDistB) return rowDistA > rowDistB;
            if (cells[a].x != cells[b].x) return cells[a].x < cells[b].x;
            return cells[a].name < cells[b].name;
        });

        vector<int> suffixWidth(order.size() + 1, 0);
        for (int i = static_cast<int>(order.size()) - 1; i >= 0; --i) {
            suffixWidth[i] = suffixWidth[i + 1] + cellWidthSites[order[i]];
        }

        size_t nodes = 0;
        function<bool(int)> search = [&](int pos) -> bool {
            if (pos == static_cast<int>(order.size())) return true;
            if (++nodes > kBacktrackNodeLimit) return false;

            int totalRemaining = 0;
            for (int cap : remaining) totalRemaining += cap;
            if (totalRemaining < suffixWidth[pos]) return false;

            int idx = order[pos];
            int widthSites = cellWidthSites[idx];
            for (int subrowId : candidates[idx]) {
                if (remaining[subrowId] < widthSites) continue;

                remaining[subrowId] -= widthSites;
                assigned[subrowId].push_back(idx);

                bool futureFeasible = true;
                for (int next = pos + 1; next < static_cast<int>(order.size()); ++next) {
                    int nextIdx = order[next];
                    bool fitsSomewhere = false;
                    for (int candidate : candidates[nextIdx]) {
                        if (remaining[candidate] >= cellWidthSites[nextIdx]) {
                            fitsSomewhere = true;
                            break;
                        }
                    }
                    if (!fitsSomewhere) {
                        futureFeasible = false;
                        break;
                    }
                }

                if (futureFeasible && search(pos + 1)) return true;

                assigned[subrowId].pop_back();
                remaining[subrowId] += widthSites;
            }
            return false;
        };

        if (!search(0)) return false;

        for (int subrowId = 0; subrowId < static_cast<int>(subrows.size()); ++subrowId) {
            AbacusSubrowState &state = subrows[subrowId];
            state.cells = std::move(assigned[subrowId]);
            state.usedSites = state.interval.r - state.interval.l - remaining[subrowId];
            vector<pair<int, int>> packed = packAbacusSubrow(state.cells, state.interval);
            state.cost = packedSubrowCost(packed, state.row, alpha);
        }
        return true;
    };

    bool exactInsertion = singleRowCells.size() <= 5000;
    bool rebuiltSingleRows = false;
    for (int idx : singleRowCells) {
        if (rebuiltSingleRows) break;
        int targetR = targetRows[idx];
        int targetS = targetSites[idx];
        int widthSites = cellWidthSites[idx];
        int maxWindow = max(0, numRows - 1);
        int currentWindow = min(maxWindow, max(8, numRows / 20));
        int bestSubrow = -1;
        double bestCost = 0.0;

        if (exactInsertion) {
            bool found = false;
            double bestDelta = numeric_limits<double>::infinity();
            double bestDisp = numeric_limits<double>::infinity();

            while (true) {
                for (int d = 0; d <= currentWindow; ++d) {
                    for (int sign : {0, -1, 1}) {
                        if (sign == 0 && d != 0) continue;
                        if (sign != 0 && d == 0) continue;
                        int row = targetR + sign * d;
                        if (row < 0 || row >= numRows) continue;

                        for (int subrowId : rowSubrows[row]) {
                            const AbacusSubrowState &state = subrows[subrowId];
                            if (state.usedSites + widthSites > state.interval.r - state.interval.l) continue;

                            vector<int> trial = state.cells;
                            trial.push_back(idx);
                            vector<pair<int, int>> packed = packAbacusSubrow(trial, state.interval);
                            double cost = packedSubrowCost(packed, row, alpha);
                            double delta = cost - state.cost;
                            double newDisp = numeric_limits<double>::infinity();
                            for (const auto &[cellIdx, site] : packed) {
                                if (cellIdx != idx) continue;
                                long long x = dieLLX + static_cast<long long>(site) * siteW;
                                long long y = dieLLY + static_cast<long long>(row) * siteH;
                                newDisp = displacementCost(idx, x, y);
                                break;
                            }

                            if (delta < bestDelta - 1e-9 ||
                                (fabs(delta - bestDelta) < 1e-9 && newDisp < bestDisp - 1e-9) ||
                                (fabs(delta - bestDelta) < 1e-9 && fabs(newDisp - bestDisp) < 1e-9 &&
                                 subrowId < bestSubrow)) {
                                found = true;
                                bestSubrow = subrowId;
                                bestDelta = delta;
                                bestCost = cost;
                                bestDisp = newDisp;
                            }
                        }
                    }
                }

                double minFutureDisp = static_cast<double>(currentWindow + 1) * siteH / dbuPerMicron;
                if ((found && minFutureDisp >= bestDisp - 1e-9) || currentWindow >= maxWindow) break;
                currentWindow = min(maxWindow, max(currentWindow + 1, currentWindow * 2));
            }
        } else {
            bool found = false;
            double bestScore = numeric_limits<double>::infinity();
            double bestDisp = numeric_limits<double>::infinity();
            int bestSiteDistance = numeric_limits<int>::max();

            while (true) {
                for (int d = 0; d <= currentWindow; ++d) {
                    for (int sign : {0, -1, 1}) {
                        if (sign == 0 && d != 0) continue;
                        if (sign != 0 && d == 0) continue;
                        int row = targetR + sign * d;
                        if (row < 0 || row >= numRows) continue;

                        for (int subrowId : rowSubrows[row]) {
                            const AbacusSubrowState &state = subrows[subrowId];
                            int capacity = state.interval.r - state.interval.l;
                            if (state.usedSites + widthSites > capacity) continue;

                            int site = max(state.interval.l, min(targetS, state.interval.r - widthSites));
                            long long x = dieLLX + static_cast<long long>(site) * siteW;
                            long long y = dieLLY + static_cast<long long>(row) * siteH;
                            double score = displacementCost(idx, x, y);
                            double disp = displacementCost(idx, x, y);
                            int siteDistance = abs(site - targetS);

                            if (score < bestScore - 1e-9 ||
                                (fabs(score - bestScore) < 1e-9 && disp < bestDisp - 1e-9) ||
                                (fabs(score - bestScore) < 1e-9 && fabs(disp - bestDisp) < 1e-9 &&
                                 siteDistance < bestSiteDistance) ||
                                (fabs(score - bestScore) < 1e-9 && fabs(disp - bestDisp) < 1e-9 &&
                                 siteDistance == bestSiteDistance && subrowId < bestSubrow)) {
                                found = true;
                                bestSubrow = subrowId;
                                bestScore = score;
                                bestDisp = disp;
                                bestSiteDistance = siteDistance;
                            }
                        }
                    }
                }

                double minFutureDisp = static_cast<double>(currentWindow + 1) * siteH / dbuPerMicron;
                if ((found && minFutureDisp >= bestDisp - 1e-9) || currentWindow >= maxWindow) break;
                currentWindow = min(maxWindow, max(currentWindow + 1, currentWindow * 2));
            }
        }

        if (bestSubrow < 0) {
            if (rebuildSingleRowByBacktracking()) {
                rebuiltSingleRows = true;
                break;
            }
            cerr << "Cannot legalize cell: " << cells[idx].name << '\n';
            return false;
        }

        AbacusSubrowState &state = subrows[bestSubrow];
        state.cells.push_back(idx);
        state.usedSites += widthSites;
        if (exactInsertion) state.cost = bestCost;
    }

    for (const AbacusSubrowState &state : subrows) {
        vector<pair<int, int>> packed = packAbacusSubrow(state.cells, state.interval);
        for (const auto &[idx, site] : packed) {
            long long x = dieLLX + static_cast<long long>(site) * siteW;
            long long y = dieLLY + static_cast<long long>(state.row) * siteH;
            placementX[idx] = x;
            placementY[idx] = y;
            placementRow[idx] = state.row;
            placementSite[idx] = site;
            setOwner(idx, state.row, site, idx);
            addCellToDensity(x, y, cells[idx].w, cells[idx].h);
        }
    }

    return true;
}

vector<pair<int, int>> Legalizer::packAbacusSubrow(const vector<int> &cellIdxs,
                                                   const Interval &interval) const {
    vector<int> ordered = cellIdxs;
    stable_sort(ordered.begin(), ordered.end(), [&](int a, int b) {
        int ta = targetSites[a];
        int tb = targetSites[b];
        if (ta != tb) return ta < tb;
        if (cells[a].x != cells[b].x) return cells[a].x < cells[b].x;
        return cells[a].name < cells[b].name;
    });

    auto placeCluster = [&](AbacusCluster &cluster) {
        int maxStart = interval.r - cluster.width;
        vector<int> ideals;
        ideals.reserve(cluster.cells.size());
        int offset = 0;
        for (int idx : cluster.cells) {
            ideals.push_back(targetSites[idx] - offset);
            offset += cellWidthSites[idx];
        }
        size_t mid = ideals.size() / 2;
        nth_element(ideals.begin(), ideals.begin() + mid, ideals.end());
        int ideal = ideals[mid];
        cluster.x = max(interval.l, min(maxStart, ideal));
    };

    vector<AbacusCluster> clusters;
    clusters.reserve(ordered.size());
    for (int idx : ordered) {
        AbacusCluster cluster;
        cluster.width = cellWidthSites[idx];
        cluster.weight = 1;
        cluster.q = targetSites[idx];
        cluster.cells.push_back(idx);
        placeCluster(cluster);
        clusters.push_back(std::move(cluster));

        while (clusters.size() >= 2) {
            AbacusCluster current = std::move(clusters.back());
            clusters.pop_back();
            AbacusCluster previous = std::move(clusters.back());
            clusters.pop_back();

            if (previous.x + previous.width <= current.x) {
                clusters.push_back(std::move(previous));
                clusters.push_back(std::move(current));
                break;
            }

            AbacusCluster merged;
            merged.width = previous.width + current.width;
            merged.weight = previous.weight + current.weight;
            merged.q = previous.q + current.q -
                       static_cast<long long>(current.weight) * previous.width;
            merged.cells = std::move(previous.cells);
            merged.cells.insert(merged.cells.end(), current.cells.begin(), current.cells.end());
            placeCluster(merged);
            clusters.push_back(std::move(merged));
        }
    }

    vector<pair<int, int>> packed;
    packed.reserve(ordered.size());
    for (const AbacusCluster &cluster : clusters) {
        int site = cluster.x;
        for (int idx : cluster.cells) {
            packed.push_back({idx, site});
            site += cellWidthSites[idx];
        }
    }
    return packed;
}

double Legalizer::packedSubrowCost(const vector<pair<int, int>> &packed,
                                   int row, double alpha) const {
    (void)alpha;
    double total = 0.0;
    long long y = dieLLY + static_cast<long long>(row) * siteH;
    for (const auto &[idx, site] : packed) {
        long long x = dieLLX + static_cast<long long>(site) * siteW;
        total += displacementCost(idx, x, y);
    }
    return total;
}

vector<int> Legalizer::candidateSitesInRow(int row, int target, int widthSites) const {
    vector<int> result;
    const auto &iv = rows[row].freeSites;
    if (iv.empty()) return result;

    auto addSite = [&](int s) {
        if (s < 0 || s + widthSites > numSites) return;
        for (int old : result) {
            if (old == s) return;
        }
        result.push_back(s);
    };

    auto addNearestForTarget = [&](int localTarget, int limit) {
        localTarget = max(0, min(numSites - widthSites, localTarget));
        auto it = lower_bound(iv.begin(), iv.end(), localTarget,
                              [](const Interval &a, int value) { return a.l < value; });
        int pos = static_cast<int>(it - iv.begin());

        int bestDx = numeric_limits<int>::max();
        for (int i = pos; i < static_cast<int>(iv.size()); ++i) {
            if (iv[i].r - iv[i].l < widthSites) continue;
            if (iv[i].l > localTarget + bestDx) break;
            int s = max(iv[i].l, min(localTarget, iv[i].r - widthSites));
            int dx = abs(s - localTarget);
            if (dx < bestDx) bestDx = dx;
            addSite(s);
            addSite(iv[i].l);
            addSite(iv[i].r - widthSites);
            if (static_cast<int>(result.size()) >= limit) break;
        }

        for (int i = pos - 1; i >= 0; --i) {
            if (iv[i].r - iv[i].l < widthSites) continue;
            int maxStart = iv[i].r - widthSites;
            if (localTarget - maxStart > bestDx && !result.empty()) break;
            int s = max(iv[i].l, min(localTarget, maxStart));
            int dx = abs(s - localTarget);
            if (dx < bestDx) bestDx = dx;
            addSite(s);
            addSite(iv[i].l);
            addSite(maxStart);
            if (static_cast<int>(result.size()) >= limit) break;
        }
    };

    int limit = 3;
    addNearestForTarget(target, limit);

    int gridSites = max(1, static_cast<int>(gridSize / siteW));
    for (int k = 1; k <= 8 && static_cast<int>(result.size()) < limit; ++k) {
        addNearestForTarget(target - k * gridSites, limit);
        addNearestForTarget(target + k * gridSites, limit);
    }

    if (result.empty()) {
        for (const auto &seg : iv) {
            if (seg.r - seg.l >= widthSites) {
                addSite(max(seg.l, min(target, seg.r - widthSites)));
                break;
            }
        }
    }

    return result;
}

bool Legalizer::canPlaceMultiRow(int row, int site, int widthSites, int heightRows) const {
    for (int rr = row; rr < row + heightRows; ++rr) {
        bool covered = false;
        for (const auto &seg : rows[rr].freeSites) {
            if (seg.l <= site && site + widthSites <= seg.r) {
                covered = true;
                break;
            }
        }
        if (!covered) return false;
    }
    return true;
}

Candidate Legalizer::findBestCandidate(int cellIdx) const {
    int widthSites = cellWidthSites[cellIdx];
    int heightRows = cellHeightRows[cellIdx];
    int tRow = targetRows[cellIdx];
    int tSite = targetSites[cellIdx];

    int baseWindow = 64;

    Candidate best;
    int maxStartRow = numRows - heightRows;
    int currentWindow = min(maxStartRow, baseWindow);
    bool found = false;

    while (!found) {
        for (int d = 0; d <= currentWindow; ++d) {
            for (int sign : {0, -1, 1}) {
                if (sign == 0 && d != 0) continue;
                if (sign != 0 && d == 0) continue;
                int row = tRow + sign * d;
                if (row < 0 || row > maxStartRow) continue;

                vector<int> sites = candidateSitesInRow(row, tSite, widthSites);
                for (int site : sites) {
                    if (!canPlaceMultiRow(row, site, widthSites, heightRows)) continue;
                    long long x = dieLLX + static_cast<long long>(site) * siteW;
                    long long y = dieLLY + static_cast<long long>(row) * siteH;
                    double disp = displacementCost(cellIdx, x, y);
                    double cost = disp;
                    if (cost < best.cost || (fabs(cost - best.cost) < 1e-9 && disp < best.displacement)) {
                        best = {row, site, cost, disp};
                        found = true;
                    }
                }
            }
        }

        double minFutureDisp = static_cast<double>(currentWindow + 1) * siteH / dbuPerMicron;
        if ((found && minFutureDisp >= best.displacement - 1e-9) || currentWindow >= maxStartRow) break;
        currentWindow = min(maxStartRow, max(currentWindow + 1, currentWindow * 2));
    }

    return best;
}

bool Legalizer::ownerCanPlace(int cellIdx, int row, int site) const {
    (void)cellIdx;
    int widthSites = cellWidthSites[cellIdx];
    int heightRows = cellHeightRows[cellIdx];
    if (row < 0 || site < 0 || row + heightRows > numRows || site + widthSites > numSites) {
        return false;
    }

    for (int rr = row; rr < row + heightRows; ++rr) {
        for (int ss = site; ss < site + widthSites; ++ss) {
            if (siteOwner[static_cast<size_t>(rr) * numSites + ss] != -1) return false;
        }
    }
    return true;
}

void Legalizer::setOwner(int cellIdx, int row, int site, int value) {
    int widthSites = cellWidthSites[cellIdx];
    int heightRows = cellHeightRows[cellIdx];
    for (int rr = row; rr < row + heightRows; ++rr) {
        for (int ss = site; ss < site + widthSites; ++ss) {
            siteOwner[static_cast<size_t>(rr) * numSites + ss] = value;
        }
    }
}

void Legalizer::removePlacedCell(int cellIdx) {
    setOwner(cellIdx, placementRow[cellIdx], placementSite[cellIdx], -1);
    addCellToDensityDelta(placementX[cellIdx], placementY[cellIdx],
                          cells[cellIdx].w, cells[cellIdx].h, -1.0);
}

void Legalizer::restorePlacedCell(int cellIdx, int row, int site) {
    long long x = dieLLX + static_cast<long long>(site) * siteW;
    long long y = dieLLY + static_cast<long long>(row) * siteH;
    placementRow[cellIdx] = row;
    placementSite[cellIdx] = site;
    placementX[cellIdx] = x;
    placementY[cellIdx] = y;
    setOwner(cellIdx, row, site, cellIdx);
    addCellToDensityDelta(x, y, cells[cellIdx].w, cells[cellIdx].h, 1.0);
}

double Legalizer::displacementCost(int cellIdx, long long x, long long y) const {
    double dx = static_cast<double>(x - cells[cellIdx].x) / dbuPerMicron;
    double dy = static_cast<double>(y - cells[cellIdx].y) / dbuPerMicron;
    return fabs(dx) + fabs(dy);
}

double Legalizer::placementCost(int cellIdx, long long x, long long y, double alpha) const {
    double a = max(0.0, min(1.0, alpha));
    double disp = displacementCost(cellIdx, x, y);
    double densityPenalty = estimateDensityPenalty(x, y, cells[cellIdx].w, cells[cellIdx].h);
    return a * kDisplacementNormFactor * disp + (1.0 - a) * densityPenalty;
}

double Legalizer::averageDisplacement() const {
    if (cells.empty()) return 0.0;
    double total = 0.0;
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        total += displacementCost(i, placementX[i], placementY[i]);
    }
    return total / static_cast<double>(cells.size());
}

double Legalizer::densityOverflowRatio() const {
    if (gridCols <= 0 || gridRows <= 0) return 0.0;

    double threshold = densityThreshold;
    vector<unsigned char> touched(static_cast<size_t>(gridCols) * gridRows, 0);
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        if (!cells[i].name.empty() && cells[i].name[0] == 'h') continue;

        long long x1 = placementX[i];
        long long y1 = placementY[i];
        long long x2 = x1 + cells[i].w;
        long long y2 = y1 + cells[i].h;
        if (x2 < gridLLX || x1 > gridURX || y2 < gridLLY || y1 > gridURY) continue;

        int c0 = max(gridLineLowerBound(x1, gridLLX, gridURX, gridSize, gridCols) - 1, 0);
        int c1 = min(gridLineUpperBound(x2, gridLLX, gridURX, gridSize, gridCols), gridCols);
        int r0 = max(gridLineLowerBound(y1, gridLLY, gridURY, gridSize, gridRows) - 1, 0);
        int r1 = min(gridLineUpperBound(y2, gridLLY, gridURY, gridSize, gridRows), gridRows);
        for (int r = r0; r < r1; ++r) {
            for (int c = c0; c < c1; ++c) {
                touched[static_cast<size_t>(r) * gridCols + c] = 1;
            }
        }
    }

    int overflow = 0;
    int total = 0;
    for (int r = 0; r < gridRows; ++r) {
        long long gy1 = gridLLY + static_cast<long long>(r) * gridSize;
        long long gy2 = min(gy1 + gridSize, gridURY);
        for (int c = 0; c < gridCols; ++c) {
            size_t idx = static_cast<size_t>(r) * gridCols + c;
            if (gridExcluded[idx]) continue;
            if (!touched[idx]) continue;

            long long gx1 = gridLLX + static_cast<long long>(c) * gridSize;
            long long gx2 = min(gx1 + gridSize, gridURX);
            double gridArea = static_cast<double>(gx2 - gx1) * static_cast<double>(gy2 - gy1);
            if (gridArea <= 0.0) continue;

            ++total;
            if (gridAreaUsed[idx] / gridArea > threshold) ++overflow;
        }
    }

    return total == 0 ? 0.0 : static_cast<double>(overflow) / total * 100.0;
}

double Legalizer::qualityScore(double alpha) const {
    double a = max(0.0, min(1.0, alpha));
    return a * kDisplacementNormFactor * averageDisplacement() + (1.0 - a) * densityOverflowRatio();
}

bool Legalizer::runtimeLimitReached() {
    if (chrono::steady_clock::now() < runtimeDeadline) return false;
    runtimeLimitHit = true;
    return true;
}

void Legalizer::saveBestPlacement(double quality) {
    bestQuality = quality;
    bestPlacementX = placementX;
    bestPlacementY = placementY;
    bestPlacementRow = placementRow;
    bestPlacementSite = placementSite;
    bestSiteOwner = siteOwner;
    bestGridAreaUsed = gridAreaUsed;
}

bool Legalizer::updateBestPlacement(double alpha) {
    double currentQuality = qualityScore(alpha);
    if (currentQuality >= bestQuality - 1e-9) return false;
    saveBestPlacement(currentQuality);
    return true;
}

void Legalizer::restoreBestPlacement() {
    if (bestPlacementX.empty()) return;
    placementX = bestPlacementX;
    placementY = bestPlacementY;
    placementRow = bestPlacementRow;
    placementSite = bestPlacementSite;
    siteOwner = bestSiteOwner;
    gridAreaUsed = bestGridAreaUsed;
}

bool Legalizer::tryMoveCell(int cellIdx, int rowRadius, int siteRadius, double alpha) {
    int widthSites = cellWidthSites[cellIdx];
    int heightRows = cellHeightRows[cellIdx];
    int oldRow = placementRow[cellIdx];
    int oldSite = placementSite[cellIdx];
    long long oldX = placementX[cellIdx];
    long long oldY = placementY[cellIdx];

    removePlacedCell(cellIdx);

    int targetR = targetRows[cellIdx];
    int targetS = targetSites[cellIdx];
    int maxStartRow = numRows - heightRows;
    int maxStartSite = numSites - widthSites;

    double bestCost = placementCost(cellIdx, oldX, oldY, alpha);
    double bestDisp = displacementCost(cellIdx, oldX, oldY);
    double oldDensityPenalty = estimateDensityPenalty(oldX, oldY, cells[cellIdx].w, cells[cellIdx].h);
    int bestRow = oldRow;
    int bestSite = oldSite;

    vector<int> rowCenters = {targetR, oldRow};
    vector<int> siteCenters = {targetS, oldSite};

    vector<int> rowsToTry;
    for (int center : rowCenters) {
        for (int d = -rowRadius; d <= rowRadius; ++d) {
            int r = center + d;
            if (r < 0 || r > maxStartRow) continue;
            if (find(rowsToTry.begin(), rowsToTry.end(), r) == rowsToTry.end()) rowsToTry.push_back(r);
        }
    }

    for (int r : rowsToTry) {
        vector<int> sitesToTry;
        for (int center : siteCenters) {
            int l = max(0, center - siteRadius);
            int rr = min(maxStartSite, center + siteRadius);
            for (int s = l; s <= rr; ++s) {
                sitesToTry.push_back(s);
            }
        }

        for (int s : sitesToTry) {
            if (!ownerCanPlace(cellIdx, r, s)) continue;
            long long x = dieLLX + static_cast<long long>(s) * siteW;
            long long y = dieLLY + static_cast<long long>(r) * siteH;
            double densityPenalty = estimateDensityPenalty(x, y, cells[cellIdx].w, cells[cellIdx].h);
            if (densityPenalty > oldDensityPenalty + 1e-9) continue;
            double disp = displacementCost(cellIdx, x, y);
            double a = max(0.0, min(1.0, alpha));
            double cost = a * kDisplacementNormFactor * disp + (1.0 - a) * densityPenalty;
            if (cost < bestCost - 1e-9 ||
                (fabs(cost - bestCost) < 1e-9 && disp < bestDisp - 1e-9)) {
                bestCost = cost;
                bestDisp = disp;
                bestRow = r;
                bestSite = s;
            }
        }
    }

    restorePlacedCell(cellIdx, bestRow, bestSite);
    return bestRow != oldRow || bestSite != oldSite;
}

void Legalizer::detailPlace(double alpha) {
    int gridSites = max(1, static_cast<int>(ceilDiv(gridSize, siteW)));
    int rowRadius = min(4, max(0, numRows - 1));
    int siteRadius = min(max(12, 2 * gridSites), max(12, numSites / 20));
    int passes = 2;
    int maxCellsToVisit = static_cast<int>(cells.size());

    if (cells.size() <= 5000) {
        siteRadius = max(0, numSites - 1);
        passes = 4;
    } else if (cells.size() > 50000) {
        rowRadius = min(rowRadius, 4);
        siteRadius = min(max(16, (11 * gridSites) / 2), max(16, numSites / 13));
        passes = 3;
        maxCellsToVisit = static_cast<int>(cells.size());
    }

    vector<int> order(cells.size());
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) order[i] = i;

    for (int pass = 0; pass < passes; ++pass) {
        if (runtimeLimitReached()) {
            updateBestPlacement(alpha);
            return;
        }

        stable_sort(order.begin(), order.end(), [&](int a, int b) {
            double da = displacementCost(a, placementX[a], placementY[a]);
            double db = displacementCost(b, placementX[b], placementY[b]);
            if (fabs(da - db) > 1e-9) return da > db;
            return cells[a].name < cells[b].name;
        });

        int moved = 0;
        for (int i = 0; i < maxCellsToVisit; ++i) {
            if (tryMoveCell(order[i], rowRadius, siteRadius, alpha)) ++moved;
            if ((i & 1023) == 1023 && runtimeLimitReached()) {
                updateBestPlacement(alpha);
                return;
            }
        }
        updateBestPlacement(alpha);
        if (moved == 0) break;
    }
}

void Legalizer::refineCells(double alpha) {
    for (int pass = 0; pass < 3; ++pass) {
        if (runtimeLimitReached()) {
            updateBestPlacement(alpha);
            return;
        }

        bool changed = false;

        for (int r = 0; r < numRows; ++r) {
            if ((r & 63) == 63 && runtimeLimitReached()) {
                updateBestPlacement(alpha);
                return;
            }

            vector<int> rowCells;
            int last = -3;
            for (int s = 0; s < numSites; ++s) {
                int owner = siteOwner[static_cast<size_t>(r) * numSites + s];
                if (owner >= 0 && owner != last && placementRow[owner] == r) {
                    rowCells.push_back(owner);
                }
                last = owner;
            }

            for (size_t i = 1; i < rowCells.size(); ++i) {
                int a = rowCells[i - 1];
                int b = rowCells[i];
                if (cellWidthSites[a] != cellWidthSites[b] ||
                    cellHeightRows[a] != cellHeightRows[b]) {
                    continue;
                }

                double oldDisp = displacementCost(a, placementX[a], placementY[a]) +
                                 displacementCost(b, placementX[b], placementY[b]);
                double newDisp = displacementCost(a, placementX[b], placementY[b]) +
                                 displacementCost(b, placementX[a], placementY[a]);
                if (newDisp >= oldDisp - 1e-9) continue;

                int aRow = placementRow[a], aSite = placementSite[a];
                int bRow = placementRow[b], bSite = placementSite[b];
                long long aX = placementX[a], aY = placementY[a];
                long long bX = placementX[b], bY = placementY[b];

                setOwner(a, aRow, aSite, -1);
                setOwner(b, bRow, bSite, -1);

                placementRow[a] = bRow;
                placementSite[a] = bSite;
                placementX[a] = bX;
                placementY[a] = bY;
                placementRow[b] = aRow;
                placementSite[b] = aSite;
                placementX[b] = aX;
                placementY[b] = aY;

                setOwner(a, placementRow[a], placementSite[a], a);
                setOwner(b, placementRow[b], placementSite[b], b);
                changed = true;
            }
        }

        updateBestPlacement(alpha);
        if (!changed) break;
    }
}

void Legalizer::refinePairSwaps(double alpha) {
    if (cells.size() > 5000) return;

    for (int pass = 0; pass < 4; ++pass) {
        bool changed = false;

        for (int a = 0; a < static_cast<int>(cells.size()); ++a) {
            for (int b = a + 1; b < static_cast<int>(cells.size()); ++b) {
                if (cellWidthSites[a] != cellWidthSites[b] ||
                    cellHeightRows[a] != cellHeightRows[b]) {
                    continue;
                }

                int aRow = placementRow[a], aSite = placementSite[a];
                int bRow = placementRow[b], bSite = placementSite[b];
                long long aX = placementX[a], aY = placementY[a];
                long long bX = placementX[b], bY = placementY[b];

                double oldCost = displacementCost(a, aX, aY) + displacementCost(b, bX, bY);
                double newCost = displacementCost(a, bX, bY) + displacementCost(b, aX, aY);
                if (newCost >= oldCost - 1e-9) continue;

                removePlacedCell(a);
                removePlacedCell(b);

                bool canSwap = ownerCanPlace(a, bRow, bSite) && ownerCanPlace(b, aRow, aSite);
                if (canSwap) {
                    restorePlacedCell(a, bRow, bSite);
                    restorePlacedCell(b, aRow, aSite);
                    changed = true;
                } else {
                    restorePlacedCell(a, aRow, aSite);
                    restorePlacedCell(b, bRow, bSite);
                }
            }
        }

        updateBestPlacement(alpha);
        if (!changed) break;
    }
}

void Legalizer::refineSmallReinsert(double alpha) {
    if (cells.size() > 128) return;

    vector<Interval> rowIntervals(numRows);
    for (int r = 0; r < numRows; ++r) {
        if (rows[r].freeSites.size() != 1) return;
        rowIntervals[r] = rows[r].freeSites.front();
    }
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        if (cellHeightRows[i] != 1) return;
    }

    vector<vector<int>> rowCells(numRows);
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        if (placementRow[i] < 0 || placementRow[i] >= numRows) return;
        rowCells[placementRow[i]].push_back(i);
    }
    auto sortByPlacement = [&](vector<int> &items) {
        stable_sort(items.begin(), items.end(), [&](int a, int b) {
            if (placementSite[a] != placementSite[b]) return placementSite[a] < placementSite[b];
            return cells[a].name < cells[b].name;
        });
    };
    for (auto &items : rowCells) sortByPlacement(items);

    auto rowPackedCost = [&](int row, const vector<int> &items, vector<pair<int, int>> *packedOut) {
        int usedSites = 0;
        for (int idx : items) usedSites += cellWidthSites[idx];
        const Interval &interval = rowIntervals[row];
        if (usedSites > interval.r - interval.l) return numeric_limits<double>::infinity();

        vector<pair<int, int>> packed = packAbacusSubrow(items, interval);
        double cost = 0.0;
        long long y = dieLLY + static_cast<long long>(row) * siteH;
        for (const auto &[idx, site] : packed) {
            long long x = dieLLX + static_cast<long long>(site) * siteW;
            cost += displacementCost(idx, x, y);
        }
        if (packedOut) *packedOut = std::move(packed);
        return cost;
    };

    auto rowCurrentCost = [&](const vector<int> &items) {
        double cost = 0.0;
        for (int idx : items) cost += displacementCost(idx, placementX[idx], placementY[idx]);
        return cost;
    };

    auto withoutCell = [](const vector<int> &items, int idx) {
        vector<int> result;
        result.reserve(items.size());
        for (int item : items) {
            if (item != idx) result.push_back(item);
        }
        return result;
    };

    for (int pass = 0; pass < 32; ++pass) {
        double bestGain = 1e-9;
        int bestCell = -1;
        int bestTargetRow = -1;
        vector<int> bestOldRowCells;
        vector<int> bestNewRowCells;

        for (int idx = 0; idx < static_cast<int>(cells.size()); ++idx) {
            int oldRow = placementRow[idx];
            for (int row = 0; row < numRows; ++row) {
                if (row == oldRow) continue;

                vector<int> oldTrial = withoutCell(rowCells[oldRow], idx);
                vector<int> newTrial = rowCells[row];
                newTrial.push_back(idx);

                double oldCost = rowCurrentCost(rowCells[oldRow]) + rowCurrentCost(rowCells[row]);
                double newCost = rowPackedCost(oldRow, oldTrial, nullptr) +
                                 rowPackedCost(row, newTrial, nullptr);
                double gain = oldCost - newCost;
                if (gain > bestGain) {
                    bestGain = gain;
                    bestCell = idx;
                    bestTargetRow = row;
                    bestOldRowCells = std::move(oldTrial);
                    bestNewRowCells = std::move(newTrial);
                }
            }
        }

        if (bestCell < 0) break;

        int oldRow = placementRow[bestCell];
        vector<pair<int, int>> oldPacked;
        vector<pair<int, int>> newPacked;
        rowPackedCost(oldRow, bestOldRowCells, &oldPacked);
        rowPackedCost(bestTargetRow, bestNewRowCells, &newPacked);

        vector<int> affected = rowCells[oldRow];
        affected.insert(affected.end(), rowCells[bestTargetRow].begin(), rowCells[bestTargetRow].end());
        sort(affected.begin(), affected.end());
        affected.erase(unique(affected.begin(), affected.end()), affected.end());
        for (int idx : affected) removePlacedCell(idx);
        for (const auto &[idx, site] : oldPacked) restorePlacedCell(idx, oldRow, site);
        for (const auto &[idx, site] : newPacked) restorePlacedCell(idx, bestTargetRow, site);

        rowCells[oldRow] = std::move(bestOldRowCells);
        rowCells[bestTargetRow] = std::move(bestNewRowCells);
        sortByPlacement(rowCells[oldRow]);
        sortByPlacement(rowCells[bestTargetRow]);
        updateBestPlacement(alpha);
    }

    for (int pass = 0; pass < 32; ++pass) {
        double bestGain = 1e-9;
        int bestA = -1;
        int bestB = -1;
        vector<int> bestARowCells;
        vector<int> bestBRowCells;

        for (int a = 0; a < static_cast<int>(cells.size()); ++a) {
            int aRow = placementRow[a];
            for (int b = a + 1; b < static_cast<int>(cells.size()); ++b) {
                int bRow = placementRow[b];
                if (aRow == bRow) continue;

                vector<int> aTrial = withoutCell(rowCells[aRow], a);
                vector<int> bTrial = withoutCell(rowCells[bRow], b);
                aTrial.push_back(b);
                bTrial.push_back(a);

                double oldCost = rowCurrentCost(rowCells[aRow]) + rowCurrentCost(rowCells[bRow]);
                double newCost = rowPackedCost(aRow, aTrial, nullptr) +
                                 rowPackedCost(bRow, bTrial, nullptr);
                double gain = oldCost - newCost;
                if (gain > bestGain) {
                    bestGain = gain;
                    bestA = a;
                    bestB = b;
                    bestARowCells = std::move(aTrial);
                    bestBRowCells = std::move(bTrial);
                }
            }
        }

        if (bestA < 0) break;

        int aRow = placementRow[bestA];
        int bRow = placementRow[bestB];
        vector<pair<int, int>> aPacked;
        vector<pair<int, int>> bPacked;
        rowPackedCost(aRow, bestARowCells, &aPacked);
        rowPackedCost(bRow, bestBRowCells, &bPacked);

        vector<int> affected = rowCells[aRow];
        affected.insert(affected.end(), rowCells[bRow].begin(), rowCells[bRow].end());
        sort(affected.begin(), affected.end());
        affected.erase(unique(affected.begin(), affected.end()), affected.end());
        for (int idx : affected) removePlacedCell(idx);
        for (const auto &[idx, site] : aPacked) restorePlacedCell(idx, aRow, site);
        for (const auto &[idx, site] : bPacked) restorePlacedCell(idx, bRow, site);

        rowCells[aRow] = std::move(bestARowCells);
        rowCells[bRow] = std::move(bestBRowCells);
        sortByPlacement(rowCells[aRow]);
        sortByPlacement(rowCells[bRow]);
        updateBestPlacement(alpha);
    }
}

void Legalizer::refineDisplacementSwaps(double alpha) {
    int rowSpan = cells.size() > 50000 ? min(30, max(0, numRows - 1)) : max(0, numRows - 1);
    int siteSpan = cells.size() > 50000 ? min(max(240, numSites / 25), max(0, numSites - 1))
                                        : max(0, numSites - 1);
    int passes = cells.size() > 50000 ? 2 : 4;
    int maxCellsToVisit = cells.size() > 50000 ? min(static_cast<int>(cells.size()), 60000)
                                               : static_cast<int>(cells.size());

    vector<int> order(cells.size());
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) order[i] = i;

    for (int pass = 0; pass < passes; ++pass) {
        if (runtimeLimitReached()) {
            updateBestPlacement(alpha);
            return;
        }

        stable_sort(order.begin(), order.end(), [&](int a, int b) {
            double da = displacementCost(a, placementX[a], placementY[a]);
            double db = displacementCost(b, placementX[b], placementY[b]);
            if (fabs(da - db) > 1e-9) return da > db;
            return cells[a].name < cells[b].name;
        });

        int swapped = 0;
        for (int i = 0; i < maxCellsToVisit; ++i) {
            int idx = order[i];
            int targetR = targetRows[idx];
            int targetS = targetSites[idx];
            double oldIdxCost = displacementCost(idx, placementX[idx], placementY[idx]);

            int bestOther = -1;
            double bestGain = 1e-9;

            int r0 = max(0, targetR - rowSpan);
            int r1 = min(numRows - cellHeightRows[idx], targetR + rowSpan);
            int s0 = max(0, targetS - siteSpan);
            int s1 = min(numSites - cellWidthSites[idx], targetS + siteSpan);

            int lastOwner = -3;
            for (int r = r0; r <= r1; ++r) {
                for (int s = s0; s <= s1; ++s) {
                    int other = siteOwner[static_cast<size_t>(r) * numSites + s];
                    if (other < 0 || other == idx || other == lastOwner) {
                        lastOwner = other;
                        continue;
                    }
                    lastOwner = other;
                    if (cellWidthSites[other] != cellWidthSites[idx] ||
                        cellHeightRows[other] != cellHeightRows[idx]) {
                        continue;
                    }

                    double oldCost = oldIdxCost +
                                     displacementCost(other, placementX[other], placementY[other]);
                    double newCost = displacementCost(idx, placementX[other], placementY[other]) +
                                     displacementCost(other, placementX[idx], placementY[idx]);
                    double gain = oldCost - newCost;
                    if (gain > bestGain) {
                        bestGain = gain;
                        bestOther = other;
                    }
                }
            }

            if (bestOther < 0) continue;

            int a = idx;
            int b = bestOther;
            int aRow = placementRow[a], aSite = placementSite[a];
            int bRow = placementRow[b], bSite = placementSite[b];
            long long aX = placementX[a], aY = placementY[a];
            long long bX = placementX[b], bY = placementY[b];

            setOwner(a, aRow, aSite, -1);
            setOwner(b, bRow, bSite, -1);

            placementRow[a] = bRow;
            placementSite[a] = bSite;
            placementX[a] = bX;
            placementY[a] = bY;
            placementRow[b] = aRow;
            placementSite[b] = aSite;
            placementX[b] = aX;
            placementY[b] = aY;

            setOwner(a, placementRow[a], placementSite[a], a);
            setOwner(b, placementRow[b], placementSite[b], b);
            ++swapped;

            if ((swapped & 1023) == 1023 && runtimeLimitReached()) {
                updateBestPlacement(alpha);
                return;
            }
        }

        updateBestPlacement(alpha);
        if (swapped == 0) break;
    }
}

void Legalizer::refineRowSliding(double alpha) {
    int passes = cells.size() > 50000 ? 5 : 8;
    int maxBlock = cells.size() > 50000 ? 128 : 512;

    auto moveBlock = [&](const vector<int> &rowCells, int lo, int hi, int deltaSites) {
        for (int p = lo; p <= hi; ++p) removePlacedCell(rowCells[p]);
        for (int p = lo; p <= hi; ++p) {
            int idx = rowCells[p];
            restorePlacedCell(idx, placementRow[idx], placementSite[idx] + deltaSites);
        }
    };

    auto blockGain = [&](const vector<int> &rowCells, int lo, int hi, int deltaSites) {
        double oldCost = 0.0;
        double newCost = 0.0;
        for (int p = lo; p <= hi; ++p) {
            int idx = rowCells[p];
            oldCost += displacementCost(idx, placementX[idx], placementY[idx]);
            long long newX = placementX[idx] + static_cast<long long>(deltaSites) * siteW;
            newCost += displacementCost(idx, newX, placementY[idx]);
        }
        return oldCost - newCost;
    };

    auto medianShift = [&](vector<int> shifts, int maxShift) {
        if (shifts.empty()) return 0;
        size_t mid = shifts.size() / 2;
        nth_element(shifts.begin(), shifts.begin() + mid, shifts.end());
        return max(0, min(maxShift, shifts[mid]));
    };

    for (int pass = 0; pass < passes; ++pass) {
        if (runtimeLimitReached()) {
            updateBestPlacement(alpha);
            return;
        }

        bool changed = false;
        for (int r = 0; r < numRows; ++r) {
            for (const Interval &seg : rows[r].freeSites) {
                vector<int> rowCells;
                int last = -3;
                for (int s = seg.l; s < seg.r; ++s) {
                    int owner = siteOwner[static_cast<size_t>(r) * numSites + s];
                    if (owner >= 0 && owner != last &&
                        placementRow[owner] == r && cellHeightRows[owner] == 1) {
                        rowCells.push_back(owner);
                    }
                    last = owner;
                }
                if (rowCells.empty()) continue;

                stable_sort(rowCells.begin(), rowCells.end(), [&](int a, int b) {
                    if (placementSite[a] != placementSite[b]) return placementSite[a] < placementSite[b];
                    return cells[a].name < cells[b].name;
                });

                int gaps = static_cast<int>(rowCells.size()) + 1;
                for (int gapIdx = 0; gapIdx < gaps; ++gapIdx) {
                    int gapL = seg.l;
                    int gapR = seg.r;
                    if (gapIdx > 0) {
                        int left = rowCells[gapIdx - 1];
                        gapL = placementSite[left] + cellWidthSites[left];
                    }
                    if (gapIdx < static_cast<int>(rowCells.size())) {
                        int right = rowCells[gapIdx];
                        gapR = placementSite[right];
                    }
                    int gap = gapR - gapL;
                    if (gap <= 0) continue;

                    int bestLo = -1;
                    int bestHi = -1;
                    int bestDelta = 0;
                    double bestGain = 1e-9;

                    if (gapIdx > 0) {
                        vector<int> shifts;
                        int hi = gapIdx - 1;
                        int loLimit = max(0, hi - maxBlock + 1);
                        for (int lo = hi; lo >= loLimit; --lo) {
                            int idx = rowCells[lo];
                            shifts.push_back(targetSites[idx] - placementSite[idx]);
                            int delta = medianShift(shifts, gap);
                            if (delta <= 0) continue;
                            double gain = blockGain(rowCells, lo, hi, delta);
                            if (gain > bestGain) {
                                bestGain = gain;
                                bestLo = lo;
                                bestHi = hi;
                                bestDelta = delta;
                            }
                        }
                    }

                    if (gapIdx < static_cast<int>(rowCells.size())) {
                        vector<int> shifts;
                        int lo = gapIdx;
                        int hiLimit = min(static_cast<int>(rowCells.size()) - 1, lo + maxBlock - 1);
                        for (int hi = lo; hi <= hiLimit; ++hi) {
                            int idx = rowCells[hi];
                            shifts.push_back(placementSite[idx] - targetSites[idx]);
                            int delta = medianShift(shifts, gap);
                            if (delta <= 0) continue;
                            double gain = blockGain(rowCells, lo, hi, -delta);
                            if (gain > bestGain) {
                                bestGain = gain;
                                bestLo = lo;
                                bestHi = hi;
                                bestDelta = -delta;
                            }
                        }
                    }

                    if (bestLo >= 0) {
                        moveBlock(rowCells, bestLo, bestHi, bestDelta);
                        changed = true;
                    }
                }
            }

            if ((r & 31) == 31 && runtimeLimitReached()) {
                updateBestPlacement(alpha);
                return;
            }
        }

        updateBestPlacement(alpha);
        if (!changed) break;
    }
}

void Legalizer::refineBoundaryTouches(double alpha) {
    if (cells.size() < 5000 || gridCols <= 0 || gridRows <= 0) return;

    double threshold = densityThreshold;
    double a = max(0.0, min(1.0, alpha));
    size_t gridCount = static_cast<size_t>(gridCols) * gridRows;

    auto gridArea = [&](int idx) {
        int r = idx / gridCols;
        int c = idx % gridCols;
        long long gx1 = gridLLX + static_cast<long long>(c) * gridSize;
        long long gx2 = min(gx1 + gridSize, gridURX);
        long long gy1 = gridLLY + static_cast<long long>(r) * gridSize;
        long long gy2 = min(gy1 + gridSize, gridURY);
        return static_cast<double>(gx2 - gx1) * static_cast<double>(gy2 - gy1);
    };

    auto overflowedWithArea = [&](int idx, double usedArea) {
        double area = gridArea(idx);
        return area > 0.0 && usedArea / area > threshold;
    };

    auto appendTouchedBins = [&](long long x, long long y, long long w, long long h,
                                 vector<int> &bins) {
        long long x2 = x + w;
        long long y2 = y + h;
        if (x2 < gridLLX || x > gridURX || y2 < gridLLY || y > gridURY) return;

        int c0 = max(gridLineLowerBound(x, gridLLX, gridURX, gridSize, gridCols) - 1, 0);
        int c1 = min(gridLineUpperBound(x2, gridLLX, gridURX, gridSize, gridCols), gridCols);
        int r0 = max(gridLineLowerBound(y, gridLLY, gridURY, gridSize, gridRows) - 1, 0);
        int r1 = min(gridLineUpperBound(y2, gridLLY, gridURY, gridSize, gridRows), gridRows);
        for (int r = r0; r < r1; ++r) {
            for (int c = c0; c < c1; ++c) {
                bins.push_back(r * gridCols + c);
            }
        }
    };

    auto appendDensityBins = [&](long long x, long long y, long long w, long long h,
                                 vector<pair<int, double>> &bins) {
        long long x1 = max(x, gridLLX);
        long long y1 = max(y, gridLLY);
        long long x2 = min(x + w, gridURX);
        long long y2 = min(y + h, gridURY);
        if (x1 >= x2 || y1 >= y2) return;

        int c0 = static_cast<int>(max(0LL, floorDiv(x1 - gridLLX, gridSize)));
        int c1 = static_cast<int>(min(static_cast<long long>(gridCols), ceilDiv(x2 - gridLLX, gridSize)));
        int r0 = static_cast<int>(max(0LL, floorDiv(y1 - gridLLY, gridSize)));
        int r1 = static_cast<int>(min(static_cast<long long>(gridRows), ceilDiv(y2 - gridLLY, gridSize)));

        for (int r = r0; r < r1; ++r) {
            long long gy1 = gridLLY + static_cast<long long>(r) * gridSize;
            long long gy2 = min(gy1 + gridSize, gridURY);
            for (int c = c0; c < c1; ++c) {
                long long gx1 = gridLLX + static_cast<long long>(c) * gridSize;
                long long gx2 = min(gx1 + gridSize, gridURX);
                long long ox = max(0LL, min(x + w, gx2) - max(x, gx1));
                long long oy = max(0LL, min(y + h, gy2) - max(y, gy1));
                if (ox <= 0 || oy <= 0) continue;
                bins.push_back({r * gridCols + c, static_cast<double>(ox) * oy});
            }
        }
    };

    auto sortUnique = [](vector<int> &items) {
        sort(items.begin(), items.end());
        items.erase(unique(items.begin(), items.end()), items.end());
    };

    vector<int> touchCounts(gridCount, 0);
    int totalTouched = 0;
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        if (!cells[i].name.empty() && cells[i].name[0] == 'h') continue;
        vector<int> bins;
        appendTouchedBins(placementX[i], placementY[i], cells[i].w, cells[i].h, bins);
        sortUnique(bins);
        for (int bin : bins) {
            if (touchCounts[bin]++ == 0) ++totalTouched;
        }
    }

    auto countedOverflow = [&](int idx) {
        return touchCounts[idx] > 0 && overflowedWithArea(idx, gridAreaUsed[idx]);
    };

    int overflowCount = 0;
    for (int idx = 0; idx < static_cast<int>(gridCount); ++idx) {
        if (countedOverflow(idx)) ++overflowCount;
    }

    double totalDisp = 0.0;
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        totalDisp += displacementCost(i, placementX[i], placementY[i]);
    }

    vector<int> order(cells.size());
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) order[i] = i;
    stable_sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        double dl = displacementCost(lhs, placementX[lhs], placementY[lhs]);
        double dr = displacementCost(rhs, placementX[rhs], placementY[rhs]);
        if (fabs(dl - dr) > 1e-9) return dl < dr;
        return cells[lhs].name < cells[rhs].name;
    });

    auto addCandidate = [](vector<int> &items, int value, int lo, int hi) {
        if (value < lo || value > hi) return;
        if (find(items.begin(), items.end(), value) == items.end()) items.push_back(value);
    };

    int moved = 0;
    int maxCellsToVisit = static_cast<int>(cells.size());
    for (int pos = 0; pos < maxCellsToVisit; ++pos) {
        if ((pos & 4095) == 4095 && runtimeLimitReached()) {
            updateBestPlacement(alpha);
            return;
        }

        int idx = order[pos];
        if (!cells[idx].name.empty() && cells[idx].name[0] == 'h') continue;

        int oldRow = placementRow[idx];
        int oldSite = placementSite[idx];
        long long oldX = placementX[idx];
        long long oldY = placementY[idx];
        double oldDisp = displacementCost(idx, oldX, oldY);
        double oldDor = totalTouched == 0 ? 0.0 : static_cast<double>(overflowCount) / totalTouched * 100.0;
        double oldScore = a * kDisplacementNormFactor * (totalDisp / static_cast<double>(cells.size())) +
                          (1.0 - a) * oldDor;

        vector<int> oldTouched;
        vector<pair<int, double>> oldDensity;
        appendTouchedBins(oldX, oldY, cells[idx].w, cells[idx].h, oldTouched);
        appendDensityBins(oldX, oldY, cells[idx].w, cells[idx].h, oldDensity);
        sortUnique(oldTouched);

        vector<int> oldAffected = oldTouched;
        for (const auto &[bin, area] : oldDensity) {
            (void)area;
            oldAffected.push_back(bin);
        }
        sortUnique(oldAffected);

        int savedTouched = totalTouched;
        int savedOverflow = overflowCount;
        double savedTotalDisp = totalDisp;

        for (int bin : oldAffected) {
            if (countedOverflow(bin)) --overflowCount;
        }
        for (int bin : oldTouched) {
            if (--touchCounts[bin] == 0) --totalTouched;
        }
        removePlacedCell(idx);
        for (int bin : oldAffected) {
            if (countedOverflow(bin)) ++overflowCount;
        }

        int maxStartRow = numRows - cellHeightRows[idx];
        int maxStartSite = numSites - cellWidthSites[idx];
        vector<int> rowCandidates;
        vector<int> siteCandidates;
        addCandidate(rowCandidates, oldRow, 0, maxStartRow);
        addCandidate(siteCandidates, oldSite, 0, maxStartSite);

        long long rowCenterLine = floorDiv(oldY - gridLLY, gridSize);
        for (long long k = rowCenterLine - 2; k <= rowCenterLine + 2; ++k) {
            long long boundary = gridLLY + k * gridSize;
            long long relStart = boundary - dieLLY;
            if (relStart % siteH == 0) {
                int row = static_cast<int>(relStart / siteH);
                if (abs(row - oldRow) <= 7) addCandidate(rowCandidates, row, 0, maxStartRow);
            }
            long long relEnd = boundary - cells[idx].h - dieLLY;
            if (relEnd % siteH == 0) {
                int row = static_cast<int>(relEnd / siteH);
                if (abs(row - oldRow) <= 7) addCandidate(rowCandidates, row, 0, maxStartRow);
            }
        }

        long long siteCenterLine = floorDiv(oldX - gridLLX, gridSize);
        for (long long k = siteCenterLine - 2; k <= siteCenterLine + 2; ++k) {
            long long boundary = gridLLX + k * gridSize;
            long long relStart = boundary - dieLLX;
            if (relStart % siteW == 0) {
                int site = static_cast<int>(relStart / siteW);
                if (abs(site - oldSite) <= 80) addCandidate(siteCandidates, site, 0, maxStartSite);
            }
            long long relEnd = boundary - cells[idx].w - dieLLX;
            if (relEnd % siteW == 0) {
                int site = static_cast<int>(relEnd / siteW);
                if (abs(site - oldSite) <= 80) addCandidate(siteCandidates, site, 0, maxStartSite);
            }
        }

        int bestRow = oldRow;
        int bestSite = oldSite;
        int bestTouched = savedTouched;
        int bestOverflow = savedOverflow;
        double bestTotalDisp = savedTotalDisp;
        double bestScore = oldScore;
        vector<int> bestTouchedBins;

        for (int row : rowCandidates) {
            for (int site : siteCandidates) {
                if (row == oldRow && site == oldSite) continue;
                if (!ownerCanPlace(idx, row, site)) continue;

                long long x = dieLLX + static_cast<long long>(site) * siteW;
                long long y = dieLLY + static_cast<long long>(row) * siteH;
                vector<int> newTouched;
                vector<pair<int, double>> newDensity;
                appendTouchedBins(x, y, cells[idx].w, cells[idx].h, newTouched);
                appendDensityBins(x, y, cells[idx].w, cells[idx].h, newDensity);
                sortUnique(newTouched);

                vector<int> affected = newTouched;
                for (const auto &[bin, area] : newDensity) {
                    (void)area;
                    affected.push_back(bin);
                }
                sortUnique(affected);

                int candidateTouched = totalTouched;
                for (int bin : newTouched) {
                    if (touchCounts[bin] == 0) ++candidateTouched;
                }
                if (candidateTouched == 0) continue;

                int candidateOverflow = overflowCount;
                for (int bin : affected) {
                    if (countedOverflow(bin)) --candidateOverflow;
                }

                for (int bin : affected) {
                    bool hasTouch = touchCounts[bin] > 0 ||
                                    find(newTouched.begin(), newTouched.end(), bin) != newTouched.end();
                    if (!hasTouch) continue;

                    double usedArea = gridAreaUsed[bin];
                    for (const auto &[areaBin, addArea] : newDensity) {
                        if (areaBin == bin) usedArea += addArea;
                    }
                    if (overflowedWithArea(bin, usedArea)) ++candidateOverflow;
                }

                double newDisp = displacementCost(idx, x, y);
                double candidateTotalDisp = totalDisp - oldDisp + newDisp;
                double candidateDor = static_cast<double>(candidateOverflow) / candidateTouched * 100.0;
                double candidateScore = a * kDisplacementNormFactor *
                                            (candidateTotalDisp / static_cast<double>(cells.size())) +
                                        (1.0 - a) * candidateDor;

                if (candidateScore < bestScore - 1e-9 ||
                    (fabs(candidateScore - bestScore) < 1e-9 && newDisp < oldDisp)) {
                    bestScore = candidateScore;
                    bestRow = row;
                    bestSite = site;
                    bestTouched = candidateTouched;
                    bestOverflow = candidateOverflow;
                    bestTotalDisp = candidateTotalDisp;
                    bestTouchedBins = std::move(newTouched);
                }
            }
        }

        if (bestRow != oldRow || bestSite != oldSite) {
            restorePlacedCell(idx, bestRow, bestSite);
            for (int bin : bestTouchedBins) ++touchCounts[bin];
            totalTouched = bestTouched;
            overflowCount = bestOverflow;
            totalDisp = bestTotalDisp;
            ++moved;
        } else {
            restorePlacedCell(idx, oldRow, oldSite);
            for (int bin : oldTouched) ++touchCounts[bin];
            totalTouched = savedTouched;
            overflowCount = savedOverflow;
            totalDisp = savedTotalDisp;
        }

        if ((moved & 4095) == 4095) updateBestPlacement(alpha);
    }

    updateBestPlacement(alpha);
}

double Legalizer::estimateDensityPenalty(long long x, long long y, long long w, long long h) const {
    if (gridCols <= 0 || gridRows <= 0) return 0.0;
    double threshold = densityThreshold;
    long long x1 = max(x, gridLLX);
    long long y1 = max(y, gridLLY);
    long long x2 = min(x + w, gridURX);
    long long y2 = min(y + h, gridURY);
    if (x1 >= x2 || y1 >= y2) return 0.0;

    int c0 = static_cast<int>(max(0LL, floorDiv(x1 - gridLLX, gridSize)));
    int c1 = static_cast<int>(min(static_cast<long long>(gridCols), ceilDiv(x2 - gridLLX, gridSize)));
    int r0 = static_cast<int>(max(0LL, floorDiv(y1 - gridLLY, gridSize)));
    int r1 = static_cast<int>(min(static_cast<long long>(gridRows), ceilDiv(y2 - gridLLY, gridSize)));
    double penalty = 0.0;
    int touched = 0;

    for (int r = r0; r < r1; ++r) {
        long long gy1 = gridLLY + static_cast<long long>(r) * gridSize;
        long long gy2 = min(gy1 + gridSize, gridURY);
        for (int c = c0; c < c1; ++c) {
            long long gx1 = gridLLX + static_cast<long long>(c) * gridSize;
            long long gx2 = min(gx1 + gridSize, gridURX);
            long long ox = max(0LL, min(x + w, gx2) - max(x, gx1));
            long long oy = max(0LL, min(y + h, gy2) - max(y, gy1));
            if (ox == 0 || oy == 0) continue;
            if (gridExcluded[static_cast<size_t>(r) * gridCols + c]) continue;
            double gridArea = static_cast<double>(gx2 - gx1) * static_cast<double>(gy2 - gy1);
            if (gridArea <= 0.0) continue;
            double addArea = static_cast<double>(ox) * oy;
            double util = (gridAreaUsed[static_cast<size_t>(r) * gridCols + c] + addArea) / gridArea;
            penalty += max(0.0, util - threshold) * 100.0;
            ++touched;
        }
    }
    return touched == 0 ? 0.0 : penalty / touched;
}

void Legalizer::occupy(int row, int site, int widthSites, int heightRows) {
    int end = site + widthSites;
    for (int rr = row; rr < row + heightRows; ++rr) {
        vector<Interval> next;
        next.reserve(rows[rr].freeSites.size() + 1);
        for (const auto &seg : rows[rr].freeSites) {
            if (end <= seg.l || seg.r <= site) {
                next.push_back(seg);
            } else {
                if (seg.l < site) next.push_back({seg.l, site});
                if (end < seg.r) next.push_back({end, seg.r});
            }
        }
        rows[rr].freeSites.swap(next);
    }
}

void Legalizer::addCellToDensity(long long x, long long y, long long w, long long h) {
    addCellToDensityDelta(x, y, w, h, 1.0);
}

void Legalizer::addCellToDensityDelta(long long x, long long y, long long w, long long h, double sign) {
    if (gridCols <= 0 || gridRows <= 0) return;
    long long x1 = max(x, gridLLX);
    long long y1 = max(y, gridLLY);
    long long x2 = min(x + w, gridURX);
    long long y2 = min(y + h, gridURY);
    if (x1 >= x2 || y1 >= y2) return;

    int c0 = static_cast<int>(max(0LL, floorDiv(x1 - gridLLX, gridSize)));
    int c1 = static_cast<int>(min(static_cast<long long>(gridCols), ceilDiv(x2 - gridLLX, gridSize)));
    int r0 = static_cast<int>(max(0LL, floorDiv(y1 - gridLLY, gridSize)));
    int r1 = static_cast<int>(min(static_cast<long long>(gridRows), ceilDiv(y2 - gridLLY, gridSize)));

    for (int r = r0; r < r1; ++r) {
        long long gy1 = gridLLY + static_cast<long long>(r) * gridSize;
        long long gy2 = min(gy1 + gridSize, gridURY);
        for (int c = c0; c < c1; ++c) {
            long long gx1 = gridLLX + static_cast<long long>(c) * gridSize;
            long long gx2 = min(gx1 + gridSize, gridURX);
            long long ox = max(0LL, min(x + w, gx2) - max(x, gx1));
            long long oy = max(0LL, min(y + h, gy2) - max(y, gy1));
            if (gridExcluded[static_cast<size_t>(r) * gridCols + c]) continue;
            double &used = gridAreaUsed[static_cast<size_t>(r) * gridCols + c];
            used += sign * static_cast<double>(ox) * oy;
            if (used < 0.0 && used > -1e-6) used = 0.0;
        }
    }
}
