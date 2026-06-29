#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

using int64 = long long;
static constexpr int64 kMinCoordinate = -1000000000LL;
static constexpr int64 kMaxCoordinate = 1000000000LL;

struct Point {
    int32_t x;
    int32_t y;
};

struct Edge {
    int64 weight;
    int u;
    int v;
};

class FastInput {
public:
    explicit FastInput(const char* path) : file_(std::fopen(path, "rb")) {
        if (file_ == nullptr) {
            throw std::runtime_error("cannot open input file");
        }
    }

    ~FastInput() {
        if (file_ != nullptr) {
            std::fclose(file_);
        }
    }

    bool readInt(int64& value) {
        int c = readChar();
        while (c != EOF && c <= ' ') {
            c = readChar();
        }
        if (c == EOF) {
            return false;
        }

        int sign = 1;
        if (c == '-') {
            sign = -1;
            c = readChar();
        }

        std::uint64_t result = 0;
        const std::uint64_t limit =
            sign < 0
                ? static_cast<std::uint64_t>(std::numeric_limits<int64>::max()) + 1ULL
                : static_cast<std::uint64_t>(std::numeric_limits<int64>::max());
        bool has_digit = false;
        while (c >= '0' && c <= '9') {
            has_digit = true;
            const std::uint64_t digit = static_cast<std::uint64_t>(c - '0');
            if (result > (limit - digit) / 10) {
                return false;
            }
            result = result * 10 + digit;
            c = readChar();
        }

        if (!has_digit) {
            return false;
        }
        if (sign < 0 && result == limit) {
            value = std::numeric_limits<int64>::min();
        } else {
            value = sign < 0 ? -static_cast<int64>(result) : static_cast<int64>(result);
        }
        return true;
    }

private:
    static constexpr std::size_t kBufferSize = 1 << 20;

    int readChar() {
        if (pos_ == size_) {
            size_ = std::fread(buffer_, 1, kBufferSize, file_);
            pos_ = 0;
            if (size_ == 0) {
                return EOF;
            }
        }
        return buffer_[pos_++];
    }

    FILE* file_;
    char buffer_[kBufferSize];
    std::size_t pos_ = 0;
    std::size_t size_ = 0;
};

class DisjointSet {
public:
    explicit DisjointSet(std::size_t n) : parent_(n), rank_(n, 0) {
        std::iota(parent_.begin(), parent_.end(), 0);
    }

    int find(int x) {
        int root = x;
        while (parent_[root] != root) {
            root = parent_[root];
        }
        while (parent_[x] != x) {
            int next = parent_[x];
            parent_[x] = root;
            x = next;
        }
        return root;
    }

    bool unite(int a, int b) {
        int root_a = find(a);
        int root_b = find(b);
        if (root_a == root_b) {
            return false;
        }
        if (rank_[root_a] < rank_[root_b]) {
            std::swap(root_a, root_b);
        }
        parent_[root_b] = root_a;
        if (rank_[root_a] == rank_[root_b]) {
            ++rank_[root_a];
        }
        return true;
    }

private:
    std::vector<int> parent_;
    std::vector<unsigned char> rank_;
};

class FenwickTree {
public:
    explicit FenwickTree(std::size_t n) : tree_(n + 1, 0), total_(0) {}

    void add(int index, int delta) {
        total_ += delta;
        for (int i = index + 1; i < static_cast<int>(tree_.size()); i += i & -i) {
            tree_[i] += delta;
        }
    }

    int nextAtOrAfter(int index) const {
        const int before = prefixCount(index);
        if (before == total_) {
            return -1;
        }
        return findByOrder(before + 1);
    }

private:
    int prefixCount(int count) const {
        int result = 0;
        for (int i = count; i > 0; i -= i & -i) {
            result += tree_[i];
        }
        return result;
    }

    int findByOrder(int order) const {
        int index = 0;
        int step = 1;
        while ((step << 1) < static_cast<int>(tree_.size())) {
            step <<= 1;
        }

        for (; step > 0; step >>= 1) {
            const int next = index + step;
            if (next < static_cast<int>(tree_.size()) && tree_[next] < order) {
                index = next;
                order -= tree_[next];
            }
        }
        return index;
    }

    std::vector<int> tree_;
    int total_;
};

static int64 transformedX(const Point& p, int dir) {
    switch (dir) {
        case 0:
            return p.x;
        case 1:
            return p.y;
        case 2:
            return -p.y;
        default:
            return p.x;
    }
}

static int64 transformedY(const Point& p, int dir) {
    switch (dir) {
        case 0:
            return p.y;
        case 1:
            return p.x;
        case 2:
            return p.x;
        default:
            return -p.y;
    }
}

static int64 manhattanDistance(const Point& a, const Point& b) {
    return std::llabs(static_cast<int64>(a.x) - static_cast<int64>(b.x)) +
           std::llabs(static_cast<int64>(a.y) - static_cast<int64>(b.y));
}

static int32_t activeKey(const Point& p, int dir) {
    return static_cast<int32_t>(-transformedY(p, dir));
}

static std::vector<Point> readUniquePoints(const char* path) {
    FastInput input(path);

    int64 raw_count = 0;
    if (!input.readInt(raw_count) || raw_count < 0) {
        throw std::runtime_error("invalid point count");
    }
    if (raw_count > std::numeric_limits<int>::max()) {
        throw std::runtime_error("too many points for this implementation");
    }

    std::vector<Point> points;
    points.reserve(static_cast<std::size_t>(raw_count));
    for (int64 i = 0; i < raw_count; ++i) {
        int64 x = 0;
        int64 y = 0;
        if (!input.readInt(x) || !input.readInt(y)) {
            throw std::runtime_error("input ended before all points were read");
        }
        if (x < kMinCoordinate || x > kMaxCoordinate || y < kMinCoordinate || y > kMaxCoordinate) {
            throw std::runtime_error("coordinate is outside assignment bounds");
        }
        points.push_back({static_cast<int32_t>(x), static_cast<int32_t>(y)});
    }

    std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) {
        if (a.x != b.x) {
            return a.x < b.x;
        }
        return a.y < b.y;
    });
    points.erase(
        std::unique(points.begin(), points.end(), [](const Point& a, const Point& b) {
            return a.x == b.x && a.y == b.y;
        }),
        points.end());
    return points;
}

static void addOctantCandidateEdges(const std::vector<Point>& points, std::vector<Edge>& edges) {
    const std::size_t n = points.size();
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);

    // These four symmetric coordinate systems enumerate the eight directed Manhattan octants.
    for (int dir = 0; dir < 4; ++dir) {
        std::vector<int32_t> keys;
        keys.reserve(n);
        for (const Point& point : points) {
            keys.push_back(activeKey(point, dir));
        }
        std::sort(keys.begin(), keys.end());
        keys.erase(std::unique(keys.begin(), keys.end()), keys.end());

        std::sort(order.begin(), order.end(), [&](int lhs, int rhs) {
            const Point& a = points[lhs];
            const Point& b = points[rhs];
            const int64 ax = transformedX(a, dir);
            const int64 ay = transformedY(a, dir);
            const int64 bx = transformedX(b, dir);
            const int64 by = transformedY(b, dir);
            const int64 asum = ax + ay;
            const int64 bsum = bx + by;
            if (asum != bsum) {
                return asum < bsum;
            }
            if (ax != bx) {
                return ax < bx;
            }
            if (ay != by) {
                return ay < by;
            }
            return lhs < rhs;
        });

        std::vector<int> active_id(keys.size(), -1);
        FenwickTree active_tree(keys.size());
        for (const int id : order) {
            const int64 x = transformedX(points[id], dir);
            const int64 y = transformedY(points[id], dir);
            const int32_t key = activeKey(points[id], dir);
            const int key_index = static_cast<int>(
                std::lower_bound(keys.begin(), keys.end(), key) - keys.begin());

            int active_index = active_tree.nextAtOrAfter(key_index);
            while (active_index != -1) {
                const int other = active_id[active_index];
                const int64 ox = transformedX(points[other], dir);
                const int64 oy = transformedY(points[other], dir);
                if (x - ox < y - oy) {
                    break;
                }

                edges.push_back({manhattanDistance(points[id], points[other]), id, other});
                active_id[active_index] = -1;
                active_tree.add(active_index, -1);
                active_index = active_tree.nextAtOrAfter(active_index);
            }

            if (active_id[key_index] == -1) {
                active_tree.add(key_index, 1);
            }
            active_id[key_index] = id;
        }
    }
}

static int64 kruskalWeight(std::size_t n, std::vector<Edge>& edges) {
    if (n <= 1) {
        return 0;
    }

    std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
        if (a.weight != b.weight) {
            return a.weight < b.weight;
        }
        if (a.u != b.u) {
            return a.u < b.u;
        }
        return a.v < b.v;
    });

    DisjointSet dsu(n);
    int64 total = 0;
    std::size_t used_edges = 0;
    for (const Edge& edge : edges) {
        if (dsu.unite(edge.u, edge.v)) {
            total += edge.weight;
            ++used_edges;
            if (used_edges + 1 == n) {
                break;
            }
        }
    }

    if (used_edges + 1 != n) {
        throw std::runtime_error("candidate graph is disconnected");
    }
    return total;
}

static void writeAnswer(const char* path, int64 answer) {
    FILE* output = std::fopen(path, "wb");
    if (output == nullptr) {
        throw std::runtime_error("cannot open output file");
    }
    std::fprintf(output, "%lld\n", answer);
    std::fclose(output);
}

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            std::fprintf(stderr, "Usage: %s <input.dat> <output.dat>\n", argv[0]);
            return 1;
        }

        std::vector<Point> points = readUniquePoints(argv[1]);
        std::vector<Edge> edges;
        if (points.size() > 1) {
            edges.reserve(4 * (points.size() - 1));
            addOctantCandidateEdges(points, edges);
        }

        const int64 answer = kruskalWeight(points.size(), edges);
        writeAnswer(argv[2], answer);
        return 0;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "RMST error: %s\n", error.what());
        return 1;
    }
}
