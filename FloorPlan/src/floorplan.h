#pragma once
#include <string>
#include <vector>
#include <map>
#include <random>
#include <climits>
#include <algorithm>

enum SymType { SYM_VERT, SYM_HORIZ };

struct Module {
    std::string name;
    int w = 0, h = 0;
    double x = 0, y = 0;
    int rot = 0; // 1 means rotated by 90 degrees from input orientation
};

struct SymPair {
    int a = 0, b = 0;
    int rep = 1; // 0=a is representative, 1=b is representative
};

struct SymGroup {
    SymType type = SYM_VERT;
    std::vector<SymPair> pairs;
    std::vector<int> selfsym; // module indices
};

// B*-tree node (used for both HB-tree and ASF-tree)
// HB:  id >= 0 → module index; id < 0 → group index encoded as -(gid+1)
// ASF: id = rep_idx in [0, np+ns)
//        0..np-1  → pair rep index
//        np..np+ns-1 → selfsym index (these must stay on rightmost branch)
struct BNode {
    int id   = -1;
    int par  = -1, left = -1, right = -1;
    double px = 0, py = 0;      // placed bottom-left (set during packing)
    double pw = 0, ph = 0;      // placed dimensions
};

struct ProfileSeg {
    double l = 0, r = 0, y = 0; // constant y on [l, r)
    ProfileSeg() {}
    ProfileSeg(double l_, double r_, double y_) : l(l_), r(r_), y(y_) {}
};

// Contour: maps x_start → height for all x in [x_start, next_key)
struct Contour {
    std::map<double,double> seg;
    Contour()  { seg[0] = 0; seg[1e100] = 0; }
    void reset(){ seg.clear(); seg[0] = 0; seg[1e100] = 0; }
    double query(double x, double w) const;
    void update(double x, double w, double h);
};

class Floorplan {
public:
    std::vector<Module>              mods;
    std::vector<SymGroup>            grps;
    std::vector<BNode>               hb;      // HB-tree nodes
    std::vector<std::vector<BNode>>  asf;     // ASF-tree nodes per group
    std::vector<std::vector<int>>    asf_pos; // rep id -> node index
    std::vector<double>              iw, ih;  // island bounding-box dims
    std::vector<double>              island_x, island_y;
    std::vector<std::vector<ProfileSeg>> island_top;
    std::vector<std::vector<ProfileSeg>> island_bottom;
    std::vector<char>                group_dirty;

    void      init(unsigned int seed = 2024);
    void      greedy_pack_large(unsigned int seed = 2024);
    void      pack();
    double    cost() const;
    void      perturb(std::mt19937& rng);
    void      save();
    void      restore();

    // Accessors
    int  hb_root()           const;
    int  asf_root(int g)     const;
    bool on_rb(const std::vector<BNode>& tr, int i) const;

private:
    void pack_asf(int g);

    // Tree surgery helpers
    void detach(std::vector<BNode>& tr, int i);
    void attach(std::vector<BNode>& tr, int i, int par, bool lft);

    // SA perturbation ops
    void op_swap_hb(std::mt19937& rng);
    void op_move_hb(std::mt19937& rng);
    void op_swap_asf_pairs(int g, std::mt19937& rng);
    void op_swap_asf_sel(int g, std::mt19937& rng);
    void op_move_asf_pair(int g, std::mt19937& rng);
    void op_change_rep(int g, std::mt19937& rng);
    void op_rotate(int g, std::mt19937& rng);
    void op_flip_axis(int g);
    void op_flip_sym(int g);

    bool is_selfsym_node(int g, int node_idx) const;
    int  rand_node(const std::vector<BNode>& tr, std::mt19937& rng) const;

    // Reused ASF working buffers to avoid per-pack allocations.
    std::vector<std::vector<double>> asf_rw_buf;
    std::vector<std::vector<double>> asf_rh_buf;
    std::vector<std::vector<double>> asf_rx_buf;
    std::vector<std::vector<double>> asf_ry_buf;

    double    packed_w = 0;
    double    packed_h = 0;
    double    packed_area = 0;

    // Saved state
    std::vector<Module>              sv_mods;
    std::vector<SymGroup>            sv_grps;
    std::vector<BNode>               sv_hb;
    std::vector<std::vector<BNode>>  sv_asf;
    std::vector<std::vector<int>>    sv_asf_pos;
    std::vector<double>              sv_iw, sv_ih;
    std::vector<double>              sv_island_x, sv_island_y;
    std::vector<std::vector<ProfileSeg>> sv_island_top;
    std::vector<std::vector<ProfileSeg>> sv_island_bottom;
    std::vector<char>                sv_group_dirty;
    double                           sv_packed_w = 0;
    double                           sv_packed_h = 0;
    double                           sv_packed_area = 0;
};
