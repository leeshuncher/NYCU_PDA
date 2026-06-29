#include "floorplan.h"
#include <cmath>
#include <limits>
#include <numeric>
#include <stack>
using namespace std;

namespace {
const double EPS = 1e-6;

void rotate_module(Module& mod) {
    swap(mod.w, mod.h);
    mod.rot ^= 1;
}

bool selfsym_type_feasible(const SymGroup& grp, const vector<Module>& mods, SymType type) {
    (void)grp;
    (void)mods;
    (void)type;
    return true;
}

void clear_tree_links(vector<BNode>& tr) {
    for (auto& nd : tr) {
        nd.par = nd.left = nd.right = -1;
        nd.px = nd.py = nd.pw = nd.ph = 0;
    }
}

void build_complete_binary_tree(vector<BNode>& tr, const vector<int>& nodes) {
    for (int pos = 1; pos < (int)nodes.size(); pos++) {
        int cur = nodes[pos];
        int par = nodes[(pos - 1) / 2];
        tr[cur].par = par;
        if (pos & 1) tr[par].left = cur;
        else         tr[par].right = cur;
    }
}

struct Rect {
    double x = 0, y = 0, w = 0, h = 0;
    Rect() {}
    Rect(double x_, double y_, double w_, double h_) : x(x_), y(y_), w(w_), h(h_) {}
};

void contour_split(Contour& c, double x) {
    auto it = c.seg.upper_bound(x);
    --it;
    if (it->first + EPS < x) c.seg[x] = it->second;
}

void contour_normalize(Contour& c) {
    for (auto it = next(c.seg.begin()); it != c.seg.end(); ) {
        if (next(it) == c.seg.end()) break; // keep the sentinel boundary
        auto prev = it;
        --prev;
        if (fabs(prev->second - it->second) <= EPS) {
            it = c.seg.erase(it);
        } else {
            ++it;
        }
    }
}

void append_profile_seg(vector<ProfileSeg>& prof, double l, double r, double y) {
    if (l + EPS >= r) return;
    if (!prof.empty() && fabs(prof.back().r - l) <= EPS && fabs(prof.back().y - y) <= EPS) {
        prof.back().r = r;
        return;
    }
    prof.push_back(ProfileSeg{l, r, y});
}

pair<vector<ProfileSeg>, vector<ProfileSeg>> extract_profiles(const vector<Rect>& rects) {
    vector<double> xs;
    xs.reserve(rects.size() * 2);
    for (const auto& rc : rects) {
        xs.push_back(rc.x);
        xs.push_back(rc.x + rc.w);
    }
    sort(xs.begin(), xs.end());
    xs.erase(unique(xs.begin(), xs.end()), xs.end());

    vector<ProfileSeg> bottom;
    vector<ProfileSeg> top;
    for (int i = 0; i + 1 < (int)xs.size(); i++) {
        double l = xs[i];
        double r = xs[i + 1];
        double min_y = numeric_limits<double>::infinity();
        double max_y = 0;
        for (const auto& rc : rects) {
            if (rc.x < r - EPS && rc.x + rc.w > l + EPS) {
                min_y = min(min_y, rc.y);
                max_y = max(max_y, rc.y + rc.h);
            }
        }
        if (!isfinite(min_y)) continue;
        append_profile_seg(bottom, l, r, min_y);
        append_profile_seg(top, l, r, max_y);
    }
    return {bottom, top};
}

double query_profile(const Contour& c, double x0, const vector<ProfileSeg>& bottom) {
    if (bottom.empty()) return 0;

    double best_y = 0;
    for (const auto& seg : bottom) {
        double cur = x0 + seg.l;
        double end = x0 + seg.r;
        auto it = c.seg.upper_bound(cur);
        --it;
        while (cur + EPS < end) {
            auto next_it = next(it);
            double next_x = min(end, next_it->first);
            best_y = max(best_y, it->second - seg.y);
            cur = next_x;
            if (cur + EPS < end) it = next_it;
        }
    }
    return max(best_y, 0.0);
}

void update_profile(Contour& c, double x0, const vector<ProfileSeg>& top, double base_y) {
    if (top.empty()) return;

    double l = x0 + top.front().l;
    double r = x0 + top.back().r;
    vector<double> cuts;
    cuts.reserve(top.size() * 2 + c.seg.size());
    cuts.push_back(l);
    cuts.push_back(r);
    for (const auto& seg : top) {
        cuts.push_back(x0 + seg.l);
        cuts.push_back(x0 + seg.r);
    }
    for (auto it = c.seg.lower_bound(l); it != c.seg.end() && it->first < r; ++it)
        cuts.push_back(it->first);
    sort(cuts.begin(), cuts.end());
    cuts.erase(unique(cuts.begin(), cuts.end()), cuts.end());

    vector<pair<double, double>> new_spans;
    new_spans.reserve(cuts.size());
    int prof_idx = 0;
    for (int i = 0; i + 1 < (int)cuts.size(); i++) {
        double start = cuts[i];
        double end = cuts[i + 1];
        if (start + EPS >= end) continue;
        double local_x = start - x0;
        while (prof_idx + 1 < (int)top.size() && top[prof_idx].r <= local_x + EPS) ++prof_idx;
        new_spans.push_back({start, base_y + top[prof_idx].y});
    }

    contour_split(c, l);
    contour_split(c, r);
    c.seg.erase(c.seg.find(l), c.seg.find(r));
    for (const auto& span : new_spans) c.seg[span.first] = span.second;
    contour_normalize(c);
}
}

// ─── Contour ────────────────────────────────────────────────────────────────

double Contour::query(double x, double w) const {
    auto it = seg.upper_bound(x); --it;
    double h = 0;
    while (it->first < x + w - EPS) { h = max(h, it->second); ++it; }
    return h;
}

void Contour::update(double x, double w, double h) {
    // split at x
    auto it = seg.upper_bound(x); --it;
    if (it->first + EPS < x) seg[x] = it->second;
    // split at x+w
    it = seg.upper_bound(x + w); --it;
    if (it->first + EPS < x + w) seg[x + w] = it->second;
    // erase [x, x+w) and set height
    seg.erase(seg.find(x), seg.find(x + w));
    seg[x] = h;
}

// ─── Tree helpers ─────────────────────────────────────────────────────────

int Floorplan::hb_root() const {
    for (int i = 0; i < (int)hb.size(); i++)
        if (hb[i].par == -1) return i;
    return 0;
}

int Floorplan::asf_root(int g) const {
    for (int i = 0; i < (int)asf[g].size(); i++)
        if (asf[g][i].par == -1) return i;
    return 0;
}

// Is node i on the rightmost branch? (reachable from root via right-child edges only)
bool Floorplan::on_rb(const vector<BNode>& tr, int i) const {
    int cur = i;
    while (tr[cur].par != -1) {
        int p = tr[cur].par;
        if (tr[p].right != cur) return false;
        cur = p;
    }
    return true;
}

bool Floorplan::is_selfsym_node(int g, int node_idx) const {
    return asf[g][node_idx].id >= (int)grps[g].pairs.size();
}

int Floorplan::rand_node(const vector<BNode>& tr, mt19937& rng) const {
    return (int)(rng() % tr.size());
}

// Detach node i from tree tr, reconnecting its children to its parent.
void Floorplan::detach(vector<BNode>& tr, int i) {
    int p = tr[i].par;
    int l = tr[i].left;
    int r = tr[i].right;

    int repl = -1;
    if (l == -1 && r == -1) {
        repl = -1;
    } else if (l == -1) {
        repl = r;
    } else if (r == -1) {
        repl = l;
    } else {
        // Both children: graft left subtree to leftmost position in right subtree
        repl = r;
        int lm = r;
        while (tr[lm].left != -1) lm = tr[lm].left;
        tr[lm].left = l;
        tr[l].par = lm;
    }

    if (repl != -1) tr[repl].par = p;
    if (p == -1) {
        // i was root; repl becomes new root
    } else if (tr[p].left == i) {
        tr[p].left = repl;
    } else {
        tr[p].right = repl;
    }

    tr[i].par = tr[i].left = tr[i].right = -1;
}

// Attach node i (must have no children) as left (lft=true) or right child of par.
// Existing child in that slot is pushed down to i's left child.
void Floorplan::attach(vector<BNode>& tr, int i, int par, bool lft) {
    int existing = lft ? tr[par].left : tr[par].right;
    if (lft) tr[par].left  = i;
    else     tr[par].right = i;
    tr[i].par = par;
    if (existing != -1) {
        tr[i].left    = existing;
        tr[existing].par = i;
    }
}

// ─── Init ─────────────────────────────────────────────────────────────────

void Floorplan::init(unsigned int seed) {
    int nm = mods.size();
    int ng = grps.size();
    iw.resize(ng, 0);
    ih.resize(ng, 0);
    island_x.assign(ng, 0);
    island_y.assign(ng, 0);
    asf.resize(ng);
    asf_pos.resize(ng);
    asf_rw_buf.resize(ng);
    asf_rh_buf.resize(ng);
    asf_rx_buf.resize(ng);
    asf_ry_buf.resize(ng);
    island_top.resize(ng);
    island_bottom.resize(ng);
    group_dirty.assign(ng, 1);

    // Mark which modules belong to a symmetry group
    vector<bool> in_grp(nm, false);
    for (int g = 0; g < ng; g++) {
        for (auto& pp : grps[g].pairs) {
            in_grp[pp.a] = in_grp[pp.b] = true;
        }
        for (int m : grps[g].selfsym) in_grp[m] = true;
    }

    auto rep_dims = [&](int g, int rep_id, SymType type) {
        int np = grps[g].pairs.size();
        if (rep_id < np) {
            int mid = (grps[g].pairs[rep_id].rep == 0) ? grps[g].pairs[rep_id].a
                                                       : grps[g].pairs[rep_id].b;
            return pair<double, double>{mods[mid].w, mods[mid].h};
        }

        int mid = grps[g].selfsym[rep_id - np];
        if (type == SYM_VERT) return pair<double, double>{mods[mid].w / 2.0, mods[mid].h};
        return pair<double, double>{mods[mid].w, mods[mid].h / 2.0};
    };

    auto rep_before = [&](int g, int lhs, int rhs, SymType type) {
        pair<double, double> ldim = rep_dims(g, lhs, type);
        pair<double, double> rdim = rep_dims(g, rhs, type);
        double lw = ldim.first;
        double lh = ldim.second;
        double rw = rdim.first;
        double rh = rdim.second;
        double la = lw * lh;
        double ra = rw * rh;
        if (la != ra) return la > ra;
        double lmax = max(lw, lh), rmax = max(rw, rh);
        if (fabs(lmax - rmax) > EPS) return lmax > rmax;
        double lmin = min(lw, lh), rmin = min(rw, rh);
        if (fabs(lmin - rmin) > EPS) return lmin > rmin;
        return lhs < rhs;
    };

    auto rebuild_asf_group = [&](int g, SymType type) {
        auto& grp = grps[g];
        grp.type = type;

        int np = grp.pairs.size();
        int ns = grp.selfsym.size();
        int nn = np + ns;

        asf[g].assign(nn, BNode());
        asf_pos[g].assign(nn, -1);
        asf_rw_buf[g].resize(nn);
        asf_rh_buf[g].resize(nn);
        asf_rx_buf[g].resize(nn);
        asf_ry_buf[g].resize(nn);
        clear_tree_links(asf[g]);

        vector<int> pair_ids(np);
        iota(pair_ids.begin(), pair_ids.end(), 0);
        sort(pair_ids.begin(), pair_ids.end(),
             [&](int lhs, int rhs) { return rep_before(g, lhs, rhs, type); });

        vector<int> self_ids(ns);
        for (int i = 0; i < ns; i++) self_ids[i] = np + i;
        sort(self_ids.begin(), self_ids.end(),
             [&](int lhs, int rhs) { return rep_before(g, lhs, rhs, type); });

        for (int i = 0; i < np; i++) {
            asf[g][i].id = pair_ids[i];
            asf_pos[g][pair_ids[i]] = i;
        }
        for (int i = 0; i < ns; i++) {
            int node_idx = np + i;
            asf[g][node_idx].id = self_ids[i];
            asf_pos[g][self_ids[i]] = node_idx;
        }

        if (np > 0) {
            vector<int> pair_nodes(np);
            iota(pair_nodes.begin(), pair_nodes.end(), 0);
            build_complete_binary_tree(asf[g], pair_nodes);
        }

        if (ns > 0) {
            int root = np;
            asf[g][root].par = -1;
            for (int i = 1; i < ns; i++) {
                int prev = np + i - 1;
                int cur = np + i;
                asf[g][prev].right = cur;
                asf[g][cur].par = prev;
            }
            if (np > 0) {
                asf[g][root].left = 0;
                asf[g][0].par = root;
            }
        }
    };

    // Build symmetry groups first so HB ordering can use their island sizes.
    for (int g = 0; g < ng; g++) {
        bool can_vert = selfsym_type_feasible(grps[g], mods, SYM_VERT);
        bool can_horiz = selfsym_type_feasible(grps[g], mods, SYM_HORIZ);

        double area_vert = numeric_limits<double>::infinity();
        double area_horiz = numeric_limits<double>::infinity();
        double span_vert = numeric_limits<double>::infinity();
        double span_horiz = numeric_limits<double>::infinity();

        if (can_vert || !can_horiz) {
            rebuild_asf_group(g, SYM_VERT);
            pack_asf(g);
            area_vert = iw[g] * ih[g];
            span_vert = max(iw[g], ih[g]);
        }

        if (can_horiz) {
            rebuild_asf_group(g, SYM_HORIZ);
            pack_asf(g);
            area_horiz = iw[g] * ih[g];
            span_horiz = max(iw[g], ih[g]);
        }

        SymType best = SYM_VERT;
        if (area_horiz < area_vert ||
            (area_horiz == area_vert && span_horiz < span_vert)) {
            best = SYM_HORIZ;
        }

        rebuild_asf_group(g, best);
        pack_asf(g);
    }

    // Build HB-tree from multiple deterministic orderings, then greedily insert
    // each item where the current packed area is smallest.
    vector<int> hb_ids;
    hb_ids.reserve(nm);
    for (int m = 0; m < nm; m++)
        if (!in_grp[m]) hb_ids.push_back(m);
    for (int g = 0; g < ng; g++)
        hb_ids.push_back(-(g + 1));

    auto dims_of = [&](int id) {
        if (id >= 0) return pair<double, double>{mods[id].w, mods[id].h};
        int g = -(id + 1);
        return pair<double, double>{iw[g], ih[g]};
    };

    auto hb_before = [&](int lhs, int rhs) {
        pair<double, double> ldim = dims_of(lhs);
        pair<double, double> rdim = dims_of(rhs);
        double lw = ldim.first;
        double lh = ldim.second;
        double rw = rdim.first;
        double rh = rdim.second;
        double la = lw * lh;
        double ra = rw * rh;
        if (fabs(la - ra) > EPS) return la > ra;
        double lmax = max(lw, lh), rmax = max(rw, rh);
        if (fabs(lmax - rmax) > EPS) return lmax > rmax;
        double lmin = min(lw, lh), rmin = min(rw, rh);
        if (fabs(lmin - rmin) > EPS) return lmin > rmin;
        return lhs < rhs;
    };

    auto eval_hb_tree = [&](const vector<BNode>& tree) {
        hb = tree;
        group_dirty.assign(ng, 1);
        pack();
        return cost();
    };

    auto make_insertion_trial = [](const vector<BNode>& base, int id, int par, bool lft) {
        vector<BNode> trial = base;
        int idx = (int)trial.size();
        trial.push_back(BNode());
        trial[idx].id = id;
        trial[idx].par = par;

        int existing = lft ? trial[par].left : trial[par].right;
        if (lft) trial[par].left = idx;
        else     trial[par].right = idx;
        if (existing != -1) {
            trial[idx].left = existing;
            trial[existing].par = idx;
        }
        return trial;
    };

    auto build_greedy_hb = [&](const vector<int>& order) {
        vector<BNode> cur;
        if (order.empty()) return cur;

        cur.push_back(BNode());
        cur[0].id = order[0];

        for (int pos = 1; pos < (int)order.size(); pos++) {
            double best_area = numeric_limits<double>::infinity();
            vector<BNode> best_tree;

            for (int p = 0; p < (int)cur.size(); p++) {
                for (int side = 0; side < 2; side++) {
                    vector<BNode> trial = make_insertion_trial(cur, order[pos], p, side == 0);
                    double area = eval_hb_tree(trial);
                    if (area < best_area) {
                        best_area = area;
                        best_tree.swap(trial);
                    }
                }
            }
            cur.swap(best_tree);
        }
        return cur;
    };

    auto build_legacy_hb = [&]() {
        vector<int> ids = hb_ids;
        vector<BNode> tree(ids.size(), BNode());
        clear_tree_links(tree);
        if (ids.size() > 24)
            sort(ids.begin(), ids.end(), hb_before);
        for (int i = 0; i < (int)ids.size(); i++) tree[i].id = ids[i];
        if (tree.size() > 24) {
            vector<int> hb_nodes(tree.size());
            iota(hb_nodes.begin(), hb_nodes.end(), 0);
            build_complete_binary_tree(tree, hb_nodes);
        } else {
            for (int i = 0; i < (int)tree.size(); i++) {
                tree[i].par = (i > 0) ? i - 1 : -1;
                tree[i].left = (i + 1 < (int)tree.size()) ? i + 1 : -1;
                tree[i].right = -1;
            }
        }
        return tree;
    };

    vector<vector<int>> orders;
    auto add_order = [&](vector<int> order) {
        for (const auto& old : orders)
            if (old == order) return;
        orders.push_back(std::move(order));
    };

    add_order(hb_ids);

    vector<int> order = hb_ids;
    sort(order.begin(), order.end(), hb_before);
    add_order(order);

    order = hb_ids;
    sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        pair<double, double> ldim = dims_of(lhs);
        pair<double, double> rdim = dims_of(rhs);
        if (ldim.first != rdim.first) return ldim.first > rdim.first;
        return hb_before(lhs, rhs);
    });
    add_order(order);

    order = hb_ids;
    sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        pair<double, double> ldim = dims_of(lhs);
        pair<double, double> rdim = dims_of(rhs);
        if (ldim.second != rdim.second) return ldim.second > rdim.second;
        return hb_before(lhs, rhs);
    });
    add_order(order);

    order = hb_ids;
    sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        pair<double, double> ldim = dims_of(lhs);
        pair<double, double> rdim = dims_of(rhs);
        int lspan = max(ldim.first, ldim.second);
        int rspan = max(rdim.first, rdim.second);
        if (lspan != rspan) return lspan > rspan;
        return hb_before(lhs, rhs);
    });
    add_order(order);

    order = hb_ids;
    sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        pair<double, double> ldim = dims_of(lhs);
        pair<double, double> rdim = dims_of(rhs);
        double lratio = max(ldim.first, ldim.second) * min(rdim.first, rdim.second);
        double rratio = max(rdim.first, rdim.second) * min(ldim.first, ldim.second);
        if (fabs(lratio - rratio) > EPS) return lratio > rratio;
        return hb_before(lhs, rhs);
    });
    add_order(order);

    for (int trial = 0; trial < 2; trial++) {
        order = hb_ids;
        mt19937 rng(seed + trial);
        shuffle(order.begin(), order.end(), rng);
        add_order(order);
    }

    double best_area = numeric_limits<double>::infinity();
    vector<BNode> best_hb;

    vector<BNode> legacy_hb = build_legacy_hb();
    best_area = eval_hb_tree(legacy_hb);
    best_hb.swap(legacy_hb);

    if (hb_ids.size() <= 40) {
        for (const auto& candidate_order : orders) {
            vector<BNode> candidate = build_greedy_hb(candidate_order);
            double area = eval_hb_tree(candidate);
            if (area < best_area) {
                best_area = area;
                best_hb.swap(candidate);
            }
        }
    }

    hb.swap(best_hb);
    group_dirty.assign(ng, 1);
}

// ─── Large-case greedy stitching ──────────────────────────────────────────

void Floorplan::greedy_pack_large(unsigned int seed) {
    int nm = mods.size();
    int ng = grps.size();

    hb.clear();
    asf.assign(ng, vector<BNode>());
    asf_pos.assign(ng, vector<int>());
    asf_rw_buf.assign(ng, vector<double>());
    asf_rh_buf.assign(ng, vector<double>());
    asf_rx_buf.assign(ng, vector<double>());
    asf_ry_buf.assign(ng, vector<double>());
    iw.assign(ng, 0);
    ih.assign(ng, 0);
    island_x.assign(ng, 0);
    island_y.assign(ng, 0);
    island_top.assign(ng, vector<ProfileSeg>());
    island_bottom.assign(ng, vector<ProfileSeg>());
    group_dirty.assign(ng, 0);

    vector<bool> in_grp(nm, false);
    for (int g = 0; g < ng; g++) {
        for (const auto& pp : grps[g].pairs) {
            in_grp[pp.a] = true;
            in_grp[pp.b] = true;
        }
        for (int m : grps[g].selfsym) in_grp[m] = true;
    }

    auto axis_dims = [&](int mid, SymType type) {
        if (type == SYM_VERT) return pair<double, double>{mods[mid].w, mods[mid].h};
        return pair<double, double>{mods[mid].h, mods[mid].w};
    };

    auto axis_dims_oriented = [&](int mid, SymType type, bool rot) {
        double w = mods[mid].w;
        double h = mods[mid].h;
        if (rot) swap(w, h);
        if (type == SYM_VERT) return pair<double, double>{w, h};
        return pair<double, double>{h, w};
    };

    auto add_actual_rect = [&](vector<Rect>& rects, int mid, double u, double v,
                               SymType type) {
        if (type == SYM_VERT) {
            mods[mid].x = u;
            mods[mid].y = v;
        } else {
            mods[mid].x = v;
            mods[mid].y = u;
        }
        rects.push_back(Rect(mods[mid].x, mods[mid].y, mods[mid].w, mods[mid].h));
    };

    auto collect_candidates = [](const Contour& c, double max_u, double width,
                                 double target_u = -1) {
        vector<double> xs;
        xs.reserve(c.seg.size() + 2);
        xs.push_back(0);
        xs.push_back(max_u);
        if (target_u >= width) xs.push_back(target_u - width);
        for (const auto& kv : c.seg) {
            if (kv.first >= 1e90) break;
            xs.push_back(kv.first);
            xs.push_back(kv.first - width);
        }
        double limit = max(max_u, target_u);
        if (limit < 0) limit = max_u;
        xs.erase(remove_if(xs.begin(), xs.end(),
                           [&](double x) { return x < -EPS || x > limit + EPS; }),
                 xs.end());
        sort(xs.begin(), xs.end());
        xs.erase(unique(xs.begin(), xs.end()), xs.end());
        return xs;
    };

    auto collect_profile_candidates = [&](const Contour& c, double max_u,
                                          const vector<ProfileSeg>& bottom,
                                          double width, double target_u = -1) {
        vector<double> xs = collect_candidates(c, max_u, width, target_u);
        double limit = max(max_u, target_u);
        if (limit < 0) limit = max_u;
        for (const auto& kv : c.seg) {
            if (kv.first >= 1e90) break;
            for (const auto& seg : bottom) {
                xs.push_back(kv.first - seg.l);
                xs.push_back(kv.first - seg.r);
            }
        }
        xs.erase(remove_if(xs.begin(), xs.end(),
                           [&](double x) { return x < -EPS || x > limit + EPS; }),
                 xs.end());
        sort(xs.begin(), xs.end());
        xs.erase(unique(xs.begin(), xs.end()), xs.end());
        return xs;
    };

    struct AxisItem {
        int idx = 0;
        double u = 0;
        double v = 0;
        double du = 0;
        double dv = 0;
        double left = 0;
        bool rot = false;
    };

    auto stitch_group = [&](int g, SymType type, int mode, double target_span) {
        SymGroup& grp = grps[g];
        grp.type = type;
        Contour half;
        vector<AxisItem> pair_place;
        vector<AxisItem> self_place;
        vector<Rect> rects;
        rects.reserve(grp.pairs.size() * 2 + grp.selfsym.size());

        double max_left = 0;
        double max_right = 0;
        double max_v = 0;

        vector<int> self_order(grp.selfsym.size());
        iota(self_order.begin(), self_order.end(), 0);
        sort(self_order.begin(), self_order.end(), [&](int lhs, int rhs) {
            auto ld = axis_dims(grp.selfsym[lhs], type);
            auto rd = axis_dims(grp.selfsym[rhs], type);
            double la = ld.first * ld.second;
            double ra = rd.first * rd.second;
            if (fabs(la - ra) > EPS) return la > ra;
            if (fabs(ld.second - rd.second) > EPS) return ld.second > rd.second;
            return ld.first > rd.first;
        });

        auto place_self = [&](int self_idx) {
            int mid = grp.selfsym[self_idx];
            AxisItem item;
            double best_score = numeric_limits<double>::infinity();
            double best_v = 0;
            for (int rot = 0; rot < 2; rot++) {
                if (rot && fabs(mods[mid].w - mods[mid].h) <= EPS) continue;
                auto dim = axis_dims_oriented(mid, type, rot != 0);
                double du = dim.first / 2.0;
                double dv = dim.second;
                double v = half.query(0, du);
                double left = dim.first - du;
                double next_left = max(max_left, left);
                double next_right = max(max_right, du);
                double next_v = max(max_v, v + dv);
                double next_span = next_left + next_right;
                double overflow = (target_span > 0) ? max(0.0, next_span - target_span) : 0.0;
                double score = (target_span > 0)
                                   ? (1e9 * overflow + 1e6 * next_v + next_span)
                                   : (next_span * next_v);
                if (score < best_score - EPS ||
                    (fabs(score - best_score) <= EPS && v < best_v)) {
                    best_score = score;
                    best_v = v;
                    item.idx = self_idx;
                    item.u = 0;
                    item.v = v;
                    item.du = du;
                    item.dv = dv;
                    item.left = left;
                    item.rot = rot != 0;
                }
            }
            if (item.rot) rotate_module(mods[mid]);
            half.update(0, item.du, item.v + item.dv);

            self_place.push_back(item);

            max_left = max(max_left, item.left);
            max_right = max(max_right, item.du);
            max_v = max(max_v, item.v + item.dv);
        };

        vector<int> pair_order(grp.pairs.size());
        iota(pair_order.begin(), pair_order.end(), 0);
        sort(pair_order.begin(), pair_order.end(), [&](int lhs, int rhs) {
            int lmid = (grp.pairs[lhs].rep == 0) ? grp.pairs[lhs].a
                                                 : grp.pairs[lhs].b;
            int rmid = (grp.pairs[rhs].rep == 0) ? grp.pairs[rhs].a
                                                 : grp.pairs[rhs].b;
            auto ld = axis_dims(lmid, type);
            auto rd = axis_dims(rmid, type);
            double la = ld.first * ld.second;
            double ra = rd.first * rd.second;
            if (fabs(la - ra) > EPS) return la > ra;
            double lspan = max(ld.first, ld.second);
            double rspan = max(rd.first, rd.second);
            if (fabs(lspan - rspan) > EPS) return lspan > rspan;
            return lhs < rhs;
        });

        auto place_pair = [&](int pair_idx) {
            SymPair& pp = grp.pairs[pair_idx];
            int mid = (pp.rep == 0) ? pp.a : pp.b;

            double best_score = numeric_limits<double>::infinity();
            double best_u = 0;
            double best_v = 0;
            AxisItem item;
            for (int rot = 0; rot < 2; rot++) {
                if (rot && fabs(mods[mid].w - mods[mid].h) <= EPS) continue;
                auto dim = axis_dims_oriented(mid, type, rot != 0);
                double du = dim.first;
                double dv = dim.second;
                double target_half = (target_span > 0) ? (target_span / 2.0) : -1;
                vector<double> candidates = collect_candidates(half, max_right, du, target_half);
                for (double cand_u : candidates) {
                    double cand_v = half.query(cand_u, du);
                    double next_left = max(max_left, cand_u + du);
                    double next_right = max(max_right, cand_u + du);
                    double next_v = max(max_v, cand_v + dv);
                    double next_span = next_left + next_right;
                    double overflow = (target_span > 0) ? max(0.0, next_span - target_span) : 0.0;
                    double score = (target_span > 0)
                                       ? (1e9 * overflow + 1e6 * next_v + next_span)
                                       : (next_span * next_v);
                    if (score < best_score - EPS ||
                        (fabs(score - best_score) <= EPS &&
                         (cand_v < best_v || (fabs(cand_v - best_v) <= EPS && cand_u < best_u)))) {
                        best_score = score;
                        best_u = cand_u;
                        best_v = cand_v;
                        item.idx = pair_idx;
                        item.u = best_u;
                        item.v = best_v;
                        item.du = du;
                        item.dv = dv;
                        item.rot = rot != 0;
                    }
                }
            }

            if (item.rot) {
                rotate_module(mods[pp.a]);
                rotate_module(mods[pp.b]);
            }
            half.update(item.u, item.du, item.v + item.dv);
            pair_place.push_back(item);

            max_left = max(max_left, item.u + item.du);
            max_right = max(max_right, item.u + item.du);
            max_v = max(max_v, item.v + item.dv);
        };

        if (mode == 0) {
            for (int self_idx : self_order) place_self(self_idx);
            for (int pair_idx : pair_order) place_pair(pair_idx);
        } else if (mode == 1) {
            for (int pair_idx : pair_order) place_pair(pair_idx);
            for (int self_idx : self_order) place_self(self_idx);
        } else {
            struct GroupItem {
                bool self = false;
                int idx = 0;
                double w = 0;
                double h = 0;
                double area = 0;
            };
            vector<GroupItem> items;
            items.reserve(grp.pairs.size() + grp.selfsym.size());
            for (int self_idx : self_order) {
                int mid = grp.selfsym[self_idx];
                auto dim = axis_dims(mid, type);
                GroupItem item;
                item.self = true;
                item.idx = self_idx;
                item.w = dim.first;
                item.h = dim.second;
                item.area = dim.first * dim.second;
                items.push_back(item);
            }
            for (int pair_idx : pair_order) {
                const SymPair& pp = grp.pairs[pair_idx];
                int mid = (pp.rep == 0) ? pp.a : pp.b;
                auto dim = axis_dims(mid, type);
                GroupItem item;
                item.self = false;
                item.idx = pair_idx;
                item.w = dim.first;
                item.h = dim.second;
                item.area = dim.first * dim.second;
                items.push_back(item);
            }

            if (mode < 6) {
                sort(items.begin(), items.end(), [&](const GroupItem& lhs,
                                                     const GroupItem& rhs) {
                    if (mode == 2) {
                        if (fabs(lhs.area - rhs.area) > EPS) return lhs.area > rhs.area;
                    } else if (mode == 3) {
                        if (fabs(lhs.h - rhs.h) > EPS) return lhs.h > rhs.h;
                    } else if (mode == 4) {
                        if (fabs(lhs.w - rhs.w) > EPS) return lhs.w > rhs.w;
                    } else {
                        double lspan = max(lhs.w, lhs.h);
                        double rspan = max(rhs.w, rhs.h);
                        if (fabs(lspan - rspan) > EPS) return lspan > rspan;
                    }
                    double lmin = min(lhs.w, lhs.h);
                    double rmin = min(rhs.w, rhs.h);
                    if (fabs(lmin - rmin) > EPS) return lmin > rmin;
                    if (lhs.self != rhs.self) return !lhs.self;
                    return lhs.idx < rhs.idx;
                });
            } else {
                mt19937 rng(seed + 100003 * (g + 1) + 9176 * (int)type + mode);
                shuffle(items.begin(), items.end(), rng);
            }

            for (const auto& item : items) {
                if (item.self) place_self(item.idx);
                else           place_pair(item.idx);
            }
        }

        double axis = max_left;
        for (const auto& item : pair_place) {
            SymPair& pp = grp.pairs[item.idx];
            int mid_rep = (pp.rep == 0) ? pp.a : pp.b;
            int mid_mir = (pp.rep == 0) ? pp.b : pp.a;
            mods[mid_mir].rot = mods[mid_rep].rot;
            mods[mid_mir].w = mods[mid_rep].w;
            mods[mid_mir].h = mods[mid_rep].h;

            add_actual_rect(rects, mid_rep, axis + item.u, item.v, type);
            add_actual_rect(rects, mid_mir, axis - item.u - item.du, item.v, type);
        }

        for (const auto& item : self_place) {
            int mid = grp.selfsym[item.idx];
            add_actual_rect(rects, mid, axis - item.left, item.v, type);
        }

        if (type == SYM_VERT) {
            iw[g] = max_left + max_right;
            ih[g] = max_v;
        } else {
            iw[g] = max_v;
            ih[g] = max_left + max_right;
        }

        pair<vector<ProfileSeg>, vector<ProfileSeg>> profiles = extract_profiles(rects);
        island_bottom[g] = std::move(profiles.first);
        island_top[g] = std::move(profiles.second);
        island_x[g] = 0;
        island_y[g] = 0;
        return iw[g] * ih[g];
    };

    // Build each symmetry group as a directly stitched island.  Try both legal
    // symmetry axes and keep the smaller island.
    for (int g = 0; g < ng; g++) {
        vector<Module> before = mods;
        vector<Module> best_mods;
        vector<ProfileSeg> best_bottom;
        vector<ProfileSeg> best_top;
        SymType best_type = SYM_VERT;
        double best_w = 0;
        double best_h = 0;
        double best_area = numeric_limits<double>::infinity();

        for (int pass = 0; pass < 2; pass++) {
            SymType type = (pass == 0) ? SYM_VERT : SYM_HORIZ;
            double group_area = 0;
            for (const auto& pp : grps[g].pairs) {
                group_area += before[pp.a].w * before[pp.a].h;
                group_area += before[pp.b].w * before[pp.b].h;
            }
            for (int m : grps[g].selfsym)
                group_area += before[m].w * before[m].h;
            double base_span = max(1.0, sqrt(group_area));
            vector<double> group_targets;
            group_targets.push_back(0);
            for (int pct : {65, 75, 85, 95, 105, 115, 130, 150, 175}) {
                double target = max(1.0, base_span * pct / 100.0);
                if (find(group_targets.begin(), group_targets.end(), target) == group_targets.end())
                    group_targets.push_back(target);
            }

            for (int mode = 0; mode < 40; mode++) {
                for (double target_span : group_targets) {
                    mods = before;
                    double area = stitch_group(g, type, mode, target_span);
                    double span = max(iw[g], ih[g]);
                    double best_span = max(best_w, best_h);
                    if (area < best_area - EPS ||
                        (fabs(area - best_area) <= EPS &&
                         (!isfinite(best_area) || span < best_span))) {
                        best_area = area;
                        best_type = type;
                        best_w = iw[g];
                        best_h = ih[g];
                        best_bottom = island_bottom[g];
                        best_top = island_top[g];
                        best_mods = mods;
                    }
                }
            }
        }

        mods = std::move(best_mods);
        grps[g].type = best_type;
        iw[g] = best_w;
        ih[g] = best_h;
        island_bottom[g] = std::move(best_bottom);
        island_top[g] = std::move(best_top);
    }

    vector<Module> base_mods = mods;
    vector<int> item_ids;
    item_ids.reserve(nm + ng);
    for (int g = 0; g < ng; g++) {
        if (iw[g] > 0 && ih[g] > 0) item_ids.push_back(-(g + 1));
    }
    for (int m = 0; m < nm; m++) {
        if (!in_grp[m]) item_ids.push_back(m);
    }
    if (item_ids.empty()) {
        packed_w = packed_h = 0;
        packed_area = 0;
        return;
    }

    auto item_dims = [&](int id) {
        if (id >= 0) return pair<double, double>{base_mods[id].w, base_mods[id].h};
        int g = -(id + 1);
        return pair<double, double>{iw[g], ih[g]};
    };

    auto item_before = [&](int lhs, int rhs) {
        auto ld = item_dims(lhs);
        auto rd = item_dims(rhs);
        double la = ld.first * ld.second;
        double ra = rd.first * rd.second;
        if (fabs(la - ra) > EPS) return la > ra;
        double lspan = max(ld.first, ld.second);
        double rspan = max(rd.first, rd.second);
        if (fabs(lspan - rspan) > EPS) return lspan > rspan;
        double lmin = min(ld.first, ld.second);
        double rmin = min(rd.first, rd.second);
        if (fabs(lmin - rmin) > EPS) return lmin > rmin;
        return lhs < rhs;
    };

    vector<vector<int>> orders;
    auto add_order = [&](vector<int> order) {
        for (const auto& old : orders)
            if (old == order) return;
        orders.push_back(std::move(order));
    };

    vector<int> order = item_ids;
    sort(order.begin(), order.end(), item_before);
    add_order(order);

    order = item_ids;
    sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        auto ld = item_dims(lhs);
        auto rd = item_dims(rhs);
        if (fabs(ld.second - rd.second) > EPS) return ld.second > rd.second;
        return item_before(lhs, rhs);
    });
    add_order(order);

    order = item_ids;
    sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        auto ld = item_dims(lhs);
        auto rd = item_dims(rhs);
        if (fabs(ld.first - rd.first) > EPS) return ld.first > rd.first;
        return item_before(lhs, rhs);
    });
    add_order(order);

    order = item_ids;
    sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        auto ld = item_dims(lhs);
        auto rd = item_dims(rhs);
        double lratio = max(ld.first, ld.second) * min(rd.first, rd.second);
        double rratio = max(rd.first, rd.second) * min(ld.first, ld.second);
        if (fabs(lratio - rratio) > EPS) return lratio > rratio;
        return item_before(lhs, rhs);
    });
    add_order(order);

    order = item_ids;
    sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        auto ld = item_dims(lhs);
        auto rd = item_dims(rhs);
        if (fabs(ld.second - rd.second) > EPS) return ld.second < rd.second;
        return item_before(lhs, rhs);
    });
    add_order(order);

    order = item_ids;
    sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        auto ld = item_dims(lhs);
        auto rd = item_dims(rhs);
        if (fabs(ld.first - rd.first) > EPS) return ld.first < rd.first;
        return item_before(lhs, rhs);
    });
    add_order(order);

    order = item_ids;
    sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        auto ld = item_dims(lhs);
        auto rd = item_dims(rhs);
        double lratio = max(ld.first, ld.second) * min(rd.first, rd.second);
        double rratio = max(rd.first, rd.second) * min(ld.first, ld.second);
        if (fabs(lratio - rratio) > EPS) return lratio < rratio;
        return item_before(lhs, rhs);
    });
    add_order(order);

    for (int trial = 0; trial < 16; trial++) {
        order = item_ids;
        mt19937 rng(seed + 1009 * trial);
        shuffle(order.begin(), order.end(), rng);
        int sorted_prefix = (int)order.size() * ((trial % 4) + 1) / 5;
        stable_sort(order.begin(), order.begin() + sorted_prefix, item_before);
        add_order(order);
    }

    double total_module_area = 0;
    double widest_item = 1;
    for (int id : item_ids) {
        auto dim = item_dims(id);
        widest_item = max(widest_item, dim.first);
        widest_item = max(widest_item, dim.second);
    }
    for (const auto& mod : mods) total_module_area += mod.w * mod.h;
    double base_target_w = max(widest_item, sqrt(total_module_area));
    vector<double> target_widths;
    for (int pct : {90, 110, 130, 150, 175, 190, 205, 210}) {
        double target = max(widest_item, base_target_w * pct / 100.0);
        for (int delta : {-3, -2, -1, 0, 1, 2, 3}) {
            double near_target = max(widest_item, target + delta);
            if (find(target_widths.begin(), target_widths.end(), near_target) == target_widths.end())
                target_widths.push_back(near_target);
        }
    }

    auto place_order = [&](const vector<int>& cur_order, double target_w) {
        mods = base_mods;
        island_x.assign(ng, 0);
        island_y.assign(ng, 0);
        Contour c;
        double max_x = 0;
        double max_y = 0;

        for (int id : cur_order) {
            struct Candidate {
                double overflow = numeric_limits<double>::infinity();
                double next_h = numeric_limits<double>::infinity();
                double next_w = numeric_limits<double>::infinity();
                double x = 0;
                double y = 0;
                double w = 0;
                double h = 0;
                bool rot = false;
            };

            Candidate best;
            auto take_candidate = [&](double x, double y, double w, double h, bool rot) {
                double next_w = max(max_x, x + w);
                double next_h = max(max_y, y + h);
                double overflow = max(0.0, x + w - target_w);
                bool better = false;
                if (fabs(overflow - best.overflow) > EPS) {
                    better = overflow < best.overflow;
                } else if (fabs(next_h - best.next_h) > EPS) {
                    better = next_h < best.next_h;
                } else if (fabs(y - best.y) > EPS) {
                    better = y < best.y;
                } else if (fabs(next_w - best.next_w) > EPS) {
                    better = next_w < best.next_w;
                } else {
                    better = x < best.x;
                }

                if (better) {
                    best.overflow = overflow;
                    best.next_h = next_h;
                    best.next_w = next_w;
                    best.x = x;
                    best.y = y;
                    best.w = w;
                    best.h = h;
                    best.rot = rot;
                }
            };

            if (id >= 0) {
                for (int rot = 0; rot < 2; rot++) {
                    double w = rot ? mods[id].h : mods[id].w;
                    double h = rot ? mods[id].w : mods[id].h;
                    if (rot && fabs(w - h) <= EPS) continue;
                    vector<double> rect_candidates = collect_candidates(c, max_x, w, target_w);
                    for (double x : rect_candidates) {
                        double y = c.query(x, w);
                        take_candidate(x, y, w, h, rot != 0);
                    }
                }

                if (best.rot) rotate_module(mods[id]);
                mods[id].x = best.x;
                mods[id].y = best.y;
                c.update(best.x, best.w, best.y + best.h);
                max_x = max(max_x, best.x + best.w);
                max_y = max(max_y, best.y + best.h);
            } else {
                int g = -(id + 1);
                double w = iw[g];
                double h = ih[g];
                vector<double> profile_candidates =
                    collect_profile_candidates(c, max_x, island_bottom[g], w, target_w);
                for (double x : profile_candidates) {
                    double y = query_profile(c, x, island_bottom[g]);
                    take_candidate(x, y, w, h, false);
                }

                double dx = best.x - island_x[g];
                double dy = best.y - island_y[g];
                for (const auto& pp : grps[g].pairs) {
                    mods[pp.a].x += dx;
                    mods[pp.a].y += dy;
                    mods[pp.b].x += dx;
                    mods[pp.b].y += dy;
                }
                for (int m : grps[g].selfsym) {
                    mods[m].x += dx;
                    mods[m].y += dy;
                }
                island_x[g] = best.x;
                island_y[g] = best.y;
                update_profile(c, best.x, island_top[g], best.y);
                max_x = max(max_x, best.x + best.w);
                max_y = max(max_y, best.y + best.h);
            }
        }

        packed_w = max_x;
        packed_h = max_y;
        packed_area = packed_w * packed_h;
        return packed_area;
    };

    double best_area = numeric_limits<double>::infinity();
    double best_w = 0;
    double best_h = 0;
    vector<Module> best_mods;
    vector<double> best_island_x;
    vector<double> best_island_y;

    for (double target_w : target_widths) {
        for (const auto& candidate_order : orders) {
            double area = place_order(candidate_order, target_w);
            double span = max(packed_w, packed_h);
            double best_span = max(best_w, best_h);
            if (area < best_area - EPS ||
                (fabs(area - best_area) <= EPS &&
                 (!isfinite(best_area) || span < best_span))) {
                best_area = area;
                best_w = packed_w;
                best_h = packed_h;
                best_mods = mods;
                best_island_x = island_x;
                best_island_y = island_y;
            }
        }
    }

    mods = std::move(best_mods);
    island_x = std::move(best_island_x);
    island_y = std::move(best_island_y);
    packed_w = best_w;
    packed_h = best_h;
    packed_area = best_area;

    auto members_of = [&](int id) {
        vector<int> members;
        if (id >= 0) {
            members.push_back(id);
        } else {
            int g = -(id + 1);
            for (const auto& pp : grps[g].pairs) {
                members.push_back(pp.a);
                members.push_back(pp.b);
            }
            for (int m : grps[g].selfsym) members.push_back(m);
        }
        return members;
    };

    auto move_item = [&](int id, double dx, double dy) {
        if (fabs(dx) <= EPS && fabs(dy) <= EPS) return;
        vector<int> members = members_of(id);
        for (int m : members) {
            mods[m].x += dx;
            mods[m].y += dy;
        }
        if (id < 0) {
            int g = -(id + 1);
            island_x[g] += dx;
            island_y[g] += dy;
        }
    };

    auto recompute_bbox = [&]() {
        double max_x = 0;
        double max_y = 0;
        for (const auto& mod : mods) {
            max_x = max(max_x, mod.x + mod.w);
            max_y = max(max_y, mod.y + mod.h);
        }
        packed_w = max_x;
        packed_h = max_y;
        packed_area = packed_w * packed_h;
    };

    auto item_min_x = [&](int id) {
        double v = numeric_limits<double>::infinity();
        for (int m : members_of(id)) v = min(v, mods[m].x);
        return v;
    };
    auto item_min_y = [&](int id) {
        double v = numeric_limits<double>::infinity();
        for (int m : members_of(id)) v = min(v, mods[m].y);
        return v;
    };

    for (int pass = 0; pass < 3; pass++) {
        vector<int> order_x = item_ids;
        sort(order_x.begin(), order_x.end(), [&](int lhs, int rhs) {
            double lx = item_min_x(lhs);
            double rx = item_min_x(rhs);
            if (fabs(lx - rx) > EPS) return lx < rx;
            return lhs < rhs;
        });

        for (int id : order_x) {
            vector<int> members = members_of(id);
            vector<char> moving(nm, 0);
            double min_x = numeric_limits<double>::infinity();
            for (int m : members) {
                moving[m] = 1;
                min_x = min(min_x, mods[m].x);
            }

            double dx = -min_x;
            for (int m : members) {
                const Module& a = mods[m];
                for (int o = 0; o < nm; o++) {
                    if (moving[o]) continue;
                    const Module& b = mods[o];
                    bool y_overlap = max(a.y, b.y) < min(a.y + a.h, b.y + b.h) - EPS;
                    if (y_overlap && b.x + b.w <= a.x + EPS) {
                        dx = max(dx, b.x + b.w - a.x);
                    }
                }
            }
            if (dx < -EPS) move_item(id, dx, 0);
        }

        vector<int> order_y = item_ids;
        sort(order_y.begin(), order_y.end(), [&](int lhs, int rhs) {
            double ly = item_min_y(lhs);
            double ry = item_min_y(rhs);
            if (fabs(ly - ry) > EPS) return ly < ry;
            return lhs < rhs;
        });

        for (int id : order_y) {
            vector<int> members = members_of(id);
            vector<char> moving(nm, 0);
            double min_y = numeric_limits<double>::infinity();
            for (int m : members) {
                moving[m] = 1;
                min_y = min(min_y, mods[m].y);
            }

            double dy = -min_y;
            for (int m : members) {
                const Module& a = mods[m];
                for (int o = 0; o < nm; o++) {
                    if (moving[o]) continue;
                    const Module& b = mods[o];
                    bool x_overlap = max(a.x, b.x) < min(a.x + a.w, b.x + b.w) - EPS;
                    if (x_overlap && b.y + b.h <= a.y + EPS) {
                        dy = max(dy, b.y + b.h - a.y);
                    }
                }
            }
            if (dy < -EPS) move_item(id, 0, dy);
        }
    }
    recompute_bbox();
}

// ─── ASF packing ──────────────────────────────────────────────────────────

void Floorplan::pack_asf(int g) {
    SymGroup& grp = grps[g];
    vector<BNode>& nd = asf[g];
    int n  = nd.size();
    int np = (int)grp.pairs.size();
    int ns = (int)grp.selfsym.size();
    if (n == 0) {
        iw[g] = ih[g] = 0;
        island_top[g].clear();
        island_bottom[g].clear();
        return;
    }

    bool horiz = (grp.type == SYM_HORIZ);

    // Rep dimensions (in packing space, before transposition for horiz)
    vector<double>& rw = asf_rw_buf[g];
    vector<double>& rh = asf_rh_buf[g];
    for (int j = 0; j < np; j++) {
        int mid = (grp.pairs[j].rep == 0) ? grp.pairs[j].a : grp.pairs[j].b;
        rw[j] = mods[mid].w;
        rh[j] = mods[mid].h;
    }
    for (int j = 0; j < ns; j++) {
        int idx = np + j;
        int mid = grp.selfsym[j];
        if (!horiz) {
            rw[idx] = mods[mid].w / 2.0;
            rh[idx] = mods[mid].h;
        } else {
            rw[idx] = mods[mid].w;
            rh[idx] = mods[mid].h / 2.0;
        }
    }

    // Transpose for horizontal packing (swap w↔h, pack as vertical, swap x↔y after)
    if (horiz)
        for (int j = 0; j < n; j++) swap(rw[j], rh[j]);

    // Pack representative B*-tree
    int root = asf_root(g);
    Contour c;
    vector<double>& rx = asf_rx_buf[g];
    vector<double>& ry = asf_ry_buf[g];

    // Preorder traversal via stack (push right before left so left is processed first)
    stack<int> stk;
    stk.push(root);
    while (!stk.empty()) {
        int i = stk.top(); stk.pop();
        int j = nd[i].id; // rep index

        double x;
        if (nd[i].par == -1) {
            x = 0;
        } else {
            int p = nd[i].par;
            int pj = nd[p].id;
            x = (nd[p].left == i) ? (rx[p] + rw[pj]) : rx[p];
        }
        double y = c.query(x, rw[j]);
        rx[i] = x; ry[i] = y;
        nd[i].pw = rw[j]; nd[i].ph = rh[j];
        c.update(x, rw[j], y + rh[j]);

        if (nd[i].right != -1) stk.push(nd[i].right);
        if (nd[i].left  != -1) stk.push(nd[i].left);
    }

    // Untranspose
    if (horiz) {
        for (int i = 0; i < n; i++) swap(rx[i], ry[i]);
        for (int j = 0; j < n; j++) swap(rw[j], rh[j]);
    }

    // Compute shift so all module positions are non-negative
    // vertical: axis at x=0; horizontal: axis at y=0
    double shift = 0;
    for (int i = 0; i < n; i++) {
        int j = nd[i].id;
        bool sel = (j >= np);
        if (!horiz) {
            double neg_x = sel ? -rw[j] : -(rx[i] + rw[j]);
            shift = max(shift, -neg_x);
        } else {
            double neg_y = sel ? -rh[j] : -(ry[i] + rh[j]);
            shift = max(shift, -neg_y);
        }
    }

    // Assign final module coordinates (relative to island origin (0,0))
    double max_x = 0, max_y = 0;
    vector<Rect> rects;
    rects.reserve(np * 2 + ns);
    for (int i = 0; i < n; i++) {
        int j = nd[i].id;
        bool sel = (j >= np);

        if (!horiz) {
            int mid_rep, mid_mir = -1;
            if (!sel) {
                mid_rep = (grp.pairs[j].rep == 0) ? grp.pairs[j].a : grp.pairs[j].b;
                mid_mir = (grp.pairs[j].rep == 0) ? grp.pairs[j].b : grp.pairs[j].a;
            } else {
                mid_rep = grp.selfsym[j - np];
            }

            if (!sel) {
                mods[mid_rep].x = rx[i] + shift;
                mods[mid_rep].y = ry[i];
                mods[mid_mir].x = -(rx[i] + rw[j]) + shift;
                mods[mid_mir].y = ry[i];
                mods[mid_mir].w = mods[mid_rep].w;
                mods[mid_mir].h = mods[mid_rep].h;
                mods[mid_mir].rot = mods[mid_rep].rot;
                rects.push_back(Rect(
                    mods[mid_mir].x, mods[mid_mir].y, mods[mid_mir].w, mods[mid_mir].h
                ));
                max_x = max(max_x, mods[mid_mir].x + mods[mid_mir].w);
                max_y = max(max_y, mods[mid_mir].y + mods[mid_mir].h);
            } else {
                // Full selfsym module: BL so that right half aligns with axis
                mods[mid_rep].x = -(mods[mid_rep].w - rw[j]) + shift;
                mods[mid_rep].y = ry[i];
            }
            rects.push_back(Rect(
                mods[mid_rep].x, mods[mid_rep].y, mods[mid_rep].w, mods[mid_rep].h
            ));
            max_x = max(max_x, mods[mid_rep].x + mods[mid_rep].w);
            max_y = max(max_y, mods[mid_rep].y + mods[mid_rep].h);
        } else {
            int mid_rep, mid_mir = -1;
            if (!sel) {
                mid_rep = (grp.pairs[j].rep == 0) ? grp.pairs[j].a : grp.pairs[j].b;
                mid_mir = (grp.pairs[j].rep == 0) ? grp.pairs[j].b : grp.pairs[j].a;
            } else {
                mid_rep = grp.selfsym[j - np];
            }

            if (!sel) {
                mods[mid_rep].x = rx[i];
                mods[mid_rep].y = ry[i] + shift;
                mods[mid_mir].x = rx[i];
                mods[mid_mir].y = -(ry[i] + rh[j]) + shift;
                mods[mid_mir].w = mods[mid_rep].w;
                mods[mid_mir].h = mods[mid_rep].h;
                mods[mid_mir].rot = mods[mid_rep].rot;
                rects.push_back(Rect(
                    mods[mid_mir].x, mods[mid_mir].y, mods[mid_mir].w, mods[mid_mir].h
                ));
                max_x = max(max_x, mods[mid_mir].x + mods[mid_mir].w);
                max_y = max(max_y, mods[mid_mir].y + mods[mid_mir].h);
            } else {
                mods[mid_rep].x = rx[i];
                mods[mid_rep].y = -(mods[mid_rep].h - rh[j]) + shift;
            }
            rects.push_back(Rect(
                mods[mid_rep].x, mods[mid_rep].y, mods[mid_rep].w, mods[mid_rep].h
            ));
            max_x = max(max_x, mods[mid_rep].x + mods[mid_rep].w);
            max_y = max(max_y, mods[mid_rep].y + mods[mid_rep].h);
        }
    }

    iw[g] = max_x;
    ih[g] = max_y;
    pair<vector<ProfileSeg>, vector<ProfileSeg>> profiles = extract_profiles(rects);
    island_bottom[g] = std::move(profiles.first);
    island_top[g] = std::move(profiles.second);
}

// ─── HB packing ──────────────────────────────────────────────────────────

void Floorplan::pack() {
    if (hb.empty()) {
        packed_w = packed_h = 0;
        packed_area = 0;
        return;
    }
    if (iw.size() != grps.size()) iw.resize(grps.size());
    if (ih.size() != grps.size()) ih.resize(grps.size());
    if (island_x.size() != grps.size()) island_x.assign(grps.size(), 0);
    if (island_y.size() != grps.size()) island_y.assign(grps.size(), 0);
    if (island_top.size() != grps.size()) island_top.resize(grps.size());
    if (island_bottom.size() != grps.size()) island_bottom.resize(grps.size());
    if (group_dirty.size() != grps.size()) group_dirty.assign(grps.size(), 1);

    Contour c;
    int root = hb_root();
    stack<int> stk;
    stk.push(root);
    double max_x = 0;
    double max_y = 0;

    while (!stk.empty()) {
        int i = stk.top(); stk.pop();
        int id = hb[i].id;

        // Compute x from parent
        double x;
        if (hb[i].par == -1) {
            x = 0;
        } else {
            int p = hb[i].par;
            x = (hb[p].left == i) ? (hb[p].px + hb[p].pw) : hb[p].px;
        }

        double w, h, y;
        if (id >= 0) {
            // Regular module
            w = mods[id].w; h = mods[id].h;
            y = c.query(x, w);
            mods[id].x = x; mods[id].y = y;
        } else {
            // Hierarchy node: pack ASF tree
            int g = -(id + 1);
            bool dirty = group_dirty[g];
            if (dirty) pack_asf(g);
            w = iw[g]; h = ih[g];
            if (w == 0) {
                group_dirty[g] = 0;
                island_x[g] = x;
                island_y[g] = 0;
                hb[i].px = x; hb[i].py = 0; hb[i].pw = 0; hb[i].ph = 0;
                goto push_children;
            }
            y = query_profile(c, x, island_bottom[g]);
            // Translate all group modules by (x, y)
            double dx = dirty ? x : (x - island_x[g]);
            double dy = dirty ? y : (y - island_y[g]);
            for (auto& pp : grps[g].pairs) {
                mods[pp.a].x += dx; mods[pp.a].y += dy;
                mods[pp.b].x += dx; mods[pp.b].y += dy;
            }
            for (int m : grps[g].selfsym) { mods[m].x += dx; mods[m].y += dy; }
            island_x[g] = x;
            island_y[g] = y;
            group_dirty[g] = 0;
        }

        hb[i].px = x; hb[i].py = y; hb[i].pw = w; hb[i].ph = h;
        if (id >= 0) c.update(x, w, y + h);
        else         update_profile(c, x, island_top[-(id + 1)], y);
        max_x = max(max_x, x + w);
        max_y = max(max_y, y + h);

        push_children:
        if (hb[i].right != -1) stk.push(hb[i].right);
        if (hb[i].left  != -1) stk.push(hb[i].left);
    }

    packed_w = max_x;
    packed_h = max_y;
    packed_area = packed_w * packed_h;
}

// ─── Cost ────────────────────────────────────────────────────────────────

double Floorplan::cost() const {
    return packed_area;
}

// ─── Save / Restore ──────────────────────────────────────────────────────

void Floorplan::save() {
    sv_mods = mods;
    sv_grps = grps;
    sv_hb   = hb;
    sv_asf  = asf;
    sv_asf_pos = asf_pos;
    sv_iw = iw;
    sv_ih = ih;
    sv_island_x = island_x;
    sv_island_y = island_y;
    sv_island_top = island_top;
    sv_island_bottom = island_bottom;
    sv_group_dirty = group_dirty;
    sv_packed_w = packed_w;
    sv_packed_h = packed_h;
    sv_packed_area = packed_area;
}

void Floorplan::restore() {
    mods.swap(sv_mods);
    grps.swap(sv_grps);
    hb.swap(sv_hb);
    asf.swap(sv_asf);
    asf_pos.swap(sv_asf_pos);
    iw.swap(sv_iw);
    ih.swap(sv_ih);
    island_x.swap(sv_island_x);
    island_y.swap(sv_island_y);
    island_top.swap(sv_island_top);
    island_bottom.swap(sv_island_bottom);
    group_dirty.swap(sv_group_dirty);
    swap(packed_w, sv_packed_w);
    swap(packed_h, sv_packed_h);
    swap(packed_area, sv_packed_area);
}

// ─── Perturbation operations ─────────────────────────────────────────────

void Floorplan::op_swap_hb(mt19937& rng) {
    if (hb.size() < 2) return;
    int n1 = rand_node(hb, rng);
    int n2 = rand_node(hb, rng);
    if (n1 == n2) return;
    swap(hb[n1].id, hb[n2].id);
}

void Floorplan::op_move_hb(mt19937& rng) {
    if (hb.size() < 2) return;
    int n = rand_node(hb, rng);
    detach(hb, n);
    int t = n;
    while (t == n) t = rand_node(hb, rng);
    bool lft = (rng() & 1);
    attach(hb, n, t, lft);
}

void Floorplan::op_swap_asf_pairs(int g, mt19937& rng) {
    int np = grps[g].pairs.size();
    if (np < 2) return;
    int aid = rng() % np;
    int bid = rng() % np;
    if (aid == bid) return;
    int a = asf_pos[g][aid];
    int b = asf_pos[g][bid];
    swap(asf[g][a].id, asf[g][b].id);
    swap(asf_pos[g][aid], asf_pos[g][bid]);
    group_dirty[g] = 1;
}

void Floorplan::op_swap_asf_sel(int g, mt19937& rng) {
    int np = grps[g].pairs.size();
    int ns = grps[g].selfsym.size();
    if (ns < 2) return;
    int aid = np + (rng() % ns);
    int bid = np + (rng() % ns);
    if (aid == bid) return;
    int a = asf_pos[g][aid];
    int b = asf_pos[g][bid];
    swap(asf[g][a].id, asf[g][b].id);
    swap(asf_pos[g][aid], asf_pos[g][bid]);
    group_dirty[g] = 1;
}

void Floorplan::op_move_asf_pair(int g, mt19937& rng) {
    if (asf[g].size() < 2) return;
    int np = grps[g].pairs.size();
    if (np == 0) return;
    int n = asf_pos[g][rng() % np];
    // Build valid targets (not n; and since we always attach as left child,
    // any node is safe — selfsym rightmost branch is unchanged)
    detach(asf[g], n);
    int t = n;
    while (t == n) t = rand_node(asf[g], rng);
    attach(asf[g], n, t, true); // always left child: keeps rightmost branch intact
    group_dirty[g] = 1;
}

void Floorplan::op_change_rep(int g, mt19937& rng) {
    if (grps[g].pairs.empty()) return;
    int idx = rng() % grps[g].pairs.size();
    grps[g].pairs[idx].rep ^= 1; // flip 0↔1
    group_dirty[g] = 1;
}

void Floorplan::op_rotate(int g, mt19937& rng) {
    // Randomly pick a pair or selfsym and rotate (swap w,h)
    auto& grp = grps[g];
    int total = grp.pairs.size() + grp.selfsym.size();
    if (total == 0) return;
    int choice = rng() % total;
    if (choice < (int)grp.pairs.size()) {
        int mid_a = grp.pairs[choice].a;
        int mid_b = grp.pairs[choice].b;
        rotate_module(mods[mid_a]);
        rotate_module(mods[mid_b]);
        group_dirty[g] = 1;
    } else {
        int mid = grp.selfsym[choice - grp.pairs.size()];
        rotate_module(mods[mid]);
        group_dirty[g] = 1;
    }
}

void Floorplan::op_flip_axis(int g) {
    SymType next = (grps[g].type == SYM_VERT) ? SYM_HORIZ : SYM_VERT;
    if (selfsym_type_feasible(grps[g], mods, next)) {
        grps[g].type = next;
        group_dirty[g] = 1;
    }
}

void Floorplan::op_flip_sym(int g) {
    grps[g].type = (grps[g].type == SYM_VERT) ? SYM_HORIZ : SYM_VERT;
    // Rotate all modules in group 90 degrees
    for (auto& pp : grps[g].pairs) {
        rotate_module(mods[pp.a]);
        rotate_module(mods[pp.b]);
    }
    for (int m : grps[g].selfsym) rotate_module(mods[m]);
    group_dirty[g] = 1;
}

// ─── Top-level perturb ────────────────────────────────────────────────────

void Floorplan::perturb(mt19937& rng) {
    int ng = grps.size();
    // Pick operation with weighted probabilities
    // Bucket: 0-29 swap_hb, 30-49 move_hb, 50-64 swap_asf_pairs,
    //         65-74 move_asf, 75-84 change_rep, 85-89 rotate, 90-94 swap_sel,
    //         95-96 flip_axis, 97-99 flip_sym
    int op = rng() % 100;

    if (op < 30) {
        op_swap_hb(rng);
    } else if (op < 50) {
        op_move_hb(rng);
    } else if (op < 65) {
        if (ng == 0) { op_swap_hb(rng); return; }
        int g = rng() % ng;
        op_swap_asf_pairs(g, rng);
    } else if (op < 75) {
        if (ng == 0) { op_swap_hb(rng); return; }
        int g = rng() % ng;
        op_move_asf_pair(g, rng);
    } else if (op < 85) {
        if (ng == 0) { op_swap_hb(rng); return; }
        int g = rng() % ng;
        op_change_rep(g, rng);
    } else if (op < 90) {
        if (ng == 0) { op_swap_hb(rng); return; }
        int g = rng() % ng;
        op_rotate(g, rng);
    } else if (op < 95) {
        if (ng == 0) { op_swap_hb(rng); return; }
        int g = rng() % ng;
        op_swap_asf_sel(g, rng);
    } else if (op < 97) {
        if (ng == 0) { op_swap_hb(rng); return; }
        int g = rng() % ng;
        op_flip_axis(g);
    } else {
        if (ng == 0) { op_swap_hb(rng); return; }
        int g = rng() % ng;
        op_flip_sym(g);
    }
}
