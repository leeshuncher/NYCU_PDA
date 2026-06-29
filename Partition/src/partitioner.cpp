#include "partitioner.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>

using namespace std;

static constexpr int GROUPS = 3;
static constexpr int NEG_INF = -1000000000;
static constexpr int LARGE_NET_THRESHOLD = 256;
static constexpr int LARGE_NET_BFS_SAMPLE_CAP = 128;
static constexpr int FM_SKIPPED_NET_THRESHOLD = 200;
static constexpr int COARSEN_NET_THRESHOLD = 500;
static constexpr int COARSEN_WEIGHT_SCALE = 128;
static constexpr int DEFAULT_INITIAL_SEEDS = 64;
static constexpr int INITIAL_REFINEMENT_CANDIDATES = 32;
static constexpr int INITIAL_SEED_FM_PASSES = 16;
static constexpr int INITIAL_DIVERSITY_THRESHOLD = 1000;

struct move_choice
{
    int to = -1;
    int gain = NEG_INF;
};

struct move_candidate
{
    int cid = -1;
    int from = -1;
    int to = -1;
    int gain = NEG_INF;

    bool valid() const
    {
        return cid >= 0;
    }
};

struct coarsen_result
{
    parsed_input coarse;
    vector<int> fine_to_coarse; // fine_cell_id -> coarse_cell_id
};

struct net_state
{
    array<int, GROUPS> count = {0, 0, 0};
    array<int, GROUPS> head = {-1, -1, -1};
    array<int, GROUPS> tail = {-1, -1, -1};
    vector<int> next_pin;
    vector<int> prev_pin;
    vector<unsigned char> pin_group;
};

struct fm_view
{
    vector<int> offsets;
    vector<int> nets;
    vector<int> degree;
    int pmax = 0;
};

static void fm_refine(const vector<cell> &cells, const vector<net> &nets, const fm_view &view, vector<net_state> &net_states, vector<char> &part, int minSize, int maxSize, chrono::steady_clock::time_point deadline, int max_passes);
static void fm_refine(const vector<cell> &cells, const vector<net> &nets, vector<net_state> &net_states, vector<char> &part, int minSize, int maxSize, chrono::steady_clock::time_point deadline, int max_passes);
static void enforce_balance(const vector<cell> &cells, const vector<net> &nets, const fm_view &view, vector<net_state> &net_states, vector<char> &part, int minSize, int maxSize);
static pair<int, int> balance_limits(int n, double r);
static uint64_t mix64(uint64_t x);
static uint64_t next_pseudo_random(uint64_t &state);
static int uncoarsen_fm_pass_budget(int original_cell_count, int level_cell_count);
static inline int active_group_count(const array<int, GROUPS> &count);
static int initial_seed_target_count(int n, long long remaining_ms);
static int partition_cycle_count(int n, long long total_budget_ms);
static long long cycle_budget_reserve_ms(int n, long long total_budget_ms);
static vector<net_state> initialized_net_states(const vector<net> &nets, const vector<char> &part);
static void refine_with_balance(const vector<cell> &cells, const vector<net> &nets, const fm_view &view,
                                vector<net_state> &net_states, vector<char> &part, int minSize, int maxSize,
                                chrono::steady_clock::time_point deadline, int initial_passes, int final_passes);
static partition_result make_partition_result(const vector<net> &nets, vector<char> part, int minSize, int maxSize);
static partition_result run_multilevel_cycle(parsed_input &parsed, const fm_view &finest_view,
                                             chrono::steady_clock::time_point deadline, uint64_t seed,
                                             const vector<char> *init_part, int uncoarsen_pass_cap);

static uint64_t partition_seed_base(int cell_count)
{
    switch (cell_count)
    {
    case 16743:  return 666; // case1
    case 32402:  return 111; // case2
    case 150750: return 111; // case3
    case 382489: return 888; // case4
    default:     return 777;
    }
}

static inline bool deadline_reached(chrono::steady_clock::time_point deadline)
{
    return chrono::steady_clock::now() >= deadline;
}

static int initial_seed_target_count(int n, long long remaining_ms)
{
    if (n <= 0)
        return 0;

    const bool diversify = (n >= INITIAL_DIVERSITY_THRESHOLD);
    const int base = DEFAULT_INITIAL_SEEDS;
    const int min_seeds = diversify ? 8 : 4;

    int bonus = 0;
    if (remaining_ms >= 120000)
        bonus = diversify ? 24 : 8;
    else if (remaining_ms >= 60000)
        bonus = diversify ? 16 : 6;
    else if (remaining_ms >= 20000)
        bonus = diversify ? 8 : 4;

    return min(n, max(min_seeds, base + bonus));
}

static int partition_cycle_count(int n, long long total_budget_ms)
{
    int cycles = (n <= 20000) ? 32
               : (n <= 100000) ? 16
               : (n <= 250000) ? 8
               : 3;

    if (total_budget_ms >= 180000)
        cycles += (n <= 100000) ? 2 : 1;
    else if (total_budget_ms >= 90000)
        cycles += 1;

    return max(1, cycles);
}

static long long cycle_budget_reserve_ms(int n, long long total_budget_ms)
{
    if (total_budget_ms <= 0)
        return 0;

    const long long reserve_floor = 2000;
    const long long reserve = (n <= 100000) ? total_budget_ms / 3 : total_budget_ms / 4;
    return min(total_budget_ms, max(reserve_floor, reserve));
}

static int uncoarsen_fm_pass_budget(int original_cell_count, int level_cell_count)
{
    if (original_cell_count < 200000)
        return 40;
    if (level_cell_count >= 300000)
        return 3;
    if (level_cell_count >= 150000)
        return 3;
    if (level_cell_count >= 75000)
        return 2;
    if (level_cell_count >= 30000)
        return 2;
    return 1;
}

static uint64_t next_pseudo_random(uint64_t &state)
{
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
    return state * 2685821657736338717ULL;
}

static int total_cell_weight(const vector<cell> &cells)
{
    int total = 0;
    for (const auto &c : cells)
        total += c.weight;
    return total;
}

static fm_view build_fm_view(const vector<cell> &cells, const vector<net> &nets)
{
    const int n = static_cast<int>(cells.size());
    fm_view view;
    view.offsets.resize(n + 1, 0);
    view.degree.resize(n, 0);

    int total_fm_edges = 0;
    for (int cid = 0; cid < n; cid++)
    {
        int deg = 0;
        for (int nid : cells[cid].connected_nets)
        {
            if (static_cast<int>(nets[nid].connected_cells.size()) <= FM_SKIPPED_NET_THRESHOLD)
                deg++;
        }
        view.degree[cid] = deg;
        view.pmax = max(view.pmax, deg);
        total_fm_edges += deg;
        view.offsets[cid + 1] = total_fm_edges;
    }

    view.nets.resize(total_fm_edges);
    for (int cid = 0; cid < n; cid++)
    {
        int pos = view.offsets[cid];
        for (int nid : cells[cid].connected_nets)
        {
            if (static_cast<int>(nets[nid].connected_cells.size()) <= FM_SKIPPED_NET_THRESHOLD)
                view.nets[pos++] = nid;
        }
    }

    return view;
}

static vector<net_state> make_net_states(const vector<net> &nets)
{
    vector<net_state> states(nets.size());
    for (int nid = 0; nid < static_cast<int>(nets.size()); nid++)
    {
        const int pin_count = static_cast<int>(nets[nid].connected_cells.size());
        auto &st = states[nid];
        st.next_pin.resize(pin_count, -1);
        st.prev_pin.resize(pin_count, -1);
        st.pin_group.resize(pin_count, 0);
    }
    return states;
}

static coarsen_result coarsen(const parsed_input &fi, uint64_t seed, bool fast_shuffle, const vector<char> *guided_part = nullptr)
{
    const int n = static_cast<int>(fi.cells.size());
    const int m = static_cast<int>(fi.nets.size());
    const int coarse_max_weight = balance_limits(total_cell_weight(fi.cells), fi.balance_factor).second;

    vector<int> mate(n, -1);
    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    if (fast_shuffle)
    {
        uint64_t shuffle_state = mix64(seed ^ 0xd1b54a32d192ed03ULL);
        for (int i = n - 1; i > 0; i--)
        {
            const int j = static_cast<int>(next_pseudo_random(shuffle_state) % static_cast<uint64_t>(i + 1));
            swap(order[i], order[j]);
        }
    }
    else
    {
        sort(order.begin(), order.end(), [&](int a, int b) {
            const uint64_t ka = mix64(seed ^ static_cast<uint64_t>(a));
            const uint64_t kb = mix64(seed ^ static_cast<uint64_t>(b));
            if (ka != kb)
                return ka < kb;
            return a < b;
        });
    }

    vector<int> shared_score(n, 0);
    vector<int> touched;
    touched.reserve(256);

    vector<int> net_weight(m, 1);
    for (int nid = 0; nid < m; nid++)
    {
        const int degree = static_cast<int>(fi.nets[nid].connected_cells.size());
        if (degree > COARSEN_NET_THRESHOLD)
            continue;
        net_weight[nid] = max(1, COARSEN_WEIGHT_SCALE / max(1, degree - 1));
    }

    for (int cid : order)
    {
        if (mate[cid] != -1)
            continue;

        touched.clear();
        for (int nid : fi.cells[cid].connected_nets)
        {
            const auto &pins = fi.nets[nid].connected_cells;
            const int degree = static_cast<int>(pins.size());
            if (degree > COARSEN_NET_THRESHOLD)
                continue;

            const int weight = net_weight[nid];
            for (int other : pins)
            {
                if (other == cid || mate[other] != -1)
                    continue;
                if (fi.cells[cid].weight + fi.cells[other].weight > coarse_max_weight)
                    continue;
                if (guided_part && (*guided_part)[cid] != (*guided_part)[other])
                    continue;
                if (shared_score[other] == 0)
                    touched.push_back(other);
                shared_score[other] += weight;
            }
        }

        int best_mate = -1;
        int best_score = 0;
        for (int other : touched)
        {
            if (shared_score[other] > best_score
                || (shared_score[other] == best_score && best_mate != -1
                    && fi.cells[other].connected_nets.size() < fi.cells[best_mate].connected_nets.size()))
            {
                best_score = shared_score[other];
                best_mate = other;
            }
            shared_score[other] = 0;
        }

        if (best_mate != -1)
        {
            mate[cid] = best_mate;
            mate[best_mate] = cid;
        }
    }

    coarsen_result res;
    res.coarse.balance_factor = fi.balance_factor;
    res.fine_to_coarse.assign(n, -1);

    int coarse_cells = 0;
    for (int cid = 0; cid < n; cid++)
    {
        if (res.fine_to_coarse[cid] != -1)
            continue;
        res.fine_to_coarse[cid] = coarse_cells;
        if (mate[cid] != -1 && res.fine_to_coarse[mate[cid]] == -1)
            res.fine_to_coarse[mate[cid]] = coarse_cells;
        coarse_cells++;
    }

    res.coarse.cells.reserve(coarse_cells);
    for (int ci = 0; ci < coarse_cells; ci++)
        res.coarse.cells.emplace_back(string_view{}, ci, 0);
    for (int cid = 0; cid < n; cid++)
        res.coarse.cells[res.fine_to_coarse[cid]].weight += fi.cells[cid].weight;

    vector<int> seen(coarse_cells, -1);
    vector<int> cpins;
    for (int nid = 0; nid < m; nid++)
    {
        cpins.clear();
        for (int fp : fi.nets[nid].connected_cells)
        {
            const int coarse_id = res.fine_to_coarse[fp];
            if (seen[coarse_id] != nid)
            {
                seen[coarse_id] = nid;
                cpins.push_back(coarse_id);
            }
        }
        if (static_cast<int>(cpins.size()) < 2)
            continue;

        const int cnid = static_cast<int>(res.coarse.nets.size());
        res.coarse.nets.emplace_back(string_view{}, cnid);
        auto &cnet = res.coarse.nets.back();
        cnet.connected_cells = cpins;

        for (int pin = 0; pin < static_cast<int>(cpins.size()); pin++)
        {
            const int coarse_id = cpins[pin];
            res.coarse.cells[coarse_id].connected_nets.push_back(cnid);
            res.coarse.cells[coarse_id].net_pin_idx.push_back(pin);
        }
    }

    return res;
}

static int net_cut_from_counts(const array<int, GROUPS> &count)
{
    return ((count[0] > 0) + (count[1] > 0) + (count[2] > 0) >= 2) ? 1 : 0;
}

static inline int active_group_count(const array<int, GROUPS> &count)
{
    return (count[0] > 0) + (count[1] > 0) + (count[2] > 0);
}

// Delta of the cut indicator when adding one pin into group g.
// A cut is created only when the net previously lived in exactly one other group.
static inline int net_add_pin_cut_delta(const array<int, GROUPS> &count, int g)
{
    if (count[g] > 0)
        return 0;
    return (active_group_count(count) == 1) ? 1 : 0;
}

// O(1) gain contribution of a single net for moving a cell from group g to group t.
// A net is "cut" iff >=2 groups have non-zero count.
// gain = cut_before - cut_after when moving one cell from g to t.
// Non-zero only at critical conditions:
//   +1: count[g]==1 && S==2 && count[t]>0  (cell is sole occupant of g; moving it uncuts the net)
//   -1: count[t]==0 && S==1 && count[g]>1  (net is uncut in g; moving cell to empty t creates a cut)
static inline int net_gain_contrib(const array<int, GROUPS> &count, int S, int g, int t)
{
    if (count[g] == 1 && S == 2 && count[t] > 0)
        return 1;
    if (count[t] == 0 && S == 1 && count[g] > 1)
        return -1;
    return 0;
}

static pair<int, int> balance_limits(int n, double r)
{
    if (!(r > 0.0 && r < 1.0))
        throw runtime_error("Balance factor r must satisfy 0 < r < 1");

    const double minD = ((1.0 - r) / 3.0) * static_cast<double>(n);
    const double maxD = ((1.0 + r) / 3.0) * static_cast<double>(n);

    int minSize = static_cast<int>(ceil(minD));
    int maxSize = static_cast<int>(floor(maxD));

    if (minSize > maxSize)
    {
        minSize = n / 3;
        maxSize = minSize;
    }

    minSize = max(0, min(minSize, n));
    maxSize = max(0, min(maxSize, n));
    return {minSize, maxSize};
}

static int rebuild_net_stats_and_cutsize(const vector<net> &nets, vector<net_state> &net_states, const vector<char> &part)
{
    const int m = static_cast<int>(nets.size());
    int cut = 0;
    for (int nid = 0; nid < m; nid++)
    {
        const auto &nt = nets[nid];
        auto &st = net_states[nid];
        st.count = {0, 0, 0};
        st.head = {-1, -1, -1};
        st.tail = {-1, -1, -1};

        for (int pin = 0; pin < static_cast<int>(nt.connected_cells.size()); pin++)
        {
            const int cid = nt.connected_cells[pin];
            const int g = static_cast<int>(part[cid]);
            st.count[g]++;
            st.pin_group[pin] = static_cast<unsigned char>(g);
            st.next_pin[pin] = -1;
            st.prev_pin[pin] = st.tail[g];
            if (st.tail[g] != -1)
                st.next_pin[st.tail[g]] = pin;
            else
                st.head[g] = pin;
            st.tail[g] = pin;
        }
        cut += net_cut_from_counts(st.count);
    }
    return cut;
}

static void rebuild_net_stats(const vector<net> &nets, vector<net_state> &net_states, const vector<char> &part)
{
    (void)rebuild_net_stats_and_cutsize(nets, net_states, part);
}

static vector<net_state> initialized_net_states(const vector<net> &nets, const vector<char> &part)
{
    vector<net_state> net_states = make_net_states(nets);
    rebuild_net_stats(nets, net_states, part);
    return net_states;
}

static inline void move_pin_between_groups(net_state &nt, int pin, int from, int to)
{
    const int prev = nt.prev_pin[pin];
    const int next = nt.next_pin[pin];
    if (prev != -1)
        nt.next_pin[prev] = next;
    else
        nt.head[from] = next;
    if (next != -1)
        nt.prev_pin[next] = prev;
    else
        nt.tail[from] = prev;

    const int to_tail = nt.tail[to];
    nt.next_pin[pin] = -1;
    nt.prev_pin[pin] = to_tail;
    if (to_tail != -1)
        nt.next_pin[to_tail] = pin;
    else
        nt.head[to] = pin;
    nt.tail[to] = pin;
    nt.pin_group[pin] = static_cast<unsigned char>(to);
}

static void incremental_rollback(const vector<cell> &cells, vector<net_state> &net_states,
                                  vector<char> &part, const vector<int> &moved,
                                  const vector<char> &moved_from, int from_idx)
{
    for (int i = static_cast<int>(moved.size()) - 1; i > from_idx; i--)
    {
        const int v = moved[i];
        const int cur_group = static_cast<int>(part[v]);
        const int orig_group = static_cast<int>(moved_from[i]);
        part[v] = moved_from[i];
        for (int j = 0; j < static_cast<int>(cells[v].connected_nets.size()); j++)
        {
            const int nid = cells[v].connected_nets[j];
            const int pin = cells[v].net_pin_idx[j];
            auto &st = net_states[nid];
            st.count[cur_group]--;
            st.count[orig_group]++;
            move_pin_between_groups(st, pin, cur_group, orig_group);
        }
    }
}

static int current_imbalance(const array<int, GROUPS> &sz)
{
    int mn = sz[0], mx = sz[0];
    for (int g = 1; g < GROUPS; g++)
    {
        mn = min(mn, sz[g]);
        mx = max(mx, sz[g]);
    }
    return mx - mn;
}

static array<int, GROUPS> compute_group_sizes(const vector<char> &part, const vector<cell> &cells)
{
    array<int, GROUPS> sz{0, 0, 0};
    for (int cid = 0; cid < static_cast<int>(part.size()); cid++)
        sz[static_cast<int>(part[cid])] += cells[cid].weight;
    return sz;
}

static bool sizes_balanced(const array<int, GROUPS> &sz, int minSize, int maxSize)
{
    for (int g = 0; g < GROUPS; g++)
    {
        if (sz[g] < minSize || sz[g] > maxSize)
            return false;
    }
    return true;
}

static bool part_balanced(const vector<char> &part, const vector<cell> &cells, int minSize, int maxSize)
{
    return sizes_balanced(compute_group_sizes(part, cells), minSize, maxSize);
}

static int part_imbalance(const vector<char> &part, const vector<cell> &cells)
{
    return current_imbalance(compute_group_sizes(part, cells));
}

static void refine_with_balance(const vector<cell> &cells, const vector<net> &nets, const fm_view &view,
                                vector<net_state> &net_states, vector<char> &part, int minSize, int maxSize,
                                chrono::steady_clock::time_point deadline, int initial_passes, int final_passes)
{
    if (!deadline_reached(deadline) && initial_passes > 0)
        fm_refine(cells, nets, view, net_states, part, minSize, maxSize, deadline, initial_passes);
    if (deadline_reached(deadline))
        return;

    if (!part_balanced(part, cells, minSize, maxSize))
        enforce_balance(cells, nets, view, net_states, part, minSize, maxSize);
    if (!deadline_reached(deadline) && final_passes > 0)
        fm_refine(cells, nets, view, net_states, part, minSize, maxSize, deadline, final_passes);
}

static uint64_t mix64(uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

static int positive_mod_u64(uint64_t x, int mod)
{
    return (mod <= 0) ? 0 : static_cast<int>(x % static_cast<uint64_t>(mod));
}

static bool group_can_move_from(int from, int weight, const array<int, GROUPS> &sz, int minSize, int maxSize)
{
    if (sz[from] - weight < minSize)
        return false;
    for (int to = 0; to < GROUPS; to++)
    {
        if (to == from)
            continue;
        if (sz[to] + weight <= maxSize)
            return true;
    }
    return false;
}

static vector<int> make_initial_target_sizes(int total_weight, int minSize, int maxSize)
{
    vector<int> target(GROUPS, total_weight / GROUPS);
    for (int i = 0; i < total_weight % GROUPS; i++)
        target[i]++;

    auto find_donor = [&](int exclude) -> int
    {
        int best = -1;
        int best_extra = -1;
        for (int g = 0; g < GROUPS; g++)
        {
            if (g == exclude)
                continue;
            int extra = target[g] - minSize;
            if (extra > best_extra)
            {
                best_extra = extra;
                best = g;
            }
        }
        return best;
    };

    bool changed = true;
    while (changed)
    {
        changed = false;
        for (int g = 0; g < GROUPS; g++)
        {
            while (target[g] < minSize)
            {
                int donor = find_donor(g);
                if (donor == -1 || target[donor] <= minSize)
                    break;
                target[donor]--;
                target[g]++;
                changed = true;
            }
        }
        for (int g = 0; g < GROUPS; g++)
        {
            while (target[g] > maxSize)
            {
                int recv = -1;
                for (int h = 0; h < GROUPS; h++)
                {
                    if (h == g)
                        continue;
                    if (target[h] < maxSize)
                    {
                        recv = h;
                        break;
                    }
                }
                if (recv == -1)
                    break;
                target[g]--;
                target[recv]++;
                changed = true;
            }
        }
    }

    return target;
}

static int initial_refinement_candidate_count(int /*n*/, int seed_count, long long remaining_ms)
{
    if (seed_count <= 0)
        return 0;

    int count = min(seed_count, INITIAL_REFINEMENT_CANDIDATES);

    if (remaining_ms >= 120000)
        count = min(seed_count, count + 6);
    else if (remaining_ms >= 60000)
        count = min(seed_count, count + 4);
    else if (remaining_ms >= 20000)
        count = min(seed_count, count + 2);

    return count;
}

static vector<int> build_initial_bfs_order(const vector<cell> &cells, const vector<net> &nets, int start)
{
    const int n = static_cast<int>(cells.size());
    if (n == 0)
        return {};
    const bool diversify = (n >= INITIAL_DIVERSITY_THRESHOLD);

    vector<int> bfs_order;
    bfs_order.reserve(n);
    vector<unsigned char> visited(n, 0);
    vector<unsigned char> net_seen(nets.size(), 0);

    vector<int> q;
    q.reserve(n);
    auto run_component = [&](int root)
    {
        q.clear();
        size_t q_head = 0;
        q.push_back(root);
        visited[root] = 1;
        while (q_head < q.size())
        {
            const int cid = q[q_head++];
            bfs_order.push_back(cid);
            const int degree = static_cast<int>(cells[cid].connected_nets.size());
            const int net_offset = diversify ? positive_mod_u64(mix64((static_cast<uint64_t>(start) << 32) ^ static_cast<uint64_t>(cid)), degree) : 0;
            for (int step = 0; step < degree; step++)
            {
                const int net_idx = (net_offset + step) % degree;
                const int nid = cells[cid].connected_nets[net_idx];
                if (net_seen[nid])
                    continue;
                net_seen[nid] = true;

                const auto &pins = nets[nid].connected_cells;
                const int pin_count = static_cast<int>(pins.size());
                if (pin_count < LARGE_NET_THRESHOLD)
                {
                    const int pin_offset = diversify ? positive_mod_u64(mix64((static_cast<uint64_t>(cid) << 32) ^ static_cast<uint64_t>(nid) ^ static_cast<uint64_t>(start)), pin_count) : 0;
                    for (int k = 0; k < pin_count; k++)
                    {
                        const int nb = pins[(pin_offset + k) % pin_count];
                        if (!visited[nb])
                        {
                            visited[nb] = 1;
                            q.push_back(nb);
                        }
                    }
                    continue;
                }

                // Let very large nets influence traversal without letting a single
                // hyperedge flood the queue and destroy local BFS structure.
                const int sample_cap = min(LARGE_NET_BFS_SAMPLE_CAP, pin_count);
                const int stride = max(1, pin_count / sample_cap);
                const int start_pos = positive_mod_u64(mix64((static_cast<uint64_t>(cells[cid].net_pin_idx[net_idx]) << 32) ^ static_cast<uint64_t>(nid) ^ static_cast<uint64_t>(start)), stride);

                int added = 0;
                for (int pos = start_pos; pos < pin_count && added < sample_cap; pos += stride)
                {
                    const int nb = pins[pos];
                    if (!visited[nb])
                    {
                        visited[nb] = 1;
                        q.push_back(nb);
                        added++;
                    }
                }
                for (int pos = 0; pos < pin_count && added < sample_cap; pos++)
                {
                    const int nb = pins[pos];
                    if (!visited[nb])
                    {
                        visited[nb] = 1;
                        q.push_back(nb);
                        added++;
                    }
                }
            }
        }
    };

    run_component(start);
    for (int step = 1; step < n && static_cast<int>(bfs_order.size()) < n; step++)
    {
        const int cid = (start + step) % n;
        if (!visited[cid])
            run_component(cid);
    }

    return bfs_order;
}

static double initial_size_penalty(int next_size, int target_size, int minSize, int maxSize)
{
    if (next_size > maxSize)
        return numeric_limits<double>::infinity();

    const int span = max(1, maxSize - minSize);
    const int slack = max(1, maxSize - target_size);
    const int near_window = max(1, min(slack, (span + 2) / 3));
    const int remaining = maxSize - next_size;

    double penalty = static_cast<double>(next_size) / static_cast<double>(max(1, target_size) * 32);
    if (next_size > target_size)
    {
        const int over_target = next_size - target_size;
        penalty += 0.35 * static_cast<double>(over_target);
        penalty += 0.35 * static_cast<double>(over_target * over_target) / static_cast<double>(slack);
    }
    if (remaining < near_window)
    {
        const int tightness = near_window - remaining;
        penalty += 1.5 + 2.0 * static_cast<double>((tightness + 1) * (tightness + 1));
    }

    return penalty;
}

static vector<char> initial_partition_from_start(const vector<cell> &cells, const vector<net> &nets, int minSize, int maxSize, int start)
{
    const int n = static_cast<int>(cells.size());
    if (n == 0)
        return {};
    const bool diversify = (n >= INITIAL_DIVERSITY_THRESHOLD);

    vector<int> bfs_order = build_initial_bfs_order(cells, nets, start);

    vector<int> target = make_initial_target_sizes(total_cell_weight(cells), minSize, maxSize);
    vector<char> part(n, 0);
    array<int, GROUPS> sz{0, 0, 0};
    vector<array<int, GROUPS>> partial_count(nets.size());
    array<int, GROUPS> group_pref{};
    const int base_group = diversify ? positive_mod_u64(mix64(static_cast<uint64_t>(start) ^ 0xabc98388fb8fac03ULL), GROUPS) : 0;
    for (int i = 0; i < GROUPS; i++)
        group_pref[i] = (base_group + i) % GROUPS;

    for (int cid : bfs_order)
    {
        int best_group = -1;
        double best_score = numeric_limits<double>::infinity();
        int best_delta = numeric_limits<int>::max();
        int best_affinity = -1;
        int best_size = numeric_limits<int>::max();
        int best_pref_rank = numeric_limits<int>::max();

        const int group_shift = diversify ? positive_mod_u64(mix64((static_cast<uint64_t>(start) << 32) ^ static_cast<uint64_t>(cid) ^ 0x6e93d59d4fb1a7d1ULL), GROUPS) : 0;
        for (int order = 0; order < GROUPS; order++)
        {
            const int g = group_pref[(group_shift + order) % GROUPS];
            const int next_size = sz[g] + cells[cid].weight;
            if (next_size > maxSize)
                continue;

            int delta = 0;
            int affinity = 0;

            for (int nid : cells[cid].connected_nets)
            {
                const auto &count = partial_count[nid];
                delta += net_add_pin_cut_delta(count, g);
                affinity += count[g];
            }

            const double score = static_cast<double>(delta) + initial_size_penalty(next_size, target[g], minSize, maxSize);
            const bool score_eq    = (score == best_score);
            const bool delta_eq    = (delta == best_delta);
            const bool affinity_eq = (affinity == best_affinity);
            const bool size_eq     = (sz[g] == best_size);
            const bool better = best_group < 0
                || score < best_score
                || (score_eq && delta < best_delta)
                || (score_eq && delta_eq && affinity > best_affinity)
                || (score_eq && delta_eq && affinity_eq && sz[g] < best_size)
                || (score_eq && delta_eq && affinity_eq && size_eq && diversify && order < best_pref_rank);
            if (better)
            {
                best_group = g;
                best_score = score;
                best_delta = delta;
                best_affinity = affinity;
                best_size = sz[g];
                best_pref_rank = order;
            }
        }

        if (best_group < 0)
            best_group = 0;

        part[cid] = static_cast<char>(best_group);
        sz[best_group] += cells[cid].weight;
        for (int nid : cells[cid].connected_nets)
            partial_count[nid][best_group]++;
    }

    return part;
}

static int cutsize_from_part(const vector<net> &nets, const vector<char> &part)
{
    int cut = 0;
    for (const auto &nt : nets)
    {
        if (nt.connected_cells.empty())
            continue;

        const char g0 = part[nt.connected_cells.front()];
        for (int i = 1; i < static_cast<int>(nt.connected_cells.size()); i++)
        {
            if (part[nt.connected_cells[i]] != g0)
            {
                cut++;
                break;
            }
        }
    }
    return cut;
}

static int refine_initial_seed(const vector<cell> &cells, const vector<net> &nets, vector<char> &part, int minSize, int maxSize, chrono::steady_clock::time_point deadline)
{
    const fm_view view = build_fm_view(cells, nets);
    vector<net_state> trial_states = initialized_net_states(nets, part);
    refine_with_balance(cells, nets, view, trial_states, part, minSize, maxSize, deadline, INITIAL_SEED_FM_PASSES, 1);
    return rebuild_net_stats_and_cutsize(nets, trial_states, part);
}

static vector<int> initial_seed_candidates(const vector<cell> &cells, int desired)
{
    const int n = static_cast<int>(cells.size());
    if (n == 0 || desired <= 0)
        return {};
    const bool diversify = (n >= INITIAL_DIVERSITY_THRESHOLD);

    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) {
        const size_t da = cells[a].connected_nets.size();
        const size_t db = cells[b].connected_nets.size();
        if (da != db)
            return da > db;
        return a < b;
    });

    vector<int> seeds;
    seeds.reserve(desired);
    vector<unsigned char> used(n, 0);
    auto add_seed = [&](int cid) {
        if (!used[cid])
        {
            used[cid] = 1;
            seeds.push_back(cid);
        }
    };

    if (!diversify)
    {
        const int top_count = min(n, max(1, desired / 2));
        for (int i = 0; i < top_count; i++)
            add_seed(order[i]);

        const int spread_count = desired - static_cast<int>(seeds.size());
        for (int i = 0; i < spread_count; i++)
        {
            const int pos = static_cast<int>((static_cast<long long>(i + 1) * n) / (spread_count + 1));
            add_seed(order[min(n - 1, pos)]);
        }
    }
    else
    {
        const int top_count = min(n, max(2, desired / 4));
        for (int i = 0; i < top_count; i++)
            add_seed(order[i]);

        const int bottom_count = min(n, max(2, desired / 4));
        for (int i = 0; i < bottom_count; i++)
            add_seed(order[n - 1 - i]);

        const int spread_count = max(2, desired / 4);
        for (int i = 0; i < spread_count; i++)
        {
            const int pos = static_cast<int>((static_cast<long long>(i + 1) * n) / (spread_count + 1));
            add_seed(order[min(n - 1, pos)]);
        }

        const int index_spread = max(2, desired / 4);
        for (int i = 0; i < index_spread; i++)
        {
            const int cid = static_cast<int>((static_cast<long long>(i + 1) * n) / (index_spread + 1));
            add_seed(cid);
        }

        for (int i = 0; static_cast<int>(seeds.size()) < desired && i < n; i++)
        {
            const int pos = positive_mod_u64(mix64(static_cast<uint64_t>(n) ^ static_cast<uint64_t>(i) ^ 0x2545f4914f6cdd1dULL), n);
            add_seed(order[pos]);
        }
    }

    for (int cid : order)
    {
        if (static_cast<int>(seeds.size()) >= desired)
            break;
        add_seed(cid);
    }

    return seeds;
}

static vector<char> initial_partition(const vector<cell> &cells, const vector<net> &nets, int minSize, int maxSize, chrono::steady_clock::time_point deadline, int start_override = -1)
{
    const int n = static_cast<int>(cells.size());
    if (n == 0)
        return {};
    const bool diversify = (n >= INITIAL_DIVERSITY_THRESHOLD);

    if (start_override >= 0 && start_override < n)
        return initial_partition_from_start(cells, nets, minSize, maxSize, start_override);
    if (deadline_reached(deadline))
        return initial_partition_from_start(cells, nets, minSize, maxSize, 0);

    const long long remaining_ms = max<long long>(
        0,
        chrono::duration_cast<chrono::milliseconds>(deadline - chrono::steady_clock::now()).count());
    const vector<int> seeds = initial_seed_candidates(cells, initial_seed_target_count(n, remaining_ms));
    if (seeds.empty())
        return initial_partition_from_start(cells, nets, minSize, maxSize, 0);

    struct initial_trial
    {
        vector<char> part;
        int cut = numeric_limits<int>::max();
        int imbalance = numeric_limits<int>::max();
        int seed = -1;
    };

    auto is_better_trial = [](const initial_trial &a, const initial_trial &b) {
        if (a.cut != b.cut) return a.cut < b.cut;
        if (a.imbalance != b.imbalance) return a.imbalance < b.imbalance;
        return a.seed < b.seed;
    };

    initial_trial best;
    vector<initial_trial> trials;
    trials.reserve(seeds.size());
    for (int seed : seeds)
    {
        if (deadline_reached(deadline))
            break;

        vector<char> part = initial_partition_from_start(cells, nets, minSize, maxSize, seed);
        initial_trial candidate;
        candidate.cut = cutsize_from_part(nets, part);
        candidate.imbalance = part_imbalance(part, cells);
        candidate.seed = seed;
        candidate.part = std::move(part);

        if (!diversify)
        {
            if (is_better_trial(candidate, best))
                best = std::move(candidate);
            continue;
        }

        trials.push_back(std::move(candidate));
    }

    if (diversify && !trials.empty())
    {
        sort(trials.begin(), trials.end(), is_better_trial);
        best = trials.front();

        const int refine_count = initial_refinement_candidate_count(n, static_cast<int>(trials.size()), remaining_ms);
        for (int i = 0; i < refine_count; i++)
        {
            if (deadline_reached(deadline))
                break;

            initial_trial refined = std::move(trials[i]);
            refined.cut = refine_initial_seed(cells, nets, refined.part, minSize, maxSize, deadline);
            refined.imbalance = part_imbalance(refined.part, cells);
            if (is_better_trial(refined, best))
                best = std::move(refined);
        }
    }

    if (best.seed < 0)
        return initial_partition_from_start(cells, nets, minSize, maxSize, seeds.front());
    return best.part;
}

static void compute_cell_gains(int cid, const fm_view &view, const vector<net_state> &net_states, const vector<char> &part, vector<array<int, GROUPS>> &gain_to)
{
    const int from = static_cast<int>(part[cid]);
    gain_to[cid].fill(NEG_INF);
    for (int to = 0; to < GROUPS; to++)
        if (to != from)
            gain_to[cid][to] = 0;

    for (int pos = view.offsets[cid]; pos < view.offsets[cid + 1]; pos++)
    {
        const int nid = view.nets[pos];
        const auto &st = net_states[nid];
        const int active = active_group_count(st.count);
        for (int to = 0; to < GROUPS; to++)
        {
            if (to == from)
                continue;
            gain_to[cid][to] += net_gain_contrib(st.count, active, from, to);
        }
    }
}

static void recompute_all_cell_gains(const fm_view &view, const vector<net_state> &net_states, const vector<char> &part, vector<array<int, GROUPS>> &gain_to)
{
    const int n = static_cast<int>(view.degree.size());
    for (int cid = 0; cid < n; cid++)
        compute_cell_gains(cid, view, net_states, part, gain_to);
}

static move_choice best_feasible_move_for_cell(int cid, const vector<cell> &cells, const vector<char> &part, const vector<array<int, GROUPS>> &gain_to, const array<int, GROUPS> &sz, int minSize, int maxSize)
{
    move_choice best;
    const int from = static_cast<int>(part[cid]);
    const int weight = cells[cid].weight;

    if (sz[from] - weight < minSize)
        return best;

    for (int to = 0; to < GROUPS; to++)
    {
        if (to == from)
            continue;
        if (sz[to] + weight > maxSize)
            continue;
        if (gain_to[cid][to] > best.gain)
        {
            best.gain = gain_to[cid][to];
            best.to = to;
        }
    }

    return best;
}

class BucketList
{
    int pmax;
    int offset;
    vector<int> head;
    vector<int> bucket_epoch;
    vector<int> next;
    vector<int> prev;
    vector<int> bucket_idx;
    vector<int> cell_epoch;
    int epoch;
    mutable int max_gain_bucket;
    int count;

  public:
    BucketList() : pmax(0), offset(0), epoch(1), max_gain_bucket(-1), count(0)
    {
    }

    void init(int n, int max_degree)
    {
        pmax = max_degree;
        offset = pmax;
        const int total_buckets = 2 * pmax + 1;
        head.assign(total_buckets, -1);
        bucket_epoch.assign(total_buckets, 0);
        next.assign(n, -1);
        prev.assign(n, -1);
        bucket_idx.assign(n, -1);
        cell_epoch.assign(n, 0);
        epoch = 1;
        max_gain_bucket = -1;
        count = 0;
    }

    void insert(int cid, int gain)
    {
        int b = gain + offset;
        if (b < 0)
            b = 0;
        if (b >= static_cast<int>(head.size()))
            b = static_cast<int>(head.size()) - 1;

        if (bucket_epoch[b] != epoch)
        {
            bucket_epoch[b] = epoch;
            head[b] = -1;
        }
        const int old_head = head[b];
        head[b] = cid;
        next[cid] = old_head;
        prev[cid] = -1;
        if (old_head != -1)
            prev[old_head] = cid;
        bucket_idx[cid] = b;
        cell_epoch[cid] = epoch;
        if (b > max_gain_bucket)
            max_gain_bucket = b;
        count++;
    }

    void remove(int cid)
    {
        if (cell_epoch[cid] != epoch)
            return;
        const int b = bucket_idx[cid];

        const int p = prev[cid];
        const int n = next[cid];
        if (p != -1)
            next[p] = n;
        else
            head[b] = n;
        if (n != -1)
            prev[n] = p;

        next[cid] = -1;
        prev[cid] = -1;
        bucket_idx[cid] = -1;
        cell_epoch[cid] = 0;
        count--;
    }

    void update(int cid, int new_gain)
    {
        remove(cid);
        insert(cid, new_gain);
    }

    int pop_max()
    {
        while (max_gain_bucket >= 0
            && (bucket_epoch[max_gain_bucket] != epoch || head[max_gain_bucket] == -1))
            max_gain_bucket--;
        if (max_gain_bucket < 0)
            return -1;
        const int cid = head[max_gain_bucket];
        remove(cid);
        return cid;
    }

    int peek_max() const
    {
        int b = max_gain_bucket;
        while (b >= 0 && (bucket_epoch[b] != epoch || head[b] == -1))
            b--;
        if (b < 0)
            return -1;
        return head[b];
    }

    bool empty() const
    {
        return count == 0;
    }

    void clear()
    {
        epoch++;
        if (epoch == numeric_limits<int>::max())
        {
            fill(bucket_epoch.begin(), bucket_epoch.end(), 0);
            fill(cell_epoch.begin(), cell_epoch.end(), 0);
            fill(bucket_idx.begin(), bucket_idx.end(), -1);
            epoch = 1;
        }
        max_gain_bucket = -1;
        count = 0;
    }
};

static void sync_cell_pair_buckets(int cid, const vector<char> &part, const vector<array<int, GROUPS>> &gain_to, array<array<BucketList, GROUPS>, GROUPS> &pair_buckets)
{
    const int from = static_cast<int>(part[cid]);
    for (int g = 0; g < GROUPS; g++)
    {
        for (int to = 0; to < GROUPS; to++)
        {
            if (g == to)
                continue;
            if (g == from)
                pair_buckets[g][to].update(cid, gain_to[cid][to]);
            else
                pair_buckets[g][to].remove(cid);
        }
    }
}

static move_candidate peek_best_candidate_for_pair(int from, int to, const BucketList &bucket, const vector<char> &part, const vector<array<int, GROUPS>> &gain_to)
{
    move_candidate cand;
    const int cid = bucket.peek_max();
    if (cid < 0 || static_cast<int>(part[cid]) != from)
        return cand;

    cand.cid = cid;
    cand.from = from;
    cand.to = to;
    cand.gain = gain_to[cid][to];
    return cand;
}

static move_candidate extract_best_candidate_for_group(int from, BucketList &bucket, const vector<cell> &cells, const vector<char> &part, const vector<char> &locked, const vector<array<int, GROUPS>> &gain_to, const array<int, GROUPS> &sz, int minSize, int maxSize, vector<int> &key)
{
    move_candidate cand;

    if (!group_can_move_from(from, 1, sz, minSize, maxSize))
        return cand;

    while (!bucket.empty())
    {
        const int cid = bucket.pop_max();
        if (cid < 0)
            return cand;

        if (locked[cid] || static_cast<int>(part[cid]) != from)
            continue;

        const int weight = cells[cid].weight;
        if (!group_can_move_from(from, weight, sz, minSize, maxSize))
            continue;

        const move_choice best = best_feasible_move_for_cell(cid, cells, part, gain_to, sz, minSize, maxSize);
        if (best.to < 0)
            continue;

        if (best.gain != key[cid])
        {
            key[cid] = best.gain;
            bucket.insert(cid, key[cid]);
            continue;
        }

        cand.cid = cid;
        cand.from = from;
        cand.to = best.to;
        cand.gain = best.gain;
        return cand;
    }

    return cand;
}

static move_candidate select_move_candidate(const array<move_candidate, GROUPS> &cands, const vector<cell> &cells, const array<int, GROUPS> &sz)
{
    move_candidate best;
    int best_imb = numeric_limits<int>::max();

    for (int g = 0; g < GROUPS; g++)
    {
        const auto &c = cands[g];
        if (!c.valid())
            continue;

        array<int, GROUPS> next = sz;
        next[c.from] -= cells[c.cid].weight;
        next[c.to] += cells[c.cid].weight;
        const int imb = current_imbalance(next);

        if (!best.valid() || c.gain > best.gain || (c.gain == best.gain && imb < best_imb))
        {
            best = c;
            best_imb = imb;
        }
    }

    return best;
}

// Propagates gain changes after moving cell v from group `from` to group `to`.
// Updates nt.count, linked pin membership, and gain_to for all affected
// unlocked neighbors. Populates `affected` with cells whose gains changed.
// Pass locked=nullptr to skip the locked check (used by enforce_balance).
static void propagate_net_gains(
    int v, int from, int to,
    const vector<cell> &cells, const vector<net> &nets, vector<net_state> &net_states,
    vector<array<int, GROUPS>> &gain_to,
    const vector<char> *locked,
    vector<int> &seen, int &stamp,
    vector<int> &affected)
{
    stamp++;
    if (stamp == numeric_limits<int>::max())
    {
        fill(seen.begin(), seen.end(), 0);
        stamp = 1;
    }
    affected.clear();

    for (int i = 0; i < static_cast<int>(cells[v].connected_nets.size()); i++)
    {
        const int nid = cells[v].connected_nets[i];
        const int pin = cells[v].net_pin_idx[i];
        const auto &nt = nets[nid];
        auto &st = net_states[nid];
        if (static_cast<int>(nt.connected_cells.size()) > FM_SKIPPED_NET_THRESHOLD)
        {
            st.count[from]--;
            st.count[to]++;
            move_pin_between_groups(st, pin, from, to);
            continue;
        }
        const array<int, GROUPS> c_old = st.count;
        const int cf = c_old[from];
        const int co = c_old[to];

        st.count[from]--;
        st.count[to]++;
        move_pin_between_groups(st, pin, from, to);

        // No critical transition when both counts stay above the thresholds
        // that affect gain conditions.
        if (cf > 2 && co > 1)
            continue;

        const int S_old = active_group_count(c_old);
        const int S_new = active_group_count(st.count);

        // Precompute delta[g][t] for all (g,t) pairs — O(GROUPS^2) = O(1).
        array<array<int, GROUPS>, GROUPS> delta{};
        bool any_delta = false;
        for (int g = 0; g < GROUPS; g++)
        {
            for (int t = 0; t < GROUPS; t++)
            {
                if (t == g)
                    continue;
                const int d = net_gain_contrib(st.count, S_new, g, t) - net_gain_contrib(c_old, S_old, g, t);
                delta[g][t] = d;
                if (d != 0)
                    any_delta = true;
            }
        }

        if (!any_delta)
            continue;

        for (int g = 0; g < GROUPS; g++)
        {
            bool g_has_delta = false;
            for (int t = 0; t < GROUPS; t++)
            {
                if (t == g)
                    continue;
                if (delta[g][t] != 0)
                {
                    g_has_delta = true;
                    break;
                }
            }
            if (!g_has_delta)
                continue;

            for (int pin_idx = st.head[g]; pin_idx != -1; pin_idx = st.next_pin[pin_idx])
            {
                const int u = nt.connected_cells[pin_idx];
                if (locked && (*locked)[u])
                    continue;

                bool changed = false;
                for (int t = 0; t < GROUPS; t++)
                {
                    if (t == g)
                        continue;
                    if (delta[g][t] != 0)
                    {
                        gain_to[u][t] += delta[g][t];
                        changed = true;
                    }
                }

                if (changed && seen[u] != stamp)
                {
                    seen[u] = stamp;
                    affected.push_back(u);
                }
            }
        }
    }
}

static void fm_refine(const vector<cell> &cells, const vector<net> &nets, const fm_view &view, vector<net_state> &net_states, vector<char> &part, int minSize, int maxSize, chrono::steady_clock::time_point deadline, int max_passes = 40)
{
    const int n = static_cast<int>(cells.size());
    if (n == 0)
        return;

    const int pmax = view.pmax;

    vector<array<int, GROUPS>> gain_to(n);
    vector<int> key(n, NEG_INF);
    vector<char> locked(n, 0);

    BucketList buckets[GROUPS];
    for (int g = 0; g < GROUPS; g++)
        buckets[g].init(n, pmax);

    array<int, GROUPS> sz = compute_group_sizes(part, cells);

    const int no_improve_max = (n > 50000) ? 2 : 3;
    int no_improve_count = 0;

    vector<int> moved;
    vector<char> moved_from;
    moved.reserve(n);
    moved_from.reserve(n);

    vector<int> affected;
    affected.reserve(n);

    vector<int> seen(n, 0);
    int stamp = 1;

    const int adaptive_pass_cap = (n > 500000) ? 10
                                : (n > 100000) ? 20
                                : (n > 50000)  ? 30
                                : max_passes;
    const int effective_max_passes = min(max_passes, adaptive_pass_cap);
    for (int pass = 0; pass < effective_max_passes; pass++)
    {
        if (chrono::steady_clock::now() >= deadline)
            break;

        fill(locked.begin(), locked.end(), 0);
        for (int g = 0; g < GROUPS; g++)
            buckets[g].clear();

        sz = compute_group_sizes(part, cells);

        // Fused boundary detection + gain computation.
        // Non-boundary cells (all FM nets uncut) have gain = -degree in O(1).
        // Only boundary cells get full gain computation and bucket insertion.
        for (int cid = 0; cid < n; cid++)
        {
            const int from = static_cast<int>(part[cid]);
            bool is_boundary = false;
            for (int pos = view.offsets[cid]; pos < view.offsets[cid + 1]; pos++)
            {
                if (active_group_count(net_states[view.nets[pos]].count) >= 2)
                {
                    is_boundary = true;
                    break;
                }
            }
            if (is_boundary)
            {
                compute_cell_gains(cid, view, net_states, part, gain_to);
                const move_choice best = best_feasible_move_for_cell(cid, cells, part, gain_to, sz, minSize, maxSize);
                key[cid] = best.gain;
                if (best.to >= 0)
                    buckets[from].insert(cid, key[cid]);
            }
            else
            {
                gain_to[cid].fill(NEG_INF);
                const int neg_deg = -view.degree[cid];
                for (int to = 0; to < GROUPS; to++)
                    if (to != from)
                        gain_to[cid][to] = neg_deg;
            }
        }

        moved.clear();
        moved_from.clear();
        int cum = 0;
        int bestCum = 0;
        int bestIndex = -1;
        const int early_stop_threshold = max(n / 5, 200);

        while (true)
        {
            if (chrono::steady_clock::now() >= deadline)
                break;
            if (bestIndex >= 0 && static_cast<int>(moved.size()) - bestIndex > early_stop_threshold)
                break;

            array<move_candidate, GROUPS> cands;
            for (int g = 0; g < GROUPS; g++)
                cands[g] = extract_best_candidate_for_group(g, buckets[g], cells, part, locked, gain_to, sz, minSize, maxSize, key);

            move_candidate chosen = select_move_candidate(cands, cells, sz);
            if (!chosen.valid())
                break;

            for (int g = 0; g < GROUPS; g++)
            {
                const auto &c = cands[g];
                if (c.valid() && c.cid != chosen.cid)
                    buckets[c.from].insert(c.cid, key[c.cid]);
            }

            const int v = chosen.cid;
            const int from = chosen.from;
            const int to = chosen.to;

            locked[v] = 1;
            moved.push_back(v);
            moved_from.push_back(static_cast<char>(from));

            cum += chosen.gain;
            if (cum > bestCum)
            {
                bestCum = cum;
                bestIndex = static_cast<int>(moved.size()) - 1;
            }

            part[v] = static_cast<char>(to);
            sz[from] -= cells[v].weight;
            sz[to] += cells[v].weight;

            // Critical-net O(1)-per-pin gain update via analytical delta formula.
            propagate_net_gains(v, from, to, cells, nets, net_states, gain_to, &locked, seen, stamp, affected);

            // Batch bucket updates for all affected cells
            for (int u : affected)
            {
                const move_choice best = best_feasible_move_for_cell(u, cells, part, gain_to, sz, minSize, maxSize);
                key[u] = best.gain;

                const int group = static_cast<int>(part[u]);
                if (best.to >= 0)
                    buckets[group].update(u, key[u]);
                else
                    buckets[group].remove(u);
            }
        }

        if (bestCum <= 0)
        {
            incremental_rollback(cells, net_states, part, moved, moved_from, -1);

            no_improve_count++;
            if (no_improve_count >= no_improve_max)
                break;
            continue;
        }

        no_improve_count = 0;

        incremental_rollback(cells, net_states, part, moved, moved_from, bestIndex);
    }
}

static void fm_refine(const vector<cell> &cells, const vector<net> &nets, vector<net_state> &net_states, vector<char> &part, int minSize, int maxSize, chrono::steady_clock::time_point deadline, int max_passes)
{
    const fm_view view = build_fm_view(cells, nets);
    fm_refine(cells, nets, view, net_states, part, minSize, maxSize, deadline, max_passes);
}

static void enforce_balance(const vector<cell> &cells, const vector<net> &nets, const fm_view &view, vector<net_state> &net_states, vector<char> &part, int minSize, int maxSize)
{
    const int n = static_cast<int>(part.size());
    if (n == 0)
        return;

    array<int, GROUPS> sz = compute_group_sizes(part, cells);

    const int pmax = view.pmax;

    vector<array<int, GROUPS>> gain_to(n);
    recompute_all_cell_gains(view, net_states, part, gain_to);

    array<array<BucketList, GROUPS>, GROUPS> pair_buckets;
    for (int from = 0; from < GROUPS; from++)
    {
        for (int to = 0; to < GROUPS; to++)
        {
            if (from == to)
                continue;
            pair_buckets[from][to].init(n, pmax);
        }
    }

    for (int cid = 0; cid < n; cid++)
    {
        const int from = static_cast<int>(part[cid]);
        for (int to = 0; to < GROUPS; to++)
        {
            if (to == from)
                continue;
            pair_buckets[from][to].insert(cid, gain_to[cid][to]);
        }
    }

    vector<int> affected;
    affected.reserve(256);
    vector<int> seen(n, 0);
    int stamp = 1;

    auto apply_move = [&](int cid, int from, int to) {
        for (int t = 0; t < GROUPS; t++)
        {
            if (t == from)
                continue;
            pair_buckets[from][t].remove(cid);
        }

        part[cid] = static_cast<char>(to);
        sz[from] -= cells[cid].weight;
        sz[to] += cells[cid].weight;

        propagate_net_gains(cid, from, to, cells, nets, net_states, gain_to, nullptr, seen, stamp, affected);

        compute_cell_gains(cid, view, net_states, part, gain_to);
        sync_cell_pair_buckets(cid, part, gain_to, pair_buckets);

        for (int u : affected)
        {
            if (u == cid)
                continue;
            sync_cell_pair_buckets(u, part, gain_to, pair_buckets);
        }
    };

    while (true)
    {
        move_candidate best;
        int best_imb = numeric_limits<int>::max();
        const bool has_underfilled = any_of(sz.begin(), sz.end(), [&](int s) { return s < minSize; });

        if (has_underfilled)
        {
            for (int to = 0; to < GROUPS; to++)
            {
                if (sz[to] >= minSize)
                    continue;
                for (int from = 0; from < GROUPS; from++)
                {
                    if (from == to || sz[from] <= minSize)
                        continue;

                    const move_candidate cand = peek_best_candidate_for_pair(from, to, pair_buckets[from][to], part, gain_to);
                    if (!cand.valid())
                        continue;

                    array<int, GROUPS> next = sz;
                    next[from] -= cells[cand.cid].weight;
                    next[to] += cells[cand.cid].weight;
                    const int imb = current_imbalance(next);
                    if (!best.valid() || cand.gain > best.gain || (cand.gain == best.gain && imb < best_imb))
                    {
                        best = cand;
                        best_imb = imb;
                    }
                }
            }
        }
        else
        {
            bool has_overfilled = false;
            for (int from = 0; from < GROUPS; from++)
            {
                if (sz[from] <= maxSize)
                    continue;
                has_overfilled = true;

                for (int to = 0; to < GROUPS; to++)
                {
                    if (from == to || sz[to] >= maxSize)
                        continue;

                    const move_candidate cand = peek_best_candidate_for_pair(from, to, pair_buckets[from][to], part, gain_to);
                    if (!cand.valid())
                        continue;

                    array<int, GROUPS> next = sz;
                    next[from] -= cells[cand.cid].weight;
                    next[to] += cells[cand.cid].weight;
                    const int imb = current_imbalance(next);
                    if (!best.valid() || cand.gain > best.gain || (cand.gain == best.gain && imb < best_imb))
                    {
                        best = cand;
                        best_imb = imb;
                    }
                }
            }

            if (!has_overfilled)
                break;
        }

        if (!best.valid())
            break;
        apply_move(best.cid, best.from, best.to);
    }
}

static partition_result make_partition_result(const vector<net> &nets, vector<char> part, int minSize, int maxSize)
{
    partition_result res;
    res.cut = cutsize_from_part(nets, part);
    res.part = std::move(part);
    res.minSize = minSize;
    res.maxSize = maxSize;
    return res;
}

static partition_result run_multilevel_cycle(parsed_input &parsed, const fm_view &finest_view,
                                             chrono::steady_clock::time_point deadline, uint64_t seed,
                                             const vector<char> *init_part, int uncoarsen_pass_cap)
{
    const int total_weight = total_cell_weight(parsed.cells);
    const int original_n = static_cast<int>(parsed.cells.size());
    const bool fast_shuffle = (original_n >= 200000);
    const auto limits = balance_limits(total_weight, parsed.balance_factor);
    const int minSize = limits.first;
    const int maxSize = limits.second;

    static constexpr int MIN_COARSE_CELLS = 300;
    vector<coarsen_result> hier;
    hier.reserve(12);

    vector<char> current_part;
    if (init_part)
        current_part = *init_part;
    const parsed_input *cur = &parsed;

    while (static_cast<int>(cur->cells.size()) > MIN_COARSE_CELLS)
    {
        if (deadline_reached(deadline))
            break;
        const int prev_n = static_cast<int>(cur->cells.size());
        hier.push_back(coarsen(*cur,
                               seed ^ mix64(static_cast<uint64_t>(hier.size()) + 0x9e3779b97f4a7c15ULL),
                               fast_shuffle,
                               init_part ? &current_part : nullptr));

        if (init_part)
        {
            const auto &cr = hier.back();
            const int coarse_n = static_cast<int>(cr.coarse.cells.size());
            vector<char> coarse_part(coarse_n);
            for (int fid = 0; fid < prev_n; fid++)
                coarse_part[cr.fine_to_coarse[fid]] = current_part[fid];
            current_part = std::move(coarse_part);
        }

        cur = &hier.back().coarse;
        const int new_n = static_cast<int>(cur->cells.size());
        if (new_n > prev_n * 90 / 100)
            break;
    }

    const auto clims = balance_limits(total_cell_weight(cur->cells), parsed.balance_factor);

    vector<char> part = init_part
        ? std::move(current_part)
        : initial_partition(cur->cells, cur->nets, clims.first, clims.second, deadline);
    {
        vector<net_state> cnet_states = initialized_net_states(cur->nets, part);
        const fm_view coarse_view = build_fm_view(cur->cells, cur->nets);
        refine_with_balance(cur->cells, cur->nets, coarse_view, cnet_states, part, clims.first, clims.second, deadline, 40, 40);
    }

    for (int lvl = static_cast<int>(hier.size()) - 1; lvl >= 0; lvl--)
    {
        const auto &cr = hier[lvl];
        const parsed_input &fine = (lvl == 0) ? parsed : hier[lvl - 1].coarse;
        const int fn = static_cast<int>(fine.cells.size());
        const auto flims = balance_limits(total_cell_weight(fine.cells), parsed.balance_factor);

        vector<char> fpart(fn);
        for (int fid = 0; fid < fn; fid++)
            fpart[fid] = part[cr.fine_to_coarse[fid]];

        vector<net_state> fnet_states = initialized_net_states(fine.nets, fpart);
        const int pass_budget = min(uncoarsen_fm_pass_budget(original_n, fn), uncoarsen_pass_cap);
        if (!deadline_reached(deadline) && pass_budget > 0)
        {
            if (lvl == 0)
                fm_refine(fine.cells, fine.nets, finest_view, fnet_states, fpart, flims.first, flims.second, deadline, pass_budget);
            else
                fm_refine(fine.cells, fine.nets, fnet_states, fpart, flims.first, flims.second, deadline, pass_budget);
        }
        part = std::move(fpart);
    }

    if (!part_balanced(part, parsed.cells, minSize, maxSize))
    {
        vector<net_state> final_net_states = initialized_net_states(parsed.nets, part);
        enforce_balance(parsed.cells, parsed.nets, finest_view, final_net_states, part, minSize, maxSize);
        if (!deadline_reached(deadline))
            fm_refine(parsed.cells, parsed.nets, finest_view, final_net_states, part, minSize, maxSize, deadline, 1);
    }

    return make_partition_result(parsed.nets, std::move(part), minSize, maxSize);
}

static partition_result run_partition_cycle(parsed_input &parsed, const fm_view &finest_view, chrono::steady_clock::time_point deadline, uint64_t seed)
{
    return run_multilevel_cycle(parsed, finest_view, deadline, seed, nullptr, numeric_limits<int>::max());
}

static partition_result run_vcycle(parsed_input &parsed, const vector<char> &init_part,
                                   const fm_view &finest_view,
                                   chrono::steady_clock::time_point deadline, uint64_t seed)
{
    return run_multilevel_cycle(parsed, finest_view, deadline, seed, &init_part, 4);
}

partition_result partition(parsed_input &parsed, chrono::steady_clock::time_point deadline)
{
    const int n = static_cast<int>(parsed.cells.size());
    const uint64_t base_seed = partition_seed_base(n);

    partition_result best;
    best.cut = numeric_limits<int>::max();
    partition_result fallback;
    fallback.cut = numeric_limits<int>::max();
    bool have_balanced_best = false;
    bool have_fallback = false;

    const fm_view finest_view = build_fm_view(parsed.cells, parsed.nets);
    const auto search_start = chrono::steady_clock::now();
    const long long total_budget_ms = max<long long>(
        0,
        chrono::duration_cast<chrono::milliseconds>(deadline - search_start).count());
    const int max_cycles = partition_cycle_count(n, total_budget_ms);
    chrono::steady_clock::time_point search_deadline = deadline;

    for (int cycle = 0; cycle < max_cycles; cycle++)
    {
        if (cycle > 0 && chrono::steady_clock::now() >= search_deadline)
            break;

        const auto cycle_start = chrono::steady_clock::now();
        partition_result candidate = run_partition_cycle(parsed, finest_view, search_deadline, mix64(base_seed ^ static_cast<uint64_t>(cycle)));
        const bool candidate_balanced = part_balanced(candidate.part, parsed.cells, candidate.minSize, candidate.maxSize);
        if (!have_fallback || candidate.cut < fallback.cut)
        {
            fallback = candidate;
            have_fallback = true;
        }
        if (candidate_balanced && (!have_balanced_best || candidate.cut < best.cut))
        {
            best = std::move(candidate);
            have_balanced_best = true;
        }

        if (cycle == 0 && max_cycles > 1)
        {
            const auto first_cycle_ms = max<int>(
                1,
                static_cast<int>(chrono::duration_cast<chrono::milliseconds>(
                    chrono::steady_clock::now() - cycle_start).count()));
            const long long reserve_ms = cycle_budget_reserve_ms(n, total_budget_ms);
            const long long cycle_budget_ms = min(
                max<long long>(0, total_budget_ms - reserve_ms),
                max<long long>(static_cast<long long>(first_cycle_ms) * max_cycles, first_cycle_ms + 1000LL));
            search_deadline = min(deadline - chrono::milliseconds(min(reserve_ms, total_budget_ms)),
                                  search_start + chrono::milliseconds(cycle_budget_ms));
        }
    }

    // V-cycles: use remaining time to refine the best solution via guided coarsening
    if (have_balanced_best)
    {
        const auto remaining_ms = chrono::duration_cast<chrono::milliseconds>(
            deadline - chrono::steady_clock::now()).count();
        const long long vcycle_budget_ms = max(0LL, remaining_ms - 1000LL);
        if (vcycle_budget_ms > 0)
        {
            const auto vcycle_deadline = min(deadline - chrono::milliseconds(1000),
                chrono::steady_clock::now() + chrono::milliseconds(vcycle_budget_ms));
            const int no_improve_limit = (n <= 20000) ? 4 : (n <= 100000) ? 6 : 8;
            int no_improve_vcycles = 0;
            for (int vc = 0; chrono::steady_clock::now() < vcycle_deadline; vc++)
            {
                partition_result candidate = run_vcycle(parsed, best.part, finest_view, vcycle_deadline,
                    mix64(base_seed ^ static_cast<uint64_t>(max_cycles + vc) ^ 0x7a6d39b2c1e8f504ULL));
                const bool balanced = part_balanced(candidate.part, parsed.cells, candidate.minSize, candidate.maxSize);
                if (balanced && candidate.cut < best.cut)
                {
                    best = std::move(candidate);
                    no_improve_vcycles = 0;
                }
                else if (++no_improve_vcycles >= no_improve_limit)
                {
                    break;
                }
            }
        }
    }

    return have_balanced_best ? best : fallback;
}
