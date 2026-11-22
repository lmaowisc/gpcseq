#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <algorithm>
#include <array>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
// Fenwick tree on ordered clinical states; a higher rank wins.
class RankCounts {
    std::vector<int> tree;
    int total = 0;
    int through(size_t rank) const {
        int count = 0;
        for (; rank; rank -= rank & -rank) count += tree[rank];
        return count;
    }
public:
    explicit RankCounts(size_t n) : tree(n + 1, 0) {}
    void add(size_t rank, int delta) {
        total += delta;
        for (; rank < tree.size(); rank += rank & -rank) tree[rank] += delta;
    }
    std::array<int,2> sides(size_t rank) const {
        return {{through(rank - 1), total - through(rank)}};
    }
};

void compare(const double* end, const int* dead, const int* arm, int n,
             const int* id, const double* time, R_xlen_t m, double* result) {
    std::vector<size_t> admissions(n,0);
    for (int i = 0; i < n; ++i)
        if (!R_FINITE(end[i]) || end[i] < 0 ||
            (arm[i] != 0 && arm[i] != 1) || (dead[i] != 0 && dead[i] != 1))
            throw std::invalid_argument("Invalid shared-window patient record.");
    std::vector<R_xlen_t> events;
    events.reserve(m);
    for (R_xlen_t e = 0; e < m; ++e) {
        if (id[e] < 1 || id[e] > n || !R_FINITE(time[e]) ||
            time[e] < 0)
            throw std::invalid_argument("Invalid shared-window admission.");
        // The public snapshot allows roundoff at the follow-up limit. The
        // comparison uses the exact window, as in the direct pair scorer.
        if (time[e] > end[id[e]-1]) continue;
        ++admissions[id[e]-1];
        events.push_back(e);
    }
    std::sort(events.begin(), events.end(), [&](R_xlen_t i, R_xlen_t j) {
        return time[i] < time[j];
    });
    std::vector<int> retire(n);
    std::iota(retire.begin(), retire.end(), 0);
    std::sort(retire.begin(), retire.end(), [&](int i, int j) {
        if (end[i] != end[j]) return end[i] < end[j];
        return arm[i] > arm[j]; // Treatment owns pairs with equal follow-up.
    });
    // Rank 1 is death, rank 2 is the largest admission count.
    const size_t best = 2 + *std::max_element(admissions.begin(), admissions.end());
    std::vector<size_t> rank(n,best);
    RankCounts active[] = {RankCounts(best), RankCounts(best)};
    RankCounts retired[] = {RankCounts(best), RankCounts(best)};
    std::vector<std::array<int,2>> score(n), start(n);
    for (int i = 0; i < n; ++i) active[arm[i]].add(best,1);
    auto collect = [&](int i) {
        const auto now = retired[1-arm[i]].sides(rank[i]);
        for (int side=0; side<2; ++side) score[i][side] += now[side]-start[i][side];
    };
    auto change = [&](int i, size_t to) {
        collect(i);
        active[arm[i]].add(rank[i],-1);
        active[arm[i]].add(to,1);
        rank[i] = to;
        start[i] = retired[1-arm[i]].sides(to);
    };
    size_t next_event = 0;
    for (int begin = 0; begin < n;) {
        const double t = end[retire[begin]];
        while (next_event < events.size() && time[events[next_event]] <= t) {
            const int i = id[events[next_event++]]-1;
            change(i,rank[i]-1);
        }
        int pos = begin;
        while (pos < n && end[retire[pos]] == t) {
            const int i = retire[pos++];
            if (dead[i]) change(i,1);
        }
        for (int j = begin; j < pos; ++j) {
            const int i = retire[j], a = arm[i];
            collect(i);
            const auto now = active[1-a].sides(rank[i]);
            for (int side=0; side<2; ++side) score[i][side] += now[side];
            active[a].add(rank[i],-1);
            retired[a].add(rank[i],1);
            result[i] = score[i][a==1 ? 0 : 1];
            result[static_cast<R_xlen_t>(n)+i] = score[i][a==1 ? 1 : 0];
        }
        begin = pos;
    }
}
}

extern "C" SEXP gpc_shared_count(SEXP end, SEXP dead, SEXP arm, SEXP id, SEXP time) {
    const R_xlen_t n=XLENGTH(end), m=XLENGTH(time);
    if (n < 1 || n > std::numeric_limits<int>::max() || TYPEOF(end)!=REALSXP ||
        TYPEOF(dead)!=LGLSXP || TYPEOF(arm)!=INTSXP || TYPEOF(id)!=INTSXP ||
        TYPEOF(time)!=REALSXP || XLENGTH(dead)!=n || XLENGTH(arm)!=n || XLENGTH(id)!=m)
        Rf_error("Invalid shared-window input types or lengths.");
    SEXP ans=PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(n), 2));
    try {
        compare(REAL(end),LOGICAL(dead),INTEGER(arm),static_cast<int>(n),
                INTEGER(id),REAL(time),m,REAL(ans));
    } catch (const std::exception& e) { Rf_error("%s",e.what()); }
    UNPROTECT(1);
    return ans;
}

extern "C" SEXP gpc_copy_channels(SEXP channels) {
    return Rf_duplicate(channels);
}

extern "C" SEXP gpc_mark_measurement_repeats(SEXP person, SEXP time, SEXP order);

static const R_CallMethodDef calls[] = {
    {"gpc_mark_measurement_repeats", (DL_FUNC) &gpc_mark_measurement_repeats, 3},
    {"gpc_copy_channels", (DL_FUNC) &gpc_copy_channels, 1},
    {"gpc_shared_count", (DL_FUNC) &gpc_shared_count, 5},
    {NULL, NULL, 0}
};

extern "C" void attribute_visible R_init_gpcseq(DllInfo* dll) {
    R_registerRoutines(dll, NULL, calls, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
