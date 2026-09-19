// Strict MRCA-node spans, matching ibd_jackknife.calculate_ibd_blocks_mrca.
// Called once per marginal tree; no segment is filtered until its span ends.
#include <cstdint>
extern "C" void scan(
    int n, const int32_t* u, const int32_t* v, const int32_t* pair_cell,
    const int32_t* parent, const double* time, double position, double cm_per_bp,
    int nb, const double* edges, int32_t* previous, double* start,
    int64_t* counts, double* lengths, int finish) {
  for (int k = 0; k < n; ++k) {
    int a = u[k], b = v[k];
    if (!finish) {
      while (a != b && a != -1 && b != -1) {
        if (time[a] <= time[b]) a = parent[a];
        else b = parent[b];
      }
      if (a == -1 || b == -1) a = -1;
    }
    if (finish || previous[k] != a) {
      if (previous[k] >= 0) {
        double len = (position - start[k]) * cm_per_bp;
        for (int j = 0; j < nb; ++j) {
          if (edges[j] < len && len <= edges[j+1]) {
            int cell = j * 9 + pair_cell[k];
            counts[cell] += 1;
            lengths[cell] += len;
            break;
          }
        }
      }
      previous[k] = a;
      start[k] = position;
    }
  }
}
