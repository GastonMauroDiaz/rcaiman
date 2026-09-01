// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// ============================================================
// HEALPix Math
// ============================================================

inline uint64_t spreadbits(uint32_t v) {
  uint64_t x = v;
  x = (x | (x << 16)) & 0x0000FFFF0000FFFFULL;
  x = (x | (x << 8))  & 0x00FF00FF00FF00FFULL;
  x = (x | (x << 4))  & 0x0F0F0F0F0F0F0F0FULL;
  x = (x | (x << 2))  & 0x3333333333333333ULL;
  x = (x | (x << 1))  & 0x5555555555555555ULL;
  return x;
}

inline int xyf2pixNest(int nside, int order, int ix, int iy, int face) {
  uint64_t pix = ((uint64_t)face << (2 * order))
  + spreadbits(ix)
  + (spreadbits(iy) << 1);

  return (int)(pix + 1);
}

inline int calcNestPosForPole(
    int nside,
    int order,
    double z,
    double z_abs,
    double scaled_phi
) {
  int ntt = (int)std::floor(scaled_phi);
  if(ntt >= 4) ntt = 3;

  double tp = scaled_phi - ntt;
  double tmp = nside * std::sqrt(3.0 * (1.0 - z_abs));

  int jp = (int)std::floor(tp * tmp);
  int jm = (int)std::floor((1.0 - tp) * tmp);

  jp = std::min(jp, nside - 1);
  jm = std::min(jm, nside - 1);

  int face, ix, iy;

  if(z >= 0.0) {
    face = ntt;
    ix = nside - jm - 1;
    iy = nside - jp - 1;
  } else {
    face = ntt + 8;
    ix = jp;
    iy = jm;
  }

  return xyf2pixNest(nside, order, ix, iy, face);
}

inline int calcNestPosForEquator(
    int nside,
    int order,
    double z,
    double scaled_phi
) {
  double temp1 = nside * (0.5 + scaled_phi);
  double temp2 = nside * (z * 0.75);

  int jp = (int)std::floor(temp1 - temp2);
  int jm = (int)std::floor(temp1 + temp2);

  int ifp = jp / nside;
  int ifm = jm / nside;

  int face;

  if(ifp == ifm) {
    face = (ifp == 4) ? 4 : ifp + 4;
  } else if(ifp < ifm) {
    face = ifp;
  } else {
    face = ifm + 8;
  }

  int mod_jm = jm % nside;
  if(mod_jm < 0) mod_jm += nside;

  int mod_jp = jp % nside;
  if(mod_jp < 0) mod_jp += nside;

  int ix = mod_jm;
  int iy = nside - mod_jp - 1;

  return xyf2pixNest(nside, order, ix, iy, face);
}

inline int ang2pixNest(
    int nside,
    int order,
    double theta,
    double phi
) {
  double z = std::cos(theta);
  double z_abs = std::abs(z);

  double scaled_phi = std::fmod(phi * 2.0 / M_PI, 4.0);
  if(scaled_phi < 0.0) scaled_phi += 4.0;

  if(z_abs > 2.0/3.0) {
    return calcNestPosForPole(
      nside,
      order,
      z,
      z_abs,
      scaled_phi
    );
  } else {
    return calcNestPosForEquator(
      nside,
      order,
      z,
      scaled_phi
    );
  }
}

// ============================================================
// Tree structures
// ============================================================

struct Node {
  int depth;
  int hp_id;
  std::vector<int> pixel_ids;
};

double compute_sd(
    const NumericVector& values,
    const std::vector<int>& ids
) {
  int n = ids.size();

  if(n < 2) return 0;

  double mean = 0;

  for(int i : ids)
    mean += values[i];

  mean /= n;

  double var = 0;

  for(int i : ids) {
    double d = values[i] - mean;
    var += d*d;
  }

  return std::sqrt(var/(n-1));
}

// ============================================================
// Functions
// ============================================================


// [[Rcpp::export]]
IntegerVector healpix_cpp(
    NumericVector zenith_deg,
    NumericVector azimuth_deg,
    int nside
)
{
  int n =
    zenith_deg.size();

  IntegerVector out(n);

  int order =
    (int)std::round(
        std::log2(
          (double)nside
        )
    );

  const double pi =
    M_PI;

  for(
    int i=0;
    i<n;
    i++
  )
  {
    if(
      NumericVector::is_na(
        zenith_deg[i]
      ) ||
        NumericVector::is_na(
          azimuth_deg[i]
        )
    )
    {
      out[i]=NA_INTEGER;
      continue;
    }

    double theta =
      zenith_deg[i]*
      pi/180.0;

    double phi =
      azimuth_deg[i]*
      pi/180.0;

    phi =
      M_PI/2.0 -
      phi;

    if(phi<0)
      phi+=2*pi;

    out[i] =
      ang2pixNest(
        nside,
        order,
        theta,
        phi
      );
  }

  return out;
}


// [[Rcpp::export]]
IntegerVector healpixtree_cpp(
    NumericVector z,
    NumericVector a,
    NumericVector values,
    double scale_parameter = 0,
    int max_splittings = 3,
    int min_size_px = 20
) {

  int n = z.size();

  IntegerVector labels(n, NA_INTEGER);

  std::vector<Node> stack;
  std::vector<std::vector<int>> base_map(12);

  // Base level

  for(int px=0; px<n; ++px) {

    if(!NumericVector::is_na(z[px]) &&
       !NumericVector::is_na(a[px])) {

       int hp_id =
       ang2pixNest(
         1,
         0,
         z[px],
          a[px]
       );

      base_map[hp_id-1].push_back(px);
    }
  }

  for(int i=0; i<12; ++i) {

    if(!base_map[i].empty()) {

      Node root;

      root.depth = 0;
      root.hp_id = i + 1;
      root.pixel_ids = base_map[i];

      stack.push_back(root);
    }
  }

  int current_label = 1;

  while(!stack.empty()) {

    Node current = stack.back();
    stack.pop_back();

    if(current.depth >= max_splittings ||
       current.pixel_ids.size() < (size_t)min_size_px) {

      for(int px : current.pixel_ids)
        labels[px] = current_label;

      current_label++;

      continue;
    }

    double sd_parent =
      compute_sd(
        values,
        current.pixel_ids
      );

    int next_order = current.depth + 1;
    int next_nside = 1 << next_order;

    std::vector<std::vector<int>> children_px(4);

    // Primer hijo real del nodo actual
    int first_child =
      (current.hp_id - 1) * 4 + 1;

    for(int px : current.pixel_ids) {

      int cid =
        ang2pixNest(
          next_nside,
          next_order,
          z[px],
           a[px]
        );

      int local_child =
        cid - first_child;

      if(local_child >= 0 &&
         local_child < 4) {

        children_px[local_child]
        .push_back(px);
      }
    }

    double sd_children_sum = 0;

    int valid_children_count = 0;

    bool valid_split = true;

    for(int i=0; i<4; ++i) {

      size_t csize =
        children_px[i].size();

      if(csize > 0) {

        if(csize < (size_t)min_size_px) {

          valid_split = false;
          break;
        }

        sd_children_sum +=
          compute_sd(
            values,
            children_px[i]
          );

        valid_children_count++;
      }
    }

    if(valid_children_count < 2)
      valid_split = false;

    if(valid_split &&
       (sd_children_sum - sd_parent)
         > scale_parameter) {

      for(int i=0; i<4; ++i) {

        if(children_px[i].size()>0) {

          Node child;

          child.depth = next_order;
          child.hp_id =
            first_child + i;

          child.pixel_ids =
            children_px[i];

          stack.push_back(child);
        }
      }

    } else {

      for(int px : current.pixel_ids)
        labels[px] = current_label;

      current_label++;
    }
  }

  return labels;
}
