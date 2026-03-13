// --------------------------------------------------------------------------
// Delaunay refinement heuristics
//
// Goal: for meshes with few vertices (nVertices ≈ k), grow toward
// targetN = min(4*k, 128) vertices by driving a global circumradius criterion.
// Angle quality is relaxed smoothly as n grows and refinement is disabled
// entirely for large meshes (n ≥ 8192).
//
// ── angleThreshDegrees ────────────────────────────────────────────────────
//
//  Minimum angle guarantee for Ruppert's algorithm.  Decreasing with nVertices
//  prevents excessive Steiner insertions on meshes that are already large.
//  Note: convergence is guaranteed only for ≤ 20.7°; 25° is empirically
//  stable (geometry-central itself warns for values > 30°).
//
//   deg ▲
//       │
//    25 ┤x x
//       │     x
//    15 │       x
//       │         x
//     0 ┤            x
//       └──┬────┬────┬──▶ n
//          96  4096  8192
//
// ── circumradiusThresh ────────────────────────────────────────────────────
//
//  Absolute face-size criterion that drives the mesh toward targetN vertices.
//  Every face with circumradius > threshold is split; the algorithm terminates
//  only when the criterion holds globally — guaranteeing spatial uniformity.
//  (maxInsertions would stop mid-queue, leaving inhomogeneities that degrade
//  eigenvector accuracy.)
//
//  Derivation: a uniform triangulation of area A with targetN equilateral
//  faces has expected circumradius ≈ sqrt(A / targetN) per face.  With
//  shapeLengthScale = sqrt(A), this simplifies to shapeLengthScale / sqrt(targetN).
//
//        ▲
//    ∞   ┤        x x x x
//        │        x
//  s/√tN ┤x x x x x            tN = min(4k, 128)
//     0  ┤
//        └────────┬───────▶ n
//               targetN
//
// --------------------------------------------------------------------------
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

namespace geometrycentral {
namespace surface {

static constexpr size_t REFINE_TARGET_FACTOR = 4;   // targetN = REFINE_TARGET_FACTOR * k
static constexpr size_t REFINE_TARGET_MAX    = 128; // cap targetN to avoid over-refinement

// Linear interpolation, clamped to [x0, x1].
inline double lerpClamped(double x, double x0, double x1, double y0, double y1)
{
    const double t = std::max(0.0, std::min(1.0, (x - x0) / (x1 - x0)));
    return y0 + (y1 - y0) * t;
}

// Whether Delaunay refinement should run for a mesh of the given size.
// Refinement is disabled for large meshes (n >= 8192).
inline bool shouldRefine(size_t nVertices)
{
    return nVertices < 8192;
}

// Angle threshold in degrees for intrinsic Delaunay refinement.
// Returns std::nullopt for n >= 8192 — the call site should skip refinement.
inline std::optional<double> computeRefinementAngleThresh(size_t nVertices)
{
    if (!shouldRefine(nVertices))
        return std::nullopt;
    if (nVertices < 4096)
        return lerpClamped(static_cast<double>(nVertices), 96.0, 4096.0, 25.0, 15.0);
    return lerpClamped(static_cast<double>(nVertices), 4096.0, 8192.0, 15.0, 0.0);
}

// Absolute circumradius threshold that drives the mesh toward targetN = min(4k, 128) vertices.
// shapeLengthScale = sqrt(totalArea), so sqrt(A / targetN) = shapeLengthScale / sqrt(targetN).
// Returns infinity when the mesh already has >= targetN vertices (only angle quality then matters).
inline double computeRefinementCircumradiusThresh(size_t nVertices, size_t k, double shapeLengthScale)
{
    const size_t targetN = std::min(REFINE_TARGET_FACTOR * k, REFINE_TARGET_MAX);
    if (nVertices >= targetN)
        return std::numeric_limits<double>::infinity();
    return shapeLengthScale / std::sqrt(static_cast<double>(targetN));
}

} // namespace surface
} // namespace geometrycentral
