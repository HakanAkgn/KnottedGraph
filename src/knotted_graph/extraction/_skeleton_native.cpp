#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

namespace py = pybind11;

namespace {

using I3 = std::array<int, 3>;
using D3 = std::array<double, 3>;

struct SparseSkeleton {
    I3 shape{0, 0, 0};
    I3 origin{0, 0, 0};
    std::vector<I3> coords;
    std::vector<std::vector<int>> adjacency;
    std::vector<std::vector<int>> components;
};

struct EdgeData {
    int u{-1};
    int v{-1};
    std::vector<D3> points;
    double weight{0.0};
};

struct GraphData {
    std::vector<std::pair<int, D3>> nodes;
    std::vector<EdgeData> edges;
};

const std::vector<I3>& neighbor_offsets() {
    static const std::vector<I3> offsets = [] {
        std::vector<I3> result;
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz)
                    if (dx || dy || dz) result.push_back({dx, dy, dz});
        return result;
    }();
    return offsets;
}

int norm_sq(const I3& value) {
    return value[0] * value[0] + value[1] * value[1] + value[2] * value[2];
}

const std::vector<std::vector<I3>>& shorter_intermediates() {
    static const std::vector<std::vector<I3>> table = [] {
        std::vector<std::vector<I3>> result;
        const auto& offsets = neighbor_offsets();
        result.resize(offsets.size());
        for (std::size_t index = 0; index < offsets.size(); ++index) {
            const I3 delta = offsets[index];
            const int edge_sq = norm_sq(delta);
            if (edge_sq <= 1) continue;
            for (const I3& step : offsets) {
                I3 remainder{
                    delta[0] - step[0],
                    delta[1] - step[1],
                    delta[2] - step[2],
                };
                if ((remainder[0] == 0 && remainder[1] == 0 && remainder[2] == 0) ||
                    std::abs(remainder[0]) > 1 ||
                    std::abs(remainder[1]) > 1 ||
                    std::abs(remainder[2]) > 1)
                    continue;
                if (norm_sq(step) < edge_sq && norm_sq(remainder) < edge_sq)
                    result[index].push_back(step);
            }
        }
        return result;
    }();
    return table;
}

std::size_t local_flat(const I3& local, const I3& shape) {
    return (
        static_cast<std::size_t>(local[0]) * static_cast<std::size_t>(shape[1])
        + static_cast<std::size_t>(local[1])
    ) * static_cast<std::size_t>(shape[2]) + static_cast<std::size_t>(local[2]);
}

bool in_bounds(const I3& local, const I3& shape) {
    return (
        local[0] >= 0 && local[0] < shape[0] &&
        local[1] >= 0 && local[1] < shape[1] &&
        local[2] >= 0 && local[2] < shape[2]
    );
}

std::vector<std::vector<int>> connected_components(
    const std::vector<std::vector<int>>& adjacency
) {
    std::vector<char> seen(adjacency.size(), 0);
    std::vector<std::vector<int>> components;
    for (int seed = 0; seed < static_cast<int>(adjacency.size()); ++seed) {
        if (seen[static_cast<std::size_t>(seed)]) continue;
        std::vector<int> stack{seed};
        seen[static_cast<std::size_t>(seed)] = 1;
        std::vector<int> component;
        while (!stack.empty()) {
            const int u = stack.back();
            stack.pop_back();
            component.push_back(u);
            const auto& row = adjacency[static_cast<std::size_t>(u)];
            for (auto it = row.rbegin(); it != row.rend(); ++it) {
                const int v = *it;
                if (!seen[static_cast<std::size_t>(v)]) {
                    seen[static_cast<std::size_t>(v)] = 1;
                    stack.push_back(v);
                }
            }
        }
        std::sort(component.begin(), component.end());
        components.push_back(std::move(component));
    }
    return components;
}

SparseSkeleton build_sparse(
    py::array_t<bool, py::array::c_style | py::array::forcecast> image
) {
    const auto view = image.unchecked<3>();
    I3 full_shape{
        static_cast<int>(view.shape(0)),
        static_cast<int>(view.shape(1)),
        static_cast<int>(view.shape(2)),
    };

    I3 lo{full_shape[0], full_shape[1], full_shape[2]};
    I3 hi{-1, -1, -1};
    for (int x = 0; x < full_shape[0]; ++x)
        for (int y = 0; y < full_shape[1]; ++y)
            for (int z = 0; z < full_shape[2]; ++z)
                if (view(x, y, z)) {
                    lo[0] = std::min(lo[0], x);
                    lo[1] = std::min(lo[1], y);
                    lo[2] = std::min(lo[2], z);
                    hi[0] = std::max(hi[0], x);
                    hi[1] = std::max(hi[1], y);
                    hi[2] = std::max(hi[2], z);
                }

    SparseSkeleton sparse;
    if (hi[0] < 0) return sparse;

    sparse.origin = lo;
    sparse.shape = {
        hi[0] - lo[0] + 1,
        hi[1] - lo[1] + 1,
        hi[2] - lo[2] + 1,
    };
    const std::size_t volume =
        static_cast<std::size_t>(sparse.shape[0]) *
        static_cast<std::size_t>(sparse.shape[1]) *
        static_cast<std::size_t>(sparse.shape[2]);
    std::vector<int> index_map(volume, -1);

    for (int lx = 0; lx < sparse.shape[0]; ++lx) {
        const int x = lo[0] + lx;
        for (int ly = 0; ly < sparse.shape[1]; ++ly) {
            const int y = lo[1] + ly;
            for (int lz = 0; lz < sparse.shape[2]; ++lz) {
                const int z = lo[2] + lz;
                if (!view(x, y, z)) continue;
                const int index = static_cast<int>(sparse.coords.size());
                sparse.coords.push_back({x, y, z});
                index_map[local_flat({lx, ly, lz}, sparse.shape)] = index;
            }
        }
    }

    sparse.adjacency.resize(sparse.coords.size());
    const auto& offsets = neighbor_offsets();
    const auto& intermediates = shorter_intermediates();
    for (int u = 0; u < static_cast<int>(sparse.coords.size()); ++u) {
        const I3 global = sparse.coords[static_cast<std::size_t>(u)];
        const I3 local{
            global[0] - lo[0],
            global[1] - lo[1],
            global[2] - lo[2],
        };
        auto& row = sparse.adjacency[static_cast<std::size_t>(u)];
        for (std::size_t offset_index = 0; offset_index < offsets.size(); ++offset_index) {
            const I3 delta = offsets[offset_index];
            const I3 target{
                local[0] + delta[0],
                local[1] + delta[1],
                local[2] + delta[2],
            };
            if (!in_bounds(target, sparse.shape)) continue;
            const int v = index_map[local_flat(target, sparse.shape)];
            if (v < 0) continue;

            bool redundant = false;
            for (const I3& step : intermediates[offset_index]) {
                const I3 middle{
                    local[0] + step[0],
                    local[1] + step[1],
                    local[2] + step[2],
                };
                if (!in_bounds(middle, sparse.shape)) continue;
                if (index_map[local_flat(middle, sparse.shape)] >= 0) {
                    redundant = true;
                    break;
                }
            }
            if (!redundant) row.push_back(v);
        }
    }

    sparse.components = connected_components(sparse.adjacency);
    return sparse;
}

D3 to_double(const I3& value) {
    return {
        static_cast<double>(value[0]),
        static_cast<double>(value[1]),
        static_cast<double>(value[2]),
    };
}

bool same_point(const D3& left, const D3& right) {
    return left[0] == right[0] && left[1] == right[1] && left[2] == right[2];
}

double distance(const D3& left, const D3& right) {
    const double dx = left[0] - right[0];
    const double dy = left[1] - right[1];
    const double dz = left[2] - right[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

std::uint64_t edge_tag(int u, int v) {
    if (u > v) std::swap(u, v);
    return (
        static_cast<std::uint64_t>(static_cast<std::uint32_t>(u)) << 32
    ) | static_cast<std::uint32_t>(v);
}

std::vector<D3> drop_consecutive_duplicates_exact(std::vector<D3> points) {
    if (points.empty()) return points;
    std::vector<D3> kept;
    kept.reserve(points.size());
    kept.push_back(points.front());
    for (std::size_t index = 1; index < points.size(); ++index) {
        if (!same_point(points[index], points[index - 1]))
            kept.push_back(points[index]);
    }
    return kept;
}

double path_weight(const std::vector<D3>& points) {
    double total = 0.0;
    for (std::size_t index = 1; index < points.size(); ++index)
        total += distance(points[index - 1], points[index]);
    return total;
}

std::vector<int> trace_cycle(
    const std::vector<int>& component,
    const std::vector<std::vector<int>>& adjacency
) {
    const int start = component.front();
    int previous = -1;
    int current = start;
    std::vector<int> order{start};
    for (std::size_t guard = 0; guard <= component.size(); ++guard) {
        std::vector<int> candidates;
        for (int v : adjacency[static_cast<std::size_t>(current)])
            if (v != previous) candidates.push_back(v);
        if (candidates.empty()) break;
        const int next = candidates.front();
        if (next == start) {
            order.push_back(start);
            return order;
        }
        order.push_back(next);
        previous = current;
        current = next;
    }
    throw std::runtime_error("pure skeleton cycle could not be traced");
}

struct LocalGraph {
    std::vector<D3> node_positions;
    std::vector<EdgeData> edges;
};

LocalGraph trace_component(
    const SparseSkeleton& sparse,
    const std::vector<int>& component,
    int hops
) {
    LocalGraph graph;
    if (component.empty()) return graph;

    std::vector<char> in_component(sparse.coords.size(), 0);
    for (int index : component) in_component[static_cast<std::size_t>(index)] = 1;

    std::vector<int> special;
    for (int index : component)
        if (sparse.adjacency[static_cast<std::size_t>(index)].size() != 2)
            special.push_back(index);

    if (special.empty()) {
        const auto order = trace_cycle(component, sparse.adjacency);
        std::vector<D3> points;
        points.reserve(order.size());
        for (int index : order) points.push_back(to_double(sparse.coords[static_cast<std::size_t>(index)]));
        graph.node_positions.push_back(points.front());
        EdgeData edge;
        edge.u = 0;
        edge.v = 0;
        edge.points = std::move(points);
        edge.weight = path_weight(edge.points);
        graph.edges.push_back(std::move(edge));
        return graph;
    }

    std::vector<int> dist(sparse.coords.size(), -1);
    std::vector<int> frontier = special;
    for (int seed : special) dist[static_cast<std::size_t>(seed)] = 0;
    for (int depth = 0; depth < hops; ++depth) {
        std::vector<int> next;
        for (int u : frontier) {
            for (int v : sparse.adjacency[static_cast<std::size_t>(u)]) {
                if (!in_component[static_cast<std::size_t>(v)]) continue;
                if (dist[static_cast<std::size_t>(v)] >= 0) continue;
                dist[static_cast<std::size_t>(v)] = depth + 1;
                next.push_back(v);
            }
        }
        frontier.swap(next);
        if (frontier.empty()) break;
    }

    std::vector<char> in_zone(sparse.coords.size(), 0);
    for (int index : component)
        if (dist[static_cast<std::size_t>(index)] >= 0)
            in_zone[static_cast<std::size_t>(index)] = 1;

    std::vector<char> seen(sparse.coords.size(), 0);
    std::vector<std::vector<int>> zones;
    for (int seed : component) {
        if (!in_zone[static_cast<std::size_t>(seed)] ||
            seen[static_cast<std::size_t>(seed)])
            continue;
        std::vector<int> stack{seed};
        seen[static_cast<std::size_t>(seed)] = 1;
        std::vector<int> zone;
        while (!stack.empty()) {
            const int u = stack.back();
            stack.pop_back();
            zone.push_back(u);
            const auto& row = sparse.adjacency[static_cast<std::size_t>(u)];
            for (auto it = row.rbegin(); it != row.rend(); ++it) {
                const int v = *it;
                if (
                    in_zone[static_cast<std::size_t>(v)] &&
                    !seen[static_cast<std::size_t>(v)]
                ) {
                    seen[static_cast<std::size_t>(v)] = 1;
                    stack.push_back(v);
                }
            }
        }
        std::sort(zone.begin(), zone.end());
        zones.push_back(std::move(zone));
    }
    std::sort(zones.begin(), zones.end(), [](const auto& left, const auto& right) {
        return left.front() < right.front();
    });

    std::vector<int> zone_of(sparse.coords.size(), -1);
    for (int zone_index = 0; zone_index < static_cast<int>(zones.size()); ++zone_index) {
        const auto& zone = zones[static_cast<std::size_t>(zone_index)];
        D3 mean{0.0, 0.0, 0.0};
        for (int voxel : zone) {
            const I3 value = sparse.coords[static_cast<std::size_t>(voxel)];
            mean[0] += value[0];
            mean[1] += value[1];
            mean[2] += value[2];
            zone_of[static_cast<std::size_t>(voxel)] = zone_index;
        }
        const double scale = 1.0 / static_cast<double>(zone.size());
        mean[0] = std::nearbyint(mean[0] * scale);
        mean[1] = std::nearbyint(mean[1] * scale);
        mean[2] = std::nearbyint(mean[2] * scale);
        graph.node_positions.push_back(mean);
    }

    std::unordered_set<std::uint64_t> visited;
    for (int source_zone = 0; source_zone < static_cast<int>(zones.size()); ++source_zone) {
        for (int source_voxel : zones[static_cast<std::size_t>(source_zone)]) {
            for (int neighbour : sparse.adjacency[static_cast<std::size_t>(source_voxel)]) {
                if (in_zone[static_cast<std::size_t>(neighbour)]) continue;
                const std::uint64_t first = edge_tag(source_voxel, neighbour);
                if (visited.count(first)) continue;

                std::vector<int> path{source_voxel, neighbour};
                visited.insert(first);
                int previous = source_voxel;
                int current = neighbour;
                while (!in_zone[static_cast<std::size_t>(current)]) {
                    std::vector<int> candidates;
                    for (int v : sparse.adjacency[static_cast<std::size_t>(current)])
                        if (v != previous) candidates.push_back(v);
                    if (candidates.empty()) break;
                    int next = candidates.front();
                    for (int candidate : candidates) {
                        if (!visited.count(edge_tag(current, candidate))) {
                            next = candidate;
                            break;
                        }
                    }
                    visited.insert(edge_tag(current, next));
                    previous = current;
                    current = next;
                    path.push_back(current);
                    if (path.size() > component.size() + 2)
                        throw std::runtime_error("skeleton path tracing did not terminate");
                }
                if (!in_zone[static_cast<std::size_t>(current)]) continue;

                const int target_zone = zone_of[static_cast<std::size_t>(current)];
                std::vector<D3> points;
                points.reserve(path.size());
                for (int voxel : path)
                    points.push_back(to_double(sparse.coords[static_cast<std::size_t>(voxel)]));
                points.front() = graph.node_positions[static_cast<std::size_t>(source_zone)];
                points.back() = graph.node_positions[static_cast<std::size_t>(target_zone)];
                points = drop_consecutive_duplicates_exact(std::move(points));
                if (points.size() < 2) continue;

                EdgeData edge;
                edge.u = source_zone;
                edge.v = target_zone;
                edge.weight = path_weight(points);
                edge.points = std::move(points);
                graph.edges.push_back(std::move(edge));
            }
        }
    }

    if (hops > 0) {
        const std::size_t max_zero_loop_points =
            static_cast<std::size_t>(std::max(4, 2 * hops + 3));
        graph.edges.erase(
            std::remove_if(
                graph.edges.begin(),
                graph.edges.end(),
                [&](const EdgeData& edge) {
                    return edge.u == edge.v && edge.points.size() <= max_zero_loop_points;
                }
            ),
            graph.edges.end()
        );
    }
    return graph;
}

GraphData trace_all(const SparseSkeleton& sparse, int hops) {
    GraphData out;
    int next_id = 0;
    for (const auto& component : sparse.components) {
        LocalGraph local = trace_component(sparse, component, hops);
        std::vector<char> active(local.node_positions.size(), 0);
        for (const auto& edge : local.edges) {
            active[static_cast<std::size_t>(edge.u)] = 1;
            active[static_cast<std::size_t>(edge.v)] = 1;
        }
        for (int local_id = 0; local_id < static_cast<int>(local.node_positions.size()); ++local_id)
            if (active[static_cast<std::size_t>(local_id)])
                out.nodes.emplace_back(
                    next_id + local_id,
                    local.node_positions[static_cast<std::size_t>(local_id)]
                );
        for (auto edge : local.edges) {
            edge.u += next_id;
            edge.v += next_id;
            out.edges.push_back(std::move(edge));
        }
        next_id += static_cast<int>(local.node_positions.size());
    }
    return out;
}

bool geometry_safe(const GraphData& graph) {
    for (const auto& edge : graph.edges) {
        if (edge.points.size() < 2) return false;
        if (edge.u == edge.v) {
            if (edge.points.size() < 3) return false;
            const D3 start = edge.points.front();
            bool excursion = false;
            for (std::size_t index = 1; index + 1 < edge.points.size(); ++index) {
                if (distance(edge.points[index], start) > 1e-10) {
                    excursion = true;
                    break;
                }
            }
            if (!excursion) return false;
        }
    }
    return true;
}

py::dict compact_candidate(const GraphData& graph, int hops) {
    py::dict result;
    result["hops"] = hops;
    py::list nodes;
    for (const auto& [node, _position] : graph.nodes) nodes.append(node);
    result["nodes"] = nodes;

    py::list edges;
    for (const auto& edge : graph.edges) {
        edges.append(py::make_tuple(
            edge.u,
            edge.v,
            edge.weight,
            static_cast<int>(edge.points.size())
        ));
    }
    result["edges"] = edges;
    result["geometry_safe"] = geometry_safe(graph);
    return result;
}

py::dict materialized_candidate(const GraphData& graph, int hops) {
    py::dict result;
    result["hops"] = hops;
    py::list nodes;
    for (const auto& [node, position] : graph.nodes)
        nodes.append(py::make_tuple(node, position));
    result["nodes"] = nodes;

    py::list edges;
    for (const auto& edge : graph.edges)
        edges.append(py::make_tuple(edge.u, edge.v, edge.points, edge.weight));
    result["edges"] = edges;
    return result;
}

py::tuple sparse_to_python(const SparseSkeleton& sparse) {
    py::list coords;
    for (const I3& point : sparse.coords) coords.append(point);
    py::list adjacency;
    for (const auto& row : sparse.adjacency) adjacency.append(row);
    return py::make_tuple(coords, adjacency);
}

}  // namespace

PYBIND11_MODULE(_skeleton_native, module) {
    module.doc() = "Native exact sparse skeleton adjacency and multi-scale tracing";

    module.def(
        "sparse_adjacency",
        [](py::array_t<bool, py::array::c_style | py::array::forcecast> image) {
            SparseSkeleton sparse = build_sparse(image);
            return sparse_to_python(sparse);
        }
    );

    module.def(
        "plan_candidates",
        [](
            py::array_t<bool, py::array::c_style | py::array::forcecast> image,
            int max_hops
        ) {
            if (max_hops < 0) throw std::invalid_argument("max_hops must be non-negative");
            SparseSkeleton sparse = build_sparse(image);
            std::vector<GraphData> candidates;
            candidates.reserve(static_cast<std::size_t>(max_hops + 1));
            for (int hops = 0; hops <= max_hops; ++hops)
                candidates.push_back(trace_all(sparse, hops));
            py::list out;
            for (int hops = 0; hops <= max_hops; ++hops)
                out.append(compact_candidate(candidates[static_cast<std::size_t>(hops)], hops));
            return out;
        }
    );

    module.def(
        "materialize_candidate",
        [](
            py::array_t<bool, py::array::c_style | py::array::forcecast> image,
            int hops
        ) {
            if (hops < 0) throw std::invalid_argument("hops must be non-negative");
            SparseSkeleton sparse = build_sparse(image);
            GraphData graph = trace_all(sparse, hops);
            return materialized_candidate(graph, hops);
        }
    );
}
