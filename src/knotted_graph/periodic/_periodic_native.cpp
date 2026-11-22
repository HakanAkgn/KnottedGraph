#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>
#include <cstdint>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

namespace py = pybind11;

namespace {

using Index = std::int64_t;
using Shift = std::array<Index, 3>;
using Cycle = std::array<Index, 4>;

Shift add(const Shift& a, const Shift& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

Shift subtract(const Shift& a, const Shift& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

class WindingPartition {
public:
    explicit WindingPartition(Index count)
        : parent_(count), size_(count, 1), to_parent_(count) {
        std::iota(parent_.begin(), parent_.end(), 0LL);
    }

    void join(Index from, Index to, const Shift& step) {
        Index from_root = find(from);
        Index to_root = find(to);
        // Deck shift of to_root relative to from_root imposed by this edge.
        Shift relative = subtract(add(to_parent_[from], step), to_parent_[to]);
        if (from_root == to_root) {
            if (relative != Shift{}) cycles_[from_root].insert(relative);
            return;
        }
        if (size_[from_root] < size_[to_root]) {
            std::swap(from_root, to_root);
            relative = subtract(Shift{}, relative);
        }
        parent_[to_root] = from_root;
        to_parent_[to_root] = relative;
        size_[from_root] += size_[to_root];
        auto moved = cycles_.extract(to_root);
        if (!moved.empty()) cycles_[from_root].merge(moved.mapped());
    }

    std::vector<Cycle> labels_and_cycles(Index* labels) {
        std::vector<Index> component(parent_.size(), -1);
        Index next_label = 0;
        for (Index vertex = 0; vertex < static_cast<Index>(parent_.size()); ++vertex) {
            auto root = find(vertex);
            auto& label = component[root];
            if (label < 0) label = next_label++;
            labels[vertex] = label;
        }
        std::vector<Cycle> result;
        for (const auto& [root, vectors] : cycles_) {
            for (const auto& w : vectors) result.push_back({component[root], w[0], w[1], w[2]});
        }
        return result;
    }

private:
    Index find(Index vertex) {
        const auto previous = parent_[vertex];
        if (previous == vertex) return vertex;
        const auto root = find(previous);
        to_parent_[vertex] = add(to_parent_[vertex], to_parent_[previous]);
        parent_[vertex] = root;
        return root;
    }

    std::vector<Index> parent_;
    std::vector<Index> size_;
    std::vector<Shift> to_parent_;
    std::map<Index, std::set<Shift>> cycles_;
};

py::tuple grid_partition(const py::array_t<bool, py::array::c_style>& occupied) {
    if (occupied.ndim() != 3) throw std::invalid_argument("expected a three-dimensional mask");
    const auto nx = occupied.shape(0);
    const auto ny = occupied.shape(1);
    const auto nz = occupied.shape(2);
    const auto* sign = occupied.data();
    py::array_t<Index> labels({nx, ny, nz});
    auto* output = labels.mutable_data();
    std::vector<Cycle> cycles;
    {
        py::gil_scoped_release release;
        WindingPartition partition(occupied.size());
        for (Index x = 0; x < nx; ++x) {
            for (Index y = 0; y < ny; ++y) {
                for (Index z = 0; z < nz; ++z) {
                    const auto start = (x * ny + y) * nz + z;
                    // The seven positive edge offsets of the Freudenthal triangulation.
                    for (int mask = 1; mask < 8; ++mask) {
                        const auto x1 = x + ((mask >> 2) & 1);
                        const auto y1 = y + ((mask >> 1) & 1);
                        const auto z1 = z + (mask & 1);
                        const auto end = (((x1 == nx ? 0 : x1) * ny + (y1 == ny ? 0 : y1)) * nz
                                          + (z1 == nz ? 0 : z1));
                        if (sign[start] == sign[end]) {
                            partition.join(start, end, {x1 == nx, y1 == ny, z1 == nz});
                        }
                    }
                }
            }
        }
        cycles = partition.labels_and_cycles(output);
    }
    return py::make_tuple(labels, cycles);
}

py::tuple graph_partition(
    Index count,
    const py::array_t<Index, py::array::c_style>& endpoints,
    const py::array_t<Index, py::array::c_style>& shifts
) {
    if (count < 0 || endpoints.ndim() != 2 || shifts.ndim() != 2 ||
        endpoints.shape(0) != shifts.shape(0) || endpoints.shape(1) != 2 || shifts.shape(1) != 3) {
        throw std::invalid_argument("invalid graph partition arrays");
    }
    const auto edges = endpoints.unchecked<2>();
    const auto deck = shifts.unchecked<2>();
    py::array_t<Index> labels(count);
    auto* output = labels.mutable_data();
    std::vector<Cycle> cycles;
    {
        py::gil_scoped_release release;
        WindingPartition partition(count);
        for (Index edge = 0; edge < edges.shape(0); ++edge) {
            const auto from = edges(edge, 0);
            const auto to = edges(edge, 1);
            if (from < 0 || from >= count || to < 0 || to >= count) {
                throw std::invalid_argument("graph edge index out of bounds");
            }
            partition.join(from, to, {deck(edge, 0), deck(edge, 1), deck(edge, 2)});
        }
        cycles = partition.labels_and_cycles(output);
    }
    return py::make_tuple(labels, cycles);
}

}  // namespace

PYBIND11_MODULE(_periodic_native, module) {
    module.def("grid_partition", &grid_partition);
    module.def("graph_partition", &graph_partition);
}
