"""Materialize the reviewed native-class change, without import-time patches.

Only NodalSkeleton.skeleton_graph and clear_cache are changed. The script is
idempotent and refuses an unexpected source body. CI commits this source edit
only after the complete suite and regression audit pass on the edited tree.
"""
from pathlib import Path

PATH = Path("src/knotted_graph/applications/nodal/skeleton.py")
NEW = '''    def skeleton_graph(
        self,
        simplify: bool = True,
        smooth_epsilon: float = 0,
        *,
        skeleton_image: Optional[NDArray] = None,
        reconstruction: str = "guarded",
        cubical_options: dict | None = None,
    ) -> nx.MultiGraph:
        """Reconstruct an index-coordinate graph with explicit evidence.

        ``guarded`` supplies source component/cycle counts and checks cleanup.
        ``cubical`` instead replays elementary collapses of the source voxel
        solid and retains the certificate as ``self.spine_certificate``. It
        permits neither an external skeleton nor geometric smoothing. Neither
        mode by itself certifies correspondence to the analytic source field.
        Existing coordinate conversion and visualization APIs are unchanged.
        """
        from ._reconstruction import reconstruct

        return reconstruct(
            self, simplify=simplify, smooth_epsilon=smooth_epsilon,
            skeleton_image=skeleton_image, reconstruction=reconstruction,
            cubical_options=cubical_options, extract=skeleton_image_to_graph,
            prune=remove_leaf_nodes, simplify_edges=simplify_edges,
            smooth_edges=smooth_edges, is_trivalent=is_trivalent,
        )

'''


def main():
    text = PATH.read_text()
    start = text.index("    def skeleton_graph(\n")
    stop = text.index("    @property\n    def total_edge_pts", start)
    old = text[start:stop]
    if "from ._reconstruction import reconstruct" not in old:
        if "G = skeleton_image_to_graph(image)" not in old:
            raise RuntimeError("unexpected NodalSkeleton implementation")
        text = text[:start] + NEW + "\n" + text[stop:]
    marker = "        self._pv_data_args = None\n\n\n    @cached_property\n    def fields_pv"
    replacement = ("        self._pv_data_args = None\n"
                   "        self._skeleton_graph_digest = None\n"
                   "        self.spine_certificate = None\n"
                   "        self.__dict__.pop(\"PDCode\", None)\n\n\n"
                   "    @cached_property\n    def fields_pv")
    if marker in text:
        text = text.replace(marker, replacement, 1)
    elif replacement not in text:
        raise RuntimeError("unexpected clear_cache implementation")
    PATH.write_text(text)


if __name__ == "__main__":
    main()
