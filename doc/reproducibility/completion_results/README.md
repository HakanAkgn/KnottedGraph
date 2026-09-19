# Full-grid results package

Results_Report.pdf explains the numerical results, figures, certificate scope and remaining gaps.
Figures/ contains vector categorical maps and full-geometry representative views. Surface artists are rasterized without deleting faces.
PhaseDIagrams.pdf retains the five vertically stacked source-notebook panels; PhaseDIagrams_Compact.pdf supplies a balanced manuscript alternative.
PorousMaterialPhaseMap.pdf restores three spatial-signature maps plus all five original representative parameter points.
PorousMaterialContinuationMap.pdf uses analytic regularity-cover components instead. These are separate scientific quantities.
No small-region filter or abstract polynomial fallback is used.

The current 42-page manuscript source is not included. manuscript_insert.tex is an insertion draft, not a reconstructed complete paper.
The all-word formulas and existing cavity implementation are unchanged.

Build the report from this directory with:
latexmk -pdf -no-shell-escape -interaction=nonstopmode -halt-on-error Results_Report.tex

All full record arrays, event certificates, voxel-link tests and example comparison witnesses are retained in data/.
The five complete source/witness/graph archives are retained under representatives/.
Input-run IDs and frozen environment are documented in the report and repository workflow.

Scientific counts:
{
  "all_word_formulas_modified": false,
  "analytic_to_voxel_correspondence_proved": false,
  "cavity_implementation_modified": false,
  "complete_source_isotopy_classification_proved": false,
  "continuum_vs_digital": {
    "retained_digital_homology_conflicts": 12,
    "retained_spine_signature_variations": 12
  },
  "events": {
    "by_family": {
      "gyroid_to_diamond": {
        "area_fractions": {
          "event": 0.10187500000000001,
          "regular": 0.89765625,
          "unknown": 0.0004687500000000007
        },
        "event_area_refers_to_cells_not_the_critical_set": true,
        "exact_area_fractions": {
          "event": [
            634758290955882667445746815409007,
            6230756230241792923652294673694720
          ],
          "regular": [
            5593077272302984411558216593938633,
            6230756230241792923652294673694720
          ],
          "unknown": [
            73016674573146116208281608677,
            155768905756044823091307366842368
          ]
        },
        "leaf_counts": {
          "event": 652,
          "regular": 1059,
          "unknown": 3
        }
      },
      "gyroid_to_schwarz_p": {
        "area_fractions": {
          "event": 0.05546874999999999,
          "regular": 0.93203125,
          "unknown": 0.012499999999999997
        },
        "event_area_refers_to_cells_not_the_critical_set": true,
        "exact_area_fractions": {
          "event": [
            172806129823112203414172947154209,
            3115378115120896461826147336847360
          ],
          "regular": [
            580725951771754612707571729243839,
            623075623024179292365229467369472
          ],
          "unknown": [
            9735556609752798718528935868489,
            778844528780224115456536834211840
          ]
        },
        "leaf_counts": {
          "event": 355,
          "regular": 775,
          "unknown": 80
        }
      },
      "schwarz_p_to_diamond": {
        "area_fractions": {
          "event": 0.056718750000000005,
          "regular": 0.93921875,
          "unknown": 0.0040625
        },
        "event_area_refers_to_cells_not_the_critical_set": true,
        "exact_area_fractions": {
          "event": [
            88350176233506677258426952821249,
            1557689057560448230913073668423680
          ],
          "regular": [
            1463010769530602231914921620152519,
            1557689057560448230913073668423680
          ],
          "unknown": [
            791013974542415217465636931239,
            194711132195056028864134208552960
          ]
        },
        "leaf_counts": {
          "event": 363,
          "regular": 704,
          "unknown": 26
        }
      }
    },
    "distinct_isotopy_classes_proved": false,
    "event_statuses": {
      "certified_level_resolved_event": 249,
      "certified_multiseed_level_event": 7,
      "certified_stationary_event_in_cell": 1063,
      "exact_degenerate_gyroid_schwarz_event": 1,
      "exact_stationary_branch_event_in_cell": 50,
      "unknown": 109
    },
    "new_event_certificates_replayed": 256,
    "regular_leaf_count": 2538
  },
  "maps": {
    "hamiltonian": {
      "baseline_statuses": {
        "evaluated": 7116,
        "existing_cavity_route_not_revised": 1,
        "fixed_diagram_evaluated": 5,
        "time_budget": 18
      },
      "by_family": {
        "hopf_to_solomon": {
          "evaluated": 1439,
          "fixed_diagram_evaluated": 1
        },
        "hopf_to_trefoil": {
          "evaluated": 960
        },
        "trefoil_to_cinquefoil": {
          "evaluated": 2338,
          "fixed_diagram_evaluated": 2
        },
        "unknot_to_solomon": {
          "evaluated": 1437,
          "existing_cavity_route_not_revised": 1,
          "fixed_diagram_evaluated": 2
        },
        "unknot_to_trefoil": {
          "evaluated": 960
        }
      },
      "complete_exact_grid": true,
      "final_statuses": {
        "evaluated": 7134,
        "existing_cavity_route_not_revised": 1,
        "fixed_diagram_evaluated": 5
      },
      "original_record_sha256": "ec3ae06a1fa7341d76a89d4d74f9f79348cdbc3ac25557a06e33c792e8610d1a",
      "retried_cells": 18,
      "retry_commit": "b4c478cb8f4a8e79af472e4bb3e2c21015d600ec"
    },
    "tpms": {
      "baseline_statuses": {
        "evaluated": 1001,
        "existing_cavity_route_not_revised": 71,
        "fixed_diagram_evaluated": 73,
        "time_budget": 178
      },
      "by_family": {
        "gyroid_to_diamond": {
          "evaluated": 407,
          "fixed_diagram_evaluated": 34
        },
        "gyroid_to_schwarz_p": {
          "evaluated": 441
        },
        "schwarz_p_to_diamond": {
          "evaluated": 320,
          "existing_cavity_route_not_revised": 71,
          "fixed_diagram_evaluated": 43,
          "time_budget": 7
        }
      },
      "complete_exact_grid": true,
      "final_statuses": {
        "evaluated": 1168,
        "existing_cavity_route_not_revised": 71,
        "fixed_diagram_evaluated": 77,
        "time_budget": 7
      },
      "original_record_sha256": "4c54591edf12b746af218509b6fd6275c16f6127086310edcde8135fc75c3d5b",
      "retried_cells": 178,
      "retry_commit": "b4c478cb8f4a8e79af472e4bb3e2c21015d600ec"
    }
  },
  "multiseed_completion": {
    "additional_exact_events": 1,
    "additional_local_events": 7,
    "freshly_replayed": true,
    "input_unknown": 117,
    "remaining_unknown": 109,
    "search_commit": "d86677ede133dba0ed038f93b51a5e4649c6506a"
  },
  "original_sources_modified": false,
  "source_commit": "3abd27718323c0800231303582be2f6d05e5902a",
  "voxel_links": {
    "hamiltonian": {
      "by_family": {
        "hopf_to_solomon": {
          "cases": 1440,
          "nonmanifold_cases": 538
        },
        "hopf_to_trefoil": {
          "cases": 960,
          "nonmanifold_cases": 434
        },
        "trefoil_to_cinquefoil": {
          "cases": 2340,
          "nonmanifold_cases": 871
        },
        "unknot_to_solomon": {
          "cases": 1440,
          "nonmanifold_cases": 534
        },
        "unknot_to_trefoil": {
          "cases": 960,
          "nonmanifold_cases": 280
        }
      },
      "cases": 7140,
      "nonmanifold_cases": 2657,
      "pl_manifold_cases": 4483,
      "seconds": 167.79884034399998
    },
    "tpms": {
      "by_family": {
        "gyroid_to_diamond": {
          "cases": 441,
          "nonmanifold_cases": 441
        },
        "gyroid_to_schwarz_p": {
          "cases": 441,
          "nonmanifold_cases": 435
        },
        "schwarz_p_to_diamond": {
          "cases": 441,
          "nonmanifold_cases": 411
        }
      },
      "cases": 1323,
      "nonmanifold_cases": 1287,
      "pl_manifold_cases": 36,
      "seconds": 31.42297706300002
    }
  }
}
