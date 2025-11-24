import os, sys
import pytest
from WatCon import generate_static_networks


def test_waters_uneven():
  """Test simple box of waters"""
  metrics, networks, _, _ = generate_static_networks.initialize_network('water_dir',
                            network_type='water-water',
                            include_hydrogens=False,
                            custom_selection=None,
                            active_region_reference=None,
                            active_region_COM=False,
                            active_region_only=False,
                            active_region_radius=8.0,
                            water_name=None,
                            multi_model_pdb=False,
                            max_distance=3.8,
                            angle_criteria=None,
                            analysis_conditions='all',
                            analysis_selection='all',
                            project_networks=False,
                            return_network=True,
                            cluster_coordinates=False,
                            clustering_method='hdbscan',
                            cluster_water_only=True,
                            min_cluster_samples=15,
                            eps=None,
                            msa_indexing=False,
                            alignment_file=None,
                            combined_fasta=None,
                            fasta_directory=None,
                            classify_water=False,
                            classification_file_base='STATIC',
                            MSA_reference_pdb=None,
                            water_reference_resids=None,
                            num_workers=1,
                            shortest_path_nodes=None,
                            max_neighbors=10)
  assert len(networks[0].water_molecules) == 8
  assert len(networks[0].connections) == 6 


  metrics, networks, _, _ = generate_static_networks.initialize_network('water_dir',
                            network_type='water-water',
                            include_hydrogens=False,
                            custom_selection=None,
                            active_region_reference=None,
                            active_region_COM=False,
                            active_region_only=False,
                            active_region_radius=8.0,
                            water_name=None,
                            multi_model_pdb=False,
                            max_distance=5.5,
                            angle_criteria=None,
                            analysis_conditions='all',
                            analysis_selection='all',
                            project_networks=False,
                            return_network=True,
                            cluster_coordinates=False,
                            clustering_method='hdbscan',
                            cluster_water_only=True,
                            min_cluster_samples=15,
                            eps=None,
                            msa_indexing=False,
                            alignment_file=None,
                            combined_fasta=None,
                            fasta_directory=None,
                            classify_water=False,
                            classification_file_base='STATIC',
                            MSA_reference_pdb=None,
                            water_reference_resids=None,
                            num_workers=1,
                            shortest_path_nodes=None,
                            max_neighbors=10)
  assert len(networks[0].connections) == 9


