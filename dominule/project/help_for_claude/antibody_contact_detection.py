"""Contact detection module using Biotite for accurate protein-protein interface analysis."""

from typing import Iterable, List, Optional, Tuple

import numpy as np
import biotite.structure.io as strucio


def detect_interface_residues_all_atoms(
    file_path: str,
    chain1_id: str,
    chain2_id: str,
    distance_threshold: float = 4.5,
) -> List[Tuple[Tuple[str, int], Tuple[str, int]]]:
    """
    Detect protein-protein interface residue pairs between two chains,
    considering all heavy atoms (non-hydrogen atoms).

    This function identifies residue pairs at protein-protein interfaces by finding
    all atom pairs within a specified distance threshold.

    Parameters:
    - file_path (str): Path to the structure file (PDB/CIF).
    - chain1_id (str): Chain ID for the first protein.
    - chain2_id (str): Chain ID for the second protein.
    - distance_threshold (float): Maximum distance (in Ångströms) between any two atoms
      to be considered in contact. Typical values: 4.5-5.0 Å for interface detection.
      Default: 4.5 Å.

    Returns:
    - List of tuples containing interacting residue pairs with residue names and IDs.
      Each tuple is formatted as: ((res1_name, res1_id), (res2_name, res2_id))
    """
    # Load the structure file (supports PDB, CIF formats via Biotite)
    structure = strucio.load_structure(file_path)

    # Filter structure to include only heavy atoms (non-hydrogen atoms)
    # Hydrogen atoms are excluded because:
    # 1. They're often not resolved in experimental structures
    # 2. Their positions are model-dependent
    # 3. Heavy atoms provide sufficient representation of molecular contacts
    is_heavy_atom = ~np.char.startswith(structure.atom_name.astype(str), "H")
    structure = structure[is_heavy_atom]

    # Filter by chain ID to isolate the two protein chains of interest
    chain1 = structure[structure.chain_id == chain1_id]
    chain2 = structure[structure.chain_id == chain2_id]

    # Ensure both chains have residues
    if len(chain1) == 0 or len(chain2) == 0:
        raise ValueError(
            f"One of the chains ({chain1_id}, {chain2_id}) has no residues in the structure."
        )

    # Compute pairwise distances between all heavy atoms in chain1 and chain2
    # This creates a distance matrix of shape (n_atoms_chain1, n_atoms_chain2)
    # Using broadcasting to compute all pairwise differences efficiently
    diff = chain1.coord[:, np.newaxis, :] - chain2.coord[np.newaxis, :, :]
    distmat = np.linalg.norm(diff, axis=2)  # Euclidean distance in 3D space

    # Apply distance threshold filter
    # This creates a boolean mask for atom pairs within the contact distance cutoff
    pair_mask = distmat <= distance_threshold

    # Extract residue pairs from atom pairs that pass both filters
    # Use a set to automatically eliminate duplicate residue pairs
    # (multiple atom-atom contacts between the same residue pair count as one contact)
    interface_pairs = set()
    for i, j in np.argwhere(pair_mask):  # Get indices of True values in contact matrix
        # Extract residue information for the contacting atoms
        res1 = (
            chain1.res_name[i],
            int(chain1.res_id[i]),
        )  # (residue name, residue ID) for chain1
        res2 = (
            chain2.res_name[j],
            int(chain2.res_id[j]),
        )  # (residue name, residue ID) for chain2
        interface_pairs.add((res1, res2))

    return list(interface_pairs)


def count_target_residue_contacts(
    file_path: str,
    antigen_chain: str,
    antibody_chain: str,
    target_residues: Iterable[int],
    distance_threshold: float = 4.5,
    antibody_cdr_residues: Optional[Iterable[int]] = None,
):
    """
    Count how many target residues on the antigen chain have contacts with the antibody chain.

    This function is useful for epitope analysis, determining which specific antigen residues
    are engaged by the antibody. It uses the same contact detection as
    detect_interface_residues_all_atoms.

    Parameters:
    - file_path (str): Path to the structure file (PDB/CIF).
    - antigen_chain (str): Chain ID containing target residues.
    - antibody_chain (str): Chain ID for the antibody.
    - target_residues (Iterable[int]): List/set of target residue IDs to check for contacts.
    - distance_threshold (float): Maximum distance (Ångströms) for atom-atom contact consideration.
      Typical values: 4.5-5.0 Å. Default: 4.5 Å.
    - antibody_cdr_residues (Optional[Iterable[int]]): If provided, only contacts involving these
      antibody residues will be counted (e.g., CDR2+CDR3 only). If None, all antibody residues are considered.

    Returns:
    - dict with the following keys:
        - 'contact_count': Number of target residues that have at least one contact
        - 'target_residues_in_contact': Sorted list of target residue IDs with contacts
        - 'all_antigen_contacts': Sorted list of all antigen residue IDs in contact (not just targets)
        - 'total_interface_pairs': Total number of unique residue-residue contacts at the interface
        - 'target_contact_pairs': List of (antigen_res_id, antibody_res_id) pairs involving targets
    """
    # Get all interface pairs using the core contact detection function
    # This identifies all residue-residue contacts between antigen and antibody
    interface_pairs = detect_interface_residues_all_atoms(
        file_path=file_path,
        chain1_id=antigen_chain,
        chain2_id=antibody_chain,
        distance_threshold=distance_threshold,
    )

    # Convert antibody_cdr_residues to a set for fast lookup if provided
    antibody_cdr_set = (
        set(int(res_id) for res_id in antibody_cdr_residues)
        if antibody_cdr_residues
        else None
    )

    # Process interface pairs to extract target residue contact information
    antigen_contacts = set()  # All antigen residues with any contact
    target_contact_pairs = []  # Specific (antigen, antibody) pairs involving targets
    target_residue_set = set(
        int(res_id) for res_id in target_residues
    )  # Convert to int set for fast lookup

    # Iterate through all detected interface pairs
    for (res1_name, res1_id), (res2_name, res2_id) in interface_pairs:
        # If antibody_cdr_residues is specified, only count contacts involving those residues
        if antibody_cdr_set is not None and res2_id not in antibody_cdr_set:
            continue

        # Track all antigen residues at the interface
        antigen_contacts.add(int(res1_id))

        # Specifically collect pairs involving our target residues of interest
        # This allows detailed analysis of which antibody residues contact each target
        if res1_id in target_residue_set:
            target_contact_pairs.append((int(res1_id), int(res2_id)))

    # Filter the target residues to only those that actually have contacts
    # This is the intersection of our target list and the residues with detected contacts
    target_contacts = [
        int(res_id) for res_id in target_residue_set if res_id in antigen_contacts
    ]

    return {
        "contact_count": len(target_contacts),
        "target_residues_in_contact": sorted(target_contacts),
        "all_antigen_contacts": sorted(antigen_contacts),
        "total_interface_pairs": len(interface_pairs),
        "target_contact_pairs": sorted(set(target_contact_pairs)),
    }


def detect_interface_residues_with_distance(
        file_path: str,
        chain1_id: str,
        chain2_id: str,
        distance_threshold: float = 4.5,
        antibody_cdr_residues: Optional[Iterable[int]] = None,
) -> List[Tuple[Tuple[str, int], Tuple[str, int], float]]:
    structure = strucio.load_structure(file_path)

    # Filter for heavy atoms
    is_heavy_atom = ~np.char.startswith(structure.atom_name.astype(str), "H")
    structure = structure[is_heavy_atom]

    chain1 = structure[structure.chain_id == chain1_id]
    chain2 = structure[structure.chain_id == chain2_id]

    if len(chain1) == 0 or len(chain2) == 0:
        raise ValueError(f"One of the chains ({chain1_id}, {chain2_id}) has no residues.")

    # --- REPAIR: Prepare the CDR filter set ---
    antibody_cdr_set = (
        set(int(res_id) for res_id in antibody_cdr_residues)
        if antibody_cdr_residues is not None
        else None
    )

    # Compute pairwise distances
    diff = chain1.coord[:, np.newaxis, :] - chain2.coord[np.newaxis, :, :]
    distmat = np.linalg.norm(diff, axis=2)

    indices = np.where(distmat <= distance_threshold)
    interface_dict = {}

    for i, j in zip(*indices):
        # --- REPAIR: Check the antibody residue ID (chain2) ---
        res2_id = int(chain2.res_id[j])

        if antibody_cdr_set is not None and res2_id not in antibody_cdr_set:
            continue

        res1 = (chain1.res_name[i], int(chain1.res_id[i]))
        res2 = (chain2.res_name[j], res2_id)
        dist = float(distmat[i, j])

        pair_key = (res1, res2)

        if pair_key not in interface_dict or dist < interface_dict[pair_key]:
            interface_dict[pair_key] = dist

    return [(k[0], k[1], v) for k, v in interface_dict.items()]


def get_cdr_residue_ids(antibody_sequence: str, cdr_sequence: str, offset: int = 1) -> list[int]:
    """
    Finds the 1-based residue IDs of a CDR sequence within a full antibody sequence.
    """
    if type(cdr_sequence) != str:
        print(cdr_sequence)
        print(type(cdr_sequence))
        raise TypeError(
            "CDR sequence must be a string."
        )
    start_idx = antibody_sequence.find(cdr_sequence)

    if start_idx == -1:
        raise ValueError(f"CDR sequence '{cdr_sequence}' not found in the full sequence.")

    # Calculate residue IDs (start_idx is 0-based, so add offset)
    first_res = start_idx + offset
    last_res = first_res + len(cdr_sequence)

    return list(range(first_res, last_res))

# Example Usage:
# antibody = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDGYSYGWYFDYWGQGTLVTVSS"
# cdr3 = "AKDGYSYGWYFDY"
# res_ids = get_cdr_residue_ids(antibody, cdr3)
# Result: [98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110]