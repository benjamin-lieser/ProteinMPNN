import os

NAMES = glob_wildcards("data/proteinmpnn_pdbs/{pdb}.pdb").pdb

rule all:
    input:
        "data/scores_all.json",
        "data/scores_cdr3.json",
        "data/scores_justcdr3.json",
rule owndir:
    input:
        "data/proteinmpnn_pdbs/{pdb}.pdb"
    output:
        "data/{pdb}/pdb.pdb"
    shell:
        "cp {input} {output}"

rule pdb2json:
    input:
        "data/{pdb}/pdb.pdb"
    output:
        "data/json/{pdb}.json"
    shell:
        "python helper_scripts/parse_multiple_chains.py --input_path=data/{wildcards.pdb} --output_path={output}"

def get_chain_list(pdb: str):
    if len(pdb) == 4: # 8ee2, the others are longer
        return 'A C'
    else:
        return 'A B'

def range_list(start: int, end: int):
    return [str(i) for i in range(start, end+1)]

def range3_list(r1, r2, r3):
    return range_list(r1[0], r1[1]) + range_list(r2[0], r2[1]) + range_list(r3[0], r3[1])

def range2str(r1, r2, r3):
    r = range3_list(r1, r2, r3)
    l = " ".join(r)
    return ", " + l

def get_fixed_pos(pdb: str):
    if pdb.startswith('8ee2'):
        return range2str((25, 31), (51, 55), (97, 108))
    elif pdb.startswith('7olz'):
        return range2str((23, 35), (50, 64), (99, 116))
    elif pdb.startswith('8q7s'):
        return range2str((22, 33), (45, 63), (97, 117))
    elif pdb.startswith('8q93'):
        return range2str((23, 35), (51, 64), (99, 119))
    else:
        raise ValueError(f"Unknown pdb: {pdb}")

def get_fixed_pos_just_cdr3(pdb: str):
    if pdb.startswith('8ee2'):
        return range2str((25, 0), (51, 0), (97, 108))
    elif pdb.startswith('7olz'):
        return range2str((23, 0), (50, 0), (99, 116))
    elif pdb.startswith('8q7s'):
        return range2str((22, 0), (45, 0), (97, 117))
    elif pdb.startswith('8q93'):
        return range2str((23, 0), (51, 0), (99, 119))
    else:
        raise ValueError(f"Unknown pdb: {pdb}")
def get_fixed_pos_without_cdr3(pdb: str):
    if pdb.startswith('8ee2'):
        return range2str((25, 31), (51, 55), (97, 0))
    elif pdb.startswith('7olz'):
        return range2str((23, 35), (50, 64), (99, 0))
    elif pdb.startswith('8q7s'):
        return range2str((22, 33), (45, 63), (97, 0))
    elif pdb.startswith('8q93'):
        return range2str((23, 35), (51, 64), (99, 0))
    else:
        raise ValueError(f"Unknown pdb: {pdb}")

rule fixed_pos:
    input:
        "data/json/{pdb}.json"
    output:
        "data/fixed_pos_all/{pdb}.json"
    run:
        chain = get_chain_list(wildcards.pdb)
        fixed_pos = get_fixed_pos(wildcards.pdb)
        shell(f"python helper_scripts/make_fixed_positions_dict.py --input_path={{input}} --output_path={{output}} --chain_list '{chain}' --position_list '{fixed_pos}' --specify_non_fixed")

rule fixed_pos_just_cdr3:
    input:
        "data/json/{pdb}.json"
    output:
        "data/fixed_pos_justcdr3/{pdb}.json"
    run:
        chain = get_chain_list(wildcards.pdb)
        fixed_pos = get_fixed_pos_just_cdr3(wildcards.pdb)
        shell(f"python helper_scripts/make_fixed_positions_dict.py --input_path={{input}} --output_path={{output}} --chain_list '{chain}' --position_list '{fixed_pos}' --specify_non_fixed")

rule fixed_pos_cdr3:
    input:
        "data/json/{pdb}.json"
    output:
        "data/fixed_pos_cdr3/{pdb}.json"
    run:
        chain = get_chain_list(wildcards.pdb)
        fixed_pos = get_fixed_pos_without_cdr3(wildcards.pdb)
        shell(f"python helper_scripts/make_fixed_positions_dict.py --input_path={{input}} --output_path={{output}} --chain_list '{chain}' --position_list '{fixed_pos}' --specify_non_fixed")

rule run_mpnn:
    input:
        fix_pos = "data/fixed_pos{cdr3}/{pdb}.json",
        pdb = "data/json/{pdb}.json"
    output:
        "data/mpnn{cdr3}/{pdb}/score_only/pdb_pdb.npz"
    shell:
        "python protein_mpnn_run.py --jsonl_path={input.pdb} --fixed_positions_jsonl={input.fix_pos} --out_folder data/mpnn{wildcards.cdr3}/{wildcards.pdb} --seed 42 --score_only 1"

rule run_mpnn_prob:
    input:
        fix_pos = "data/fixed_pos{cdr3}/{pdb}.json",
        pdb = "data/json/{pdb}.json"
    output:
        "data/mpnn{cdr3}/{pdb}/unconditional_probs_only/pdb.npz"
    shell:
        "python protein_mpnn_run.py --jsonl_path={input.pdb} --fixed_positions_jsonl={input.fix_pos} --out_folder data/mpnn{wildcards.cdr3}/{wildcards.pdb} --seed 42 --unconditional_probs_only 1"

rule run_mpnn_uncon_prob:
    input:
        pdb = "data/json/{pdb}.json"
    output:
        "data/mpnn/{pdb}/unconditional_probs_only/pdb.npz"
    shell:
        "python protein_mpnn_run.py --jsonl_path={input.pdb} --out_folder data/mpnn/{wildcards.pdb} --seed 42 --unconditional_probs_only 1"

rule extract_cdr:
    input:
        "data/mpnn/{pdb}/unconditional_probs_only/pdb.npz"
    output:
        "data/cdr/{pdb}.txt"
    run:
        import numpy as np
        data = np.load(input[0])
        S = data['S']
        mask = data['design_mask']
        cdr = S[mask == 1]

        alphabet = 'ACDEFGHIKLMNPQRSTVWYX'
        cdr = [alphabet[i] for i in cdr]
        cdr = ''.join(cdr)
        with open(output[0], 'w') as f:
            f.write(cdr)

rule extract_scores:
    input:
        expand("data/mpnn{{cdr3}}/{pdb}/score_only/pdb_pdb.npz", pdb=NAMES)
    output:
        "data/scores{cdr3}.json"
    run:
        import json
        import numpy as np
        scores = {}
        for pdb in NAMES:
            data = np.load(f"data/mpnn{wildcards.cdr3}/{pdb}/score_only/pdb_pdb.npz")
            scores[pdb] = data['score'].item()
        with open(output[0], 'w') as f:
            json.dump(scores, f)