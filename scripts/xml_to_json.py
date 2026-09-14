#!/usr/bin/env python3
"""One-off converter: lactoferrin.xml (old tinyxml config format) -> lactoferrin.json (new schema).

This is a throwaway migration script, not part of the shipped project. It emits JSON-with-comments
(JSONC: // and /* */ comments), which the new C++ loader (nlohmann::json with ignore_comments=true)
accepts directly.
"""
import sys
import json
import re
import xml.etree.ElementTree as ET


def get_value(parent, tag, default=None):
    el = parent.find(tag) if parent is not None else None
    if el is None:
        return default
    return el.attrib.get("value", default)


def num(s):
    """Normalize a numeric literal (handles quirky source values like '0.') into valid JSON."""
    return json.dumps(float(s))


def convert(xml_path, json_path):
    tree = ET.parse(xml_path)
    root = tree.getroot()
    assert root.tag == "configuration"

    lines = []
    lines.append("// Configuration for the protein-cutting simulator (JSON format, replaces the old XML format).")
    lines.append("// Comments (// and /* */) are supported by the loader and stripped before parsing.")
    lines.append("{")

    # ---------------- parameters ----------------
    params = root.find("parameters")
    lines.append("\t// === Simulation parameters ===")
    lines.append("\t\"parameters\": {")
    lines.append("\t\t// set to a number to get reproducible runs; null = seed from current time")
    lines.append("\t\t\"randomSeed\": null,")
    lines.append("\t\t// stop condition: stop after maximum number of iterations reached")
    lines.append(f'\t\t"maxTime": {get_value(params, "maxTime", "1200000")},')
    lines.append("\t\t// stop condition: stop after reaching the given degree of hydrolysis (DH)")
    lines.append(f'\t\t"maxDH": {num(get_value(params, "maxDH", "0.10"))},')
    lines.append("\t\t// how many attempts to cut should the simulated pepsin perform, per iteration")
    lines.append(f'\t\t"maxAttemptsPerTime": {get_value(params, "maxAttemptsPerTime", "1")},')
    lines.append("\t\t// stop condition: stop after failing this many consecutive attempts to cut")
    lines.append(f'\t\t"maxAttempts": {get_value(params, "maxAttempts", "1000")},')
    lines.append("")
    lines.append("\t\t// THE FOLLOWING THREE ARE EXPERIMENTAL, ADVISE NOT TO MODIFY THEM")
    lines.append("\t\t// pepsin 'dies out' during the process, so it is more effective at the beginning; we thought it could be")
    lines.append("\t\t// interesting to have pepsin reduce its probability of activation during time, but the code is not ready")
    lines.append(f'\t\t"initialPepsin": {num(get_value(params, "initialPepsin", "1.0"))},')

    always_dying = get_value(params, "pepsinAlwaysDying", "false")
    always_dying_json = "true" if always_dying.strip().lower() in ("true", "1") else "false"
    lines.append(f'\t\t"pepsinAlwaysDying": {always_dying_json},')
    lines.append(f'\t\t"pepsinDyingRatio": {num(get_value(params, "pepsinDyingRatio", "1.0"))}')
    lines.append("\t},")
    lines.append("")

    # ---------------- proteins ----------------
    lines.append('\t// here you can define multiple proteins, each with a number of copies ("quantity").')
    lines.append('\t// "disulfideBonds" lists 1-indexed positions in the sequence that pepsin finds harder to cut')
    lines.append("\t// (not all of them are disulfide bonds, some are glycosylations)")
    lines.append('\t"proteins": [')

    protein_blocks = []
    for protein in root.find("proteins").findall("protein"):
        name = protein.attrib.get("name", "protein")
        quantity = protein.attrib.get("quantity", "1")
        bonds_str = protein.attrib.get("disulfideBonds", "")
        bonds = [int(b) for b in bonds_str.split()] if bonds_str.strip() else []
        seq = re.sub(r"\s+", "", protein.text or "")

        block = []
        block.append("\t\t{")
        block.append(f'\t\t\t"name": {json.dumps(name)},')
        block.append(f'\t\t\t"quantity": {quantity},')
        block.append(f'\t\t\t"disulfideBonds": {json.dumps(bonds)},')
        block.append(f'\t\t\t"sequence": {json.dumps(seq)}')
        block.append("\t\t}")
        protein_blocks.append("\n".join(block))

    lines.append(",\n".join(protein_blocks))
    lines.append("\t],")
    lines.append("")

    # ---------------- cuts ----------------
    lines.append("\t// basic probabilities of cutting a bond, depending on the aminoacids in position P1 (left) / P1' (right);")
    lines.append("\t// probabilities here are taken from Hamuro et al., 2008.")
    lines.append("\t// NOTE: entries with probability 0 are omitted (a missing right-aminoacid means probability 0),")
    lines.append("\t// matching the original loader, which only stored cuts with probability > 0.")
    lines.append('\t"cuts": {')

    cut_groups = []  # [(left, [(right, prob_str), ...])]
    skipped_zero_cuts = 0
    for cut in root.find("cuts").findall("cut"):
        left = cut.attrib["left"]
        right = cut.attrib["right"]
        prob = cut.attrib.get("probability", "0.0")
        if float(prob) <= 0.0:
            skipped_zero_cuts += 1
            continue
        if not cut_groups or cut_groups[-1][0] != left:
            cut_groups.append((left, []))
        cut_groups[-1][1].append((right, prob))

    cut_lines = []
    for left, entries in cut_groups:
        inner = ", ".join(f'"{r}": {num(p)}' for r, p in entries)
        cut_lines.append(f'\t\t"{left}": {{ {inner} }}')
    lines.append(",\n".join(cut_lines))
    lines.append("\t},")
    lines.append("")

    # ---------------- alterations ----------------
    lines.append("\t// position-dependent adjustments to the base cut probability, for aminoacids found in")
    lines.append("\t// positions P2-P4 (left) / P2'-P4' (right); data taken from Powers et al., 1977")
    lines.append('\t"alterations": {')

    alterations_root = root.find("alterations")
    alt_groups = []  # [(aminoacid, {"left": [(pos, prob)], "right": [(pos, prob)]})]
    for alt in alterations_root.findall("alteration"):
        aminoacid = alt.attrib["aminoacid"]
        side = alt.attrib["side"]
        position = alt.attrib["position"]
        prob = alt.attrib.get("probability", "0.0")
        if not alt_groups or alt_groups[-1][0] != aminoacid:
            alt_groups.append((aminoacid, {"left": [], "right": []}))
        alt_groups[-1][1][side].append((position, prob))

    alt_lines = []
    for aminoacid, sides in alt_groups:
        left_inner = ", ".join(f'"{p}": {num(v)}' for p, v in sides["left"])
        right_inner = ", ".join(f'"{p}": {num(v)}' for p, v in sides["right"])
        alt_lines.append(f'\t\t"{aminoacid}": {{ "left": {{ {left_inner} }}, "right": {{ {right_inner} }} }}')
    lines.append(",\n".join(alt_lines))
    lines.append("\t},")
    lines.append("")

    # ---------------- terminalAlterations ----------------
    lines.append("\t// multipliers applied to the cut probability when close to either end of a peptide chain;")
    lines.append("\t// data taken from Powers et al., 1977")
    lines.append('\t"terminalAlterations": {')

    term_left, term_right = [], []
    for term in alterations_root.findall("terminal"):
        side = term.attrib["side"]
        position = term.attrib["position"]
        multiplier = term.attrib.get("multiplier", "1.0")
        (term_left if side == "left" else term_right).append((position, multiplier))

    left_inner = ", ".join(f'"{p}": {num(v)}' for p, v in term_left)
    right_inner = ", ".join(f'"{p}": {num(v)}' for p, v in term_right)
    lines.append(f'\t\t"left": {{ {left_inner} }},')
    lines.append(f'\t\t"right": {{ {right_inner} }}')
    lines.append("\t}")

    lines.append("}")
    lines.append("")

    text = "\n".join(lines)

    with open(json_path, "w", encoding="utf-8") as f:
        f.write(text)

    # sanity check: strip comments and verify it parses as valid JSON
    stripped = re.sub(r"//[^\n]*", "", text)
    parsed = json.loads(stripped)
    print(f"Wrote {json_path} ({len(text)} bytes). Sanity-parsed OK: "
          f"{len(parsed['proteins'])} protein(s), {len(parsed['cuts'])} cut groups "
          f"({skipped_zero_cuts} zero-probability cuts omitted), "
          f"{len(parsed['alterations'])} alteration groups.")


if __name__ == "__main__":
    convert(sys.argv[1], sys.argv[2])
