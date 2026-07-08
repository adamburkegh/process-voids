#Usage: python bpmn_colour.py model.bpmn model_coloured.bpmn --var-path var

import argparse
import pickle
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, Optional

from processtree import ProcessTree, Activity, Tau


# BPMN, BPMNDI, DC, DI, BIOC Namespaces
NS = {
    "bpmn": "http://www.omg.org/spec/BPMN/20100524/MODEL",
    "bpmndi": "http://www.omg.org/spec/BPMN/20100524/DI",
    "dc": "http://www.omg.org/spec/DD/20100524/DC",
    "di": "http://www.omg.org/spec/DD/20100524/DI",
    "bioc": "http://bpmn.io/schema/bpmn/biocolor/1.0",
}

for _prefix, _uri in NS.items():
    ET.register_namespace(_prefix, _uri)


def qn(prefix: str, tag: str) -> str:
    #Builds namespace tag for ElementTree
    return f"{{{NS[prefix]}}}{tag}"


# Loading the pvoid.py - DerivationnPipeline Output
def load_pickle(path: str, name: str):
    file_path = Path(path) / name
    with open(file_path, "rb") as f:
        return pickle.load(f)


def load_skip_probs(path: str) -> Dict[ProcessTree, float]:

    # Loads the `skip_probs` dict (ProcessTree node -> P(skip)) pickled by DerivationPipeline.compute() to `<path>/skip_probs`.

    return load_pickle(path, "skip_probs")


def load_tree(path: str) -> ProcessTree:
    #Loads the process tree pickled by DerivationPipeline.compute() to "path/tree"
    return load_pickle(path, "tree")


def activity_skip_probs_by_label(skip_probs: Dict[ProcessTree, float]) -> Dict[str, float]:

    # Maps activity *labels* to their skip probability, so they can be matched against the `name`/`id` of BPMN task elements.

    sums: Dict[str, float] = {}
    counts: Dict[str, int] = {}
    for node, prob in skip_probs.items():
        if not isinstance(node, Activity):
            continue
        label = node.name
        sums[label] = sums.get(label, 0.0) + prob
        counts[label] = counts.get(label, 0) + 1
    return {label: sums[label] / counts[label] for label in sums}


# Colour Codes
GREEN = "#00A651"
ORANGE = "#F7941E"
RED = "#ED1C24"
MAROON = "#800000"


def prob_to_colour(p: float) -> str:
    """
    Skip probability- four bands to return a #RRGGBB" hex string:
        0%    <= p <= 14% -> green
        15%   <= p <= 20% -> orange
        21%   <= p <= 30% -> red
        31%   <= p <= 100% -> maroon
    """
    p = max(0.0, min(1.0, p))
    if p <= 0.14:
        return GREEN
    if p <= 0.20:
        return ORANGE
    if p <= 0.30:
        return RED
    return MAROON


# BPMN colouring
TASK_TAGS = [
    "task", "userTask", "serviceTask", "manualTask", "scriptTask",
    "sendTask", "receiveTask", "businessRuleTask", "callActivity",
]


def find_task_elements(root: ET.Element) -> Dict[str, ET.Element]:
    #Returns a dict id -> element for every task like element in every <bpmn:process>
    tasks = {}
    for process in root.iter(qn("bpmn", "process")):
        for tag in TASK_TAGS:
            for el in process.iter(qn("bpmn", tag)):
                el_id = el.get("id")
                if el_id is not None:
                    tasks[el_id] = el
    return tasks


def find_shape_for_element(root: ET.Element, element_id: str) -> Optional[ET.Element]:
    for shape in root.iter(qn("bpmndi", "BPMNShape")):
        if shape.get("bpmnElement") == element_id:
            return shape
    return None


def colour_bpmn(
    in_path: str,
    out_path: str,
    label_to_prob: Dict[str, float],
    match_by: str = "name",
) -> None:

    # Reads the .bpmn file, colours every task whose name (or id, if match_by="id") matches a key in `label_to_prob`, and writes the result to out_path

    tree_xml = ET.parse(in_path)
    root = tree_xml.getroot()

    tasks = find_task_elements(root)
    coloured = 0
    unmatched = []

    for task_id, task_el in tasks.items():
        key = task_el.get("name") if match_by == "name" else task_id
        if key not in label_to_prob:
            unmatched.append(key)
            continue

        colour = prob_to_colour(label_to_prob[key])
        shape = find_shape_for_element(root, task_id)
        if shape is None:
            continue

        shape.set(qn("bioc", "fill"), colour)
        shape.set(qn("bioc", "stroke"), "#000000")
        coloured += 1

    tree_xml.write(out_path, xml_declaration=True, encoding="UTF-8")

    print(f"Coloured {coloured} task(s).")
    if unmatched:
        print(f"{len(unmatched)} task(s) had no matching skip probability: {unmatched}")


def main():
    parser = argparse.ArgumentParser(
        description="Colour BPMN task activities by their skip probability "
                    "(as computed by pvoid.py / DerivationPipeline)."
    )
    parser.add_argument("bpmn_in", help="Path to the input .bpmn file")
    parser.add_argument("bpmn_out", help="Path to write the coloured .bpmn file to")
    parser.add_argument(
        "--var-path", default="var",
        help="Directory in which pvoid.py's DerivationPipeline.compute(path=...) "
             "wrote 'tree' and 'skip_probs' (default: 'var')",
    )
    parser.add_argument(
        "--match-by", choices=["name", "id"], default="name",
        help="Match BPMN tasks to activities by their visible 'name' or by "
             "their 'id' (default: name)",
    )
    args = parser.parse_args()

    skip_probs = load_skip_probs(args.var_path)
    label_to_prob = activity_skip_probs_by_label(skip_probs)

    colour_bpmn(args.bpmn_in, args.bpmn_out, label_to_prob, match_by=args.match_by)


if __name__ == "__main__":
    main()