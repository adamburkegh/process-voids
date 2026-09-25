#Usage: python -m process_voids.bpmn_colour log.xes model.ptml model.bpmn model_coloured.bpmn [--metric voidsalign|voidsat|voidmass_process]

import argparse
import xml.etree.ElementTree as ET
from typing import Dict, Optional

import pm4py_config as pm4py
from skipalignments import ProcessTree, Activity

from process_voids import pvoid
from process_voids.tree import from_pm4py


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


def values_by_label(values: Dict[ProcessTree, float]) -> Dict[str, float]:

    # Maps activity *labels* to their value (the mean, where a label appears more than once), so they can be matched against the `name`/`id` of BPMN task elements.

    sums: Dict[str, float] = {}
    counts: Dict[str, int] = {}
    for node, value in values.items():
        if not isinstance(node, Activity):
            continue
        label = node.name
        sums[label] = sums.get(label, 0.0) + value
        counts[label] = counts.get(label, 0) + 1
    return {label: sums[label] / counts[label] for label in sums}


# Colour Codes
RED = "#D72000FF"
ORANGE = "#EE6100FF"
AMBER = "#FFAD0AFF"
TEAL = "#1BB6AFFF"
GREY = "#9093A2FF"
NAVY = "#777F9FFF"


def value_to_colour(p: float) -> str:
    """
    A void metric's value - six even bands to return a #RRGGBBAA hex string:
        0%      <= p < 16.67%  -> NAVY
        16.67%  <= p < 33.33%  -> GREY
        33.33%  <= p < 50%     -> TEAL
        50%     <= p < 66.67%  -> AMBER
        66.67%  <= p < 83.33%  -> ORANGE
        83.33%  <= p <= 100%   -> RED
    """
    p = max(0.0, min(1.0, p))
    if p < 1 / 6:
        return NAVY
    if p < 2 / 6:
        return GREY
    if p < 3 / 6:
        return TEAL
    if p < 4 / 6:
        return AMBER
    if p < 5 / 6:
        return ORANGE
    return RED


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
    label_to_value: Dict[str, float],
    match_by: str = "name",
) -> None:

    # Reads the .bpmn file, colours every task whose name (or id, if match_by="id") matches a key in `label_to_value`, and writes the result to out_path

    tree_xml = ET.parse(in_path)
    root = tree_xml.getroot()

    tasks = find_task_elements(root)
    coloured = 0
    unmatched = []

    for task_id, task_el in tasks.items():
        key = task_el.get("name") if match_by == "name" else task_id
        if key not in label_to_value:
            unmatched.append(key)
            continue

        colour = value_to_colour(label_to_value[key])
        shape = find_shape_for_element(root, task_id)
        if shape is None:
            continue

        shape.set(qn("bioc", "fill"), colour)
        shape.set(qn("bioc", "stroke"), "#000000")
        coloured += 1

    tree_xml.write(out_path, xml_declaration=True, encoding="UTF-8")

    print(f"Coloured {coloured} task(s).")
    if unmatched:
        print(f"{len(unmatched)} task(s) had no matching value: {unmatched}")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="python -m process_voids.bpmn_colour",
        description="Colour BPMN task activities by a void metric."
    )
    parser.add_argument("log", help="XES event log")
    parser.add_argument("model", help="PTML process tree the metric is computed against")
    parser.add_argument("bpmn_in", help="Path to the input .bpmn file")
    parser.add_argument("bpmn_out", help="Path to write the coloured .bpmn file to")
    parser.add_argument(
        "--metric", choices=tuple(pvoid.METRICS), default="voidsalign",
        help="void metric to colour by (default: voidsalign)",
    )
    parser.add_argument(
        "--match-by", choices=["name", "id"], default="name",
        help="Match BPMN tasks to activities by their visible 'name' or by "
             "their 'id' (default: name)",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    log = pm4py.read_xes(args.log)
    tree = from_pm4py(pm4py.read_ptml(args.model))
    values = pvoid.METRICS[args.metric](log, tree)
    colour_bpmn(args.bpmn_in, args.bpmn_out, values_by_label(values), match_by=args.match_by)


if __name__ == "__main__":
    main()