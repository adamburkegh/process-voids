'''
process_voids.bpmn_colour: paints a void metric's per-activity values
onto a BPMN diagram's tasks.
'''

import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from skipalignments import Activity, Sequence, Tau

from process_voids import bpmn_colour
from process_voids.bpmn_colour import NAVY, RED, TEAL, qn

BPMN = '''<?xml version="1.0" encoding="UTF-8"?>
<bpmn:definitions xmlns:bpmn="http://www.omg.org/spec/BPMN/20100524/MODEL"
                  xmlns:bpmndi="http://www.omg.org/spec/BPMN/20100524/DI"
                  xmlns:dc="http://www.omg.org/spec/DD/20100524/DC">
  <bpmn:process id="P">
    <bpmn:task id="T1" name="o"/>
    <bpmn:task id="T2" name="s"/>
    <bpmn:task id="T3" name="unmodelled"/>
  </bpmn:process>
  <bpmndi:BPMNDiagram id="D">
    <bpmndi:BPMNPlane id="PL" bpmnElement="P">
      <bpmndi:BPMNShape id="S1" bpmnElement="T1"><dc:Bounds x="0" y="0" width="1" height="1"/></bpmndi:BPMNShape>
      <bpmndi:BPMNShape id="S2" bpmnElement="T2"><dc:Bounds x="0" y="0" width="1" height="1"/></bpmndi:BPMNShape>
      <bpmndi:BPMNShape id="S3" bpmnElement="T3"><dc:Bounds x="0" y="0" width="1" height="1"/></bpmndi:BPMNShape>
    </bpmndi:BPMNPlane>
  </bpmndi:BPMNDiagram>
</bpmn:definitions>
'''


class ValuesByLabelTest(unittest.TestCase):

    def test_activities_keyed_by_label_averaging_duplicates(self):
        a1, a2, t = Activity(None, 'a', 1), Activity(None, 'a', 1), Tau(None, 'skip', 0)
        seq = Sequence(None, [a1, a2, t])
        values = {seq: 0.9, a1: 0.2, a2: 0.4, t: 1.0}
        by_label = bpmn_colour.values_by_label(values)
        self.assertEqual(set(by_label), {'a'})
        self.assertAlmostEqual(by_label['a'], 0.3)


class ColourBpmnTest(unittest.TestCase):

    def test_matched_tasks_are_filled_by_value_band(self):
        with tempfile.TemporaryDirectory() as tmp:
            bpmn_in, bpmn_out = Path(tmp, 'in.bpmn'), Path(tmp, 'out.bpmn')
            bpmn_in.write_text(BPMN, encoding='utf-8')
            bpmn_colour.colour_bpmn(str(bpmn_in), str(bpmn_out), {'o': 0.0, 's': 0.4})
            shapes = {s.get('bpmnElement'): s.get(qn('bioc', 'fill'))
                      for s in ET.parse(bpmn_out).getroot().iter(qn('bpmndi', 'BPMNShape'))}
        self.assertEqual(shapes, {'T1': NAVY, 'T2': TEAL, 'T3': None})

    def test_a_wholly_void_task_is_red(self):
        self.assertEqual(bpmn_colour.value_to_colour(1.0), RED)


class CliTest(unittest.TestCase):

    def test_takes_log_model_and_two_diagrams_with_a_metric(self):
        args = bpmn_colour.parse_args(['l.xes', 'm.ptml', 'in.bpmn', 'out.bpmn',
                                       '--metric', 'voidsat'])
        self.assertEqual((args.log, args.model, args.bpmn_in, args.bpmn_out, args.metric),
                         ('l.xes', 'm.ptml', 'in.bpmn', 'out.bpmn', 'voidsat'))

    def test_default_metric_is_voidsalign(self):
        self.assertEqual(bpmn_colour.parse_args(['l', 'm', 'i', 'o']).metric, 'voidsalign')


if __name__ == '__main__':
    unittest.main()
