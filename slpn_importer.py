
'''
Python translation of 

https://github.com/promworkbench/StochasticLabelledPetriNets/blob/main/src/org/processmining/stochasticlabelledpetrinets/plugins/StochasticLabelledPetriNetImportPlugin.java

ChatGPT assisted

'''


import io
from fractions import Fraction

from slpn import StochasticLabelledPetriNet


class StochasticLabelledPetriNetImporter:
    """
    A Python translation of the Java class
    org.processmining.stochasticlabelledpetrinets.plugins.StochasticLabelledPetriNetImporter.
    """

    @staticmethod
    def import_from_stream(stream: io.TextIOBase, filename: str = None, file_size: int = None):
        """Equivalent to importFromStream in Java."""
        return StochasticLabelledPetriNetImporter.read(stream)

    @staticmethod
    def read(stream: io.TextIOBase) -> StochasticLabelledPetriNet:
        """Parse an .slpn file from an input stream."""
        result = StochasticLabelledPetriNet()
        reader = stream if hasattr(stream, "readline") else io.TextIOWrapper(stream)

        # Read header (skip first line)
        StochasticLabelledPetriNetImporter.get_next_line(reader)

        # Places
        number_of_places = int(StochasticLabelledPetriNetImporter.get_next_line(reader))
        for place in range(number_of_places):
            result.addPlace()

            in_initial_marking = int(StochasticLabelledPetriNetImporter.get_next_line(reader))
            if in_initial_marking > 0:
                result.addPlaceToInitialMarking(place, in_initial_marking)

        # Transitions
        number_of_transitions = int(StochasticLabelledPetriNetImporter.get_next_line(reader))
        for transition in range(number_of_transitions):
            line = StochasticLabelledPetriNetImporter.get_next_line(reader)
            weight_str = StochasticLabelledPetriNetImporter.get_next_line(reader)
            weight = StochasticLabelledPetriNetImporter.parse_number(weight_str)

            if line.startswith("silent"):
                result.addSilentTransition(weight)
            elif line.startswith("label "):
                label = line[6:]
                result.addTransition(label, weight)
            else:
                raise RuntimeError("Invalid transition")

            # incoming places
            num_incoming = int(StochasticLabelledPetriNetImporter.get_next_line(reader))
            for _ in range(num_incoming):
                place = int(StochasticLabelledPetriNetImporter.get_next_line(reader))
                result.addPlaceTransitionArc(place, transition)

            # outgoing places
            num_outgoing = int(StochasticLabelledPetriNetImporter.get_next_line(reader))
            for _ in range(num_outgoing):
                place = int(StochasticLabelledPetriNetImporter.get_next_line(reader))
                result.addTransitionPlaceArc(transition, place)

        return result

    @staticmethod
    def parse_number(s: str) -> float:
        """Parse a number or a rational string like '3/4'."""
        s = s.strip()
        try:
            return float(s)
        except ValueError:
            if "/" in s:
                parts = s.split("/")
                if len(parts) == 2:
                    try:
                        f = Fraction(int(parts[0]), int(parts[1]))
                        return float(f)
                    except (ValueError, ZeroDivisionError):
                        return None
            return None

    @staticmethod
    def get_next_line(reader: io.TextIOBase) -> str:
        """Read next non-comment line (skipping those starting with '#')."""
        line = reader.readline()
        while line and line.startswith("#"):
            line = reader.readline()
        if not line:
            raise EOFError("Unexpected end of file")
        return line.strip()


def read_slpn(slpn_path):
    with open(slpn_path, "r", encoding="utf-8") as f:
        return StochasticLabelledPetriNetImporter.read(f)


# Example usage:
if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python slpn_importer.py <path_to_slpn_file>")
        sys.exit(1)
    path = sys.argv[1]
    net = read_slpn(path)
    print(net)


