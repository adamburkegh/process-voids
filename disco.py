
import sys

import pm4py
from pm4py.convert import convert_to_petri_net



debug = print



def discover_pn(name,event_log,noise_threshold):
    pt = pm4py.discover_process_tree_inductive(event_log,
                                               noise_threshold=noise_threshold)
    net, initial_marking, final_marking = convert_to_petri_net(pt)
    cover = 100*(1-noise_threshold)
    fname = f"{name}_ind_n{int(cover):d}"
    pm4py.write.write_ptml(pt,f"var/{fname}")
    pm4py.write.write_pnml(net, initial_marking, final_marking, 
                           f"var/{fname}")
    pm4py.save_vis_petri_net(net, initial_marking, final_marking, 
                             f"var/{fname}_pn.png")
    debug(f"Discovered net written to {fname}.[ptree|pnml|_pn.png]")


def disco_main():
    logname = sys.argv[1]
    outprefix = sys.argv[2]
    log = pm4py.read_xes(logname)
    for noise in [0,0.05,0.1,0.2,0.3,0.4]:
        discover_pn(outprefix,log,noise)


if __name__ == '__main__':
    disco_main()


