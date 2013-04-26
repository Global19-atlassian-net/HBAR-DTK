#!/usr/bin/env python
"""
A simple script to convert Celera(R) Assembler's "best.edges" to a GML which can be used to feed into Gephi to
check the topology of the best overlapping graph.

Usage:
    python CA_best_edge_to_GML.py asm.gkp_store asm.tigStore best.edge output.gml 
"""

import networkx as nx
import os
import shlex
import sys
import subprocess

gkp_store = sys.argv[1]
tig_store = sys.argv[2]
best_edge_file = sys.argv[3]
output_file = sys.argv[4]



G=nx.DiGraph()
frg_to_tig = {}
args = shlex.split("tigStore -g %s -t %s 1 -D unitiglist" % (gkp_store, tig_store ))
out = subprocess.check_output(args)
out = out.split("\n")
for l in out:
    l = l.strip().split()
    if len(l) == 0: continue
    if l[0] == "maID": continue
    unitig_id = int(l[0])

    os.system("tigStore -g %s -t %s 1 -d fr -u %d > frag_list" % ( gkp_store, tig_store, unitig_id) )

    args = shlex.split( "tigStore -g %s -t %s 1 -d fr -u %d" % ( gkp_store, tig_store, unitig_id) )
    f_out = subprocess.check_output(args)
    f_out = f_out.split("\n")
    for l in f_out:
        """FRG    1453 179419,182165"""
        l = l.replace(",", " ")
        l = l.strip().split()
        if len(l) == 0: continue
        frg_id = l[1]
        frg_to_tig[frg_id] = unitig_id

with open(best_edge_file) as f:
    for l in f:
        if l[0] == "#": continue
        l = l.strip().split()
        id1, lib_id, best5, o1, best3, o3 = l
        G.add_node(id1, unitig="utg%s" % frg_to_tig[id1])
        if best5 != "0":
            G.add_edge(best5, id1)
        if best3 != "0":
            G.add_edge(id1, best3)
        #G[id1]["unitig"] = frg_to_tig[id1]

nx.write_gml(G, output_file)
