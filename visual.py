import matplotlib.pyplot as plt
import forgi.visual.mplotlib as fvm
import forgi.graph.bulge_graph as fgb

dbseq = open("output.txt", "r").read()

rna = fgb.BulgeGraph.from_dotbracket(dbseq)

fig = plt.figure()
ax = fig.add_subplot(111)
fvm.plot_rna(rna, text_kwargs={"fontsize":2})
plt.show()