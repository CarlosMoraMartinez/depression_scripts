# source plot_env/bin/activate
import sys
import pandas as pd
import plotly.graph_objects as go
import numpy as np

C_NS = "#A5ABBD"

name = sys.argv[1]
nodes_file = sys.argv[2]
links_file = sys.argv[3]

nodes = pd.read_csv(nodes_file, sep="\t")
links = pd.read_csv(links_file, sep="\t")

links["color"] = links["color"].fillna(C_NS)
nodes["color"] = nodes["color"].fillna(C_NS)

links["source"] = links["source"].astype(int)
links["target"] = links["target"].astype(int)
links["Size"] = links["Size"].astype(float)


#num_conds = len([i for i in nodes["xpos"] if i ==0 ])
#num_sps = len([i for i in nodes["xpos"] if i ==1 ])
#num_procs = nodes.shape[0] - num_conds - num_sps

#nodes["ypos"] = [i for i in np.linspace(0, 1, num_conds)] + [i for i in np.linspace(0, 1, num_procs)] + [i for i in np.linspace(0, 1, num_sps)]


# Build Sankey figure
fig = go.Figure(data=[go.Sankey(
    arrangement="fixed",  # Allows custom x positions
    node=dict(
        pad=15,
        thickness=15,
        line=dict(color="black", width=0.5),
        label=nodes["value"],
        color=nodes["color"]#,
        #x=nodes["xpos"].astype(float)  # Custom horizontal positions
        #y= nodes["ypos"]
        # Optional: y=nodes["ypos"] if you want vertical control
    ),
    link=dict(
        source=links["source"],
        target=links["target"],
        value=links["Size"],
        color=links["color"]
    )
)])
#fig.write_html(f"{name}_py.html")
#fig.show()

fig.write_image(f"{name}_py.pdf", format="pdf")

