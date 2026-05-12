import requests
import pandas as pd

protein = "TP53"
species = 9606  # Human

url = f"https://string-db.org/api/tsv/interaction_partners"
params = {"identifiers": protein, "species": species, "limit": 20}
response = requests.get(url, params=params)

lines = response.text.strip().split("\n")
rows = [line.split("\t") for line in lines[1:]]
df = pd.DataFrame(rows, columns=lines[0].split("\t"))
print(df[["preferredName_A", "preferredName_B", "score"]].head())

import networkx as nx

G = nx.Graph()

for _, row in df.iterrows():
    score = float(row["score"])
    if score > 0.7:  # only high-confidence interactions
        G.add_edge(row["preferredName_A"], row["preferredName_B"], weight=score)

print(f"Nodes: {G.number_of_nodes()}")
print(f"Edges: {G.number_of_edges()}")

degree_centrality = nx.degree_centrality(G)

# Sort to find top hubs
top_hubs = sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)
print("Top 5 hub proteins:")
for protein, score in top_hubs[:5]:
    print(f"  {protein}: {score:.3f}")

import matplotlib.pyplot as plt

plt.figure(figsize=(12, 10))

# Node size = degree (bigger = more connections)
node_sizes = [3000 * degree_centrality[n] + 300 for n in G.nodes()]

# Node colour = degree tier
degrees = dict(G.degree())
node_colors = []
for n in G.nodes():
    d = degrees[n]
    if d >= 5:
        node_colors.append("#1D9E75")   # hub
    elif d >= 3:
        node_colors.append("#378ADD")   # medium
    else:
        node_colors.append("#AFA9EC")   # low

pos = nx.spring_layout(G, seed=42)  # consistent layout

nx.draw_networkx(
    G, pos,
    node_size=node_sizes,
    node_color=node_colors,
    font_size=9,
    font_weight="bold",
    edge_color="#cccccc",
    width=1.2,
    with_labels=True
)

plt.title("TP53 Protein Interaction Network (STRING DB)", fontsize=14)
plt.axis("off")
plt.tight_layout()
plt.savefig("protein_network.png", dpi=150, bbox_inches="tight")
plt.show()
