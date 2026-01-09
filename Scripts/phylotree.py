import sys
import toytree
import toyplot.png
from Bio import SeqIO
from pathlib import Path

# ============================
# Arguments
# ============================
if len(sys.argv) != 6:
    print("Usage: python plot_mitophylo.py <db_fasta> <query_fasta> <tree_nwk> <species_name> <out_png>")
    sys.exit(1)

db_fasta = Path(sys.argv[1])
query_fasta = Path(sys.argv[2])
tree_file = Path(sys.argv[3])
species_name = sys.argv[4]
out_png = Path(sys.argv[5])

# ============================
# Load database headers
# ============================
id2header = {}
for record in SeqIO.parse(db_fasta, "fasta"):
    id2header[record.id] = record.description

# ============================
# Add query sequence
# ============================
for record in SeqIO.parse(query_fasta, "fasta"):
    query_id = record.id
    id2header[query_id] = record.description

# ============================
# Load tree
# ============================
tree = toytree.tree(tree_file)
tip_ids = tree.get_tip_labels()

# ============================
# Map full labels and color query
# ============================
tip_labels = [id2header.get(tip, tip) for tip in tip_ids]
colors = ["red" if tip == query_id else "black" for tip in tip_ids]

# ============================
# Draw tree
# ============================
canvas = tree.draw(
    width=1000,
    height=600,
    tip_labels=tip_labels,
    tip_labels_colors=colors,
)[0]

# ============================
# Render to PNG
# ============================
toyplot.png.render(canvas, out_png)
print(f"✅ Tree plotted: {out_png}")
