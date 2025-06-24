import toytree
import toyplot.png
from Bio import SeqIO

# 1. Make an ID -> full header dictionary for the database
db_fasta = "genome.fna"
id2header = {}
for record in SeqIO.parse(db_fasta, "fasta"):
    # Map only the first word after ">" to the full description
    id2header[record.id] = record.description

# 2. Add the query ID and header from mito.fasta
query_fasta = "mito.fasta"
for record in SeqIO.parse(query_fasta, "fasta"):
    query_id = record.id
    query_header = record.description
    id2header[query_id] = query_header

# 3. Load the tree and get its tip labels (IDs only)
tree = toytree.tree("tree.nwk")
tip_ids = tree.get_tip_labels()  # These are just the IDs

# 4. For each tip, use the full header as the label and color the query tip red
tip_labels = [id2header.get(tip, tip) for tip in tip_ids]
colors = ["red" if tip == query_id else "black" for tip in tip_ids]

canvas = tree.draw(
    width=1000,
    height=600,
    tip_labels=tip_labels,
    tip_labels_colors=colors,
)[0]

toyplot.png.render(canvas, "mito_tree.png")
