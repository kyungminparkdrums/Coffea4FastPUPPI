import torch

# path to one of pt files
path = "graphs_may26/graphs_10.pt"

torch.set_printoptions(
    precision=6,
    sci_mode=False,
)

payload = torch.load(
    path,
    map_location="cpu",
    #weights_only=False, # doesn't work in old version? boh
)

# -------------------------------------------------
# New dataset format
# -------------------------------------------------

graphs = payload["graphs"]
metadata = payload["metadata"]

feature_names = metadata["feature_names"]

print("Type:", type(graphs))
print("Number of graphs:", len(graphs))

print("\nMetadata:")
for k, v in metadata.items():
    print(f"{k}: {v}")

g = graphs[0]

print("\nFirst graph:")
print(g)

print("\nAvailable attributes:")
print(g.keys)

print("\nx shape:", g.x.shape)
print("y:", g.y)
print("event_idx:", g.event_idx)
print("center_idx:", g.center_idx)

print("\nFeature layout:")
for i, name in enumerate(feature_names):
    print(f"{i:2d}: {name}")

print("\nFirst few node features:")
print(g.x[:5])

print("\nFirst few node features with names:")

n_print = min(5, g.x.shape[0])

for inode in range(n_print):

    print(f"\nNode {inode}")

    for i, name in enumerate(feature_names):

        print(
            f"  {i:2d} {name:20s}: "
            f"{g.x[inode, i].item(): .6f}"
        )

# -------------------------------------------------
# Center node
# -------------------------------------------------

is_center_idx = feature_names.index("is_center")

print("\nCenter node:")

center_mask = g.x[:, is_center_idx] > 0.5

center_positions = torch.where(center_mask)[0]

print("center positions:", center_positions.tolist())

for pos in center_positions.tolist():

    print(f"\nCenter node position {pos}")

    for i, name in enumerate(feature_names):

        print(
            f"  {i:2d} {name:20s}: "
            f"{g.x[pos, i].item(): .6f}"
        )
