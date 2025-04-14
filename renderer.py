import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from shapely.geometry import LineString, Polygon, MultiPoint
from shapely.geometry.polygon import orient
from shapely.ops import polygonize
from scipy.spatial import Delaunay

# Load the shape data
with open("test.json", "r") as file:
    data = json.load(file)

# Flip coordinates along the Y-axis
def flip_y(points):
    return [(x, -y) for x, y in points]

# Alpha shape function with fixed alpha
def alpha_shape_fixed(points, alpha=0.7):
    if len(points) < 4:
        return None
    coords = np.array(points)
    try:
        tri = Delaunay(coords)
        triangles = coords[tri.simplices]

        def triangle_area(tri_pts):
            a = np.linalg.norm(tri_pts[0] - tri_pts[1])
            b = np.linalg.norm(tri_pts[1] - tri_pts[2])
            c = np.linalg.norm(tri_pts[2] - tri_pts[0])
            s = (a + b + c) / 2
            return np.sqrt(s * (s - a) * (s - b) * (s - c))

        def alpha_filter(tri_pts):
            area = triangle_area(tri_pts)
            max_edge = max(
                np.linalg.norm(tri_pts[0] - tri_pts[1]),
                np.linalg.norm(tri_pts[1] - tri_pts[2]),
                np.linalg.norm(tri_pts[2] - tri_pts[0])
            )
            return area / (max_edge ** 2) < alpha

        filtered_triangles = [tri for tri in triangles if alpha_filter(tri)]
        edges = set()
        for tri_pts in filtered_triangles:
            for i in range(3):
                edge = tuple(sorted((tuple(tri_pts[i]), tuple(tri_pts[(i + 1) % 3]))))
                if edge in edges:
                    edges.remove(edge)
                else:
                    edges.add(edge)

        line_strings = [LineString(edge) for edge in edges]
        polygons = list(polygonize(line_strings))
        if polygons:
            return polygons[0]
    except Exception:
        return None
    return None

# Create the plot for alpha = 0.7
fig, axs = plt.subplots(1, 2, figsize=(20, 12))
axs[0].set_title("TOP Mounting Side")
axs[1].set_title("BOTTOM Mounting Side")

for ax in axs:
    ax.set_aspect('equal')
    ax.axis('off')

for shape_entry in data[0]["shapes"]:
    outline = shape_entry.get("outline", [])
    components = shape_entry.get("components", [])

    if not components:
        components = [{"position": {"x": 0, "y": 0}, "mountingSide": "TOP"}]

    for comp in components:
        pos = comp["position"]
        side = comp.get("mountingSide", "TOP").upper()
        ax = axs[0] if side == "TOP" else axs[1]

        segments = []
        raw_points = []
        for item in outline:
            if item["type"] == "Line":
                p1 = (item["points"][0]["x"], item["points"][0]["y"])
                p2 = (item["points"][1]["x"], item["points"][1]["y"])
                p1, p2 = flip_y([p1, p2])
                x1, y1 = p1[0] + pos["x"], p1[1] + pos["y"]
                x2, y2 = p2[0] + pos["x"], p2[1] + pos["y"]
                segments.append(LineString([(x1, y1), (x2, y2)]))
                raw_points.extend([(x1, y1), (x2, y2)])

        polygon = alpha_shape_fixed(raw_points, alpha=0.7)

        if polygon and isinstance(polygon, Polygon):
            polygon = orient(polygon)
            coords = list(polygon.exterior.coords)
            patch = patches.Polygon(coords, closed=True, fill=True,
                                     facecolor='green', edgecolor='black', alpha=0.4, linewidth=0.7)
            ax.add_patch(patch)

        for seg in segments:
            x1, y1 = seg.coords[0]
            x2, y2 = seg.coords[1]
            ax.plot([x1, x2], [y1, y2], 'b-', linewidth=0.5)

        for item in outline:
            if item["type"] == "Circle":
                center = item["location"]
                radius = item["radius"]
                cx, cy = flip_y([(center["x"], center["y"])])[0]
                cx += pos["x"]
                cy += pos["y"]
                circle = patches.Circle((cx, cy), radius, fill=False, edgecolor='green', linewidth=0.5)
                ax.add_patch(circle)

# Save the final image
output_path = "output.png"
plt.tight_layout()
plt.savefig(output_path, dpi=300)
plt.close()

print("Saved to:", output_path)
