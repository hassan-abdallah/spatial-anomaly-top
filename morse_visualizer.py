import sys
import numpy as np
import plotly.graph_objs as go
from scipy.interpolate import griddata
from scipy.spatial import Delaunay
import ngl
from hdff import *
import hdtopology as hdt
from plotly.offline import plot
import colorsys
import argparse

# load data
def main(args):
    try:
        dataset = np.loadtxt(args.input)
    except:
        print('Can not load:', args.input)
        exit()
    print(dataset) 
    f = dataset[:, -1]
    sample = dataset[:, :-1].astype('float32')
    data = np.concatenate((sample, np.matrix(f).T), axis=1).astype('f')
    if args.edges == 1:
        edges = np.loadtxt(args.edges_file).astype(np.uint32)
    else: 
        edges = ngl.getSymmetricNeighborGraph('BSkeleton', sample, 10,1)
    # initialize extremum graph
    eg = hdt.ExtremumGraphExt()
    eg.initialize(data, np.array([0], dtype=np.uint8), edges, True, 300, 0, 1)

    # generate color palette
    def generate_distinct_colors(n):
        colors = []
        for i in range(n):
            hue = i / n
            saturation = 0.9
            value = 0.9
            r, g, b = colorsys.hsv_to_rgb(hue, saturation, value)
            colors.append(f'rgb({int(r*255)}, {int(g*255)}, {int(b*255)})')
        return colors

    tab10 = generate_distinct_colors(12)

    # obtain top twelve persistence values and normalize  
    raw_persistences = eg.persistences()[1:12]
    persistences = [float(p) for p in raw_persistences if float(p) < 1e10]
    persistences.reverse()
    min_p, max_p = min(persistences), max(persistences)
    slider_positions = [(p - min_p) / (max_p - min_p) for p in persistences]

    frames = []
    xmin = sample[:, 0].min()
    xmax = sample[:, 0].max()

    # build a grid to later interpolate function values over 
    grid_x, grid_y = np.mgrid[
        sample[:, 0].min():sample[:, 0].max():50j,
        sample[:, 1].min():sample[:, 1].max():50j
    ]

    # calculate convex hull of set of points
    grid_points = np.vstack((grid_x.ravel(), grid_y.ravel())).T
    tri = Delaunay(sample[:, :2])
    mask = tri.find_simplex(grid_points) >= 0

    for n in persistences:
        segs = eg.segmentation(eg.countForPersistence(n))
        seg_colors = np.full(len(sample), -1, dtype=int)
        for i, seg in enumerate(segs):
            seg_colors[seg] = i
        seg_max = seg_colors.max()
        if seg_max < 0:
            continue

        colors = [tab10[c] if c >= 0 else 'rgb(200,200,200)' for c in seg_colors]

        # interpolate and mask 
        grid_z = griddata(sample[:, :2], f, (grid_x, grid_y), method='cubic')
        grid_seg = griddata(sample[:, :2], seg_colors, (grid_x, grid_y), method='nearest')
        grid_z_flat = grid_z.ravel()
        grid_seg_flat = grid_seg.ravel().astype(float)
        grid_z_flat[~mask] = np.nan
        grid_seg_flat[~mask] = np.nan
        grid_z = grid_z_flat.reshape(grid_x.shape)
        grid_seg = grid_seg_flat.reshape(grid_x.shape)

        # generate surface plot
        surface = go.Surface(
            x=grid_x, y=grid_y, z=grid_z,
            surfacecolor=grid_seg,
            colorscale=[
                [i / max(seg_max, 1), color]
                for i, color in enumerate(tab10[:seg_max + 1])
            ],
            showscale=False,
            name="surface",
            showlegend=True,
            opacity=0.85,
            contours=dict(
                x=dict(show=True, color="black", width=1),
                y=dict(show=True, color="black", width=1),
                z=dict(show=False)
            ),
            lighting=dict(ambient=0.6, diffuse=0.8, specular=0.4, roughness=0.9, fresnel=0.1),
            lightposition=dict(x=100, y=200, z=0)
        )

        # generate scatter plot
        scatter = go.Scatter3d(
            x=sample[:, 0], y=sample[:, 1], z=f,
            mode='markers',
            marker=dict(size=2, color=colors),
            name="points",
        )

        # include edges 
        edge_x, edge_y, edge_z = [], [], []
        for i, j in edges:
            edge_x += [sample[i, 0], sample[j, 0], None]
            edge_y += [sample[i, 1], sample[j, 1], None]
            edge_z += [f[i], f[j], None]
        lines = go.Scatter3d(
            x=edge_x, y=edge_y, z=edge_z,
            mode='lines',
            line=dict(color='gray', width=1),
            name="edges",
            showlegend=True
        )

        frames.append(go.Frame(
            data=[surface, lines, scatter],
            name=f"{n}",
            layout=go.Layout(annotations=[
                dict(
                    showarrow=False,
                    text=f"<b>Extrema = {eg.countForPersistence(n)}</b>",
                    xref='paper', yref='paper',
                    x=0.05, y=0.95,
                    font=dict(size=24, color='black'),
                )
            ])
        ))

    layout = go.Layout(
        scene=dict(
            zaxis=dict(range=[f.min(), f.max() + f.max() * 0.05]),
            #xaxis_title='Latitude',
            #yaxis_title='Longitude',
            #zaxis_title='Mortality Rate (Cancer)',
            aspectmode='manual',
            aspectratio=dict(x=1, y=1, z=1)
        ),
        sliders=[{
            "steps": [
                {
                    "args": [[f"{n}"], {"frame": {"duration": 0}, "mode": "immediate"}],
                    "label": f"{n:.3f}",
                    "method": "animate"
                } for n in persistences
            ],
            "currentvalue": {"prefix": "Persistence= ", "font": {"size": 18}}
        }],
        legend=dict(font=dict(size=20))
    )

    fig = go.Figure(data=frames[0].data, frames=frames, layout=layout)
    plot(fig, filename=args.output, auto_open=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input', default=None, help="Specify data input filename", type=str,required=True)
    parser.add_argument('--edges', default=0, help="Specify if edges are provided in file (=1) or if they should be computed by program (=0)", type=str)
    parser.add_argument('--edges_file', default=None, help="Specify edges input filename", type=str)
    parser.add_argument('--output', help="Specify output filename", type=str, required=True)
    


    args = parser.parse_args()
    sys.exit(main(args))

