import numpy as np 
import matplotlib.pyplot as plt
import os 
import sys
import pandas as pd
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle

plt.rcParams['text.usetex'] = True

betaJ_plus_values = np.arange(-10.0, 11.0, 1.0)
betaJ_min_values = np.arange(10.0, -11.0, -1.0)
betaJ_values = betaJ_plus_values #+ betaJ_min_values

betaJ_max = 10.0

def load_dataset(betaJ, f):
    offset = int(f * 5e4)
    file_path = f'results/SIM_250217_0001/N_y=5/N_x=5/betaJ={betaJ}/rep1/equ_lattice_step{offset}.vertex'
    print(file_path)
    # Load the dataset from the CSV file
    dataset = pd.read_csv(file_path, sep='\s+', index_col=False)
    return dataset

def update(frame):
    plt.clf()
    ax = plt.gca()
    # There are $N_\beta$ values of $\beta J$, and $N_f$ frames for each value of $\beta_J$. 
    #If we construct the frame table in the correct order, we can extract the value of $\beta J$ as well as the frame index for that subset.
    j = frame // N_f 
    f = frame % N_f 
    
    betaJ = betaJ_values[j]
    
    data = load_dataset(betaJ, f)

    plot_vertices_from_raw_data(data, betaJ)
    x = data['x'].to_numpy()
    y = data['y'].to_numpy()
    
    draw_large_magnet_state(x, y, betaJ)


def draw_large_magnet_state(x, y, betaJ):
    angle = np.arccos(betaJ/betaJ_max)
    
    # start point
    y0 = np.median(y)
    x0 = np.max(x) + 0.5 
    # end point
    L = np.max(y) / 2.0 
    x1 = L*np.sin(angle)
    y1 = L*np.cos(angle)
    # left magnet
    plt.arrow(x0 + 0.5, y0, x1, y1, head_width=0.1, head_length=0.1, length_includes_head = True)
    # right magnet
    plt.arrow(np.min(x) - 1.0, y0, 0, L, head_width=0.1, head_length=0.1, length_includes_head = True)   

# Takes a vertex type and coordinates, and draws the corresponding
# vertex using the provided figure/axis objects.
def plot_vertex_for_movie(vert_type, x, y):
    dx = dy = 0.5
    out_kwargs = {
        "head_width": 0.1,
        "head_length": 0.1,
        "length_includes_head": True,
        "color": "red"
    }
    in_kwargs = {
        "head_width": 0.1,
        "head_length": 0.1,
        "length_includes_head": True,
        "color": "blue"
    }

    #print(f"vertex at coords ({x},{y}) has type {vert_type}")
    match vert_type:
        case 1:
            plt.arrow(x - dx, y, dx, 0, **in_kwargs)
            plt.arrow(x, y, dx, 0, **out_kwargs)
            plt.arrow(x, y, 0, dy, **out_kwargs)
            plt.arrow(x, y - dy, 0, dy, **in_kwargs)
        case 2:
            plt.arrow(x, y, -dx, 0, **out_kwargs)
            plt.arrow(x + dx, y, -dx, 0, **in_kwargs)
            plt.arrow(x, y + dy, 0, -dy, **in_kwargs)
            plt.arrow(x, y, 0, -dy, **out_kwargs)
        case 3:
            plt.arrow(x - dx, y, dx, 0, **in_kwargs)
            plt.arrow(x, y, dx, 0, **out_kwargs)
            plt.arrow(x, y + dy, 0, -dy, **in_kwargs)
            plt.arrow(x, y, 0, -dy, **out_kwargs)
        case 4:
            plt.arrow(x, y, -dx, 0, **out_kwargs)
            plt.arrow(x + dx, y, -dx, 0, **in_kwargs)
            plt.arrow(x, y, 0, dy, **out_kwargs)
            plt.arrow(x, y - dy, 0, dy, **in_kwargs)
        case 6:
            plt.arrow(x, y, -dx, 0, **out_kwargs)
            plt.arrow(x, y, dx, 0, **out_kwargs)
            plt.arrow(x, y + dy, 0, -dy, **in_kwargs)
            plt.arrow(x, y - dy, 0, dy, **in_kwargs)
        case 5:
            plt.arrow(x - dx, y, dx, 0, **in_kwargs)
            plt.arrow(x + dx, y, -dx, 0, **in_kwargs)
            plt.arrow(x, y, 0, dy, **out_kwargs)
            plt.arrow(x, y, 0, -dy, **out_kwargs)
        case _:
            print("Invalid vertex type " + str(vert_type) + ", vertex type must be from 1-6")

# Takes a vertex type and coordinates, and draws the corresponding
# vertex using the provided figure/axis objects.
def plot_vertex(x, y, vert_type, colour, fig, ax):
    print(str(vert_type))
    dx = dy = 0.5

    kwargs = {
        "head_width": 0.1,
        "head_length": 0.1,
        "length_includes_head": True,
        "color": colour
    }
    
    x = int(x)
    y = int(y)
    vert_type = int(vert_type)

    #print(f"vertex at coords ({x},{y}) has type {vert_type}")
    match vert_type:
        case 1:
            ax.arrow(x - dx, y, dx, 0, **kwargs)
            ax.arrow(x, y, dx, 0, **kwargs)
            ax.arrow(x, y, 0, dy, **kwargs)
            ax.arrow(x, y - dy, 0, dy, **kwargs)
        case 2:
            ax.arrow(x, y, -dx, 0, **kwargs)
            ax.arrow(x + dx, y, -dx, 0, **kwargs)
            ax.arrow(x, y + dy, 0, -dy, **kwargs)
            ax.arrow(x, y, 0, -dy, **kwargs)
        case 3:
            ax.arrow(x - dx, y, dx, 0, **kwargs)
            ax.arrow(x, y, dx, 0, **kwargs)
            ax.arrow(x, y + dy, 0, -dy, **kwargs)
            ax.arrow(x, y, 0, -dy, **kwargs)
        case 4:
            ax.arrow(x, y, -dx, 0, **kwargs)
            ax.arrow(x + dx, y, -dx, 0, **kwargs)
            ax.arrow(x, y, 0, dy, **kwargs)
            ax.arrow(x, y - dy, 0, dy, **kwargs)
        case 6:
            ax.arrow(x, y, -dx, 0, **kwargs)
            ax.arrow(x, y, dx, 0, **kwargs)
            ax.arrow(x, y + dy, 0, -dy, **kwargs)
            ax.arrow(x, y - dy, 0, dy, **kwargs)
        case 5:
            ax.arrow(x - dx, y, dx, 0, **kwargs)
            ax.arrow(x + dx, y, -dx, 0, **kwargs)
            ax.arrow(x, y, 0, dy, **kwargs)
            ax.arrow(x, y, 0, -dy, **kwargs)
        case _:
            print("Invalid vertex type " + str(vert_type) + ", vertex type must be from 1-6")
    #ax.text(x - 0.6, y + 0.6, f"Type {vert_type}", fontsize=8)


dataset = {
    "1": {
        "1": np.array([[0, 0, 2, "red"], [0, 1, 2, "red"], [1, 0, 2, "blue"], [1, 1, 2, "blue"]]),
        "2": np.array([[0, 0, 2, "red"], [0, 1, 2, "red"], [1, 0, 4, "blue"], [1, 1, 4, "blue"]]),
        "3": np.array([[0, 0, 2, "red"], [0, 1, 2, "red"], [1, 0, 6, "blue"], [1, 1, 2, "blue"]]),
        "4": np.array([[0, 0, 2, "red"], [0, 1, 2, "red"], [1, 0, 5, "blue"], [1, 1, 6, "blue"]])
    },
    "2": {
        "1": np.array([[0, 0, 4, "red"], [0, 1, 4, "red"], [1, 0, 2, "blue"], [1, 1, 2, "blue"]]),
        "2": np.array([[0, 0, 4, "red"], [0, 1, 4, "red"], [1, 0, 4, "blue"], [1, 1, 4, "blue"]]),
        "3": np.array([[0, 0, 4, "red"], [0, 1, 4, "red"], [1, 0, 6, "blue"], [1, 1, 2, "blue"]]),
        "4": np.array([[0, 0, 4, "red"], [0, 1, 4, "red"], [1, 0, 5, "blue"], [1, 1, 6, "blue"]]) 
    },
    "3": {
        "1": np.array([[0, 0, 4, "red"], [0, 1, 6, "red"], [1, 0, 6, "blue"], [1, 1, 3, "blue"]]),
        "2": np.array([[0, 0, 4, "red"], [0, 1, 6, "red"], [1, 0, 2, "blue"], [1, 1, 5, "blue"]]),
        "3": np.array([[0, 0, 4, "red"], [0, 1, 6, "red"], [1, 0, 4, "blue"], [1, 1, 1, "blue"]]),
        "4": np.array([[0, 0, 4, "red"], [0, 1, 6, "red"], [1, 0, 2, "blue"], [1, 1, 3, "blue"]]),
        "5": np.array([[0, 0, 4, "red"], [0, 1, 6, "red"], [1, 0, 6, "blue"], [1, 1, 5, "blue"]])
    },
    "4": {
        "1": np.array([[0, 0, 6, "red"], [0, 1, 2, "red"], [1, 0, 5, "blue"], [1, 1, 4, "blue"]]),
        "2": np.array([[0, 0, 6, "red"], [0, 1, 2, "red"], [1, 0, 1, "blue"], [1, 1, 4, "blue"]]),
        "3": np.array([[0, 0, 6, "red"], [0, 1, 2, "red"], [1, 0, 5, "blue"], [1, 1, 6, "blue"]]),
        "4": np.array([[0, 0, 6, "red"], [0, 1, 2, "red"], [1, 0, 1, "blue"], [1, 1, 6, "blue"]]),
        "5": np.array([[0, 0, 6, "red"], [0, 1, 2, "red"], [1, 0, 3, "blue"], [1, 1, 2, "blue"]])
    }
}


def plot_vertices_from_dict():
    fig, axs = plt.subplots(nrows = 4, ncols = 5)
    for ax in axs.flatten():
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
    for i in range(1, 5):
        dic = dataset[str(i)]
        for j in range(1, len(dic) + 1):
            ax = axs[i - 1][j - 1]
            state_data = dic[str(j)]
            for state in state_data:
                print(state)
                [x, y, vert_type, colour] = state
                plot_vertex(x, y, vert_type, colour, fig, ax)

    fig.savefig("all_2x2_analytic_states.pdf", dpi=300.0)



def plot_vertices_from_raw_data(data, betaJ):
    betaJ_max = 10.0
    x = data['x'].to_numpy()
    y = data['y'].to_numpy()
    states = data['state'].to_numpy() + np.ones_like(x)
    count = x.size
    

    for i in range(count):
        plot_vertex_for_movie(states[i], x[i], y[i])


    #plt.xlabel('x')
    #plt.ylabel('y')
    #plt.title(f'betaJ={betaJ}')
    #plt.legend()


def plot_vertex_different_colours(vert_type, x, y, fig, ax, bond_tracker):
    dx = dy = 0.5
    gap = 2*dx
    
    kwargs_dict = {
        "1": {
            "head_width": 0.1,
            "head_length": 0.1,
            "length_includes_head": True,
            "color": "red"
        },
        "2": {
            "head_width": 0.1,
            "head_length": 0.1,
            "length_includes_head": True,
            "color": "blue"
        },
        "3": {
            "head_width": 0.1,
            "head_length": 0.1,
            "length_includes_head": True,
            "color": "green"
        },
        "4": {
            "head_width": 0.1,
            "head_length": 0.1,
            "length_includes_head": True,
            "color": "black"
        },
        "5": {
            "head_width": 0.1,
            "head_length": 0.1,
            "length_includes_head": True,
            "color": "magenta"
        },
        "6": {
            "head_width": 0.1,
            "head_length": 0.1,
            "length_includes_head": True,
            "color": "orange"
        }
    }
    
    kwargs = kwargs_dict.get(str(vert_type), None)
    
    if not kwargs:
        print(f"Invalid vertex type {vert_type}, must be 1-6")
        return   
    
    def bond_key(x1, y1, x2, y2):
        return (int(x1), int(y1), int(x2), int(y2))

    def draw_arrow(x1, y1, arrow_coords, idx):
        match idx:
            case 0: print("LEFT\n")
            case 1: print("RIGHT\n")
            case 2: print("ABOVE\n")
            case 3: print("BELOW\n")
        print(f"(x, y) = ({x1}, {y1})\n")
        # All possible bonds that can surround the vertex
        # In order to have a consistent convention for the bonds, we need to define a direction.
        # We say that bonds go from left to right and from down to up.
        # example: (x1, y1) = (0,0) ; RIGHT = (0, 0, 1, 0)
        #          (x1, y1) = (1,0) ; LEFT = (1 - 1, 0, 1, 0) = (0, 0, 1, 0)
        # such that these two cases refer to the same bond, as they should
        bonds = [(x1-gap, y1, x1, y1),
                    (x1, y1, x1 + gap, y1),
                    (x1, y1, x1, y1+gap),
                    (x1, y1 - gap, x1, y1)]
        key = bonds[idx]
        
        (x, y, dx, dy) = arrow_coords

        if key in bond_tracker:
            print(f"I entered this if statement from key: {key}\n")
            ax.arrow(x, y, dx, dy, color=kwargs["color"], head_width=0)  # No head
        else:
            ax.arrow(x, y, dx, dy, **kwargs)  # With head
            print(key)
            bond_tracker.add(key)

    # Arrows are drawn in the following order:
    # 0 - left, 1, right, 2, above, 3, below.
    match vert_type:
        case 1:
            draw_arrow(x, y, (x - dx, y, dx, 0), 0)          
            draw_arrow(x, y, (x, y, dx, 0), 1)
            draw_arrow(x, y, (x, y, 0, dy), 2)          
            draw_arrow(x, y, (x, y - dy, 0, dy), 3)
        case 2:
            draw_arrow(x, y, (x, y, -dx, 0), 0)
            draw_arrow(x, y, (x + dx, y, -dx, 0), 1)
            draw_arrow(x, y, (x, y + dy, 0, -dy), 2)
            draw_arrow(x, y, (x, y, 0, -dy), 3)
        case 3:
            draw_arrow(x, y, (x - dx, y, dx, 0), 0)
            draw_arrow(x, y, (x, y, dx, 0), 1)
            draw_arrow(x, y, (x, y + dy, 0, -dy), 2)
            draw_arrow(x, y, (x, y, 0, -dy), 3)
        case 4:
            draw_arrow(x, y, (x, y, -dx, 0), 0)
            draw_arrow(x, y, (x + dx, y, -dx, 0), 1)
            draw_arrow(x, y, (x, y, 0, dy), 2)
            draw_arrow(x, y, (x, y - dy, 0, dy), 3)
        case 6:
            draw_arrow(x, y, (x, y, -dx, 0), 0)
            draw_arrow(x, y, (x, y, dx, 0), 1)
            draw_arrow(x, y, (x, y + dy, 0, -dy), 2)
            draw_arrow(x, y, (x, y - dy, 0, dy), 3)
        case 5:
            draw_arrow(x, y, (x - dx, y, dx, 0), 0)
            draw_arrow(x, y, (x + dx, y, -dx, 0), 1)
            draw_arrow(x, y, (x, y, 0, dy), 2)
            draw_arrow(x, y, (x, y, 0, -dy), 3)


def plot_vertices(filename, out_file):
    fig, ax = plt.subplots()
    data = pd.read_csv(filename, sep='\s+', index_col=False)
    x = data['x'].to_numpy()
    y = data['y'].to_numpy()
    states = data['state'].to_numpy() + np.ones_like(x)
    
    bond_tracker = set()

    for i in range(x.size):
        plot_vertex_different_colours(states[i], x[i], y[i], fig, ax, bond_tracker)

    print(bond_tracker)    
    x0 = np.min(x)
    y0 = np.min(y)
    xmax = np.max(x)
    ymax = np.max(y)
    midpoint = np.median(y)

    print(x0, y0, xmax, ymax)

    r1 = Rectangle((x0 - 3, y0), width=2.0, height=float(ymax), facecolor='red', edgecolor='black')
    r2 = Rectangle((xmax + 1, y0), width=2.0, height=float(ymax), facecolor='deepskyblue', edgecolor='black')
    
    ax.add_patch(r1)
    ax.add_patch(r2)

    ax.arrow(x0 - 4, y0 - 1, 1.0, 0, head_width=0.1, head_length=0.1, length_includes_head = True)
    ax.arrow(x0-4, y0 - 1, 0, 1.0, head_width=0.1, head_length=0.1, length_includes_head=True)
    ax.text(x0-4.3, y0, r'$y$', fontsize='large')
    ax.text(x0-2.9, y0 - 1.1, r'$x$', fontsize='large')


    ax.arrow(x0 - 1.5, midpoint, -1.0, 0, head_width=0.1, head_length=0.1, length_includes_head=True)
    ax.arrow(xmax + 1.5, midpoint, 1.0, 0, head_width=0.1, head_length=0.1, length_includes_head=True)
    ax.arrow(xmax + 0.75, midpoint, 0.5, 0, head_width=0.1, head_length=0.1, length_includes_head=True)
    ax.arrow(xmax + 1.25, midpoint, -0.5, 0, head_width=0.1, head_length=0.1, length_includes_head=True)
    ax.text(xmax + 0.4, midpoint - 0.5, r'$J_R$', fontsize='x-large')

    ax.set_box_aspect(5/11)
    ax.axis('off')

    fig.savefig(out_file, dpi=300.0, bbox_inches='tight', pad_inches=0.01)
    plt.close(fig)

plot_vertices("test_data_20240618.vertex", "example_configuration.pdf")

#
#fig, ax = plt.subplots()
#ax.set_aspect(aspect='equal', adjustable='box')
#ax.spines['top'].set_visible(False)
#ax.spines['bottom'].set_visible(False)
#ax.spines['left'].set_visible(False)
#ax.spines['right'].set_visible(False)
#
#ax.set_xticks([])
#ax.set_yticks([])
#
#N_f = 10
#
#ani = FuncAnimation(fig, update, frames=np.arange(0, len(betaJ_values)*N_f, 1), interval=10000)
#
#ani.save('vertex_states_vs_betaJ.mp4', writer='ffmpeg', fps=10)
