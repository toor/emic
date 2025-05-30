import numpy as np
np.set_printoptions(threshold = np.inf)

max_y = 7
y_step = 1
y_positions = np.arange(2, max_y + y_step, y_step)
max_x = 20
x_step = 1

reps = np.array([1, 2, 3])

betaJ_vals = np.arange(-10.0, 10.1, 0.1)
betaJ_vals = np.where(np.isclose(betaJ_vals, -0.0, atol=1e-9), 0.0, betaJ_vals)

def generate_x_positions(y, max_x, step):
    return list(range(y, max_x + step, step))

# Generate all (x, y) combinations
combinations = [(x, y) for y in y_positions for x in generate_x_positions(y, max_x, x_step)]

rule all:
  input: 
         expand("results/SIM_NAME/N_y={y}/N_x={x}/betaJ={betaJ:.1f}/rep{r}/data.csv", y=[y for x, y in combinations], x=[x for x,y in combinations], betaJ=betaJ_vals, r=reps),
         "figs/SIM_NAME/entropy_vs_Nx.pdf"

rule run:
  output: "results/SIM_NAME/N_y={y}/N_x={x}/betaJ={betaJ}/rep{r}/data.csv",
          "results/SIM_NAME/N_y={y}/N_x={x}/betaJ={betaJ}/rep{r}/equ_lattice.vertex",
          "results/SIM_NAME/N_y={y}/N_x={x}/betaJ={betaJ}/rep{r}/seed.txt"
  params: in_file = "none"
  shell: "./main {wildcards.betaJ} {wildcards.x} {wildcards.y} {output[0]} {output[1]} {output[2]}" 

rule analyse_data:
  output: "figs/SIM_NAME/entropy_vs_Nx.pdf"
  params: min_y = np.min(y_positions),
          max_y = np.max(y_positions),
          y_step = 1,
          max_x = max_x,
          x_step = x_step,
          prefix = "SIM_241021_0001",
          Jmax = float(np.ceil(np.max(betaJ_vals))),
          Jstep = 0.1,
          reps = 1
  script: "analyse_data.py"

  # Rule to combine the data into a shared states_visited.csv
rule combine_states_visited:
    input:
          expand("results/SIM_240923_0004/N_y={{y}}/N_x={{x}}/rep{{r}}/states_visited_{betaJ}.csv", betaJ=betaJ_vals)
    output:
        "results/SIM_240923_0004/N_y={y}/N_x={x}/rep{r}/states_visited.csv"
    shell:
        "cat {input} > {output} && rm {input}"  # Combine the state files into one

rule check_equilibration:
  output: "results/SIM_240822_0002/N_y={y}/N_x={x}/beta={beta}/J={J}/rep{r}/data.csv",
          "results/SIM_240822_0002/N_y={y}/N_x={x}/beta={beta}/J={J}/rep{r}/equ_lattice.vertex"
  params: in_file = "none"
  shell: "./sixvertex {wildcards.betaJ} {wildcards.x} {wildcards.y} {output[0]} {output[1]} {params.in_file}"

rule plot_autocorrelation:
  input: expand("results/SIM_240607/N_y={{y}}/N_x={x}/beta={{beta}}/J={{J}}/rep1/data.csv", x=generate_x_positions(20, 80, 5))
  output: "figs/autocorrelation/N_y={y}/beta={beta}/J={J}/tau_vs_N_x.pdf"
  script: "plot_autocorrelation.py"
