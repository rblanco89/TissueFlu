import random
width = 300
height = 200
spacing = 1.0  # cell spacing
z = 0.0        # flat layer

output_file = "cell_positions.csv"

with open(output_file, "w") as f:
    for y in range(height):
        for x in range(width):
            #if random.uniform(0, 1) > 0.2:
            xpos = (x + 0.5) * spacing
            ypos = (y + 0.5) * spacing
            f.write(f"{xpos},{ypos},{z}\n")
