width = 30
height = 20
spacing = 1.0  # cell spacing
z = 0.0        # flat layer

output_file = "cell_positions.csv"

with open(output_file, "w") as f:
    for y in range(height):
        for x in range(width):
            xpos = (x + 0.5) * spacing
            ypos = (y + 0.5) * spacing
            f.write(f"{xpos},{ypos},{z}\n")
