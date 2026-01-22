import random
width = 300
height = 300
depth = 0
spacing = 1.0  # cell spacing
zpos = 0.0        # flat layer

output_file = "rectangle.csv"

with open(output_file, "w") as f:
    f.write("X,Y,Z\n")
    #for z in range(depth):
    for y in range(height):
        for x in range(width):
            #if random.uniform(0, 1) > 0.1:
            xpos = (x + 0.5) * spacing
            ypos = (y + 0.5) * spacing
            #zpos = (z + 0.5) * spacing
            f.write(f"{xpos},{ypos},{zpos}\n")
