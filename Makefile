# === Configuration ===
HOST_COMPILER = gcc
NVCC = nvcc -ccbin $(HOST_COMPILER)
NVCCFLAGS = -O2 -arch=sm_86
LIBRARIES = -lm -lcurand

# === File Structure ===
SRC_DIR = src
SRCS_C = $(wildcard $(SRC_DIR)/*.c)
SRCS_CU = $(wildcard $(SRC_DIR)/*.cu)

OBJS_C = $(SRCS_C:.c=.o)
OBJS_CU = $(SRCS_CU:.cu=.o)

TARGET = aeroflu

# === Build Rules ===
all: $(TARGET)

$(TARGET): $(OBJS_C) $(OBJS_CU)
	$(NVCC) $(NVCCFLAGS) $^ -o $@ $(LIBRARIES)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# === Clean ===
clean:
	rm -f $(SRC_DIR)/*.o $(TARGET)

# === Run ===
run: all
	mkdir -p results
	./$(TARGET) --config=config.conf --structure=structures/rectangle.csv
