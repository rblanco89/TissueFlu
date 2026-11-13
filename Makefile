# TissueFlu: a GPU agent-based model of influenza spread in lung tissue.
#
#   make                          build ./tissueFlu
#   make ARCH=sm_75               build for a different GPU
#   make HOST_COMPILER=g++-12     force the host compiler nvcc calls
#   make run STRUCTURE=cells.csv  build, then run lung.conf on that tissue
#   make clean
#   make help                     show the detected settings
#
# Needs the CUDA toolkit (nvcc) and a host C++ compiler that nvcc accepts.

# === GPU architecture ===
# sm_86 covers the cards this model was developed on (RTX 30xx, A2000). Other
# common values: sm_60 P100, sm_70 V100, sm_75 T4 / RTX 20xx, sm_80 A100,
# sm_89 RTX 40xx, sm_90 H100, sm_120 RTX 50xx. The build also embeds PTX, so a
# binary runs on newer cards by just-in-time compilation, but not on older
# ones. "native" targets the card in the build machine (CUDA 11.5 or newer);
# it needs that card to be visible while compiling, so it is no use on a
# cluster login node.
ARCH ?= sm_86

# === Host compiler ===
# nvcc rejects host compilers newer than the ones its own release was tested
# against, and distributions disagree on how those are named, so the choice
# cannot be hard-coded. Unless HOST_COMPILER is set, try nvcc's default first
# and then the versioned compilers on PATH, newest first, and keep the first
# one nvcc actually accepts. The probe preprocesses an empty translation unit,
# which is enough to trip the version check in CUDA's host_config.h.
CUDA_PROBE = nvcc $(1) -E -x cu /dev/null -o /dev/null >/dev/null 2>&1
HOST_CANDIDATES = $(foreach v,20 19 18 17 16 15 14 13 12 11 10 9 8,g++-$(v) gcc-$(v))

ifeq ($(origin HOST_COMPILER),undefined)
  HOST_COMPILER := $(shell \
    $(call CUDA_PROBE,) && echo nvcc-default || \
    for c in $(HOST_CANDIDATES); do \
      command -v $$c >/dev/null 2>&1 || continue; \
      $(call CUDA_PROBE,-ccbin $$c) && { echo $$c; break; }; \
    done)
endif

ifeq ($(HOST_COMPILER),)
  $(warning No host compiler that nvcc accepts was found. Install a supported)
  $(warning gcc/g++ and pass it with `make HOST_COMPILER=gcc-13`.)
  HOST_COMPILER := nvcc-default
endif
ifneq ($(HOST_COMPILER),nvcc-default)
  CCBIN = -ccbin $(HOST_COMPILER)
endif

NVCC = nvcc $(CCBIN)
NVCCFLAGS = -O3 -arch=$(ARCH) -lineinfo
LIBRARIES = -lm -lstdc++

# === File structure ===
SRC_DIR = src
SRCS_CU = $(wildcard $(SRC_DIR)/*.cu)
OBJS_CU = $(SRCS_CU:.cu=.o)
DEPS    = $(SRC_DIR)/headers.h $(SRC_DIR)/ranNumbers.h

TARGET = tissueFlu

# === Build rules ===
all: $(TARGET)

$(TARGET): $(OBJS_CU)
	$(NVCC) $(NVCCFLAGS) $^ -o $@ $(LIBRARIES)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cu $(DEPS)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# === Clean ===
clean:
	rm -f $(SRC_DIR)/*.o $(TARGET)

# === Run ===
# The tissue geometry is not part of this repository: point STRUCTURE at a CSV
# with a header line and then x,y,z per cell.
CONF      ?= lung.conf
STRUCTURE ?= structure.csv
OUTPUT    ?= results

run: $(TARGET)
	@mkdir -p $(OUTPUT)
	./$(TARGET) --config $(CONF) --structure $(STRUCTURE) --output $(OUTPUT)

help:
	@echo "ARCH          = $(ARCH)"
	@echo "HOST_COMPILER = $(HOST_COMPILER)   (nvcc-default = no -ccbin)"
	@echo "NVCC          = $(NVCC) $(NVCCFLAGS)"
	@echo "CONF          = $(CONF)"
	@echo "STRUCTURE     = $(STRUCTURE)"
	@echo "OUTPUT        = $(OUTPUT)"

.PHONY: all clean run help
