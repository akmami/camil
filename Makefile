CC ?= gcc
AR ?= ar

TARGET := camil
SRC_DIR := src

# dependencies
LCP_DIR := deps/lcptools
LCP_PREFIX := deps/lcptools/
LCP_INC := $(LCP_PREFIX)/include
LCP_LIB := $(LCP_PREFIX)/lib/liblcptools.a
KLIB_DIR := deps/klib

CFLAGS := -O3 -Wall -Wextra -Wpedantic -Wshadow
CPPFLAGS += -std=c99 -I$(LCP_INC) -I$(KLIB_DIR)
LDLIBS += $(LCP_LIB) -lz -lpthread -lm

SRCS := $(wildcard $(SRC_DIR)/*.c)
OBJS := $(patsubst $(SRC_DIR)/%.c,$(SRC_DIR)/%.o,$(SRCS))

.PHONY: all install debug clean deep-clean

all: $(TARGET)

$(TARGET): $(OBJS) 
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDLIBS)
	@rm -f $(OBJS)

install: 
	@if [ ! -f $(LCP_DIR)/Makefile ]; then \
		echo "initializing the lcptools submodule"; \
		git submodule update --init --recursive; \
	fi
	@mkdir -p $(LCP_PREFIX)
	$(MAKE) -C $(LCP_DIR) install PREFIX=$(abspath $(LCP_PREFIX))

debug:
	$(MAKE) CFLAGS="-O0 -g -Wall -Wextra -Wpedantic -Wshadow -fsanitize=address,undefined" LDLIBS="$(LCP_LIB) -lz -lpthread -lm -fsanitize=address,undefined" all

clean:
	rm -rf $(TARGET)

deep-clean: clean
	rm -rf $(LCP_PREFIX)
	$(MAKE) -C $(LCP_DIR) clean PREFIX=$(abspath $(LCP_PREFIX)) 2>/dev/null || true