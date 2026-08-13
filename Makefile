ifeq ($(origin CC),default)
	CC := gcc
endif

AR ?= ar

TARGET := camil
SRC_DIR := src

# dependencies
LCP_DIR := deps/lcptools
LCP_PREFIX := deps/lcptools
LCP_INC := $(LCP_PREFIX)/include
LCP_LIB := $(LCP_PREFIX)/lib/liblcptools.a
KLIB_DIR := deps/klib

# why not decide like this :)
MACRO_VAL := $(shell gcc -dM -E -x c /dev/null 2>/dev/null | grep __STDC_VERSION__ | awk '{print $$3}')
$(shell echo $$MACRO_VAL)
ifeq ($(MACRO_VAL),202311L)
    CFLAGS += -std=c23
else ifeq ($(MACRO_VAL),201710L)
    CFLAGS += -std=c17
else ifeq ($(MACRO_VAL),201112L)
    CFLAGS += -std=c11
else
    CFLAGS += -std=c99
endif

CFLAGS += -Wall -Wextra -Wpedantic -Wshadow
CPPFLAGS += -I$(LCP_INC) -I$(KLIB_DIR)
LDLIBS += $(LCP_LIB) -lz -lpthread -lm

ifdef DEBUG
	CFLAGS += -Og -g3 -fno-omit-frame-pointer -fsanitize=address,undefined
	LDFLAGS += -fsanitize=address,undefined
else
	CFLAGS += -O3 -flto=auto
	LDFLAGS += -flto=auto
endif

ifdef NATIVE
	CFLAGS += -march=native
endif

SRCS := $(wildcard $(SRC_DIR)/*.c)
OBJS := $(patsubst $(SRC_DIR)/%.c,$(SRC_DIR)/%.o,$(SRCS))

.PHONY: all install debug clean deep-clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)
	@rm -f $(OBJS)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

install:
	@if [ ! -f $(LCP_DIR)/Makefile ]; then \
		echo "initializing the lcptools submodule"; \
		git submodule update --init --recursive; \
	fi
	@mkdir -p $(LCP_PREFIX)
	$(MAKE) -C $(LCP_DIR) install PREFIX=$(abspath $(LCP_PREFIX))

debug:
	@$(MAKE) DEBUG=1

clean:
	rm -rf $(TARGET) $(OBJS)

deep-clean: clean
	rm -rf $(LCP_PREFIX)
	$(MAKE) -C $(LCP_DIR) clean PREFIX=$(abspath $(LCP_PREFIX)) 2>/dev/null || true