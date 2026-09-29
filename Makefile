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

# lcptools build variant. The core key is the 64 bit label lcptools computes,
# so LABEL must be 64; the key derivation and the reverse complement parse
# need positions and the DNA alphabet. Pinned here so that an upstream change
# of the defaults cannot silently pair camil with a differently shaped
# struct core.
LCP_VARIANT := LABEL=64 POS=32 DCT=1 CORE=var ALPHABET=dna

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

.PHONY: all install debug test clean deep-clean

all: $(TARGET)

$(TARGET): $(LCP_LIB) $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)
	@rm -f $(OBJS)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c $(LCP_LIB)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# lcptools is built into its private prefix the first time it is needed
$(LCP_LIB):
	$(MAKE) install

install:
	@if [ ! -f $(LCP_DIR)/Makefile ]; then \
		echo "initializing the lcptools submodule"; \
		git submodule update --init --recursive; \
	fi
	@mkdir -p $(LCP_PREFIX)
	$(MAKE) -C $(LCP_DIR) install PREFIX=$(abspath $(LCP_PREFIX)) $(LCP_VARIANT)

debug:
	@$(MAKE) DEBUG=1

test: $(TARGET)
	sh tests/run_tests.sh

clean:
	rm -rf $(TARGET) $(OBJS)

deep-clean: clean
	rm -rf $(LCP_PREFIX)
	$(MAKE) -C $(LCP_DIR) clean PREFIX=$(abspath $(LCP_PREFIX)) $(LCP_VARIANT) 2>/dev/null || true