#
# REQS must be defined. This is normally handled by the parent make process.
#

REQS := BIN_DIR BOOST_ROOT TABIX_ROOT

assert-defined-indirect = $(if $($1),,$(error Variable '$1' must be defined))
$(foreach I,$(REQS), $(call assert-defined-indirect,$I))

BOOST_INCLUDE_DIR = $(BOOST_ROOT)/include
BOOST_LIB_DIR = $(BOOST_ROOT)/lib

TABIX_INCLUDE_DIR = $(TABIX_ROOT)
TABIX_LIB_DIR = $(TABIX_ROOT)

# debug:
#export XFLAGS= -g3 -fno-inline -O0 -Wall -Wextra
XFLAGS= -O3 -Wall -Wextra

LIBTRIO_DIR = $(CURDIR)/libtrio
LIBBLOCK_DIR = $(CURDIR)/libblock
export LIBUTIL_DIR = $(CURDIR)/libutil

CC := $(CXX) # this gets the linker to include libc++

export CXXFLAGS = $(XFLAGS) -I$(LIBUTIL_DIR) -I$(BOOST_INCLUDE_DIR) -I$(TABIX_INCLUDE_DIR)
export LDFLAGS = -L$(BOOST_LIB_DIR) -L$(TABIX_LIB_DIR)
export LDLIBS = -ltabix_and_faidx -lz -lm

PROGS = $(TRIOPROGS) $(BLOCKPROGS) 
TRIOPROGS := trio twins
BLOCKPROGS := break_blocks extract_variants gatk_to_gvcf set_haploid_region
OBJS = $(PROGS:%=%.o)

export LIBTRIO := libtrio.a
LIBTRIO_PATH = $(LIBTRIO_DIR)/$(LIBTRIO)
export LIBBLOCK := libblock.a
LIBBLOCK_PATH = $(LIBBLOCK_DIR)/$(LIBBLOCK)
export LIBUTIL := libutil.a
LIBUTIL_PATH = $(LIBUTIL_DIR)/$(LIBUTIL)


.PHONY: all build clean install test


all: build

build: $(PROGS)


$(TRIOPROGS): %: %.o $(LIBTRIO_PATH) $(LIBUTIL_PATH)
$(TRIOPROGS): CXXFLAGS += -I$(LIBUTIL_DIR) -I$(LIBTRIO_DIR) 
$(TRIOPROGS): LDLIBS += -lboost_program_options

$(BLOCKPROGS): %: %.o $(LIBBLOCK_PATH) $(LIBUTIL_PATH)
$(BLOCKPROGS): CXXFLAGS += -I$(LIBUTIL_DIR) -I$(LIBBLOCK_DIR) 
$(BLOCKPROGS): LDLIBS += -lboost_program_options

$(LIBTRIO_PATH):
	$(MAKE) -C $(LIBTRIO_DIR) $(LIBTRIO)

$(LIBBLOCK_PATH):
	$(MAKE) -C $(LIBBLOCK_DIR) $(LIBBLOCK)

$(LIBUTIL_PATH):
	$(MAKE) -C $(LIBUTIL_DIR) $(LIBUTIL)

include devel.mk

install: build
	cp $(PROGS) *.pl $(BIN_DIR)

test:
	$(MAKE) -C $(LIBUTIL_DIR) $@
	#$(MAKE) -C $(LIBBLOCK_DIR) $@
	#$(MAKE) -C $(LIBUTIL_DIR) $@

clean:
	$(MAKE) -C $(LIBTRIO_DIR) $@
	$(MAKE) -C $(LIBBLOCK_DIR) $@
	$(MAKE) -C $(LIBUTIL_DIR) $@
	rm -f $(PROGS) $(OBJS) $(VERSION)
