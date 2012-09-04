#
# src dir contains gvcftools code
#
# redist dir contains boost, tabix and samtools dependencies
#

SRC_DIR := $(CURDIR)/src
export BIN_DIR := $(CURDIR)/bin
REDIST_DIR := $(CURDIR)/redist
export TABIX_ROOT := $(REDIST_DIR)/tabix
export BOOST_ROOT := $(REDIST_DIR)/boost/stage

all: install

build:
	$(MAKE) -C $(REDIST_DIR) && \
	$(MAKE) -C $(SRC_DIR)

install: build
	mkdir -p $(BIN_DIR) && $(MAKE) -C $(SRC_DIR) $@ 

clean: srcclean
	$(MAKE) -C $(REDIST_DIR) clean


###### developer targets

# Cleans only src but leaves redist in place:
srcclean:
	$(MAKE) -C $(SRC_DIR) clean
	rm -rf $(BIN_DIR)

# Create emacs tag files
etags:
	cd src && find . -type f -iname "*.[ch]" -or -iname "*.cpp" -or -iname "*.hh" | sed "s/\.\///" | etags -
