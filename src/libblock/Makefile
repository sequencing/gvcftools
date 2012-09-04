#
# IDS must be defined:
#
IDS := BOOST_INCLUDE_DIR TABIX_INCLUDE_DIR XFLAGS LIBBLOCK LIBUTIL_DIR

assert-defined-indirect = $(if $($1),,$(error Variable '$1' must be defined))
$(foreach I,$(IDS), $(call assert-defined-indirect,$I))

CXXFLAGS = $(XFLAGS) -I$(LIBUTIL_DIR) -I$(BOOST_INCLUDE_DIR) -I$(TABIX_INCLUDE_DIR)

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:%.cpp=%.o)

$(LIBBLOCK): $(OBJS)
	$(AR) -csru $@ $(OBJS)

clean:
	rm -f $(LIBBLOCK) $(OBJS)

