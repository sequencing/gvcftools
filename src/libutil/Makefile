#
# IDS must be defined:
#
IDS := BOOST_INCLUDE_DIR XFLAGS LIBUTIL

assert-defined-indirect = $(if $($1),,$(error Variable '$1' must be defined))
$(foreach I,$(IDS), $(call assert-defined-indirect,$I))


CXXFLAGS = $(XFLAGS) -I$(BOOST_INCLUDE_DIR)

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:%.cpp=%.o)

$(LIBUTIL): $(OBJS)
	$(AR) -csru $@ $(OBJS)

clean:
	rm -f $(LIBUTIL) $(OBJS)

